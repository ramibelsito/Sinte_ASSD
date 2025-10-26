// tp2_allinone.cpp
// Trabajo práctico 2024 - Implementación todo en 1 archivo usando Qt Widgets
// Funcionalidades: FFT (Cooley-Tukey), parser MIDI simple, síntesis aditiva, FM,
// Karplus-Strong, sample playback (WAV), delay/reverb simples, espectrograma, export WAV.
// Compilar con Qt (widgets, multimedia). C++17.

#include <QtWidgets>
#include <QtMultimedia>
#include <complex>
#include <vector>
#include <cmath>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <cstdint>
#include <memory>
#include <functional>

using namespace std;
using cfloat = complex<float>;

// --------------------------- Utilidades de audio / WAV -------------------------
struct WAVHeader {
    char riff[4]; uint32_t chunkSize; char wave[4];
    char fmt[4]; uint32_t subchunk1; uint16_t audioFormat; uint16_t numChannels;
    uint32_t sampleRate; uint32_t byteRate; uint16_t blockAlign; uint16_t bitsPerSample;
    char dataTag[4]; uint32_t dataSize;
};

bool writeWAV(const QString &filename, const vector<float>& samples, uint32_t sampleRate) {
    WAVHeader h;
    memcpy(h.riff, "RIFF",4);
    memcpy(h.wave,"WAVE",4);
    memcpy(h.fmt,"fmt ",4);
    h.subchunk1 = 16;
    h.audioFormat = 1; // PCM
    h.numChannels = 1;
    h.sampleRate = sampleRate;
    h.bitsPerSample = 16;
    h.blockAlign = h.numChannels * h.bitsPerSample/8;
    h.byteRate = h.sampleRate * h.blockAlign;
    memcpy(h.dataTag,"data",4);
    h.dataSize = samples.size() * h.blockAlign;
    h.chunkSize = 36 + h.dataSize;

    ofstream f(filename.toStdString(), ios::binary);
    if(!f) return false;
    f.write((char*)&h, sizeof(h));
    for(float s: samples){
        float cl = max(-1.0f, min(1.0f, s));
        int16_t v = (int16_t) (cl * 32767.0f);
        f.write((char*)&v, sizeof(v));
    }
    f.close();
    return true;
}

bool readWAV(const QString &filename, vector<float>& out, uint32_t& sampleRate) {
    ifstream f(filename.toStdString(), ios::binary);
    if(!f) return false;
    WAVHeader h;
    f.read((char*)&h, sizeof(h));
    if(strncmp(h.riff,"RIFF",4)!=0) return false;
    if(h.audioFormat!=1) return false;
    sampleRate = h.sampleRate;
    size_t samples = h.dataSize / (h.bitsPerSample/8);
    out.resize(samples);
    for(size_t i=0;i<samples;i++){
        if(h.bitsPerSample==16){
            int16_t v; f.read((char*)&v,sizeof(v));
            out[i] = v/32768.0f;
        } else if(h.bitsPerSample==8){
            uint8_t v; f.read((char*)&v,sizeof(v));
            out[i] = (v-128)/128.0f;
        } else return false;
    }
    return true;
}

// --------------------------- FFT Cooley-Tukey (radix-2) -------------------------
void fft_recursive(const vector<cfloat>& in, vector<cfloat>& out, size_t n, size_t stride) {
    if (n==1) {
        out[0] = in[0];
        return;
    }
    size_t m = n/2;
    vector<cfloat> even(m), odd(m);
    for(size_t i=0;i<m;i++){
        even[i] = in[i*2*stride];
        odd[i]  = in[(i*2+1)*stride];
    }
    vector<cfloat> Fe(m), Fo(m);
    fft_recursive(even, Fe, m, 1);
    fft_recursive(odd,  Fo, m, 1);
    for(size_t k=0;k<m;k++){
        cfloat t = std::polar(1.0f, float(-2.0f*M_PI*k/n)) * Fo[k];
        out[k] = Fe[k] + t;
        out[k+m] = Fe[k] - t;
    }
}

// Public API: in -> out, n must be power of two
void fft(complex<float>* in, complex<float>* out, size_t n) {
    // Convert to vector and call recursive (simple)
    vector<cfloat> vin(n);
    for(size_t i=0;i<n;i++) vin[i] = in[i];
    vector<cfloat> vout(n);
    // iterative bit-reversal + iterative would be more efficient; recursive here for clarity
    // We'll implement an iterative radix-2 decimation-in-time for robustness
    size_t N = n;
    // bit-reversal permutation
    size_t j=0;
    for(size_t i=0;i<N;i++){
        if(i<j) vout[j] = vin[i], vout[i] = vin[j];
        else vout[i] = vin[i];
        size_t m = N/2;
        while(m>=1 && j>=m){ j-=m; m/=2; }
        j+=m;
    }
    // Cooley-Tukey iterative
    for(size_t len=2; len<=N; len<<=1){
        float ang = -2.0f * M_PI / len;
        cfloat wlen = polar(1.0f, ang);
        for(size_t i=0;i<N;i+=len){
            cfloat w = 1;
            for(size_t k=0;k<len/2;k++){
                cfloat u = vout[i+k];
                cfloat v = vout[i+k+len/2] * w;
                vout[i+k] = u + v;
                vout[i+k+len/2] = u - v;
                w *= wlen;
            }
        }
    }
    for(size_t i=0;i<N;i++) out[i] = vout[i];
}

// Helper: next power of two
size_t nextPow2(size_t v){
    size_t p=1;
    while(p<v) p<<=1;
    return p;
}

// Window function (Hann)
void applyHann(vector<float>& buf){
    size_t N = buf.size();
    for(size_t n=0;n<N;n++){
        float w = 0.5f*(1.0f - cos(2.0f*M_PI*n/(N-1)));
        buf[n] *= w;
    }
}

// --------------------------- MIDI parser (simple) ------------------------------
// Minimal MIDI parser to extract tracks and Note On/Off events with delta times.
// Supports format 0 and 1, common tick division, running status naive handling.
// Not a full-featured parser; extend as needed.

struct MidiEvent {
    uint32_t time; // absolute ticks
    uint8_t type; // 0x9 = note on, 0x8 = note off
    uint8_t channel;
    uint8_t note;
    uint8_t vel;
};

struct MidiTrack {
    vector<MidiEvent> events;
    string name;
};

struct MidiFile {
    uint16_t format;
    uint16_t ntrks;
    uint16_t division; // ticks per quarter
    vector<MidiTrack> tracks;
};

static uint32_t read_be32(ifstream &f){
    uint8_t b[4]; f.read((char*)b,4);
    return (b[0]<<24)|(b[1]<<16)|(b[2]<<8)|b[3];
}
static uint16_t read_be16(ifstream &f){
    uint8_t b[2]; f.read((char*)b,2);
    return (b[0]<<8)|b[1];
}
static uint8_t read8(ifstream &f){
    char c; f.get(c); return (uint8_t)c;
}
static uint32_t read_vlq(ifstream &f){
    uint32_t v=0;
    while(true){
        uint8_t b = read8(f);
        v = (v<<7) | (b & 0x7F);
        if(!(b & 0x80)) break;
    }
    return v;
}

bool loadMidi(const QString &path, MidiFile &out) {
    ifstream f(path.toStdString(), ios::binary);
    if(!f) return false;
    char id[4]; f.read(id,4);
    if(strncmp(id,"MThd",4)!=0) return false;
    uint32_t headerLen = read_be32(f);
    out.format = read_be16(f);
    out.ntrks = read_be16(f);
    out.division = read_be16(f);
    // skip extra header bytes
    if(headerLen>6) f.seekg(headerLen-6, ios::cur);
    out.tracks.clear();
    for(int t=0;t<out.ntrks;t++){
        char tid[4]; f.read(tid,4);
        if(strncmp(tid,"MTrk",4)!=0) {
            return false;
        }
        uint32_t len = read_be32(f);
        streampos trackStart = f.tellg();
        MidiTrack track;
        uint32_t absTime = 0;
        uint8_t runningStatus = 0;
        while((uint32_t)(f.tellg() - trackStart) < len){
            uint32_t dt = read_vlq(f);
            absTime += dt;
            uint8_t b = read8(f);
            if(b==0xFF){ // meta
                uint8_t meta = read8(f);
                uint32_t l = read_vlq(f);
                if(meta==0x03){ // track name
                    string name(l,' ');
                    f.read(&name[0], l);
                    track.name = name;
                } else {
                    f.seekg(l, ios::cur);
                }
            } else if(b==0xF0 || b==0xF7){ // sysex
                uint32_t l = read_vlq(f);
                f.seekg(l, ios::cur);
            } else {
                uint8_t status;
                if(b & 0x80){
                    status = b;
                    runningStatus = status;
                } else {
                    // running status: b is first data byte
                    if(runningStatus==0) return false;
                    status = runningStatus;
                    // put back one byte
                    f.seekg(-1, ios::cur);
                }
                uint8_t cmd = (status>>4)&0x0F;
                uint8_t ch = status & 0x0F;
                if(cmd==0x9 || cmd==0x8){ // note on/off
                    uint8_t note = read8(f);
                    uint8_t vel = read8(f);
                    if(cmd==0x9 && vel==0){
                        // treat as note off
                        MidiEvent e{absTime,0x8,ch,note,vel};
                        track.events.push_back(e);
                    } else {
                        MidiEvent e{absTime,(uint8_t)(cmd),ch,note,vel};
                        track.events.push_back(e);
                    }
                } else if(cmd==0xC){ // program change - 1 byte
                    uint8_t prog = read8(f);
                    // ignore for now
                } else {
                    // other channel messages: assume 2 data bytes
                    uint8_t d1 = read8(f);
                    uint8_t d2 = (cmd==0xC||cmd==0xD)?0:read8(f);
                    (void)d1; (void)d2;
                }
            }
        }
        out.tracks.push_back(track);
    }
    return true;
}

// --------------------------- Synthesis engines --------------------------------

struct ADSR {
    float A=0.01f, D=0.05f, S=0.8f, R=0.1f;
    float sampleRate=44100;
    float env(float t, float noteLen) const {
        if(t<0) return 0;
        if(t < A) return (t/A);
        else if(t < A+D) return 1.0f - (t-A)/D*(1.0f - S);
        else if(t < noteLen) return S;
        else {
            float rt = t - noteLen;
            if(rt > R) return 0;
            return S * (1.0f - rt/R);
        }
    }
};

// Instrument presets for additive synthesis
struct InstrumentPreset {
    QString name;
    vector<float> partialAmps;
    vector<ADSR> envelopes;

    static InstrumentPreset piano() {
        InstrumentPreset p;
        p.name = "Piano";
        p.partialAmps = {1.0f, 0.7f, 0.4f, 0.25f, 0.15f, 0.08f, 0.05f};
        p.envelopes.resize(p.partialAmps.size());
        for(auto& env : p.envelopes) {
            env.A = 0.005f;
            env.D = 0.1f;
            env.S = 0.7f;
            env.R = 0.3f;
        }
        return p;
    }

    static InstrumentPreset organ() {
        InstrumentPreset p;
        p.name = "Organ";
        p.partialAmps = {1.0f, 0.5f, 0.33f, 0.25f, 0.2f, 0.17f, 0.14f};
        p.envelopes.resize(p.partialAmps.size());
        for(auto& env : p.envelopes) {
            env.A = 0.02f;
            env.D = 0.1f;
            env.S = 0.9f;
            env.R = 0.1f;
        }
        return p;
    }

    static InstrumentPreset trumpet() {
        InstrumentPreset p;
        p.name = "Trumpet";
        p.partialAmps = {1.0f, 0.9f, 0.8f, 0.6f, 0.4f, 0.2f, 0.1f};
        p.envelopes.resize(p.partialAmps.size());
        for(auto& env : p.envelopes) {
            env.A = 0.05f;
            env.D = 0.1f;
            env.S = 0.8f;
            env.R = 0.1f;
        }
        return p;
    }

    static vector<InstrumentPreset> getAllPresets() {
        return {piano(), organ(), trumpet()};
    }
};

float midiNoteFreq(int note) {
    // A4 = MIDI 69 = 440Hz
    return 440.0f * pow(2.0f, (note - 69)/12.0f);
}

// Additive synth: generate note (mono) with partials and per-partial envelopes
vector<float> synthAdditive(int midiNote, float dur, float sr, const vector<float>& partialAmps, const vector<ADSR>& partialEnvs){
    float f0 = midiNoteFreq(midiNote);
    size_t N = (size_t)round(dur*sr);
    vector<float> out(N,0.0f);
    for(size_t n=0;n<N;n++){
        float t = float(n)/sr;
        float s=0;
        for(size_t p=0;p<partialAmps.size();p++){
            float amp = partialAmps[p];
            float phase = 2.0f*M_PI * f0 * (p+1) * t;
            float env = (p < partialEnvs.size()) ? partialEnvs[p].env(t, dur) : 1.0f;
            s += amp * env * sin(phase);
        }
        out[n]=s;
    }
    return out;
}

// FM synth: carrier + modulator simple implementation
vector<float> synthFM(int midiNote, float dur, float sr, float carrierRatio, float modRatio, float modIndex, const ADSR& env){
    float fc = midiNoteFreq(midiNote) * carrierRatio;
    float fm = midiNoteFreq(midiNote) * modRatio;
    size_t N = (size_t)round(dur*sr);
    vector<float> out(N);
    for(size_t n=0;n<N;n++){
        float t = float(n)/sr;
        float a = env.env(t,dur);
        float mod = modIndex * sin(2.0f*M_PI*fm*t);
        float s = a * sin(2.0f*M_PI*fc*t + mod);
        out[n]=s;
    }
    return out;
}

// Karplus-Strong pluck (monophonic)
vector<float> karplusStrong(int midiNote, float dur, float sr, float decay = 0.996f){
    float f = midiNoteFreq(midiNote);
    int L = max(2, (int)round(sr / f));
    int N = (int)round(dur*sr);
    vector<float> buffer(L);
    // initial noise
    for(int i=0;i<L;i++) buffer[i] = ((float)rand() / RAND_MAX)*2.0f - 1.0f;
    vector<float> out(N);
    int idx = 0;
    for(int n=0;n<N;n++){
        float s = buffer[idx];
        out[n] = s;
        float next = decay * 0.5f * (s + buffer[(idx+1)%L]);
        buffer[idx] = next;
        idx = (idx+1)%L;
    }
    return out;
}

// Simple delay (feedback)
void applyDelay(vector<float>& buf, float sr, float delaySec, float fb, float mix){
    int delaySamples = (int)round(delaySec*sr);
    if(delaySamples<=0) return;
    vector<float> ring(delaySamples,0.0f);
    int idx=0;
    for(size_t i=0;i<buf.size();i++){
        float in = buf[i];
        float d = ring[idx];
        float out = in*(1.0f-mix) + d*mix;
        ring[idx] = in + d*fb;
        buf[i] = out;
        idx = (idx+1)%delaySamples;
    }
}

// Flanger effect implementation
void applyFlanger(vector<float>& buf, float sr, float rate = 1.0f, float depth = 0.001f, float feedback = 0.5f, float mix = 0.5f) {
    // Parameters:
    // rate: LFO frequency in Hz (how fast the delay changes)
    // depth: maximum delay time in seconds
    // feedback: amount of feedback (0-1)
    // mix: wet/dry mix (0-1)

    // Calculate parameters
    float maxDelaySamples = depth * sr;
    vector<float> delayLine(size_t(maxDelaySamples) + 1, 0.0f);
    vector<float> dry = buf; // Store original signal

    // LFO for delay modulation
    float phase = 0.0f;
    float phaseInc = 2.0f * M_PI * rate / sr;

    size_t writePos = 0;
    for(size_t i = 0; i < buf.size(); i++) {
        // Calculate current delay from LFO
        float lfo = (sin(phase) + 1.0f) * 0.5f; // 0 to 1
        float currentDelay = lfo * maxDelaySamples;
        
        // Calculate read position
        float readPos = float(writePos) - currentDelay;
        if(readPos < 0.0f) readPos += delayLine.size();
        
        // Linear interpolation
        size_t readPos1 = size_t(readPos);
        size_t readPos2 = (readPos1 + 1) % delayLine.size();
        float frac = readPos - float(readPos1);
        float delayed = delayLine[readPos1] * (1.0f - frac) + delayLine[readPos2] * frac;
        
        // Apply feedback and mix
        float input = buf[i];
        delayLine[writePos] = input + delayed * feedback;
        buf[i] = input * (1.0f - mix) + delayed * mix;
        
        // Update positions
        writePos = (writePos + 1) % delayLine.size();
        phase += phaseInc;
        if(phase >= 2.0f * M_PI) phase -= 2.0f * M_PI;
    }
}

// Simple Schroeder reverberator (comb + allpass minimal)
void applyReverb(vector<float>& buf, float sr, float roomSize = 1.0f, float damping = 0.5f, float wetLevel = 0.5f){
    vector<int> combMs = {29, 37, 41}; // ms-ish values
    // Adjust delay times based on room size
    for(auto& ms : combMs) {
        ms = int(ms * roomSize);
    }
    
    // Store original signal for wet/dry mix
    vector<float> dry = buf;
    
    for(int ms: combMs){
        int L = int(sr * (ms/1000.0f));
        if(L<=0) continue;
        vector<float> ring(L,0.0f);
        int idx=0;
        float fb = 0.78f * (1.0f - damping * 0.3f); // Damping affects feedback
        for(size_t i=0;i<buf.size();i++){
            float in = buf[i];
            float d = ring[idx];
            float out = in + d*0.5f;
            ring[idx] = in + d*fb;
            buf[i] = out;
            idx = (idx+1)%L;
        }
    }
    
    // Apply wet/dry mix
    for(size_t i=0; i<buf.size(); i++) {
        buf[i] = dry[i] * (1.0f - wetLevel) + buf[i] * wetLevel;
    }
}

// --------------------------- Mixer: render MIDI to audio -----------------------
struct TrackSynthConfig {
    enum Type {ADDITIVE, FM, KARPLUS, SAMPLE} type = ADDITIVE;
    // additive params
    vector<float> partialAmps = {1.0f, 0.5f, 0.25f, 0.1f};
    vector<ADSR> partialEnvs;
    // FM params
    float carrierRatio=1.0f, modRatio=2.0f, modIndex=5.0f;
    ADSR env;
    // Karplus params
    float ksDecay = 0.996f;
    // sample
    vector<float> sample;
    uint32_t sampleRate=44100;
};

// Callback function type for progress reporting
using ProgressCallback = std::function<void(float)>;

vector<float> renderMidi(const MidiFile& mf, const vector<TrackSynthConfig>& configs, 
                        float secondsPerBeat=0.5f, uint32_t sampleRate=44100,
                        ProgressCallback progressCb = nullptr){
    // Convert events to absolute time in seconds using division (ticks per quarter) and secondsPerBeat
    float tickSec = secondsPerBeat / mf.division;
    // estimate total length
    uint32_t maxTick=0;
    for(const auto& tr: mf.tracks){
        if(!tr.events.empty()){
            maxTick = max(maxTick, tr.events.back().time);
        }
    }
    size_t N = (size_t)ceil((maxTick * tickSec + 1.0f) * sampleRate);
    vector<float> mix(N,0.0f);
    
    // Count total notes for progress tracking
    size_t totalNotes = 0;
    size_t processedNotes = 0;
    for(const auto& tr : mf.tracks) {
        for(const auto& ev : tr.events) {
            if(ev.type == 0x9 && ev.vel > 0) totalNotes++;
        }
    }

    for(size_t ti=0; ti<mf.tracks.size(); ++ti){
        const auto& tr = mf.tracks[ti];
        if(ti >= configs.size()) break;
        const auto& cfg = configs[ti];
        // Build list of note-on -> note-off pairs (very naive: search next off)
        for(size_t eidx=0; eidx<tr.events.size(); ++eidx){
            const auto& ev = tr.events[eidx];
            if(ev.type==0x9 && ev.vel>0){
                uint32_t startTick = ev.time;
                // find matching note-off
                uint32_t endTick = startTick + mf.division; // default 1 quarter if not found
                for(size_t j=eidx+1;j<tr.events.size();++j){
                    if((tr.events[j].note==ev.note) && (tr.events[j].type==0x8 || (tr.events[j].type==0x9 && tr.events[j].vel==0))){
                        endTick = tr.events[j].time;
                        break;
                    }
                }
                float startSec = startTick * tickSec;
                float endSec = endTick * tickSec;
                float dur = max(0.01f, endSec - startSec);
                // generate audio for this note
                vector<float> noteBuf;
                if(cfg.type == TrackSynthConfig::ADDITIVE){
                    noteBuf = synthAdditive(ev.note, dur, sampleRate, cfg.partialAmps, cfg.partialEnvs);
                } else if(cfg.type == TrackSynthConfig::FM){
                    noteBuf = synthFM(ev.note, dur, sampleRate, cfg.carrierRatio, cfg.modRatio, cfg.modIndex, cfg.env);
                } else if(cfg.type == TrackSynthConfig::KARPLUS){
                    noteBuf = karplusStrong(ev.note, dur, sampleRate, cfg.ksDecay);
                } else { // SAMPLE
                    // naive: copy sample, or stretch if needed (no resampling implemented)
                    if(cfg.sample.empty()) continue;
                    noteBuf = cfg.sample;
                }

                // Update progress after each note
                processedNotes++;
                if(progressCb) {
                    float progress = float(processedNotes) / totalNotes;
                    progressCb(progress);
                }
                // apply velocity scaling
                float velGain = ev.vel / 127.0f;
                for(size_t n=0;n<noteBuf.size();++n){
                    size_t pos = (size_t)round((startSec*sampleRate)) + n;
                    if(pos < mix.size()) mix[pos] += noteBuf[n] * velGain;
                }
            }
        }
    }
    // Apply soft clipping normalize
    float maxv=0;
    for(float s:mix) maxv = max(maxv, fabs(s));
    if(maxv>0.95f){
        float g = 0.95f / maxv;
        for(auto &x:mix) x *= g;
    }
    return mix;
}

// --------------------------- Spectrogram rendering (simple) --------------------
QImage renderSpectrogram(const vector<float>& audio, uint32_t sr, int winSize=1024, int hop=512){
    int Nwin = winSize;
    int Nhop = hop;
    size_t frames = (audio.size() < Nwin) ? 1 : 1 + (audio.size()-Nwin)/Nhop;
    int height = Nwin/2;
    QImage img(frames, height, QImage::Format_RGB32);
    img.fill(Qt::black);
    vector<complex<float>> in(Nwin), out(Nwin);
    vector<float> window(Nwin);
    // Hann
    for(int i=0;i<Nwin;i++) window[i] = 0.5f*(1 - cos(2*M_PI*i/(Nwin-1)));
    for(size_t f=0; f<frames; ++f){
        size_t start = f * Nhop;
        for(int n=0;n<Nwin;n++){
            float s = (start + n < audio.size()) ? audio[start+n] : 0.0f;
            in[n] = complex<float>(s * window[n], 0.0f);
        }
        // perform fft via our API
        vector<complex<float>> tin(Nwin), tout(Nwin);
        for(int i=0;i<Nwin;i++) tin[i]=in[i];
        fft(&tin[0], &tout[0], Nwin);
        for(int k=0;k<height;k++){
            float mag = abs(tout[k]) / Nwin;
            float db = 20*log10f(max(1e-6f, mag));
            // map db range -100..0 to 0..255
            float v = (db + 100.0f) / 100.0f;
            int iv = max(0, min(255, (int)round(v*255)));
            QColor c = QColor::fromHsv((int)(iv*0.7), 255, iv);
            img.setPixel(f, height-1-k, c.rgb());
        }
    }
    return img.mirrored();
}

// --------------------------- Qt GUI -------------------------------------------
// FFT Visualization widget
class FFTWidget : public QWidget {
    Q_OBJECT
public:
    FFTWidget(QWidget* parent = nullptr) : QWidget(parent) {
        setMinimumSize(200, 100);
        setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    }

    void updateFFT(const vector<float>& magnitudes) {
        fftData = magnitudes;
        update();
    }

protected:
    void paintEvent(QPaintEvent*) override {
        QPainter p(this);
        p.fillRect(rect(), Qt::black);
        
        if(fftData.empty()) return;
        
        p.setPen(Qt::green);
        float w = width() / float(fftData.size());
        float h = height();
        
        for(size_t i = 0; i < fftData.size(); i++) {
            float mag = fftData[i];
            // Convert to dB and normalize
            float db = 20 * log10(max(1e-6f, mag));
            float normDb = (db + 100.0f) / 100.0f; // -100dB to 0dB range
            float barHeight = normDb * h;
            p.drawLine(QPointF(i*w, h), QPointF(i*w, h - barHeight));
        }
    }

private:
    vector<float> fftData;
};

class WAVPlayer : public QObject {
    Q_OBJECT
public:
    WAVPlayer(QObject* parent = nullptr) : QObject(parent) {
        QAudioFormat format;
        format.setSampleRate(44100); // default, will be updated when loading file
        format.setChannelCount(1);
        format.setSampleSize(16);
        format.setCodec("audio/pcm");
        format.setByteOrder(QAudioFormat::LittleEndian);
        format.setSampleType(QAudioFormat::SignedInt);

        QAudioDeviceInfo info(QAudioDeviceInfo::defaultOutputDevice());
        if (!info.isFormatSupported(format)) {
            format = info.nearestFormat(format);
        }

        audio = new QAudioOutput(format, this);
        buffer = new QBuffer(this);
        connect(audio, &QAudioOutput::stateChanged, this, &WAVPlayer::handleStateChange);

        // Setup timer for FFT updates
        fftTimer = new QTimer(this);
        connect(fftTimer, &QTimer::timeout, this, &WAVPlayer::updateFFT);
        fftTimer->setInterval(50); // 20 fps
    }

    bool loadFile(const QString& path) {
        vector<float> samples;
        uint32_t sr;
        if(!readWAV(path, samples, sr)) return false;
        
        // Convert to PCM
        QByteArray pcm;
        pcm.resize(samples.size() * 2);
        int16_t* dst = (int16_t*)pcm.data();
        for(size_t i = 0; i < samples.size(); i++) {
            float v = max(-1.0f, min(1.0f, samples[i]));
            dst[i] = int16_t(v * 32767);
        }
        
        buffer->setData(pcm);
        sampleRate = sr;

        // Update audio format with new sample rate
        QAudioFormat format = audio->format();
        format.setSampleRate(sampleRate);
        
        QAudioDeviceInfo info(QAudioDeviceInfo::defaultOutputDevice());
        if (!info.isFormatSupported(format)) {
            format = info.nearestFormat(format);
        }

        delete audio;
        audio = new QAudioOutput(format, this);
        connect(audio, &QAudioOutput::stateChanged, this, &WAVPlayer::handleStateChange);
        
        return true;
    }

    void play() {
        if(!buffer->isOpen())
            buffer->open(QIODevice::ReadOnly);
        else
            buffer->seek(0);

        audio->stop();
        audio->start(buffer);
        fftTimer->start();
    }

    void stop() {
        audio->stop();
        buffer->seek(0);
        fftTimer->stop();
    }

signals:
    void fftUpdated(const vector<float>& magnitudes);

private slots:
    void handleStateChange(QAudio::State state) {
        if(state == QAudio::IdleState || state == QAudio::StoppedState) {
            fftTimer->stop();
            buffer->seek(0);
        }
    }

    void updateFFT() {
        // Get current playing position
        qint64 pos = buffer->pos();
        int samplesPerFFT = 1024;
        
        // Read samples for FFT
        vector<complex<float>> fftBuf(samplesPerFFT);
        QByteArray raw = buffer->read(samplesPerFFT * 2); // 16-bit samples
        const int16_t* samples = (const int16_t*)raw.constData();
        int numSamples = raw.size() / 2;
        
        // Fill FFT buffer
        for(int i = 0; i < samplesPerFFT; i++) {
            float sample = (i < numSamples) ? samples[i] / 32768.0f : 0.0f;
            fftBuf[i] = complex<float>(sample, 0);
        }
        
        // Apply window
        for(int i = 0; i < samplesPerFFT; i++) {
            float w = 0.5f * (1.0f - cos(2.0f * M_PI * i / (samplesPerFFT-1))); // Hann window
            fftBuf[i] *= w;
        }
        
        // Perform FFT
        vector<complex<float>> fftOut(samplesPerFFT);
        fft(&fftBuf[0], &fftOut[0], samplesPerFFT);
        
        // Calculate magnitudes (only first half, rest is mirror)
        vector<float> mags(samplesPerFFT/2);
        for(int i = 0; i < samplesPerFFT/2; i++) {
            mags[i] = abs(fftOut[i]) / samplesPerFFT;
        }
        
        emit fftUpdated(mags);
        
        // Seek back if needed
        if(pos >= buffer->size())
            buffer->seek(0);
    }

private:
    QAudioOutput* audio;
    QBuffer* buffer;
    QTimer* fftTimer;
    int sampleRate = 44100;
};

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    MainWindow(){
        sampleRate = 44100;
        central = new QWidget;
        setCentralWidget(central);
        QVBoxLayout *l = new QVBoxLayout(central);

        // Top controls
        QHBoxLayout *top = new QHBoxLayout;
        btnLoadMidi = new QPushButton("Load MIDI");
        lblMidi = new QLabel("No MIDI loaded");
        btnRender = new QPushButton("Render WAV");
        btnSave = new QPushButton("Save WAV");
        top->addWidget(btnLoadMidi);
        top->addWidget(lblMidi);
        top->addWidget(btnRender);
        top->addWidget(btnSave);
        l->addLayout(top);

        // Effects panel
        QGroupBox* effectsGroup = new QGroupBox("Effects");
        QFormLayout* effectsLayout = new QFormLayout(effectsGroup);
        
        // Delay controls
        chkDelay = new QCheckBox("Enable");
        sldDelayTime = new QSlider(Qt::Horizontal);
        sldDelayTime->setRange(1, 1000); // 1-1000ms
        sldDelayTime->setValue(200);
        sldDelayFeedback = new QSlider(Qt::Horizontal);
        sldDelayFeedback->setRange(0, 100); // 0-100%
        sldDelayFeedback->setValue(30);
        sldDelayMix = new QSlider(Qt::Horizontal);
        sldDelayMix->setRange(0, 100); // 0-100%
        sldDelayMix->setValue(50);

        QHBoxLayout* delayLayout = new QHBoxLayout;
        delayLayout->addWidget(chkDelay);
        QFormLayout* delayParamsLayout = new QFormLayout;
        delayParamsLayout->addRow("Time (ms):", sldDelayTime);
        delayParamsLayout->addRow("Feedback (%):", sldDelayFeedback);
        delayParamsLayout->addRow("Mix (%):", sldDelayMix);
        delayLayout->addLayout(delayParamsLayout);
        effectsLayout->addRow("Delay:", delayLayout);

        // Reverb controls
        chkReverb = new QCheckBox("Enable");
        sldReverbRoom = new QSlider(Qt::Horizontal);
        sldReverbRoom->setRange(50, 200); // 50-200% room size
        sldReverbRoom->setValue(100);
        sldReverbDamping = new QSlider(Qt::Horizontal);
        sldReverbDamping->setRange(0, 100); // 0-100% damping
        sldReverbDamping->setValue(50);
        sldReverbWet = new QSlider(Qt::Horizontal);
        sldReverbWet->setRange(0, 100); // 0-100% wet level
        sldReverbWet->setValue(50);

        QHBoxLayout* reverbLayout = new QHBoxLayout;
        reverbLayout->addWidget(chkReverb);
        QFormLayout* reverbParamsLayout = new QFormLayout;
        reverbParamsLayout->addRow("Room Size (%):", sldReverbRoom);
        reverbParamsLayout->addRow("Damping (%):", sldReverbDamping);
        reverbParamsLayout->addRow("Wet Level (%):", sldReverbWet);
        reverbLayout->addLayout(reverbParamsLayout);
        effectsLayout->addRow("Reverb:", reverbLayout);

        // Flanger controls
        chkFlanger = new QCheckBox("Enable");
        sldFlangerRate = new QSlider(Qt::Horizontal);
        sldFlangerRate->setRange(1, 1000); // 0.1-10 Hz
        sldFlangerRate->setValue(100); // 1 Hz default
        sldFlangerDepth = new QSlider(Qt::Horizontal);
        sldFlangerDepth->setRange(1, 200); // 0.1-20 ms
        sldFlangerDepth->setValue(10); // 1 ms default
        sldFlangerFeedback = new QSlider(Qt::Horizontal);
        sldFlangerFeedback->setRange(0, 95); // 0-95%
        sldFlangerFeedback->setValue(50);
        sldFlangerMix = new QSlider(Qt::Horizontal);
        sldFlangerMix->setRange(0, 100); // 0-100%
        sldFlangerMix->setValue(50);

        QHBoxLayout* flangerLayout = new QHBoxLayout;
        flangerLayout->addWidget(chkFlanger);
        QFormLayout* flangerParamsLayout = new QFormLayout;
        flangerParamsLayout->addRow("Rate (0.1-10 Hz):", sldFlangerRate);
        flangerParamsLayout->addRow("Depth (0.1-20 ms):", sldFlangerDepth);
        flangerParamsLayout->addRow("Feedback (%):", sldFlangerFeedback);
        flangerParamsLayout->addRow("Mix (%):", sldFlangerMix);
        flangerLayout->addLayout(flangerParamsLayout);
        effectsLayout->addRow("Flanger:", flangerLayout);

        l->addWidget(effectsGroup);

        // Instrument assignment and simple params
        QHBoxLayout *midPanel = new QHBoxLayout;
        lstTracks = new QListWidget;
        midPanel->addWidget(lstTracks);

        // Right panel: synth params
        QWidget *paramWidget = new QWidget;
        QFormLayout *pf = new QFormLayout(paramWidget);
        cmbSynth = new QComboBox;
        cmbSynth->addItems({"Additive","FM","Karplus-Strong"});
        pf->addRow("Synth type:", cmbSynth);

        // Instrument presets for additive synthesis
        cmbInstrument = new QComboBox;
        for(const auto& preset : InstrumentPreset::getAllPresets()) {
            cmbInstrument->addItem(preset.name);
        }
        pf->addRow("Instrument:", cmbInstrument);
        
        spinSecondsPerBeat = new QDoubleSpinBox; spinSecondsPerBeat->setRange(0.1,2.0); spinSecondsPerBeat->setValue(0.5);
        pf->addRow("sec/beat:", spinSecondsPerBeat);
        btnPreview = new QPushButton("Preview note (A4)");
        pf->addRow(btnPreview);

        // Show/hide instrument selector based on synth type
        auto updateInstrumentVisibility = [this]() {
            bool isAdditive = cmbSynth->currentText() == "Additive";
            cmbInstrument->setVisible(isAdditive);
            QWidget* label = qobject_cast<QFormLayout*>(cmbInstrument->parentWidget()->layout())->labelForField(cmbInstrument);
            if (label) label->setVisible(isAdditive);
        };
        connect(cmbSynth, &QComboBox::currentTextChanged, updateInstrumentVisibility);
        updateInstrumentVisibility();
        midPanel->addWidget(paramWidget);
        l->addLayout(midPanel);

        // Progress bar
        progressBar = new QProgressBar;
        progressBar->setRange(0, 10000); // Para precisión de 0.01%
        progressBar->setTextVisible(true);
        progressBar->setFormat("Renderizando: %p%");
        progressBar->hide();
        l->addWidget(progressBar);

        // Spectrogram display
        specLabel = new QLabel;
        specLabel->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        specLabel->setMinimumHeight(256);
        specLabel->setMinimumWidth(512);
        specLabel->setStyleSheet("background: black;");
        l->addWidget(specLabel);

        // Connections
        connect(btnLoadMidi, &QPushButton::clicked, this, &MainWindow::onLoadMidi);
        connect(btnRender, &QPushButton::clicked, this, &MainWindow::onRender);
        connect(btnSave, &QPushButton::clicked, this, &MainWindow::onSave);
        connect(btnPreview, &QPushButton::clicked, this, &MainWindow::onPreview);

        setWindowTitle("TP2 - Synthesis All-in-One");
        resize(900,700);
    }

private slots:
    void onLoadMidi(){
        QString file = QFileDialog::getOpenFileName(this, "Open MIDI", ".", "MIDI files (*.mid *.midi)");
        if(file.isEmpty()) return;
        MidiFile mf;
        if(!loadMidi(file, mf)){
            QMessageBox::warning(this, "Error", "Failed to parse MIDI");
            return;
        }
        midi = mf;
        lblMidi->setText(QString("MIDI loaded: %1 tracks").arg(mf.tracks.size()));
        lstTracks->clear();
        trackConfigs.clear();
        for(size_t i=0;i<mf.tracks.size();++i){
            QListWidgetItem *it = new QListWidgetItem(QString("Track %1: %2").arg((int)i).arg(QString::fromStdString(mf.tracks[i].name)));
            lstTracks->addItem(it);
            TrackSynthConfig cfg;
            // default partial envs
            cfg.partialEnvs.resize(cfg.partialAmps.size());
            cfg.env.sampleRate = sampleRate;
            for(auto &pe: cfg.partialEnvs) pe.sampleRate = sampleRate;
            trackConfigs.push_back(cfg);
        }
    }

    void onRender(){
        if(midi.tracks.empty()){
            QMessageBox::warning(this,"No MIDI", "Load a MIDI file first.");
            return;
        }

        // Mostrar y resetear barra de progreso
        progressBar->setValue(0);
        progressBar->show();
        QApplication::processEvents();

        // Set synth types according to ui (simple: same for all tracks)
        TrackSynthConfig::Type type = TrackSynthConfig::ADDITIVE;
        if(cmbSynth->currentText()=="FM") type = TrackSynthConfig::FM;
        else if(cmbSynth->currentText()=="Karplus-Strong") type = TrackSynthConfig::KARPLUS;
        
        // Update track configs
        for(auto &cfg: trackConfigs) {
            cfg.type = type;
            if(type == TrackSynthConfig::ADDITIVE) {
                // Apply selected instrument preset
                auto presets = InstrumentPreset::getAllPresets();
                const auto& preset = presets[cmbInstrument->currentIndex()];
                cfg.partialAmps = preset.partialAmps;
                cfg.partialEnvs = preset.envelopes;
                // Ensure sample rate is set for all envelopes
                for(auto &env : cfg.partialEnvs) {
                    env.sampleRate = sampleRate;
                }
            }
        }

        float spb = spinSecondsPerBeat->value();
        rendered = renderMidi(midi, trackConfigs, spb, sampleRate,
            [this](float progress) {
                // Update progress bar (0-80% for MIDI rendering)
                progressBar->setValue(int(progress * 8000));
                QApplication::processEvents();
            });

        // Update progress to 80%
        progressBar->setValue(8000);
        QApplication::processEvents();

        // Apply effects based on controls
        if (chkDelay->isChecked()) {
            float delayTime = sldDelayTime->value() / 1000.0f; // convert ms to seconds
            float feedback = sldDelayFeedback->value() / 100.0f; // convert percentage to 0-1
            float mix = sldDelayMix->value() / 100.0f; // convert percentage to 0-1
            applyDelay(rendered, sampleRate, delayTime, feedback, mix);
        }
        
        if (chkReverb->isChecked()) {
            float roomSize = sldReverbRoom->value() / 100.0f;
            float damping = sldReverbDamping->value() / 100.0f;
            float wetLevel = sldReverbWet->value() / 100.0f;
            applyReverb(rendered, sampleRate, roomSize, damping, wetLevel);
        }

        // Apply flanger if enabled
        if (chkFlanger->isChecked()) {
            float rate = sldFlangerRate->value() / 100.0f; // 0.1-10 Hz
            float depth = sldFlangerDepth->value() / 10000.0f; // convert to seconds
            float feedback = sldFlangerFeedback->value() / 100.0f;
            float mix = sldFlangerMix->value() / 100.0f;
            applyFlanger(rendered, sampleRate, rate, depth, feedback, mix);
        }

        // Update progress to 90%
        progressBar->setValue(9000);
        QApplication::processEvents();
        
        // spectrogram
        QImage spec = renderSpectrogram(rendered, sampleRate, 1024, 512);
        specLabel->setPixmap(QPixmap::fromImage(spec).scaled(specLabel->size(), Qt::IgnoreAspectRatio, Qt::SmoothTransformation));
        // Update progress to 100% and hide
        progressBar->setValue(10000);
        QApplication::processEvents();
        progressBar->hide();
        QMessageBox::information(this,"Rendered", QString("Rendered %1 samples (%2 s)").arg(rendered.size()).arg((double)rendered.size()/sampleRate));
    }

    void onSave(){
        if(rendered.empty()){
            QMessageBox::warning(this,"No audio", "Render audio before saving.");
            return;
        }
        QString file = QFileDialog::getSaveFileName(this, "Save WAV", "out.wav", "WAV file (*.wav)");
        if(file.isEmpty()) return;
        if(!writeWAV(file, rendered, sampleRate)){
            QMessageBox::warning(this,"Error","Failed to write WAV");
        } else {
            QMessageBox::information(this,"Saved", "WAV saved.");
        }
    }

    void onPreview(){
        // simple preview: synth A4 for 2s with chosen synth type
        int midiNote = 69;
        float dur = 2.0f;
        vector<float> buf;
        TrackSynthConfig cfg;
        cfg.partialEnvs.resize(cfg.partialAmps.size());
        if(cmbSynth->currentText()=="Additive"){
            cfg.type = TrackSynthConfig::ADDITIVE;
            auto presets = InstrumentPreset::getAllPresets();
            const auto& preset = presets[cmbInstrument->currentIndex()];
            buf = synthAdditive(midiNote, dur, sampleRate, preset.partialAmps, preset.envelopes);
        } else if(cmbSynth->currentText()=="FM"){
            cfg.type = TrackSynthConfig::FM;
            buf = synthFM(midiNote, dur, sampleRate, cfg.carrierRatio, cfg.modRatio, cfg.modIndex, cfg.env);
        } else {
            buf = karplusStrong(midiNote, dur, sampleRate, cfg.ksDecay);
        }
        // store rendered and show spectrogram
        rendered = buf;
        QImage spec = renderSpectrogram(rendered, sampleRate, 1024, 512);
        specLabel->setPixmap(QPixmap::fromImage(spec).scaled(specLabel->size(), Qt::KeepAspectRatio));
        // play using QAudioOutput (simple)
        playBuffer(rendered);
    }

    void playBuffer(const vector<float>& buf){
        QAudioFormat format;
        format.setSampleRate(sampleRate);
        format.setChannelCount(1);
        format.setSampleSize(16);
        format.setCodec("audio/pcm");
        format.setByteOrder(QAudioFormat::LittleEndian);
        format.setSampleType(QAudioFormat::SignedInt);
        QAudioDeviceInfo info(QAudioDeviceInfo::defaultOutputDevice());
        if(!info.isFormatSupported(format)){
            format = info.nearestFormat(format);
        }
        QAudioOutput *out = new QAudioOutput(format, this);
        QByteArray pcm; pcm.resize(buf.size()*2);
        int16_t *dst = (int16_t*)pcm.data();
        for(size_t i=0;i<buf.size();i++){
            float v = max(-1.0f, min(1.0f, buf[i]));
            dst[i] = int16_t(v*32767);
        }
        QBuffer *buffer = new QBuffer(this);
        buffer->setData(pcm);
        buffer->open(QIODevice::ReadOnly);
        out->start(buffer);
        // cleanup when done
        connect(out, &QAudioOutput::stateChanged, this, [out,buffer](QAudio::State s){
            if(s==QAudio::IdleState || s==QAudio::StoppedState){
                out->deleteLater();
                buffer->close();
                buffer->deleteLater();
            }
        });
    }

private:
    QWidget *central;
    QPushButton *btnLoadMidi, *btnRender, *btnSave, *btnPreview;
    QLabel *lblMidi;
    QListWidget *lstTracks;
    QComboBox *cmbSynth;
    QComboBox *cmbInstrument;
    QDoubleSpinBox *spinSecondsPerBeat;
    QLabel *specLabel;
    QProgressBar *progressBar;

    // Effects controls
    QCheckBox *chkDelay;
    QSlider *sldDelayTime;
    QSlider *sldDelayFeedback;
    QSlider *sldDelayMix;
    QCheckBox *chkReverb;
    QSlider *sldReverbRoom;
    QSlider *sldReverbDamping;
    QSlider *sldReverbWet;
    
    // Flanger controls
    QCheckBox *chkFlanger;
    QSlider *sldFlangerRate;
    QSlider *sldFlangerDepth;
    QSlider *sldFlangerFeedback;
    QSlider *sldFlangerMix;

    // WAV Player and FFT display
    WAVPlayer *wavPlayer;
    QPushButton *btnLoadWav;
    QPushButton *btnPlayWav;
    QPushButton *btnStopWav;
    QLabel *lblWavFile;
    FFTWidget *fftWidget;

    MidiFile midi;
    vector<TrackSynthConfig> trackConfigs;
    vector<float> rendered;
    uint32_t sampleRate;
};

// --------------------------- main --------------------------------------------
int main(int argc, char **argv){
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();
}

#include "tp2_allinone.moc"