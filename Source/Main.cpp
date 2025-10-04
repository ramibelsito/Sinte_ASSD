#include <JuceHeader.h>

class MainComponent : public juce::Component {
public:
    MainComponent() { setSize (400, 200); }

    void paint (juce::Graphics& g) override {
        g.fillAll(juce::Colours::black);
        g.setColour(juce::Colours::white);
        g.setFont(20.0f);
        g.drawFittedText("Hello JUCE!", getLocalBounds(), juce::Justification::centred, 1);
    }
};

class MainApplication : public juce::JUCEApplication {
public:
    const juce::String getApplicationName() override { return "MyProject"; }
    const juce::String getApplicationVersion() override { return "1.0"; }
    void initialise (const juce::String&) override {
        mainWindow.reset(new MainWindow("MyProject", new MainComponent(), *this));
    }
    void shutdown() override { mainWindow = nullptr; }

private:
    class MainWindow : public juce::DocumentWindow {
    public:
        MainWindow(juce::String name, juce::Component* c, JUCEApplication& a)
            : DocumentWindow(name, juce::Colours::lightgrey, DocumentWindow::allButtons),
              app(a) {
            setUsingNativeTitleBar(true);
            setContentOwned(c, true);
            centreWithSize(getWidth(), getHeight());
            setVisible(true);
        }
        void closeButtonPressed() override { app.systemRequestedQuit(); }
    private:
        JUCEApplication& app;
    };
    std::unique_ptr<MainWindow> mainWindow;
};

START_JUCE_APPLICATION(MainApplication)
