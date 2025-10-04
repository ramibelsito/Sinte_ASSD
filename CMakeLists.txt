cmake_minimum_required(VERSION 3.15)

project(Sinte_ASSD VERSION 0.1.0)

# Ruta a JUCE
add_subdirectory(JUCE/JUCE)

juce_add_gui_app(Sinte_ASSD
    PRODUCT_NAME "Sinte_ASSD"
    VERSION ${PROJECT_VERSION}
    COMPANY_NAME "Grupo 5 ASSD"
    ICON_BIG "icon.png"
    ICON_SMALL "icon_small.png"
    SOURCES
        Source/Main.cpp
        Source/MainComponent.cpp
        #Agregar el resto de archivos fuente aquí
        )

target_link_libraries(Sinte_ASSD PRIVATE
    juce::juce_audio_utils
    juce::juce_dsp
    juce::juce_gui_extra
    juce::juce_core
    juce::juce_data_structures
    juce::juce_events
    juce::juce_graphics
    juce::juce_gui_basics
)

