cmake_minimum_required(VERSION 3.5)
project(tp2_allinone LANGUAGES CXX)

# Usar C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Buscar Qt5 (Widgets y Multimedia)
find_package(Qt5 REQUIRED COMPONENTS Widgets Multimedia)

# Habilitar moc/uic/rcc automáticos (para Q_OBJECT)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

# Archivo fuente
set(SOURCES
    tp2_allinone.cpp
)

# Crear ejecutable
add_executable(${PROJECT_NAME} ${SOURCES})

# Linkear con las librerías de Qt5
target_link_libraries(${PROJECT_NAME}
    Qt5::Widgets
    Qt5::Multimedia
)