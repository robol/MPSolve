# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)

find_package(Qt5Widgets)

if(Qt5Widgets_FOUND)

set(XMPSOLVE_SOURCES
  polynomialsolver.cpp
  polynomialsolver.h
  rootsmodel.cpp
  rootsmodel.h
  polsyntaxhighlighter.cpp
  polsyntaxhighlighter.h
  mpsolveworker.cpp
  mpsolveworker.h
  mainwindow.ui
  mainwindow.cpp
  mainwindow.h
  polfileeditorwindow.cpp
  polfileeditorwindow.ui
  polfileeditorwindow.h
  polfileeditor.cpp
  polfileeditor.ui
  polfileeditor.h
  qrootsrenderer.cpp
  main.cpp
  rootsrenderer.cpp
  rootsrenderer.h
  root.cpp
  root.h
  resources.qrc
)

add_executable(xmpsolve ${XMPSOLVE_SOURCES})

include_directories(${CMAKE_SOURCE_DIR}/include
  ${CMAKE_BINARY_DIR}/include
)

add_definitions("-include ${CMAKE_BINARY_DIR}/config.h")

target_link_libraries(xmpsolve Qt5::Widgets mps)

endif(Qt5Widgets_FOUND)
