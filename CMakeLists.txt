cmake_minimum_required(VERSION 3.0...3.31.6 FATAL_ERROR)

project(CalPhenix)

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

set(CMAKE_CXX_FLAGS "-Wall -Werror -Wpedantic -pipe -O2")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

find_package(ROOT QUIET)

if(${ROOT_FOUND})
   message(STATUS "ROOT package found")
   link_libraries(${ROOT_LIBRARIES})
elseif(EXISTS "${ROOT_PATH}/bin/root")
   message(STATUS "ROOT path found")
   set(ROOT_INCLUDE_DIRS "${ROOT_PATH}/include")
   set(ROOT_root_CMD "${ROOT_PATH}/bin/root")
   execute_process(COMMAND ${ROOT_root_CMD}-config --glibs OUTPUT_VARIABLE ROOT_LIB_FLAGS)
   string(STRIP ${ROOT_LIB_FLAGS} ROOT_LIB_FLAGS)
   set(CMAKE_SHARED_LINKER_FLAGS ${ROOT_LIB_FLAGS})
   set(CMAKE_EXE_LINKER_FLAGS ${ROOT_LIB_FLAGS})
else()
   message(FATAL_ERROR "ROOT not found: install it via package manager or set ROOT_PATH if you installed it manually")
endif()

# json was changed for yaml; if json will be used this field will be uncommented
#if (EXISTS /usr/include/json)
#   message(STATUS "jsoncpp package found")
#elseif (EXISTS /usr/include/jsoncpp)
#   message(STATUS "libjsoncpp package found")
#   include_directories(${JSON_CPP_PATH}/include/jsoncpp)
#elseif (EXISTS ${JSON_CPP_PATH}/include/json)
#   message(STATUS "jsoncpp path found")
#   link_directories(${JSON_CPP_PATH}/build/lib)
#   include_directories(${JSON_CPP_PATH}/include)
#else()
#   message(FATAL_ERROR "jsoncpp/libjsoncpp not found: install it via package manager or set environment variable JSON_CPP_PATH if you installed it manually")
#endif()

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/CppTools/include ${CMAKE_SOURCE_DIR}/ROOTTools/include ${CMAKE_SOURCE_DIR}/ProgressBar/include)

add_subdirectory(yaml-cpp EXCLUDE_FROM_ALL)
add_subdirectory(CppTools EXCLUDE_FROM_ALL)
add_subdirectory(ROOTTools EXCLUDE_FROM_ALL)
add_subdirectory(ProgressBar EXCLUDE_FROM_ALL)

#link_libraries(jsoncpp)
link_libraries(yaml-cpp)
link_libraries(Box IOTools StrTools ErrorHandler)
link_libraries(TCanvasTools TF1Tools)
link_libraries(PBar)

add_library(InputYAMLReader ${CMAKE_SOURCE_DIR}/src/InputYAMLReader.cpp)

link_libraries(InputYAMLReader)

add_executable(SigmalizedResiduals ${CMAKE_SOURCE_DIR}/src/SigmalizedResiduals.cpp)
add_executable(CheckSigmalizedResiduals ${CMAKE_SOURCE_DIR}/src/CheckSigmalizedResiduals.cpp)
add_executable(EMCTTowerOffset ${CMAKE_SOURCE_DIR}/src/EMCTTowerOffset.cpp)
add_executable(EMCTRunByRunOffset ${CMAKE_SOURCE_DIR}/src/EMCTRunByRunOffset.cpp)
