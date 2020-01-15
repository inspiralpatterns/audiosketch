cmake_minimum_required(VERSION 3.8)
project(audiosketch)

set(CMAKE_CXX_STANDARD 11)

include_directories(/Users/Boggiz/audiosketch/files)
set(SOURCE_FILES main.cpp WavFile.h files/Vox.wav files/out.wav)
add_executable(audiosketch ${SOURCE_FILES})