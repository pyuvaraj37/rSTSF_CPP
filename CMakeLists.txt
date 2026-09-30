cmake_minimum_required(VERSION 3.3)

project(main)

include_directories(include)
add_executable(main src/main.cpp 
                    src/XTransformations.cpp
                    src/intBasedT.cpp
                    src/treeBasedPredict.cpp
                    src/loader.cpp)
# add everytime you create cpp 
target_link_libraries(main PRIVATE fftw3)