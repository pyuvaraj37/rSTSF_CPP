cmake_minimum_required(VERSION 3.3)

project(main)

include_directories(include)
add_executable(main src/main.cpp 
                    src/autoreg.cpp
                    src/intBasedT.cpp)
# add everytime you create cpp 
