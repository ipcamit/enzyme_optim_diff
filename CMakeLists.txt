cmake_minimum_required(VERSION 3.16)
project(libdescriptor)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fuse-ld=lld -flto")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fuse-ld=lld -flto")
add_executable(o3_trial  ${SYMFUNFILES} o3_err_example.cpp SymFun.cpp SymFunData.cpp)
target_link_options(o3_trial PRIVATE "LINKER:--lto-legacy-pass-manager")
target_link_options(o3_trial PRIVATE "LINKER:-mllvm=-load=/opt/enzyme/enzyme/build/Enzyme/LLDEnzyme-12.so")
