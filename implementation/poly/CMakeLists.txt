 cmake_minimum_required (VERSION 2.8.11)
project (vml)
set(NAME "vml")
include_directories(${CMAKE_SOURCE_DIR})

include_directories("./ate-pairing/include/")
link_directories("./ate-pairing/lib/")

include_directories(/usr/local/include/)
link_directories("/usr/local/lib/")



set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} cmake/)


set (CMAKE_C_FLAGS "-pthread  -march=native -O3 -maes ")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg -ggdb")
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -std=c++11")



set(basics
)

add_executable(inputvpd inputvpd.cpp)
target_link_libraries(inputvpd gmp gmpxx zm)


#add_executable(vpd vpd_test.cpp)
#target_link_libraries(vpd gmp gmpxx zm)

#add_executable(mult multivariate_test.cpp)
#target_link_libraries(mult gmp gmpxx zm)