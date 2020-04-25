find_package (Boost COMPONENTS system filesystem REQUIRED)
include_directories (${CMAKE_SOURCE_DIR}/include ${Boost_INCLUDE_DIRS})

file(GLOB SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
file(GLOB HEADERS RELATIVE ${CMAKE_SOURCE_DIR}/include *.?pp *.h)

list(REMOVE_ITEM SRCS "paths.cpp")

include(CheckIncludeFiles)

add_library(silent STATIC ${SRCS} ${HEADERS})
add_executable(silent_project paths.cpp)
target_link_libraries(silent_project silent ${OpenCV_LIBS} ${Boost_LIBRARIES} python2.7)

