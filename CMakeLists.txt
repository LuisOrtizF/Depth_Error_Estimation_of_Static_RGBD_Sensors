cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(depth_error)

find_package(OpenCV 3.4.1 REQUIRED)
find_package(PCL 1.8.0 REQUIRED)

include_directories ("${PROJECT_SOURCE_DIR}/fce")
add_subdirectory (fce)

include_directories(${PCL_INCLUDE_DIRS} )
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS} ${PCL_LIBRARIES} fce)