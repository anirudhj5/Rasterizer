cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(project1F)
SET(VTK_DIR C:/Users/Anirudh/Desktop/CIS441/1F/src)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(project1F project1F)

if(VTK_LIBRARIES)
  target_link_libraries(project1F ${VTK_LIBRARIES})
else()
  target_link_libraries(project1F vtkHybrid)
endif()

