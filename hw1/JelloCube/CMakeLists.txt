cmake_minimum_required (VERSION 3.1)
project (Assignment1)

#########################################################
# SET COMPILER 
#########################################################
set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -std=c++11)

#########################################################
# FIND GLUT
#########################################################
find_package(GLUT REQUIRED)
include_directories("C:/Users/Xiangshun Bei/Dropbox/Study/USC/520/hw/hw1/OpenGLUT-0.6.3")
link_directories("C:/Users/Xiangshun Bei/Dropbox/Study/USC/520/hw/hw1/OpenGLUT-0.6.3")
add_definitions(${GLUT_DEFINITIONS})
if(NOT GLUT_FOUND)
    message(ERROR " GLUT not found!")
endif(NOT GLUT_FOUND)

#########################################################
# FIND OPENGL
#########################################################
find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})
link_directories(${OpenGL_LIBRARY_DIRS})
add_definitions(${OpenGL_DEFINITIONS})
if(NOT OPENGL_FOUND)
    message(ERROR " OPENGL not found!")
endif(NOT OPENGL_FOUND)


#########################################################
# CREATE WORLD FILE VARIABLES
#########################################################
SET(CW_SOURCE
  createWorld.cpp
)
SET(CW_HEADERS
  #createWorld.h #createWorld.cpp doesn't have header
)

#########################################################
# JELLO EXE VARIABLE
#########################################################
SET(JELLO_SOURCE
  jello.cpp
  input.cpp
  physics.cpp
  ppm.cpp
  showCube.cpp
  pic.cpp
)
SET(JELLO_HEADERS
  jello.h
  input.h
  physics.h
  #ppm.h #ppm doesn't have a header
  showCube.h
  pic.h
)

#########################################################
# ADD SUB DIRECTORIES
#########################################################
#add_subdirectory(<name>)

#########################################################
# ADD EXECUTABLES
#########################################################
add_executable(jello ${JELLO_SOURCE} ${JELLO_HEADERS})
add_executable(createWorld ${CW_SOURCE} ${CW_HEADERS})

#########################################################
# LINK LIBRARIES
#########################################################
target_link_libraries(jello ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})