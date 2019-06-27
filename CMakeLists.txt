cmake_minimum_required(VERSION 2.8.9)
project(compliance_opt)

 INCLUDE_DIRECTORIES(
  include
  ${CMAKE_SOURCE_DIR}/include
  )

ADD_LIBRARY(compliance_opt ${CMAKE_SOURCE_DIR}/source/compliance_opt.cc)

ADD_EXECUTABLE(run ${CMAKE_SOURCE_DIR}/source/run.cc)
TARGET_LINK_LIBRARIES(run compliance_opt)

# Define a distclean target to remove every generated file:
ADD_CUSTOM_TARGET(distclean
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
  COMMAND ${CMAKE_COMMAND} -E remove_directory CMakeFiles
  COMMAND ${CMAKE_COMMAND} -E remove CMakeCache.txt cmake_install.cmake Makefile
  COMMENT "distclean invoked"
  )