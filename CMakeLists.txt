cmake_minimum_required(VERSION 3.0)
project(compliance_opt)

 INCLUDE_DIRECTORIES(
  include
  ${CMAKE_SOURCE_DIR}/include
  )

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})

find_package (Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR}/unsupported/)
include_directories(${EIGEN3_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR}/supported/)

FIND_PACKAGE(Umfpack REQUIRED) #requires you to copy FindUmfpack.cmake into the source directory
FIND_PACKAGE(BLAS REQUIRED)

find_package(GSL REQUIRED)    # See below (2)

include_directories(${UMFPACK_INCLUDES})

ADD_LIBRARY(compliance_opt ${CMAKE_SOURCE_DIR}/source/compliance_opt.cc)
TARGET_LINK_LIBRARIES(compliance_opt umfpack)
target_link_libraries(compliance_opt GSL::gsl GSL::gslcblas)

ADD_EXECUTABLE(run ${CMAKE_SOURCE_DIR}/source/run.cc)
TARGET_LINK_LIBRARIES(run compliance_opt)

ADD_EXECUTABLE(dakota_runner ${CMAKE_SOURCE_DIR}/source/dakota_runner.cc)
TARGET_LINK_LIBRARIES(dakota_runner compliance_opt)

ADD_EXECUTABLE(dakota_compliance_op ${CMAKE_SOURCE_DIR}/source/dakota_compliance_op.cc)
TARGET_LINK_LIBRARIES(dakota_compliance_op compliance_opt)


# Define a distclean target to remove every generated file:
ADD_CUSTOM_TARGET(distclean
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
  COMMAND ${CMAKE_COMMAND} -E remove_directory CMakeFiles
  COMMAND ${CMAKE_COMMAND} -E remove CMakeCache.txt cmake_install.cmake Makefile
  COMMENT "distclean invoked"
  )
