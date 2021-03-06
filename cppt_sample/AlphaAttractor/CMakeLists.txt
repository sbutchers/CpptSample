CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(attractor_model)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "~/.cpptransport201701/share/cmake/")

INCLUDE(CompilerFlags)

SET_COMPILER_FLAGS()

FIND_PACKAGE(CppTransport REQUIRED)

INCLUDE_DIRECTORIES(${CPPTRANSPORT_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR} ${COSMOSIS_SRC_DIR}/config/compilers.mk ${COSMOSIS_SRC_DIR})

ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/attractor_model_core.h ${CMAKE_CURRENT_BINARY_DIR}/attractor_model_mpi.h
  COMMAND CppTransport --verbose --fast ${CMAKE_CURRENT_SOURCE_DIR}/attractor_model.model
  DEPENDS attractor_model.model
)

SET(HEADERS ${CMAKE_CURRENT_BINARY_DIR}/attractor_model_core.h ${CMAKE_CURRENT_BINARY_DIR}/attractor_model_mpi.h)

ADD_CUSTOM_TARGET(Generator DEPENDS ${HEADERS})

ADD_LIBRARY(cppt_sample SHARED cosmosis_module_attractor.cpp)
ADD_DEPENDENCIES(cppt_sample Generator)
TARGET_LINK_LIBRARIES(cppt_sample ${CPPTRANSPORT_LIBRARIES} ${COSMOSIS_SRC_DIR}/cosmosis/datablock/libcosmosis.so)
SET_TARGET_PROPERTIES(cppt_sample PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED YES
  CXX_EXTENSIONS NO
  )
TARGET_COMPILE_OPTIONS(cppt_sample PRIVATE -Wno-ignored-attributes -Wno-deprecated-declarations)