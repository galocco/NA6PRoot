set(THIS_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(THIS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(MOD_NAME "base")
set(LIB_NAME "${MOD_NAME}Lib")

# Include directoroes
include_directories(${THIS_INCLUDE_DIR})
include_directories(${UTILS_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/sim/include)
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

file(GLOB THIS_SOURCES ${THIS_SRC_DIR}/*.cxx)
file(GLOB THIS_HEADERS ${THIS_INCLUDE_DIR}/*.h)
# NA6PGeometryManager uses the configured MWPC grid in its legacy fallback.
# Compile the parameter singleton in baseLib to avoid a baseLib -> simLib cycle.
list(APPEND THIS_SOURCES ${CMAKE_SOURCE_DIR}/sim/src/NA6PMWPCParam.cxx)

# Add the library
add_library(${LIB_NAME} SHARED ${THIS_SOURCES} ${THIS_HEADERS} )

set(DICTIONARY_HEADERS
  NA6PLayoutParam.h
  MagneticFieldRegion.h
  MagneticField.h
  NA6PBeamParam.h
  NA6PBeam.h
  NA6PGeometryManager.h
  NA6PMCComposedLabel.h
  NA6PMCTruthContainer.h
  ) # Add all class headers requiring dictionary
  
set(DICTIONARY_LINKDEF ${THIS_SRC_DIR}/${MOD_NAME}LinkDef.h)
set(DICTIONARY_OUTPUT ${LIB_NAME}Dict)

ROOT_GENERATE_DICTIONARY(
    ${DICTIONARY_OUTPUT}
    ${DICTIONARY_HEADERS}
    LINKDEF ${DICTIONARY_LINKDEF}
    MODULE ${LIB_NAME}
)

# Set include directories for consumers
target_include_directories(${LIB_NAME} INTERFACE ${THIS_INCLUDE_DIR})

# !!!!!!!!!!!!!!!!!!!!!!
#target_compile_options(${LIB_NAME} PRIVATE -ggdb -g -O0)

# Link against required libraries
target_link_libraries(${LIB_NAME}
    ${ROOT_LIBRARIES}
    fmt::fmt
    FairLogger::FairLogger
    utilsLib
)

# Installation instructions
install(TARGETS ${LIB_NAME} LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
install(FILES ${THIS_HEADERS} DESTINATION include)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/lib${LIB_NAME}_rdict.pcm
  ${CMAKE_CURRENT_BINARY_DIR}/lib${LIB_NAME}.rootmap DESTINATION lib)
