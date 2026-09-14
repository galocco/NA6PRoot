cmake_minimum_required(VERSION 3.20)
project(NA6PMWPCParticleGun LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

get_filename_component(NA6P_SOURCE "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
if(DEFINED ENV{NA6PROOT_ROOT} AND NOT "$ENV{NA6PROOT_ROOT}" STREQUAL "")
  set(_default_na6p_install "$ENV{NA6PROOT_ROOT}")
else()
  set(_default_na6p_install "${NA6P_SOURCE}/install")
endif()
set(NA6P_INSTALL "${_default_na6p_install}" CACHE PATH "Installed NA6PRoot tree to test")

if(NOT EXISTS "${NA6P_INSTALL}/include/NA6PMWPCChamber.h")
  message(FATAL_ERROR
    "${NA6P_INSTALL}/include/NA6PMWPCChamber.h was not found. "
    "Build and install the current feature/mwpc-dice branch first, then source install/init.sh.")
endif()

find_package(ROOT REQUIRED COMPONENTS Core RIO Tree Geom Physics EG)
find_package(Geant4 REQUIRED)
find_package(VMC REQUIRED)
find_package(VGM REQUIRED)
find_package(Geant4VMC REQUIRED)
find_package(fmt REQUIRED)
find_package(FairLogger REQUIRED)

find_library(NA6P_SIM simLib PATHS "${NA6P_INSTALL}/lib" NO_DEFAULT_PATH REQUIRED)
find_library(NA6P_BASE baseLib PATHS "${NA6P_INSTALL}/lib" NO_DEFAULT_PATH REQUIRED)
find_library(NA6P_UTILS utilsLib PATHS "${NA6P_INSTALL}/lib" NO_DEFAULT_PATH REQUIRED)

add_executable(mwpc_particle_gun
  mwpc_particle_gun.cxx
  MWPCTestModule.cxx
)

target_include_directories(mwpc_particle_gun PRIVATE
  "${CMAKE_CURRENT_LIST_DIR}"
  "${NA6P_INSTALL}/include"
  ${Geant4_INCLUDE_DIRS}
)

target_compile_definitions(mwpc_particle_gun PRIVATE
  MWPC_TEST_SOURCE_DIR="${CMAKE_CURRENT_LIST_DIR}"
)

target_link_libraries(mwpc_particle_gun PRIVATE
  "${NA6P_SIM}"
  "${NA6P_BASE}"
  "${NA6P_UTILS}"
  geant4vmc
  VMCLibrary
  ROOT::Core
  ROOT::RIO
  ROOT::Tree
  ROOT::Geom
  ROOT::Physics
  ROOT::EG
  FairLogger::FairLogger
  fmt::fmt
  ${Geant4_LIBRARIES}
)
