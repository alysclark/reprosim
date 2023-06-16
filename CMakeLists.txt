
cmake_minimum_required(VERSION 3.18.0)

set(CMAKE_OSX_DEPLOYMENT_TARGET 10.15 CACHE STRING "Minimum OS X deployment version.")

project(UnnamedProject VERSION 0.1.0 LANGUAGES C)

get_property(IS_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)

if(IS_MULTI_CONFIG)
  set(CONFIG_DIR_SUFFIX "/$<CONFIG>")
endif()

add_subdirectory(src)

