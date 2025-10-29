#
#  cmake file for "source coding library"
#
set( MY_PROJECT_NAME SourceCoding )


#----- minimum required cmake version -----
cmake_minimum_required( VERSION 3.10 FATAL_ERROR )

#----- project name -----
project( ${MY_PROJECT_NAME} )

#----- build type -----
set( default_build_type "Release" )
if( NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES )
  message( STATUS "Setting build type to '${default_build_type}' as none was specified." )
  set( CMAKE_BUILD_TYPE "${default_build_type}" CACHE
    STRING "Choose the type of build." FORCE )
  # Set the possible values of build type for cmake-gui
  set_property( CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Release" "Debug" )  
endif()
if( NOT CMAKE_CONFIGURATION_TYPES AND
    NOT CMAKE_BUILD_TYPE STREQUAL "Debug" AND
    NOT CMAKE_BUILD_TYPE STREQUAL "Release" )
  message( SEND_ERROR
    " Build type '${CMAKE_BUILD_TYPE}' is not supported.\n"
    " Choose one of: 'Debug' 'Release'"
    )
endif()

#----- set c++11 -----
set( CMAKE_CXX_STANDARD          20 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

#----- compiler options -----
if( MSVC )
  add_compile_options( /W4 /WX )
  add_compile_options( $<$<CONFIG:Release>:/Oi> )
  add_compile_options( $<$<CONFIG:Release>:/Ot> )
elseif( CMAKE_CXX_COMPILER_ID STREQUAL GNU )
  add_compile_options( -Wall -Wextra -pedantic -Werror )
  add_compile_options( $<$<CONFIG:Release>:-O3> )
elseif( CMAKE_CXX_COMPILER_ID STREQUAL Clang )
  add_compile_options( -Wall -Wextra -pedantic -Werror )
  add_compile_options( $<$<CONFIG:Release>:-O3> )
endif()

#----- use folders in IDEs for projects -----
set_property( GLOBAL PROPERTY USE_FOLDERS ON )

#----- configuration identification -----
set( MY_COMPILER_VERSION "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}" )
if( ${CMAKE_VERSION} VERSION_LESS "3.12" )
  message( "WARNING: For enabling all features, at least CMake 3.12 is required." )
  if( MSVC )
    add_compile_options(   $<$<CONFIG:Release>:/D> $<$<CONFIG:Release>:MY_BUILD_MODE="release"> )
    add_compile_options(     $<$<CONFIG:Debug>:/D> $<$<CONFIG:Debug>:MY_BUILD_MODE="debug"> )
  else()
    add_compile_options(   $<$<CONFIG:Release>:-DMY_BUILD_MODE="release"> )
    add_compile_options(     $<$<CONFIG:Debug>:-DMY_BUILD_MODE="debug"> )
  endif()
else()
  add_compile_definitions( $<$<CONFIG:Release>:MY_BUILD_MODE="release"> )
  add_compile_definitions(   $<$<CONFIG:Debug>:MY_BUILD_MODE="debug"> )
endif()

#----- use ccache -----
find_program( CCACHE_FOUND ccache )
if( CCACHE_FOUND )
  message( STATUS "ccache found. using it." )
  set_property( GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache )
  set_property( GLOBAL PROPERTY RULE_LAUNCH_LINK    ccache )
endif()



#===== library =====
file( GLOB LIB_SOURCES "src/lib/*.cpp" )
file( GLOB LIB_HEADERS "src/lib/*.h" )
set( LIB_NAME ${MY_PROJECT_NAME} )
add_library( ${LIB_NAME} ${LIB_SOURCES} ${LIB_HEADERS} )



#===== executables =====
file( GLOB APP_SOURCES "src/app/*.cpp" )

foreach( testSourceFile ${APP_SOURCES} )

    get_filename_component(       EXE_NAME  ${testSourceFile} NAME_WE ) 
    add_executable(             ${EXE_NAME} ${testSourceFile} )
    target_include_directories( ${EXE_NAME} PUBLIC "src/lib" )
    target_link_libraries(      ${EXE_NAME} ${LIB_NAME} )
    if( MSVC )
      set_target_properties(    ${EXE_NAME}
        PROPERTIES
		FOLDER                   "tests"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}" )
    else()
      if( CMAKE_BUILD_TYPE STREQUAL "Debug" )
        set_target_properties(    ${EXE_NAME}
          PROPERTIES
          RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}/debug" )
      else()
        set_target_properties(    ${EXE_NAME}
          PROPERTIES
          RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}/release" )
      endif()
    endif()

endforeach( testSourceFile ${APP_SOURCES} )



