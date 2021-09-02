cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

# GoogleTest
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()
# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()


# ---- Project ----

# Note: update this to your new project's name and version
project(IsoLasso_Parallel
        VERSION 1.0.1
        LANGUAGES CXX)

# ---- Include guards ----

message(${PROJECT_SOURCE_DIR})
message(${PROJECT_BINARY_DIR})

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.")
endif()

# --- Import tools ----

include(cmake/tools.cmake)

# ---- Add dependencies via CPM ----
# see https://github.com/TheLartians/CPM.cmake for more info

include(cmake/CPM.cmake)

# PackageProject.cmake will be used to make our target installable
CPMAddPackage(
  NAME PackageProject.cmake
  GITHUB_REPOSITORY TheLartians/PackageProject.cmake
  VERSION 1.0
)

set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS} -O3 -DDEBUG -std=c++20 -lpthread")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -std=c++20 -lpthread")

aux_source_directory(. DIR_SRCS)
include_directories("${PROJECT_SOURCE_DIR}/lib/include")

add_subdirectory(src/utils)
add_subdirectory(src/EM_module)
add_subdirectory(src/SAM_module)
link_directories(lib/include/SAM_module)
link_directories(lib/include/EM_module)
link_directories(lib/include/utils)

add_executable(IsoLasso_Parallel src/main.cpp)

#find_package(Boost 1.73.0 REQUIRED COMPONENTS iostreams)
#if(Boost_FOUND)
#    INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )
#    message(${Boost_INCLUDE_DIRS})
#    target_link_libraries(IsoLasso_Parallel ${Boost_LIBRARY} )
#    message("Boost iostream found!")
#endif()

#Pthread
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

target_link_libraries(IsoLasso_Parallel SAM_functions utils_functions EM_functions Threads::Threads )

## UnitTest
find_package(GTest REQUIRED)
enable_testing()

add_executable(unittest Unittest/src/SAM_IO_test.cpp Unittest/src/SAMProcess_test.cpp)
target_link_libraries(unittest SAM_functions utils_functions EM_functions gtest_main Threads::Threads )
