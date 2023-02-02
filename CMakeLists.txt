cmake_minimum_required(VERSION 3.14) #3.14 needed for GoogleTest

# define project name, version
project(PSEMolDyn_GroupG)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# set Release as the default build type if it is not yet set.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# collect all cpp files
file(GLOB_RECURSE MY_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    # header don't need to be included but this might be necessary for some IDEs
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
)

# create make target
add_executable(MolSim ${MY_SRC})

# set cxx standard. You may raise this if you want.
target_compile_features(MolSim
        PRIVATE
            cxx_std_17
)

target_include_directories(MolSim
        PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}/libs/libxsd
            #${CMAKE_CURRENT_SOURCE_DIR}/libs/spdlog
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src
)

find_package(PkgConfig REQUIRED)
pkg_check_modules(XERCESC REQUIRED IMPORTED_TARGET xerces-c)

target_link_libraries(MolSim
        # stuff that is used in headers and source files
        PUBLIC
            PkgConfig::XERCESC
)

set(MP_OPT_GCC "")
set(MP_OPT_INTEL "")

find_package(OpenMP)
#set(OpenMP_CXX_FOUND 0) #for testing if it works without OpenMP
if(OpenMP_CXX_FOUND)
    target_link_libraries(MolSim PUBLIC OpenMP::OpenMP_CXX)
    set(MP_OPT_GCC -fopenmp)
    set(MP_OPT_INTEL -qopenmp)
endif()

# activate all compiler warnings. Clean up your code :P
# depending on the compiler different flags are used
target_compile_options(MolSim
        PRIVATE
            $<$<CXX_COMPILER_ID:GNU>: -Wall ${MP_OPT_GCC}> 
            # disable some spam warnings for icpc...
            $<$<CXX_COMPILER_ID:Intel>: -w3 -wd383,981,1418,1572,225 ${MP_OPT_INTEL}>
)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

include(doxygen)
include(googletest)