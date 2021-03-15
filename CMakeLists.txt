cmake_minimum_required(VERSION 3.16)
project(mesh_smooth_for_finite_elements)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(DATA_SHARED_PATH ${PROJECT_SOURCE_DIR}/data CACHE PATH "location of shared data resources")
set(GTEST_PATH ${PROJECT_SOURCE_DIR}/external/googletest/googletest)

if (MSVC)
    # Enable parallel compilation
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj /openmp")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR})
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR})
else ()
    # Libigl requires a modern C++ compiler that supports c++11
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
endif ()

if (APPLE)
    set(INCLUDE /usr/local/include)
    set(SURFACE_MESH /usr/local/include/SurfaceMesh)
    set(SURFACE_MESH_LIB /usr/local/lib/libSurfaceMesh.dylib)
    set(SPD_LOG_LIB /usr/local/lib/libspdlog.dylib)
    set(FMT_LIB /usr/local/Cellar/fmt/7.0.3/lib/libfmt.dylib)
elseif (UNIX)
    set(SURFACE_MESH /usr/local/include/SurfaceMesh)
    set(SURFACE_MESH_LIB /usr/local/lib/libSurfaceMesh.so)
endif ()

INCLUDE_DIRECTORIES(
        $ENV{HOME}/usr/include/
        ${PROJECT_SOURCE_DIR}/3rd/include
)

link_directories(
        $ENV{HOME}/usr/lib
        $ENV{HOME}/usr/Linux/64/gcc4/lib/
)

# SurfaceMesh
if (APPLE)
    include_directories(${SURFACE_MESH} ${INCLUDE})
elseif (MSVC)
    include_directories("$ENV{SURFACE_MESH_INCLUDE_DIR}")
    message("SurfaceMesh path = $ENV{SURFACE_MESH_INCLUDE_DIR}")
else ()
    include_directories(${SURFACE_MESH} ${INCLUDE})
endif ()
add_definitions("-DHEADERONLY")

# Eigen
if (MSVC)
    set(EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}")
    if (NOT EIGEN3_INCLUDE_DIR)
        MESSAGE(FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
    else ()
        MESSAGE("Eigen path = $ENV{EIGEN3_INCLUDE_DIR}")
    endif ()
    INCLUDE_DIRECTORIES("${EIGEN3_INCLUDE_DIR}")
else ()
    find_package(Eigen3 3.3.0)
    if (Eigen3_FOUND)
        include_directories(${EIGEN3_INCLUDE_DIR})
    else ()
        message(can't find eigen3)
        return()
    endif ()
endif ()

find_package(Boost COMPONENTS system thread program_options REQUIRED)
if (Boost_FOUND)
    message("find boost by find_package(), " PATH=${Boost_INCLUDE_DIRS})
    include_directories(${Boost_INCLUDE_DIRS})
else (Boost_FOUND)
    if (UNIX)
        message("find boost by UNIX, " PATH=$ENV{HOME}/boost/)
        link_directories($ENV{HOME}/boost/)
    elseif (WIN32)
        message("find boost by WIN32, " PATH=$ENV{HOME})
        include_directories($ENV{HOME}/boost/boost_1_58_0)
        link_directories($ENV{HOME}/boost/boost_1_58_0/lib64-msvc-12.0)
        include_directories($ENV{HOME}/boost_1_64_0)
        link_directories($ENV{HOME}/boost_1_64_0/lib64-msvc-14.1)
    elseif (APPLE)
        message("find boost by APPLE, " PATH=$ENV{HOME})

    endif ()
endif ()

if (NOT TARGET spdlog)
    find_package(spdlog REQUIRED)
endif ()

aux_source_directory(. DIR_SRC)

if (APPLE)
    INCLUDE_DIRECTORIES(
            /usr/local/include
    )

else ()
    INCLUDE_DIRECTORIES(
            $ENV{HOME}/usr/include/
    )
endif ()


if (APPLE)
    link_directories(
            /usr/local/lib
    )

else ()
    link_directories(
            $ENV{HOME}/usr/lib
            $ENV{HOME}/usr/Linux/64/gcc4/lib/
    )
endif ()

#
add_subdirectory(io)
add_subdirectory(common)
add_subdirectory(algorithm)


# google test
add_subdirectory(test)

add_executable(${PROJECT_NAME} ${DIR_SRC} algorithm/MeshSmooth.cpp algorithm/MeshSmooth.h algorithm/energy/quality_metric.cpp algorithm/energy/quality_metric.h algorithm/common.h algorithm/common.cpp)

target_link_libraries(${PROJECT_NAME}
        PUBLIC
        Boost::thread
        Boost::system
        Boost::program_options
        spdlog::spdlog
        common
        io
        alg
        )
if (MSVC)

else ()

endif ()


target_compile_definitions(${PROJECT_NAME} PUBLIC "-DDATA_SHARED_PATH=\"${DATA_SHARED_PATH}\"")
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})