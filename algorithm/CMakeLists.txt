file(GLOB HEADERS *.h *.hpp energy/*.h energy/*.hpp)
file(GLOB SOURCES *.cpp energy/*.cpp)

include_directories(../)

if (UNIX)
    add_library(alg STATIC ${SOURCES} ${HEADERS})
elseif (WIN32)
    add_library(alg STATIC ${SOURCES} ${HEADERS})
elseif (APPLE)
    add_library(alg SHARED ${SOURCES} ${HEADERS} common.h)
endif ()

target_link_libraries(alg
        ${SURFACE_MESH_LIB}
        spdlog::spdlog
        common
        io
        )

if (MSVC)
    target_link_libraries(alg
            "$ENV{SURFACE_MESH_LIB}"
            )
else ()
    target_link_libraries(alg
            ${SURFACE_MESH_LIB}
            )
endif ()

target_compile_definitions(alg PUBLIC "-DDATA_SHARED_PATH=\"${DATA_SHARED_PATH}\"")