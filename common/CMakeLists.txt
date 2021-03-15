file(GLOB HEADERS *.h *.hpp)
file(GLOB SOURCES *.cpp)

if (UNIX)
    add_library(common STATIC ${SOURCES} ${HEADERS})
elseif (WIN32)
    add_library(common STATIC ${SOURCES} ${HEADERS})
elseif (APPLE)
    add_library(common SHARED ${SOURCES} ${HEADERS} common.h eigen_utils.cpp eigen_utils.h)
endif ()

target_link_libraries(common
        ${SURFACE_MESH_LIB}
        spdlog::spdlog
        )

if (MSVC)
    target_link_libraries(common
            "$ENV{SURFACE_MESH_LIB}"
            )
else ()
    target_link_libraries(common
            ${SURFACE_MESH_LIB}
            )
endif ()

target_compile_definitions(common PUBLIC "-DDATA_SHARED_PATH=\"${DATA_SHARED_PATH}\"")
