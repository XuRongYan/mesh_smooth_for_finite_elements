file(GLOB HEADERS *.h *.hpp)
file(GLOB SOURCES *.cpp)

if(UNIX)
  add_library(io STATIC ${SOURCES} ${HEADERS})
elseif(WIN32)
  add_library(io STATIC ${SOURCES} ${HEADERS})
elseif(APPLE)
  add_library(io SHARED ${SOURCES} ${HEADERS})
endif()

target_link_libraries(io
        Boost::thread
        Boost::system
        Boost::program_options
        spdlog::spdlog
        common
)

if (MSVC)
  target_link_libraries(io
          "$ENV{SURFACE_MESH_LIB}"
          )
else ()
  target_link_libraries(io
          ${SURFACE_MESH_LIB}
          )
endif ()

target_compile_definitions(io PUBLIC "-DDATA_SHARED_PATH=\"${DATA_SHARED_PATH}\"")