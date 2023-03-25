set(mx3d_FOUND TRUE)

get_filename_component(MX3d_DIR ${MX3d_DIRS} DIRECTORY)
set(MX3d_INCLUDE_DIRS "${MX3d_DIR}/include")
include_directories(${MX3d_INCLUDE_DIRS})
