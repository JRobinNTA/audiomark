include_directories(${PORT_DIR})

set(PORT_SOURCE
    ${PORT_DIR}/th_api.c
    ${PORT_DIR}/th_cfft.c
    ${PORT_DIR}/th_cfft_tables.c
    # Include additional source files
)

add_definitions(-DUSE_CMSIS_DSP)
