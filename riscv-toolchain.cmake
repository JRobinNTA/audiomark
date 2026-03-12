set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR riscv64)

set(TOOLCHAIN_PATH "/home/johnr/Devel/lfx/toolchain-llvm")
set(GCC_TOOLCHAIN "/home/johnr/Devel/lfx/toolchain-gcc")

set(CMAKE_C_COMPILER "${TOOLCHAIN_PATH}/bin/clang")
set(CMAKE_CXX_COMPILER "${TOOLCHAIN_PATH}/bin/clang++")

# Set targets for BOTH C and C++
set(TRIPLE riscv64-unknown-linux-gnu)
set(CMAKE_C_COMPILER_TARGET ${TRIPLE})
set(CMAKE_CXX_COMPILER_TARGET ${TRIPLE})

set(CMAKE_SYSROOT "${GCC_TOOLCHAIN}/sysroot")

# Hint to Clang where the GCC headers/binaries are
set(CMAKE_C_FLAGS_INIT "--gcc-toolchain=${GCC_TOOLCHAIN} -menable-experimental-extensions -march=rv64gcvbp0p19 -mabi=lp64d")
set(CMAKE_CXX_FLAGS_INIT "--gcc-toolchain=${GCC_TOOLCHAIN} -menable-experimental-extensions -march=rv64gcvbp0p19 -mabi=lp64d")

set(CMAKE_EXE_LINKER_FLAGS_INIT "-fuse-ld=lld --gcc-toolchain=${GCC_TOOLCHAIN} -L${GCC_TOOLCHAIN}/sysroot/lib64/lp64d")
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)
