cmake_minimum_required(VERSION 3.16)
project(cactus_md VERSION 1.0.0 LANGUAGES CXX)

# ── Compiler settings ──────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Optional: enable warnings
if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# ── Library: md ───────────────────────────────────────────────────────────────
# PeriodicTable is header-only, so only PDBReader.cpp needs compiling.
add_library(md
    src/PDBReader.cpp)

target_include_directories(md
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# ── Executable ────────────────────────────────────────────────────────────────
add_executable(cactus_resid
    src/main.cpp)

target_link_libraries(cactus_resid PRIVATE md)
