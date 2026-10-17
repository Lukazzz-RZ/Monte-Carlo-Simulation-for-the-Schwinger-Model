cmake_minimum_required(VERSION 3.20)

project(SchwingerWorm LANGUAGES C)

# The current implementation uses Windows APIs (_mkdir, Sleep and rand_s).
# A port to other operating systems requires changes to the C sources.
if(NOT WIN32)
    message(FATAL_ERROR "SchwingerWorm currently requires Windows (MSVC).")
endif()

if(NOT MSVC)
    message(FATAL_ERROR "The current build is supported with Microsoft Visual C only.")
endif()

add_executable(SchwingerWorm
    src/Main.c
    src/Initializer.c
    src/UpdCore.c
    src/Plaq_mod.c
    src/Worm_ferm.c
    src/Measures.c
    src/Visualizer.c
    src/Other.c
)

target_include_directories(SchwingerWorm PRIVATE include)
set_target_properties(SchwingerWorm PROPERTIES
    C_STANDARD 11
    C_STANDARD_REQUIRED ON
    C_EXTENSIONS OFF
)

# Microsoft's rand_s declaration is enabled by this preprocessor definition.
target_compile_definitions(SchwingerWorm PRIVATE _CRT_RAND_S)
target_compile_options(SchwingerWorm PRIVATE /W3)

