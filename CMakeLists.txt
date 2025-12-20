cmake_minimum_required(VERSION 3.24)
project(Prometheus VERSION 0.1 LANGUAGES CXX)

# High-performance compilation option
option(HIGH_PERFORMANCE "Enable high-performance compilation flags for 100+ FPS" ON)

# C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Make FetchContent available
include(FetchContent)

# Options (user can turn vendoring off)
option(USE_VENDORED_EMBREE "Download & build Embree from source when not found" ON)
option(USE_VENDORED_SDL3 "Download & build SDL3 from source when not found" ON)
option(USE_VENDORED_GLM    "Download & use GLM from source when not found" ON)


# ------------------------
# Try system packages first (QUIET), fallback to vendored via FetchContent
# ------------------------

# Embree
find_package(embree 4 QUIET)

if(NOT embree_FOUND AND USE_VENDORED_EMBREE)
    message(STATUS "Embree (v4) not found on system: fetching Embree v4.4.0 into")
    FetchContent_Declare(
            embree
            GIT_REPOSITORY https://github.com/RenderKit/embree.git
            GIT_TAG        v4.4.0
            GIT_SHALLOW    TRUE
    )

    # Simple vendored build: disable ISPC, tutorials, SYCL and install step
    set(EMBREE_ISPC_SUPPORT OFF CACHE BOOL "Disable ISPC" FORCE)
    set(EMBREE_TUTORIALS    OFF CACHE BOOL "Don't build examples" FORCE)
    set(EMBREE_SYCL_SUPPORT OFF CACHE BOOL "Disable SYCL" FORCE)
    set(EMBREE_INSTALL      OFF CACHE BOOL "Don't install embree" FORCE)
    set(EMBREE_TASKING_SYSTEM "INTERNAL" CACHE STRING "Use internal tasking" FORCE)

    # Restrict ISAs: enable only SSE4.2; disable AVX/AVX2/AVX512 to avoid building extra kernels:
    set(EMBREE_ISA_SSE2     ON  CACHE BOOL "Enable SSE2" FORCE)
    set(EMBREE_ISA_SSE42    ON  CACHE BOOL "Enable SSE4.2" FORCE)
    set(EMBREE_BACKFACE_CULLING ON CACHE BOOL "Enable backface culling" FORCE)
    set(EMBREE_BACKFACE_CULLING_SPHERES ON CACHE BOOL "Enable backface culling spheres" FORCE)
    set(EMBREE_ISA_AVX      OFF CACHE BOOL "Disable AVX" FORCE)
    set(EMBREE_ISA_AVX2     OFF CACHE BOOL "Disable AVX2" FORCE)
    set(EMBREE_ISA_AVX512SKX OFF CACHE BOOL "Disable AVX512 SKX" FORCE)
    set(EMBREE_ISA_AVX512KNL OFF CACHE BOOL "Disable AVX512 KNL" FORCE)

    # Avoid building shared libraries that can add Unix only linker flags which MinGW doesn't understand.
    if(MINGW OR (WIN32 AND CMAKE_C_COMPILER_ID STREQUAL "GNU"))
        message(STATUS "Detected MinGW/GNU on Windows => forcing static only Embree build to avoid linking issues")
        set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries" FORCE)
        set(EMBREE_BUILD_SHARED OFF CACHE BOOL "Build embree shared libs" FORCE)
    endif()

    FetchContent_MakeAvailable(embree)
endif()

# SDL3
find_package(SDL3 CONFIG QUIET)
if(NOT SDL3_FOUND)
    find_package(SDL3 QUIET) # plain mode fallback
endif()

if(NOT SDL3_FOUND AND USE_VENDORED_SDL3)
    message(STATUS "SDL3 not found on system: fetching SDL (release-3.2.20) into the build tree")
    include(FetchContent)
    FetchContent_Declare(
            sdl3
            GIT_REPOSITORY https://github.com/libsdl-org/SDL.git
            GIT_TAG        release-3.2.24
            GIT_SHALLOW    TRUE
    )

    # Configure vendored SDL3 options (tweak if you want different features)
    set(BUILD_SHARED_LIBS ON CACHE BOOL "Build SDL3 as shared library" FORCE)
    set(SDL_TEST OFF CACHE BOOL "Disable SDL tests" FORCE)
    # keep Vulkan off if you don't need it
    set(SDL_VULKAN OFF CACHE BOOL "Disable Vulkan support" FORCE)

    FetchContent_MakeAvailable(sdl3)
endif()

# GLM: header-only
find_package(glm CONFIG QUIET)

if(NOT glm_FOUND AND USE_VENDORED_GLM)
    message(STATUS "glm not found on system: fetching glm (v1.0.1) into the build tree")
    FetchContent_Declare(
            glm
            GIT_REPOSITORY https://github.com/g-truc/glm.git
            GIT_TAG        1.0.1
            GIT_SHALLOW    TRUE
    )
    # Don't let GLM's CMakeLists try to build libraries
    set(GLM_BUILD_LIBRARY OFF CACHE BOOL "Don't build GLM as a library" FORCE)
    set(GLM_BUILD_TESTS OFF CACHE BOOL "Don't build GLM tests" FORCE)
    FetchContent_MakeAvailable(glm)
endif()

# ------------------------
# Project targets
# ------------------------

add_executable(Prometheus
        src/main.cpp
        src/geometries/geometry.cpp
        src/geometries/geometry.h
        src/geometries/triangle.cpp
        src/geometries/triangle.h
        src/init/embree_scene.cpp
        src/init/embree_scene.h
        src/init/embree_device.cpp
        src/init/embree_device.h
        src/structs/vertex.h
        src/ray/ray.h
        src/ray/hit_result.h
        src/ray/ray_tracer.cpp
        src/ray/ray_tracer.h
        src/geometries/sphere.cpp
        src/geometries/sphere.h
        src/geometries/plane.cpp
        src/geometries/plane.h
        src/camera/camera.cpp
        src/camera/camera.h
        src/render/renderer.cpp
        src/render/renderer.h
        src/structs/material.h
        src/render/scene_manager.cpp
        src/render/scene_manager.h
        src/render/timer.cpp
        src/render/timer.h
        src/structs/light.h
        src/structs/global_indices.h
        external/tiny_obj_loader.h
        src/parser/obj_parser.cpp
        src/geometries/mesh.cpp
        src/geometries/mesh.h
        src/parser/obj_parser.h
        src/structs/area_light.cpp
        src/structs/area_light.h
)

# Copy resources to output folder
set(RESOURCES_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/resources")
file(GLOB_RECURSE RESOURCE_FILES
        "${RESOURCES_SOURCE_DIR}/*.jpg"
        "${RESOURCES_SOURCE_DIR}/*.png"
        "${RESOURCES_SOURCE_DIR}/*.obj"
)
set(RESOURCES_OUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/resources/")
file(MAKE_DIRECTORY ${RESOURCES_OUT_DIR})
foreach(RESOURCE ${RESOURCE_FILES})
    add_custom_command(TARGET Prometheus POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy ${RESOURCE}
            ${RESOURCES_OUT_DIR})
endforeach(RESOURCE)

# High-performance flags - TARGET SPECIFIC
if(HIGH_PERFORMANCE)
    if(MSVC)
        # For MSVC: Only apply to Release builds, not Debug
        target_compile_options(Prometheus PRIVATE
                $<$<CONFIG:Release>:/O2>
                $<$<CONFIG:RelWithDebInfo>:/O2>
                $<$<CONFIG:Release>:/arch:AVX2>
                $<$<CONFIG:RelWithDebInfo>:/arch:AVX2>
        )
        target_link_options(Prometheus PRIVATE
                $<$<CONFIG:Release>:/LTCG>
                $<$<CONFIG:RelWithDebInfo>:/LTCG>
        )
    else()
        # For GCC/MinGW: Use simpler optimizations that work
        target_compile_options(Prometheus PRIVATE
                $<$<CONFIG:Release>:-O3>
                $<$<CONFIG:RelWithDebInfo>:-O3>
                $<$<CONFIG:Release>:-march=native>
                $<$<CONFIG:RelWithDebInfo>:-march=native>
                $<$<CONFIG:Release>:-ffast-math>
                $<$<CONFIG:RelWithDebInfo>:-ffast-math>
        )
        # Only enable LTO if compiler supports it
        include(CheckCXXCompilerFlag)
        check_cxx_compiler_flag(-flto COMPILER_SUPPORTS_LTO)
        if(COMPILER_SUPPORTS_LTO)
            target_compile_options(Prometheus PRIVATE
                    $<$<CONFIG:Release>:-flto>
                    $<$<CONFIG:RelWithDebInfo>:-flto>
            )
            set_property(TARGET Prometheus PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
            set_property(TARGET Prometheus PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO TRUE)
        else()
            message(STATUS "LTO not supported by compiler, skipping...")
        endif()
    endif()

    target_compile_definitions(Prometheus PRIVATE PARALLEL_EXECUTION)
    message(STATUS "High-performance compilation flags enabled for Prometheus target")
endif()


target_include_directories(Prometheus PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/camera
        ${CMAKE_CURRENT_SOURCE_DIR}/src/geometries
        ${CMAKE_CURRENT_SOURCE_DIR}/src/init
        ${CMAKE_CURRENT_SOURCE_DIR}/src/ray
        ${CMAKE_CURRENT_SOURCE_DIR}/src/structs
        ${CMAKE_CURRENT_SOURCE_DIR}/external
)

# Ensure GLM experimental features allowed for all compilation units of this target
target_compile_definitions(Prometheus PRIVATE GLM_ENABLE_EXPERIMENTAL SDL_MAIN_HANDLED)

# Link Embree
if(TARGET embree)
    target_link_libraries(Prometheus PRIVATE embree)
elseif(TARGET embree::embree)
    target_link_libraries(Prometheus PRIVATE embree::embree)
else()
    if(embree_FOUND)
        if(TARGET ${embree_LIBRARIES})
            target_link_libraries(Prometheus PRIVATE ${embree_LIBRARIES})
        else()
            if(DEFINED EMBREE_LIBRARIES)
                target_link_libraries(Prometheus PRIVATE ${EMBREE_LIBRARIES})
            endif()
        endif()
    else()
        message(WARNING "Embree not found. Prometheus will compile but without Embree features.")
    endif()
endif()

# Link SDL3
if(TARGET SDL3::SDL3)
    target_link_libraries(Prometheus PRIVATE SDL3::SDL3)
else()
    if(SDL3_FOUND)
        if(DEFINED SDL3_LIBRARIES)
            target_link_libraries(Prometheus PRIVATE ${SDL3_LIBRARIES})
        endif()
        if(DEFINED SDL3_INCLUDE_DIRS)
            target_include_directories(Prometheus PRIVATE ${SDL3_INCLUDE_DIRS})
        endif()
    else()
        message(WARNING "SDL3 not found. Window/display functionality will be unavailable.")
    endif()
endif()

# GLM (header-only - just add include directories, never link)
if(TARGET glm::glm)
    # Modern GLM installation with proper interface target
    target_link_libraries(Prometheus PRIVATE glm::glm)
elseif(glm_FOUND)
    # Found via find_package - add include directories only
    if(DEFINED GLM_INCLUDE_DIRS)
        target_include_directories(Prometheus PRIVATE ${GLM_INCLUDE_DIRS})
    elseif(DEFINED GLM_INCLUDE_DIR)
        target_include_directories(Prometheus PRIVATE ${GLM_INCLUDE_DIR})
    endif()
else()
    # FetchContent case - manually add include directory
    if(DEFINED glm_SOURCE_DIR)
        target_include_directories(Prometheus PRIVATE ${glm_SOURCE_DIR})
    endif()
endif()

# Compiler warnings
if(MSVC)
    target_compile_options(Prometheus PRIVATE /W4 /permissive-)
else()
    target_compile_options(Prometheus PRIVATE -Wall -Wextra -Wpedantic)
endif()

# If SDL3 is a shared target, copy runtime DLL next to exe (Windows)
if(TARGET SDL3::SDL3)
    add_custom_command(TARGET Prometheus POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:SDL3::SDL3>
            $<TARGET_FILE_DIR:Prometheus>
            COMMENT "Copying SDL3 runtime to executable output directory"
    )
endif()

if(WIN32)
    set(EMBREE_TARGET_NAME "")

    if(TARGET embree::embree)
        set(EMBREE_TARGET_NAME embree::embree)
    elseif(TARGET embree)
        set(EMBREE_TARGET_NAME embree)
    endif()

    # If we have embree.dll, copy it.
    if(EMBREE_TARGET_NAME)
        add_custom_command(TARGET Prometheus POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:${EMBREE_TARGET_NAME}>
                $<TARGET_FILE_DIR:Prometheus>
                COMMENT "Copying Embree runtime DLL to executable output directory"
        )

        # Ensure embree is built before Prometheus
        add_dependencies(Prometheus ${EMBREE_TARGET_NAME})
    else()
        # If non worked - try the typical FetchContent in _deps path used by ExternalProject/FetchContent
        set(_possible_embree_dll "${CMAKE_BINARY_DIR}/_deps/embree-build/embree4.dll")
        if(EXISTS "${_possible_embree_dll}")
            add_custom_command(TARGET Prometheus POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    "${_possible_embree_dll}"
                    $<TARGET_FILE_DIR:Prometheus>
                    COMMENT "Copying vendored embree4.dll to executable output directory (fallback)"
            )
        endif()
    endif()
endif()

# Install target
install(TARGETS Prometheus RUNTIME DESTINATION bin)


message(STATUS "CMake configuration done. Configure & build: cmake -S . -B build && cmake --build build --config Release")