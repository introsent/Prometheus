cmake_minimum_required(VERSION 3.24)
project(Prometheus VERSION 0.1 LANGUAGES CXX)

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
    set(EMBREE_ISPC_SUPPORT OFF CACHE BOOL "Disable ISPC to avoid ISPC tool requirement" FORCE)
    set(EMBREE_TUTORIALS    OFF CACHE BOOL "Don't build embree tutorials/examples" FORCE)
    set(EMBREE_SYCL_SUPPORT OFF CACHE BOOL "Disable SYCL support for vendored build" FORCE)
    set(EMBREE_INSTALL      OFF CACHE BOOL "Don't install embree during vendored build" FORCE)

    # Use internal tasking to avoid requiring oneTBB
    set(EMBREE_TASKING_SYSTEM "INTERNAL" CACHE STRING "Embree tasking system (TBB | INTERNAL)" FORCE)

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
)

target_include_directories(Prometheus PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/camera
        ${CMAKE_CURRENT_SOURCE_DIR}/src/geometries
        ${CMAKE_CURRENT_SOURCE_DIR}/src/init
        ${CMAKE_CURRENT_SOURCE_DIR}/src/ray
        ${CMAKE_CURRENT_SOURCE_DIR}/src/structs
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

# GLM (header only)
if(TARGET glm::glm)
    target_link_libraries(Prometheus PRIVATE glm::glm)
else()
    if(glm_FOUND)
        if(DEFINED GLM_INCLUDE_DIR)
            target_include_directories(Prometheus PRIVATE ${GLM_INCLUDE_DIR})
        elseif(DEFINED GLM_INCLUDE_DIRS)
            target_include_directories(Prometheus PRIVATE ${GLM_INCLUDE_DIRS})
        endif()
    else()
        message(WARNING "glm not found - math helpers will be unavailable.")
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
