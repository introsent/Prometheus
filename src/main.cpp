#include <iostream>
#include <chrono>
#include <thread>

// SDL2 (SDL_MAIN_HANDLED is defined via CMake so SDL won't override main)
#if defined(__has_include)
#  if __has_include(<SDL.h>)
#    include <SDL.h>
#  elif __has_include(<SDL3/SDL.h>)
#    include <SDL3/SDL.h>
#    include <SDL3/SDL_main.h>
#  else
#    error "SDL2 header not found"
#  endif
#else
#  include <SDL.h>
#endif

// GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>


// Embree (rtc api) - safe-guard: include only if available
#if defined(__has_include)
#  if __has_include(<embree4/rtcore.h>)
#    include <embree4/rtcore.h>
#    define HAVE_EMBREE 1
#  elif __has_include(<embree/rtcore.h>)
#    include <embree/rtcore.h>
#    define HAVE_EMBREE 1
#  else
#    define HAVE_EMBREE 0
#  endif
#endif

int main(int argc, char** argv)
{
    (void)argc;
    (void)argv;

    std::cout << "=== Prometheus small lib test ===\n";

    // GLM test
    {
        glm::vec3 v1(1.0f, 2.0f, 3.0f);
        glm::vec3 v2(4.0f, -1.0f, 0.5f);
        glm::vec3 v3 = v1 + v2;
        glm::mat4 m = glm::translate(glm::mat4(1.0f), glm::vec3(2.0f, 0.0f, 0.0f));
        glm::vec4 p = m * glm::vec4(v3, 1.0f);

        std::cout << "[GLM] v1 = " << glm::to_string(v1) << "\n";
        std::cout << "[GLM] v2 = " << glm::to_string(v2) << "\n";
        std::cout << "[GLM] v1+v2 = " << glm::to_string(v3) << "\n";
        std::cout << "[GLM] translated v3 = " << glm::to_string(p) << "\n";
    }

    // Embree test
#if HAVE_EMBREE
    {
        std::cout << "[Embree] initializing device...\n";
        RTCDevice device = rtcNewDevice(nullptr);
        if (!device) {
            std::cerr << "[Embree] Failed to create device (rtcNewDevice returned null)\n";
        } else {
            std::cout << "[Embree] device created\n";

            RTCScene scene = rtcNewScene(device);
            if (!scene) {
                std::cerr << "[Embree] rtcNewScene failed\n";
            } else {
                RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

                struct Vertex { float x, y, z, w; };
                struct Triangle { int v0, v1, v2; };

                Vertex* verts = (Vertex*)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, sizeof(Vertex), 3);
                verts[0] = {0.0f, 0.0f, 0.0f, 0.0f};
                verts[1] = {1.0f, 0.0f, 0.0f, 0.0f};
                verts[2] = {0.0f, 1.0f, 0.0f, 0.0f};

                Triangle* tris = (Triangle*)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, sizeof(Triangle), 1);
                tris[0] = {0, 1, 2};

                rtcCommitGeometry(geom);
                unsigned int geomID = rtcAttachGeometry(scene, geom);
                rtcReleaseGeometry(geom);

                rtcCommitScene(scene);
                std::cout << "[Embree] scene committed, geomID = " << geomID << "\n";

                rtcReleaseScene(scene);
            }
            rtcReleaseDevice(device);
            std::cout << "[Embree] device released\n";
        }
    }
#else
    std::cout << "[Embree] Embree headers not available at compile-time; embedding test skipped.\n";
#endif

    SDL_Window *win = NULL;
    SDL_Renderer *renderer = NULL;
    int width = 640, height = 420;
    bool loopShouldStop = false;

    if (SDL_Init(SDL_INIT_VIDEO) == 1)
    {
        std::cout << "[SDL3] SDL_Init() succeeded\n";
    }

    win = SDL_CreateWindow("Prometheus", width, height, 0);

    renderer = SDL_CreateRenderer(win, NULL);

    while (!loopShouldStop)
    {
        SDL_Event event;
        while (SDL_PollEvent(&event))
        {
            switch (event.type)
            {
                case SDL_EVENT_QUIT:
                    loopShouldStop = true;
                    break;
                default: ;
            }
        }

        SDL_RenderClear(renderer);
        SDL_RenderPresent(renderer);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(win);

    SDL_Quit();

    return 0;
}
