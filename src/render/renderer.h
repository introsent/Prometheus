//
// Created by minaj on 10/27/2025.
//

#ifndef RENDERER_H
#define RENDERER_H

#include <SDL3/SDL.h>
#include <glm/glm.hpp>
#include <vector>
#include "camera.h"
#include "scene_manager.h"
#include "ray_tracer.h"

// Define this to enable parallel execution
#define PARALLEL_EXECUTION

class Renderer {
public:
    Renderer(int width, int height);
    ~Renderer();

    bool initialize();
    void render(const Camera& camera, const SceneManager& scene);
    void present();
    bool shouldQuit();

private:
    void renderPixel(const Camera& camera, const SceneManager& scene,
                     const RayTracer& tracer, uint32_t pixelIndex);

    static glm::vec3 traceRay(const Ray& ray, const RayTracer& tracer, const SceneManager& scene);

    static SDL_Color toSDLColor(const glm::vec3& color);

    int m_width;
    int m_height;
    SDL_Window* m_window;
    SDL_Renderer* m_renderer;
    SDL_Texture* m_texture;
    std::vector<Uint32> m_pixels;
    bool m_quit;
};
#endif //RENDERER_H
