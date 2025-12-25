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

    void setTestMode(bool enabled);

    bool isTestComplete() const;

    void saveScreenshot(const std::string &filename) const;

    bool initialize();
    void render(const Camera& camera, const SceneManager& scene);
    void present();
    bool shouldQuit();

private:
    static SDL_Color toSDLColor(const glm::vec3& color);

    int m_width;
    int m_height;
    SDL_Window* m_window;
    SDL_Renderer* m_renderer;
    SDL_Texture* m_texture;
    std::vector<Uint32> m_pixels;
    bool m_quit;
    int m_testMode;
    int m_currentSamples;
    int m_maxSamples;
    std::vector<int> m_sampleCounts;

    std::string m_testFolder;
    void createTestFolder();
    [[nodiscard]] std::string getTestFolderPath() const;
    void createTestSummary() const;
};
#endif //RENDERER_H
