//
// Created by minaj on 10/27/2025.
//

#ifndef RENDERER_H
#define RENDERER_H

#include <chrono>
#include <string>
#include <SDL3/SDL.h>
#include <glm/glm.hpp>
#include <vector>
#include <cstdint>
#include "camera.h"
#include "scene_manager.h"
#include "ray_tracer.h"
#include "sampler_base.h"

class Renderer {
public:
    Renderer(int width, int height);
    ~Renderer();

    void setTestMode(bool enabled);
    bool isTestComplete() const;

    // Save current framebuffer (m_pixels) to BMP with strategy prefix.
    void saveScreenshot(const std::string &filename) const;

    // Overload: save arbitrary buffer
    void saveScreenshot(const std::vector<uint32_t>& pixels, const std::string& filename) const;

    bool initialize();
    void render(const Camera& camera, const SceneManager& scene);
    void present();
    bool shouldQuit();

    void setAreaLightStrategy(SceneManager& scene, SamplingStrategy strategy);

    // helper: render into an arbitrary pixel buffer (used for ground truth rendering)
    void renderSceneIntoBuffer(const Camera& camera, const SceneManager& scene,
                               std::vector<uint32_t>& outPixels, int areaLightSamples,
                               std::chrono::duration<double>* outDuration = nullptr);

    // MSE helper
    double computeMSE(const std::vector<uint32_t>& reference,
                      const std::vector<uint32_t>& test) const;

    // convert sRGB channel to linear (declared inline)
    static inline float srgbToLinear(float c);

    // Save / load raw ground truth buffer (binary)
    void saveRawBuffer(const std::vector<uint32_t>& pixels,
                       const std::string& path) const;
    bool loadRawBuffer(std::vector<uint32_t>& pixels,
                       const std::string& path) const;

    // Generate a canonical ground truth file (writes gt.bin + gt.bmp)
    // - samplerStrategy: which sampling strategy to use for GT
    // - samples: how many area-light samples to use for GT
    // - outDir: directory where gt.bin / gt.bmp will be written
    void generateGroundTruth(const Camera& camera, SceneManager& scene,
                             SamplingStrategy samplerStrategy, int samples,
                             const std::string& outDir);


    void setSceneName(const std::string& sceneName) {
        m_sceneName = sceneName;
    }

private:
    int m_width;
    int m_height;
    SDL_Window* m_window;
    SDL_Renderer* m_renderer;
    SDL_Texture* m_texture;

    // framebuffer stored as 32-bit ARGB pixels.
    std::vector<uint32_t> m_pixels;

    bool m_quit;

    // test mode flag
    bool m_testMode;
    std::string m_sceneName;

    int m_currentSamples;
    int m_maxSamples;
    std::vector<int> m_sampleCounts;

    SamplingStrategy m_currentStrategy;
    std::string m_strategyName;

    std::string m_testFolder;
    void createTestFolder();
    [[nodiscard]] std::string getTestFolderPath() const;
    void createTestSummary() const;

    // ground truth storage (persisted in memory; we also save/load raw file on disk)
    std::vector<uint32_t> m_groundTruthPixels;
    std::string m_groundTruthDir = "tests/ground_truth";
    std::string m_groundTruthImage = "gt.bmp";
    std::string m_groundTruthRaw = "gt.bin";
    bool m_groundTruthComputed = false;
    int m_groundTruthSamples = 20000;
};
#endif //RENDERER_H
