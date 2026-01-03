//
// Created by minaj on 10/27/2025.
//

#include "renderer.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <random>
#include <chrono>
#include <filesystem>
#include "sampler_base.h"
#include "lights/mesh_area_light.h"
#include "lights/triangle_area_light.h"
#include "mis/bsdf_sampler.h"
#include "samplers/visibility_aware_sampler.h"
#include <cmath>

#if defined(PARALLEL_EXECUTION)
#include <execution>
#endif

class RandomGenerator {
public:
    RandomGenerator() : m_gen(std::random_device{}()), m_dist(0.0f, 1.0f) {}

    float get() { return m_dist(m_gen); }

private:
    std::mt19937 m_gen;
    std::uniform_real_distribution<float> m_dist;
};

inline float Renderer::srgbToLinear(float c) {
    if (c <= 0.04045f) return c / 12.92f;
    return std::pow((c + 0.055f) / 1.055f, 2.4f);
}

Renderer::Renderer(const int width, const int height)
    : m_width(width), m_height(height), m_window(nullptr),
      m_renderer(nullptr), m_texture(nullptr), m_quit(false),
      m_testMode(false), m_currentSamples(1), m_maxSamples(256), m_currentStrategy() {
    m_pixels.resize(width * height);

    m_sampleCounts = {1, 2, 4, 8, 16, 32, 64, 128, 256};
}

Renderer::~Renderer() {
    if (m_texture) SDL_DestroyTexture(m_texture);
    if (m_renderer) SDL_DestroyRenderer(m_renderer);
    if (m_window) SDL_DestroyWindow(m_window);
    SDL_Quit();
}

void Renderer::setTestMode(bool enabled) {
    m_testMode = enabled;
    if (m_testMode) {
        m_currentSamples = 1;
        createTestFolder();  // create folder when test starts

        std::cout << "\n=== Starting Area Light Sampling Test ===" << std::endl;
        std::cout << "Test folder: " << m_testFolder << std::endl;
        std::cout << "Testing " << m_sampleCounts.size() << " sample counts:" << std::endl;
        for (size_t i = 0; i < m_sampleCounts.size(); ++i) {
            std::cout << "  " << i+1 << ". " << m_sampleCounts[i] << " samples" << std::endl;
        }
        std::cout << "Results will be saved in: " << m_testFolder << std::endl;
        std::cout << "==========================================\n" << std::endl;
    }
}

bool Renderer::isTestComplete() const {
    return m_testMode && m_currentSamples > m_maxSamples;
}

void Renderer::saveScreenshot(const std::string& filename) const {
    std::string fullFilename = m_strategyName + "_" + filename;
    const std::string filepath = m_testFolder + "/" + fullFilename;

    // create an SDL_Surface from our pixel data
    SDL_Surface* surface = SDL_CreateSurface(m_width, m_height, SDL_PIXELFORMAT_ARGB8888);
    if (!surface) {
        std::cerr << "Failed to create surface: " << SDL_GetError() << std::endl;
        return;
    }

    // copy pixel data from m_pixels to surface
    for (int y = 0; y < m_height; ++y) {
        for (int x = 0; x < m_width; ++x) {
            uint32_t pixel = m_pixels[y * m_width + x];

            // get color components
            uint8_t a = (pixel >> 24) & 0xFF;
            uint8_t r = (pixel >> 16) & 0xFF;
            uint8_t g = (pixel >> 8) & 0xFF;
            uint8_t b = pixel & 0xFF;

            // write to surface (adjust for endianness)
            auto* target_pixel = reinterpret_cast<uint32_t *>(static_cast<uint8_t *>(surface->pixels) + y * surface->pitch + x * 4);

            // convert to surface format (ARGB)
            *target_pixel = (a << 24) | (r << 16) | (g << 8) | b;
        }
    }

    // Save the surface as a BMP file
    if (SDL_SaveBMP(surface, filepath.c_str()) != 0) {
        std::cerr << "Failed to save screenshot " << filepath << ": " << SDL_GetError() << std::endl;
    } else {
        std::cout << "Saved: " << filepath << std::endl;
    }

    SDL_DestroySurface(surface);
}

void Renderer::saveScreenshot(const std::vector<uint32_t>& pixels, const std::string& filename) const {
    std::string fullFilename = m_strategyName + "_" + filename;
    const std::string filepath = m_testFolder + "/" + fullFilename;

    SDL_Surface* surface = SDL_CreateSurface(m_width, m_height, SDL_PIXELFORMAT_ARGB8888);
    if (!surface) {
        std::cerr << "Failed to create surface: " << SDL_GetError() << std::endl;
        return;
    }

    for (int y = 0; y < m_height; ++y) {
        for (int x = 0; x < m_width; ++x) {
            uint32_t pixel = pixels[y * m_width + x];
            uint8_t a = (pixel >> 24) & 0xFF;
            uint8_t r = (pixel >> 16) & 0xFF;
            uint8_t g = (pixel >> 8) & 0xFF;
            uint8_t b = pixel & 0xFF;
            auto* target_pixel = reinterpret_cast<uint32_t *>(static_cast<uint8_t *>(surface->pixels) + y * surface->pitch + x * 4);
            *target_pixel = (a << 24) | (r << 16) | (g << 8) | b;
        }
    }

    if (SDL_SaveBMP(surface, filepath.c_str()) != 0) {
        std::cerr << "Failed to save screenshot " << filepath << ": " << SDL_GetError() << std::endl;
    } else {
        std::cout << "Saved: " << filepath << std::endl;
    }
    SDL_DestroySurface(surface);
}


bool Renderer::initialize() {
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        std::cerr << "SDL initialization failed: " << SDL_GetError() << std::endl;
        return false;
    }

    m_window = SDL_CreateWindow("Prometheus - Area Light Sampling Test",
                                 m_width, m_height, SDL_WINDOW_RESIZABLE);
    if (!m_window) {
        std::cerr << "Window creation failed: " << SDL_GetError() << std::endl;
        return false;
    }

    m_renderer = SDL_CreateRenderer(m_window, NULL);
    if (!m_renderer) {
        std::cerr << "Renderer creation failed: " << SDL_GetError() << std::endl;
        return false;
    }

    m_texture = SDL_CreateTexture(m_renderer, SDL_PIXELFORMAT_ARGB8888,
                                   SDL_TEXTUREACCESS_STREAMING, m_width, m_height);
    if (!m_texture) {
        std::cerr << "Texture creation failed: " << SDL_GetError() << std::endl;
        return false;
    }

    return true;
}


void Renderer::render(const Camera& camera, const SceneManager& scene) {
    auto renderStart = std::chrono::high_resolution_clock::now();

    /// GROUND TRUTH HANDLING
    // if test mode and ground truth not computed yet -> try to load it from disk first
    if (m_testMode && !m_groundTruthComputed) {
        // ensure ground-truth dir exists
        try {
            if (!std::filesystem::exists(m_groundTruthDir)) {
                std::filesystem::create_directories(m_groundTruthDir);
            }
        } catch (...) {}

        const std::string rawPath = m_groundTruthDir + "/" + m_groundTruthRaw;

        if (std::filesystem::exists(rawPath)) {
            if (loadRawBuffer(m_groundTruthPixels, rawPath)) {
                m_groundTruthComputed = true;
                std::cout << "Loaded ground truth from: " << rawPath << std::endl;
            } else {
                std::cerr << "Found ground-truth file but failed to load: " << rawPath << std::endl;
            }
        }
    }

    // if still not computed and test-mode -> generate ground truth
    if (m_testMode && !m_groundTruthComputed) {
        std::cout << "Rendering ground-truth with " << m_groundTruthSamples
                  << " samples (this may take a while)..." << std::endl;

        std::chrono::duration<double> gtDuration{};
        renderSceneIntoBuffer(camera, scene, m_groundTruthPixels, m_groundTruthSamples, &gtDuration);

        // save raw + bmp to disk so subsequent runs can reuse it
        try {
            if (!std::filesystem::exists(m_groundTruthDir)) {
                std::filesystem::create_directories(m_groundTruthDir);
            }
        } catch (...) {}

        const std::string rawPath = m_groundTruthDir + "/" + m_groundTruthRaw;
        saveRawBuffer(m_groundTruthPixels, rawPath);

        // temporarily set m_testFolder so saveScreenshot writes to the groundtruth dir
        const std::string oldTestFolder = m_testFolder;
        m_testFolder = m_groundTruthDir;
        saveScreenshot(m_groundTruthPixels, m_groundTruthImage);
        m_testFolder = oldTestFolder;

        // save ground truth timing
        std::ofstream timingFile(getTestFolderPath() + "/timing_results.txt", std::ios::app);
        if (timingFile.is_open()) {
            timingFile << "GT\t" << std::fixed << std::setprecision(3)
                      << gtDuration.count() << "\t-" << std::endl;
            timingFile.close();
        }

        m_groundTruthComputed = true;
        std::cout << "Ground-truth saved. Proceeding with tests." << std::endl;
    }

    /// MAIN RENDER
    int areaLightSamples = m_testMode ? m_currentSamples : 32;

    if (m_testMode) {
        MeshAreaLight::clearVisibilityCache();
        std::cout << "Rendering with " << areaLightSamples << " samples..." << std::endl;
    }

    std::chrono::duration<double> renderDuration{};
    renderSceneIntoBuffer(camera, scene, m_pixels, areaLightSamples, &renderDuration);

    // update texture for display
    SDL_UpdateTexture(m_texture, nullptr, m_pixels.data(), m_width * sizeof(Uint32));

    /// TEST MODE: SAVE RESULTS
    if (m_testMode) {
        // Save screenshot
        std::ostringstream filename;
        filename << "samples_" << std::setw(5) << std::setfill('0')
                 << areaLightSamples << ".bmp";
        saveScreenshot(filename.str());

        // calculate MSE against ground truth
        double mse = -1.0;
        if (m_groundTruthComputed && !m_groundTruthPixels.empty()) {
            mse = computeMSE(m_groundTruthPixels, m_pixels);
            std::cout << "  MSE: " << std::scientific << std::setprecision(6) << mse << std::endl;
        }

        // save timing and MSE data
        std::string resultsFilePath = getTestFolderPath() + "/results.txt";

        // write header if file doesn't exist
        bool writeHeader = !std::filesystem::exists(resultsFilePath);

        std::ofstream resultsFile(resultsFilePath, std::ios::app);
        if (resultsFile.is_open()) {
            if (writeHeader) {
                resultsFile << "Strategy\tSamples\tTime(s)\tMSE" << std::endl;
                resultsFile << "----------------------------------------" << std::endl;
            }

            resultsFile << m_strategyName << "\t"
                       << areaLightSamples << "\t"
                       << std::fixed << std::setprecision(3) << renderDuration.count() << "\t"
                       << std::scientific << std::setprecision(6) << mse << std::endl;
            resultsFile.close();
        }

        std::cout << "  Time: " << std::fixed << std::setprecision(3)
                 << renderDuration.count() << " seconds" << std::endl;

        // move to next sample count
        auto it = std::ranges::find(m_sampleCounts, m_currentSamples);
        if (it != m_sampleCounts.end() && ++it != m_sampleCounts.end()) {
            m_currentSamples = *it;
        } else {
            m_currentSamples = m_maxSamples + 1; // Mark as complete
            std::cout << "\n=== Area Light Sampling Test Complete ===" << std::endl;
            std::cout << "All screenshots saved to: " << m_testFolder << std::endl;
            std::cout << "Results saved to: " << resultsFilePath << std::endl;

            // create summary file with MSE analysis
            createTestSummary();
        }
    }
}


void Renderer::createTestSummary() const {
    const std::string summaryPath = getTestFolderPath() + "/test_summary.txt";

    if (std::ofstream summary(summaryPath); summary.is_open()) {
        summary << "=== Area Light Sampling Test Summary ===" << std::endl;
        summary << "Strategy: " << m_strategyName << std::endl;
        summary << "Test conducted: " << m_testFolder << std::endl;
        summary << "Resolution: " << m_width << "x" << m_height << std::endl;
        summary << "Ground truth samples: " << m_groundTruthSamples << std::endl;
        summary << "Total tests: " << m_sampleCounts.size() << std::endl;
        summary << "Sample counts tested: ";
        for (size_t i = 0; i < m_sampleCounts.size(); ++i) {
            summary << m_sampleCounts[i];
            if (i < m_sampleCounts.size() - 1) summary << ", ";
        }
        summary << std::endl << std::endl;

        // read and include results
        std::string resultsPath = getTestFolderPath() + "/results.txt";
        if (std::filesystem::exists(resultsPath)) {
            summary << "=== Results ===" << std::endl;
            std::ifstream resultsFile(resultsPath);
            std::string line;
            while (std::getline(resultsFile, line)) {
                summary << line << std::endl;
            }
            summary << std::endl;
        }

        summary << "=== Analysis Guide ===" << std::endl;
        summary << "1. Lower MSE = better quality (closer to ground truth)" << std::endl;
        summary << "2. Compare MSE vs Time for efficiency analysis" << std::endl;
        summary << "3. Look for diminishing returns at higher sample counts" << std::endl;
        summary << "4. Compare across strategies at equal sample counts" << std::endl;

        summary.close();
        std::cout << "Test summary saved to: " << summaryPath << std::endl;
    }
}

void Renderer::present() {
    SDL_RenderClear(m_renderer);
    SDL_RenderTexture(m_renderer, m_texture, nullptr, nullptr);

    // Display current sample count on window title
    if (m_testMode && m_currentSamples <= m_maxSamples) {
        std::ostringstream title;
        title << "Prometheus - Testing: " << m_currentSamples << " samples";
        SDL_SetWindowTitle(m_window, title.str().c_str());
    }

    SDL_RenderPresent(m_renderer);
}

bool Renderer::shouldQuit() {
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
        if (event.type == SDL_EVENT_QUIT) {
            m_quit = true;
        }
        if (event.type == SDL_EVENT_KEY_DOWN && event.key.key == SDLK_ESCAPE) {
            m_quit = true;
            std::cout << "Test interrupted by user." << std::endl;
        }
        // space bar to pause/resume test
        if (event.type == SDL_EVENT_KEY_DOWN && event.key.key == SDLK_SPACE) {
            m_testMode = !m_testMode;
            std::cout << "Test mode: " << (m_testMode ? "ON" : "OFF") << std::endl;
        }
    }
    return m_quit || (m_testMode && isTestComplete());
}

void Renderer::setAreaLightStrategy(SceneManager& scene, SamplingStrategy strategy) {
    m_currentStrategy = strategy;

    // set strategy name for display
    switch (strategy) {
        case SamplingStrategy::Uniform:
            m_strategyName = "Uniform";
            break;
        case SamplingStrategy::AreaImportance:
            m_strategyName = "AreaImportance";
            break;
        case SamplingStrategy::HierarchicalFlux:
            m_strategyName = "HierarchicalFlux";
            break;
        case SamplingStrategy::VisibilityAwareHierarchical:
            m_strategyName = "VisibilityAware";
            break;
    }

    std::cout << "Setting all area lights to: " << m_strategyName << std::endl;

    auto& lights = scene.getLights();
    int count = 0;

    for (auto& light : lights) {
        if (light->type == LightType::MeshArea && light->meshAreaLight) {
            light->meshAreaLight->setSamplingStrategy(strategy);
            count++;
        } else if (light->type == LightType::TriangleArea && light->triangleAreaLight) {
            light->triangleAreaLight->setSamplingStrategy(strategy);
            count++;
        }
    }

    std::cout << "  Applied to " << count << " area lights" << std::endl;
}

void Renderer::createTestFolder() {
    // Create a timestamp for the folder name
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    std::tm tm_buf;

#ifdef _WIN32
    localtime_s(&tm_buf, &in_time_t);
#else
    localtime_r(&in_time_t, &tm_buf);
#endif

    std::ostringstream folderName;
    folderName << std::put_time(&tm_buf, "%Y%m%d_%H%M%S");

    // Create scene-based directory structure
    m_testFolder = "tests/" + m_sceneName + "/" + m_strategyName + "/" + folderName.str();

    // create the directory
    try {
        std::filesystem::create_directories(m_testFolder);
        std::cout << "Created test folder: " << m_testFolder << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Failed to create test folder: " << e.what() << std::endl;
        m_testFolder = ".";  // fall back to current directory
    }
}

std::string Renderer::getTestFolderPath() const {
    return m_testFolder.empty() ? "." : m_testFolder;
}

/// Compute MSE between reference and test
// converts to linear RGB before computing squared error.
double Renderer::computeMSE(const std::vector<uint32_t>& reference,
                            const std::vector<uint32_t>& test) const {
    if (reference.size() != test.size() || reference.empty()) return -1.0;

    const size_t N = reference.size();
    double sumSq = 0.0;
    for (size_t i = 0; i < N; ++i) {
        uint32_t refPx = reference[i];
        uint32_t tPx = test[i];

        float rr = ((refPx >> 16) & 0xFF) / 255.0f;
        float rg = ((refPx >> 8) & 0xFF) / 255.0f;
        float rb = (refPx & 0xFF) / 255.0f;

        float tr = ((tPx >> 16) & 0xFF) / 255.0f;
        float tg = ((tPx >> 8) & 0xFF) / 255.0f;
        float tb = (tPx & 0xFF) / 255.0f;

        // linearize
        rr = srgbToLinear(rr);
        rg = srgbToLinear(rg);
        rb = srgbToLinear(rb);

        tr = srgbToLinear(tr);
        tg = srgbToLinear(tg);
        tb = srgbToLinear(tb);

        const auto dr = static_cast<double>(rr - tr);
        const auto dg = static_cast<double>(rg - tg);
        const auto db = static_cast<double>(rb - tb);

        sumSq += dr * dr + dg * dg + db * db;
    }

    // MSE per channel (average over pixels and over channels)
    const double mse = sumSq / (static_cast<double>(N) * 3.0);
    return mse;
}

/// Render the scene into an arbitrary pixel buffer (very similar to existing render loop).
// if outDuration is non-null, we measure render time and write it back.
void Renderer::renderSceneIntoBuffer(const Camera& camera, const SceneManager& scene,
                                     std::vector<uint32_t>& outPixels, int areaLightSamples,
                                     std::chrono::duration<double>* outDuration) {
    // Prepare output buffer
    outPixels.assign(m_width * m_height, 0);

    auto start = std::chrono::high_resolution_clock::now();

    const RayTracer tracer(&scene);
    const glm::vec3 cameraPos = camera.getPosition();
    const auto& lights = scene.getLights();

    const uint32_t tilesX = (m_width + 7) / 8;
    const uint32_t tilesY = (m_height + 7) / 8;
    const uint32_t totalTiles = tilesX * tilesY;

    std::vector<uint32_t> tileIndices(totalTiles);
    for (uint32_t i = 0; i < totalTiles; ++i) tileIndices[i] = i;

#if defined(PARALLEL_EXECUTION)
    std::for_each(std::execution::par, tileIndices.begin(), tileIndices.end(),
                  [&](uint32_t tileIdx) {
        RandomGenerator rng;
#else
    RandomGenerator rng;
    for (uint32_t tileIdx : tileIndices) {
#endif
        const uint32_t tileX = tileIdx % tilesX;
        const uint32_t tileY = tileIdx / tilesX;
        const uint32_t startX = tileX * 8;
        const uint32_t startY = tileY * 8;
        const uint32_t endX = std::min(startX + 8, static_cast<uint32_t>(m_width));
        const uint32_t endY = std::min(startY + 8, static_cast<uint32_t>(m_height));

        for (uint32_t py = startY; py < endY; ++py) {
            for (uint32_t px = startX; px < endX; ++px) {
                const float u = (static_cast<float>(px) + 0.5f) / static_cast<float>(m_width);
                const float v = (static_cast<float>(py) + 0.5f) / static_cast<float>(m_height);

                const Ray ray = camera.generateRay(u, v);
                const HitResult hit = tracer.intersect(ray);

                glm::vec3 color(0.f);

                if (hit.didHit) {
                    const unsigned char matId = scene.getGeometryMaterial(hit.geomID);
                    const Material* mat = scene.getMaterial(matId);
                    const glm::vec3 viewDir = glm::normalize(cameraPos - hit.origin);

                    if (mat->getType() == MaterialType::Emissive) {
                        color = dynamic_cast<const Material_Emissive*>(mat)->getEmission();
                    }

                    for (const auto& pLight : lights) {
                        if (!pLight->isAreaLight()) {
                            glm::vec3 lightDir = pLight->origin - hit.origin;
                            const float distSq = glm::dot(lightDir, lightDir);
                            const float dist = std::sqrt(distSq);
                            lightDir *= (1.0f / dist);

                            const float cosAngle = glm::dot(hit.normal, lightDir);
                            if (cosAngle <= 0.0f) continue;

                            const Ray shadowRay{hit.origin, lightDir, 0.0001f, dist};
                            if (tracer.isOccluded(shadowRay)) continue;

                            glm::vec3 radiance;
                            if (pLight->type == LightType::Point) {
                                radiance = pLight->color * (pLight->intensity / distSq);
                            } else {
                                radiance = pLight->color * pLight->intensity;
                            }

                            const glm::vec3 brdf = mat->shade(hit.origin, hit.normal, viewDir, lightDir);
                            color += brdf * radiance * cosAngle;
                        }
                        else if (pLight->type == LightType::MeshArea) {
                            MeshAreaLight* meshLight = pLight->meshAreaLight;
                            if (!meshLight) continue;

                            glm::vec3 bsdfSamplingContrib(0.0f);
                            glm::vec3 lightSamplingContrib(0.0f);

                            // Strategy 1: light sampling
                            for (int s = 0; s < areaLightSamples; ++s) {
                                const float u1 = rng.get();
                                const float u2 = rng.get();
                                const float u3 = rng.get();

                                const AreaLightSample lightSample = meshLight->sample(
                                    hit.origin,
                                    hit.normal,
                                    u1, u2, u3
                                );


                                if (lightSample.pdf <= 0.0f) continue;

                                glm::vec3 lightDir = lightSample.position - hit.origin;
                                const float distSq = glm::dot(lightDir, lightDir);
                                const float dist = std::sqrt(distSq);
                                lightDir *= (1.0f / dist);

                                const float cosTheta = glm::dot(hit.normal, lightDir);
                                if (cosTheta <= 0.0f) continue;

                                const float cosTheta_light = glm::dot(lightSample.normal, -lightDir);
                                if (cosTheta_light <= 0.0f) continue;

                                const Ray shadowRay{hit.origin, lightDir, 0.0001f, dist - 0.0001f};
                                bool wasVisible = !tracer.isOccluded(shadowRay);

                                if (meshLight->getVisibilitySampler() && !lightSample.traversalPath.empty()) {
                                    meshLight->getVisibilitySampler()->recordTraversalVisibility(
                                        hit.origin,
                                        lightSample.traversalPath,
                                        wasVisible
                                    );
                                }

                                if (!wasVisible) continue;
                                const glm::vec3 brdf = mat->shade(hit.origin, hit.normal, viewDir, lightDir);

                                // compute geometric term: G(x,y) = cos(θ_x) * cos(θ_y) / r^2
                                const float geometricTerm = (cosTheta * cosTheta_light) / distSq;

                                // base contribution (works for ALL sampling strategies)
                                glm::vec3 contribution = lightSample.radiance * brdf * geometricTerm;

                                // compute MIS weight if BSDF sampling is enabled
                                float misWeight = 1.0f;
                                if (mat->getType() == MaterialType::CookTorrence) {
                                    float lightPdfSolidAngle = lightSample.pdf * distSq / cosTheta_light;
                                    float bsdfPdf = BSDFSampler::pdfDiffuse(hit.normal, lightDir);

                                    misWeight = MISWeightCalculator::calculateWeight(
                                        lightPdfSolidAngle, bsdfPdf, MISHeuristic::Balance);
                                }

                                // apply MIS weight and divide by PDF
                                contribution *= misWeight / lightSample.pdf;

                                lightSamplingContrib += contribution;
                            }

                            // Strategy 2: BSDF sampling (if supported)
                            if (mat->getType() == MaterialType::CookTorrence) {
                                for (int s = 0; s < areaLightSamples; ++s) {
                                    const float u1 = rng.get();
                                    const float u2 = rng.get();

                                    BSDFSample bsdfSample = BSDFSampler::sampleDiffuse(hit.normal, u1, u2);
                                    if (bsdfSample.pdf <= 0.0f) continue;

                                    const Ray bsdfRay{hit.origin, bsdfSample.direction, 0.0001f, 1000.0f};
                                    const HitResult bsdfHit = tracer.intersect(bsdfRay);

                                    if (bsdfHit.didHit && meshLight->containsPoint(bsdfHit.origin)) {
                                        const float dist = glm::length(bsdfHit.origin - hit.origin);
                                        const float distSq = dist * dist;

                                        const float cosTheta = glm::dot(hit.normal, bsdfSample.direction);
                                        const float cosTheta_light = glm::dot(bsdfHit.normal, -bsdfSample.direction);

                                        if (cosTheta > 0.0f && cosTheta_light > 0.0f) {
                                            // get light PDF for this point (in area measure)
                                            const float lightPdfArea = meshLight->pdf(hit.origin, bsdfHit.origin);

                                            // pdf_solidAngle = pdf_area * cos(theta_light) / distance^2
                                            float lightPdfSolidAngle = lightPdfArea * cosTheta_light / distSq;

                                            // both PDFs now in solid angle measure
                                            const float misWeight = MISWeightCalculator::calculateWeight(
                                                bsdfSample.pdf,        // solid angle (cosine-weighted)
                                                lightPdfSolidAngle,    // solid angle
                                                MISHeuristic::Balance
                                            );

                                            const glm::vec3 brdf = mat->shade(hit.origin, hit.normal,
                                                                             viewDir, bsdfSample.direction);
                                            const glm::vec3 radiance = meshLight->getEmissionAt(bsdfHit.origin);

                                            // L = L_e * f_r * cos(theta) / pdf
                                            glm::vec3 contribution = radiance * brdf * cosTheta;

                                            if (misWeight > 0.0f && bsdfSample.pdf > 0.0f) {
                                                contribution *= misWeight / bsdfSample.pdf;
                                            }

                                            bsdfSamplingContrib += contribution;
                                        }
                                    }
                                }
                            }

                            // average contributions correctly: each strategy already provides an estimate
                            // averaged over its own samples
                            if (areaLightSamples > 0) {
                                lightSamplingContrib /= static_cast<float>(areaLightSamples);
                                if (mat->getType() == MaterialType::CookTorrence) {
                                    bsdfSamplingContrib /= static_cast<float>(areaLightSamples);
                                }
                            }

                            // combine the two strategies (they are already weighted by MIS)
                            color += lightSamplingContrib + bsdfSamplingContrib;
                        }
                    }

                    color = glm::clamp(color, 0.f, 1.f);
                }

                const auto r = static_cast<Uint8>(color.r * 255.f);
                const auto g = static_cast<Uint8>(color.g * 255.f);
                const auto b = static_cast<Uint8>(color.b * 255.f);

                const uint32_t pixelIndex = py * m_width + px;
                outPixels[pixelIndex] = (0xFFu << 24) | (r << 16) | (g << 8) | b;
            }
        }
#if defined(PARALLEL_EXECUTION)
    });
#else
    }
#endif

    auto end = std::chrono::high_resolution_clock::now();
    if (outDuration) {
        *outDuration = end - start;
    }
}

void Renderer::saveRawBuffer(const std::vector<uint32_t>& pixels,
                             const std::string& path) const
{
    std::ofstream out(path, std::ios::binary);
    uint32_t w = m_width, h = m_height;
    out.write((char*)&w, sizeof(uint32_t));
    out.write((char*)&h, sizeof(uint32_t));
    out.write((char*)pixels.data(), pixels.size() * sizeof(uint32_t));
}

bool Renderer::loadRawBuffer(std::vector<uint32_t>& pixels,
                             const std::string& path) const
{
    std::ifstream in(path, std::ios::binary);
    if (!in) return false;

    uint32_t w, h;
    in.read((char*)&w, sizeof(uint32_t));
    in.read((char*)&h, sizeof(uint32_t));

    if (static_cast<int>(w) != m_width || static_cast<int>(h) != m_height) return false;

    pixels.resize(w * h);
    in.read((char*)pixels.data(), pixels.size() * sizeof(uint32_t));
    return true;
}

void Renderer::generateGroundTruth(const Camera& camera, SceneManager& scene,
                                   SamplingStrategy samplerStrategy, int samples,
                                   const std::string& outDir) {
    try {
        // Create outDir if necessary
        std::filesystem::create_directories(outDir);
    } catch (const std::exception& e) {
        std::cerr << "Failed to create ground truth directory: " << e.what() << std::endl;
        return;
    }

    // Update m_groundTruthDir to be scene-specific
    m_groundTruthDir = "tests/" + m_sceneName + "/ground_truth";

    std::cout << "Generating ground truth (" << samples << " samples) using strategy: ";
    switch (samplerStrategy) {
        case SamplingStrategy::Uniform: std::cout << "Uniform\n"; break;
        case SamplingStrategy::AreaImportance: std::cout << "AreaImportance\n"; break;
        case SamplingStrategy::HierarchicalFlux: std::cout << "HierarchicalFlux\n"; break;
        case SamplingStrategy::VisibilityAwareHierarchical: std::cout << "VisibilityAware\n"; break;
        default: std::cout << "Unknown\n"; break;
    }
    // save current strategy to restore later
    SamplingStrategy oldStrategy = m_currentStrategy;

    // force strategy on all area lights for GT
    setAreaLightStrategy(scene, samplerStrategy);

    // render GT into local buffer
    std::chrono::duration<double> gtDuration;
    renderSceneIntoBuffer(camera, scene, m_groundTruthPixels, samples, &gtDuration);

    // save raw binary buffer and BMP
    const std::string rawPath = outDir + "/" + m_groundTruthRaw;
    const std::string bmpPath = outDir + "/" + m_groundTruthImage;

    saveRawBuffer(m_groundTruthPixels, rawPath);

    // temporarily set m_testFolder so saveScreenshot writes to the outDir without changing semantics
    const std::string oldTestFolder = m_testFolder;
    m_testFolder = outDir;
    saveScreenshot(m_groundTruthPixels, m_groundTruthImage);
    m_testFolder = oldTestFolder;

    // mark ground truth as computed in memory
    m_groundTruthComputed = true;

    std::cout << "Ground truth saved to: " << outDir << " (render time " << gtDuration.count() << " s)\n";

    // restore previous sampling strategy
    setAreaLightStrategy(scene, oldStrategy);
}



