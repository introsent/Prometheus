//
// Created by minaj on 10/27/2025.
//

#include "renderer.h"
#include "area_light.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <random>
#include <chrono>
#include <filesystem>

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

Renderer::Renderer(const int width, const int height)
    : m_width(width), m_height(height), m_window(nullptr),
      m_renderer(nullptr), m_texture(nullptr), m_quit(false),
      m_testMode(false), m_currentSamples(1), m_maxSamples(16384) {
    m_pixels.resize(width * height);

    // Generate all sample counts (powers of two from 1 to 16384)
    int samples = 1;
    while (samples <= m_maxSamples) {
        m_sampleCounts.push_back(samples);
        samples <<= 1;  // Multiply by 2
    }
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
        createTestFolder();  // Create folder when test starts

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
    // prepend folder path to filename
    const std::string filepath = m_testFolder + "/" + filename;

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


bool Renderer::initialize() {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
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
    // start timing
    auto renderStart = std::chrono::high_resolution_clock::now();

    const RayTracer tracer(&scene);
    const glm::vec3 cameraPos = camera.getPosition();
    const auto& lights = scene.getLights();

    // determine number of samples to use
    int areaLightSamples = m_testMode ? m_currentSamples : 1028;

    // display current sample count
    if (m_testMode) {
        std::cout << "Rendering with " << areaLightSamples << " samples..." << std::endl;
    }

    // tile-based rendering for parallelization
    constexpr uint32_t tileSize = 64;
    const uint32_t tilesX = (m_width + 7) / 8;
    const uint32_t tilesY = (m_height + 7) / 8;
    const uint32_t totalTiles = tilesX * tilesY;

    std::vector<uint32_t> tileIndices(totalTiles);
    for (uint32_t i = 0; i < totalTiles; ++i) {
        tileIndices[i] = i;
    }

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

        // Process tile
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

                    // Direct lighting from all light sources
                    for (const auto& pLight : lights) {
                        if (!pLight->isAreaLight()) {
                            // Point or directional light
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
                        else if (pLight->type == LightType::TriangleArea) {
                            // Triangle Area Light (Monte Carlo)
                            TriangleAreaLight* areaLight = pLight->triangleAreaLight;
                            if (!areaLight) continue;

                            glm::vec3 areaLightContrib(0.0f);

                            for (int s = 0; s < areaLightSamples; ++s) {
                                const float u1 = rng.get();
                                const float u2 = rng.get();
                                const AreaLightSample lightSample = areaLight->sample(hit.origin, u1, u2);

                                glm::vec3 lightDir = lightSample.position - hit.origin;
                                const float distSq = glm::dot(lightDir, lightDir);
                                const float dist = std::sqrt(distSq);
                                lightDir *= (1.0f / dist);

                                const float cosTheta = glm::dot(hit.normal, lightDir);
                                if (cosTheta <= 0.0f) continue;

                                const float cosTheta_light = glm::dot(lightSample.normal, -lightDir);
                                if (cosTheta_light <= 0.0f) continue;

                                const Ray shadowRay{hit.origin, lightDir, 0.0001f, dist - 0.0001f};
                                if (tracer.isOccluded(shadowRay)) continue;

                                const glm::vec3 brdf = mat->shade(hit.origin, hit.normal, viewDir, lightDir);
                                const float geometricTerm = cosTheta_light / distSq;
                                const glm::vec3 contribution =
                                    lightSample.radiance * brdf * cosTheta * geometricTerm / lightSample.pdf;

                                areaLightContrib += contribution;
                            }

                            color += areaLightContrib / static_cast<float>(areaLightSamples);
                        }
                        else if (pLight->type == LightType::MeshArea) {
                            // Mesh Area Light (Monte Carlo)
                            MeshAreaLight* meshLight = pLight->meshAreaLight;
                            if (!meshLight) continue;

                            glm::vec3 areaLightContrib(0.0f);

                            for (int s = 0; s < areaLightSamples; ++s) {
                                const float u1 = rng.get();
                                const float u2 = rng.get();
                                const float u3 = rng.get();

                                const AreaLightSample lightSample = meshLight->sample(hit.origin, u1, u2, u3);

                                glm::vec3 lightDir = lightSample.position - hit.origin;
                                const float distSq = glm::dot(lightDir, lightDir);
                                const float dist = std::sqrt(distSq);
                                lightDir *= (1.0f / dist);

                                const float cosTheta = glm::dot(hit.normal, lightDir);
                                if (cosTheta <= 0.0f) continue;

                                const float cosTheta_light = glm::dot(lightSample.normal, -lightDir);
                                if (cosTheta_light <= 0.0f) continue;

                                const Ray shadowRay{hit.origin, lightDir, 0.0001f, dist - 0.0001f};
                                if (tracer.isOccluded(shadowRay)) continue;

                                const glm::vec3 brdf = mat->shade(hit.origin, hit.normal, viewDir, lightDir);
                                const float geometricTerm = cosTheta_light / distSq;
                                const glm::vec3 contribution =
                                    lightSample.radiance * brdf * cosTheta * geometricTerm / lightSample.pdf;

                                areaLightContrib += contribution;
                            }

                            color += areaLightContrib / static_cast<float>(areaLightSamples);
                        }
                    }

                    // Clamp color
                    color = glm::clamp(color, 0.f, 1.f);
                }

                // Convert to pixel color
                const auto r = static_cast<Uint8>(color.r * 255.f);
                const auto g = static_cast<Uint8>(color.g * 255.f);
                const auto b = static_cast<Uint8>(color.b * 255.f);

                const uint32_t pixelIndex = py * m_width + px;
                m_pixels[pixelIndex] = (0xFFu << 24) | (r << 16) | (g << 8) | b;
            }
        }
#if defined(PARALLEL_EXECUTION)
    });
#else
    }
#endif

    // end timing
    auto renderEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> renderDuration = renderEnd - renderStart;

    // update texture for display
    SDL_UpdateTexture(m_texture, nullptr, m_pixels.data(), m_width * sizeof(Uint32));

    // handle test mode operations
    if (m_testMode) {
        // save screenshot
        std::ostringstream filename;
        filename << "samples_" << std::setw(5) << std::setfill('0')
                 << areaLightSamples << ".bmp";
        saveScreenshot(filename.str());

        // save timing data - use folder path
        std::string timingFilePath = getTestFolderPath() + "/timing_results.txt";
        std::ofstream timingFile(timingFilePath, std::ios::app);
        if (timingFile.is_open()) {
            timingFile << areaLightSamples << "\t"
                      << std::fixed << std::setprecision(3)
                      << renderDuration.count() << std::endl;
            timingFile.close();
            std::cout << "  Time: " << std::fixed << std::setprecision(3)
                     << renderDuration.count() << " seconds" << std::endl;
        }

        // Move to next sample count
        if (auto it = std::ranges::find(m_sampleCounts, m_currentSamples); it != m_sampleCounts.end() && ++it != m_sampleCounts.end()) {
            m_currentSamples = *it;
        } else {
            m_currentSamples = m_maxSamples + 1; // Mark as complete
            std::cout << "\n=== Area Light Sampling Test Complete ===" << std::endl;
            std::cout << "All screenshots saved to: " << m_testFolder << std::endl;
            std::cout << "Timing data saved to: " << timingFilePath << std::endl;

            // create a summary file
            createTestSummary();
        }
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
        // Space bar to pause/resume test
        if (event.type == SDL_EVENT_KEY_DOWN && event.key.key == SDLK_SPACE) {
            m_testMode = !m_testMode;
            std::cout << "Test mode: " << (m_testMode ? "ON" : "OFF") << std::endl;
        }
    }
    return m_quit || (m_testMode && isTestComplete());
}

void Renderer::createTestFolder() {
    // create a timestamp for the folder name
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    std::tm tm_buf;

#ifdef _WIN32
    localtime_s(&tm_buf, &in_time_t);
#else
    localtime_r(&in_time_t, &tm_buf);
#endif

    std::ostringstream folderName;
    folderName << "area_light_test_";
    folderName << std::put_time(&tm_buf, "%Y%m%d_%H%M%S");

    m_testFolder = folderName.str();

    // create the directory
    try {
        if (!std::filesystem::exists(m_testFolder)) {
            std::filesystem::create_directory(m_testFolder);
            std::cout << "Created test folder: " << m_testFolder << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Failed to create test folder: " << e.what() << std::endl;
        m_testFolder = ".";  // fall back to current directory
    }
}

std::string Renderer::getTestFolderPath() const {
    return m_testFolder.empty() ? "." : m_testFolder;
}

void Renderer::createTestSummary() const {
    const std::string summaryPath = getTestFolderPath() + "/test_summary.txt";

    if (std::ofstream summary(summaryPath); summary.is_open()) {
        summary << "=== Area Light Sampling Test Summary ===" << std::endl;
        summary << "Test conducted: " << m_testFolder << std::endl;
        summary << "Resolution: " << m_width << "x" << m_height << std::endl;
        summary << "Total tests: " << m_sampleCounts.size() << std::endl;
        summary << "Sample counts tested: ";
        for (size_t i = 0; i < m_sampleCounts.size(); ++i) {
            summary << m_sampleCounts[i];
            if (i < m_sampleCounts.size() - 1) summary << ", ";
        }
        summary << std::endl << std::endl;

        summary << "To analyze the results:" << std::endl;
        summary << "1. Check timing_results.txt for render times" << std::endl;
        summary << "2. Compare image quality between sample counts" << std::endl;
        summary << "3. Look for diminishing returns after certain sample count" << std::endl;
        summary << "4. The optimal sample count balances quality vs performance" << std::endl;

        summary.close();
        std::cout << "Test summary saved to: " << summaryPath << std::endl;
    }
}
