//
// Optimized renderer with better parallelization
//

#include "renderer.h"
#include "area_light.h"
#include <algorithm>
#include <iostream>
#include <random>

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
      m_renderer(nullptr), m_texture(nullptr), m_quit(false) {
    m_pixels.resize(width * height);
}

Renderer::~Renderer() {
    if (m_texture) SDL_DestroyTexture(m_texture);
    if (m_renderer) SDL_DestroyRenderer(m_renderer);
    if (m_window) SDL_DestroyWindow(m_window);
    SDL_Quit();
}

bool Renderer::initialize() {
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        std::cerr << "SDL initialization failed: " << SDL_GetError() << std::endl;
        return false;
    }

    m_window = SDL_CreateWindow("Prometheus", m_width, m_height, 0);
    if (!m_window) {
        std::cerr << "Window creation failed: " << SDL_GetError() << std::endl;
        return false;
    }

    m_renderer = SDL_CreateRenderer(m_window, nullptr);
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
    const RayTracer tracer(&scene);

    // pre-cache camera position and lights
    const glm::vec3 cameraPos = camera.getPosition();
    const auto& lights = scene.getLights();

    // area light settings
    constexpr int AREA_LIGHT_SAMPLES = 8; // number of samples per area light

    // use larger tiles for better cache coherency
    constexpr uint32_t tileSize = 64; // 8x8 tiles
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
        // thread-local random generator
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

        // process tile
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

                    // direct lighting from all light sources
                    for (const auto& pLight : lights) {
                        if (!pLight->isAreaLight()) {
                            // for point or directional light
                            glm::vec3 lightDir = pLight->origin - hit.origin;
                            const float distSq = glm::dot(lightDir, lightDir);
                            const float dist = std::sqrt(distSq);
                            lightDir *= (1.0f / dist);

                            const float cosAngle = glm::dot(hit.normal, lightDir);
                            if (cosAngle <= 0.0f) continue;

                            const Ray shadowRay{hit.origin, lightDir, 0.0001f, dist};
                            if (tracer.isOccluded(shadowRay)) {
                                continue;
                            }

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
                            // for Triangle Area Light (Monte Carlo)
                            TriangleAreaLight* areaLight = pLight->triangleAreaLight;
                            if (!areaLight) continue;

                            glm::vec3 areaLightContrib(0.0f);

                            // Monte Carlo integration with multiple samples
                            for (int s = 0; s < AREA_LIGHT_SAMPLES; ++s) {
                                // generate random samples
                                const float u1 = rng.get();
                                const float u2 = rng.get();

                                // sample point on area light
                                const AreaLightSample lightSample = areaLight->sample(hit.origin, u1, u2);

                                // direction from hit point to light sample
                                glm::vec3 lightDir = lightSample.position - hit.origin;
                                const float distSq = glm::dot(lightDir, lightDir);
                                const float dist = std::sqrt(distSq);
                                lightDir *= (1.0f / dist);

                                // check if light is above surface
                                const float cosTheta = glm::dot(hit.normal, lightDir);
                                if (cosTheta <= 0.0f) continue;

                                // check if we are hitting the back of the light
                                const float cosTheta_light = glm::dot(lightSample.normal, -lightDir);
                                if (cosTheta_light <= 0.0f) continue;

                                // visibility test
                                const Ray shadowRay{hit.origin, lightDir, 0.0001f, dist - 0.0001f};
                                if (tracer.isOccluded(shadowRay)) {
                                    continue;
                                }

                                // compute rendering equation components
                                // note: L = Le * BRDF * cosTheta * cosTheta_light / (distance^2 * pdf)

                                const glm::vec3 brdf = mat->shade(hit.origin, hit.normal, viewDir, lightDir);

                                // geometric term: cosTheta_light / distance^2
                                const float geometricTerm = cosTheta_light / distSq;

                                // Monte Carlo estimator
                                const glm::vec3 contribution =
                                    lightSample.radiance * brdf * cosTheta * geometricTerm / lightSample.pdf;

                                areaLightContrib += contribution;
                            }

                            // average over all samples
                            color += areaLightContrib / static_cast<float>(AREA_LIGHT_SAMPLES);
                        }
                        else if (pLight->type == LightType::MeshArea) {
                            // for Mesh Area Light (Monte Carlo)
                            MeshAreaLight* meshLight = pLight->meshAreaLight;
                            if (!meshLight) continue;

                            glm::vec3 areaLightContrib(0.0f);

                            for (int s = 0; s < AREA_LIGHT_SAMPLES; ++s) {
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
                                if (tracer.isOccluded(shadowRay)) {
                                    continue;
                                }

                                const glm::vec3 brdf = mat->shade(hit.origin, hit.normal, viewDir, lightDir);
                                const float geometricTerm = cosTheta_light / distSq;
                                const glm::vec3 contribution =
                                    lightSample.radiance * brdf * cosTheta * geometricTerm / lightSample.pdf;

                                areaLightContrib += contribution;
                            }

                            color += areaLightContrib / static_cast<float>(AREA_LIGHT_SAMPLES);
                        }
                    }

                    // clamp after all lights
                    color = glm::clamp(color, 0.f, 1.f);
                }

                // convert to pixel color
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

    SDL_UpdateTexture(m_texture, nullptr, m_pixels.data(), m_width * sizeof(Uint32));
}

void Renderer::present() {
    SDL_RenderClear(m_renderer);
    SDL_RenderTexture(m_renderer, m_texture, nullptr, nullptr);
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
        }
    }
    return m_quit;
}

SDL_Color Renderer::toSDLColor(const glm::vec3& color) {
    return SDL_Color{
        static_cast<Uint8>(std::clamp(color.r * 255.f, 0.f, 255.f)),
        static_cast<Uint8>(std::clamp(color.g * 255.f, 0.f, 255.f)),
        static_cast<Uint8>(std::clamp(color.b * 255.f, 0.f, 255.f)),
        255
    };
}