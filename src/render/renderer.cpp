//
// Optimized renderer with better parallelization
//

#include "renderer.h"
#include <algorithm>
#include <iostream>

#if defined(PARALLEL_EXECUTION)
#include <execution>
#endif

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
    const auto totalPixels = static_cast<uint32_t>(m_width * m_height);

    // Pre-cache camera position and lights to avoid repeated calls
    const glm::vec3 cameraPos = camera.getPosition();
    const auto& lights = scene.getLights();

    // Use larger tiles for better cache coherency
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
#else
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

                    // Accumulate lighting from all sources
                    for (const auto& pLight : lights) {
                        glm::vec3 lightDir = pLight->origin - hit.origin;
                        const float distSq = glm::dot(lightDir, lightDir);
                        const float dist = std::sqrt(distSq);
                        lightDir *= (1.0f / dist); // Normalize

                        // Early exit if light is behind surface
                        const float cosAngle = glm::dot(hit.normal, lightDir);
                        if (cosAngle <= 0.0f) continue;

                        // Shadow test
                        const Ray shadowRay{hit.origin, lightDir, 0.0001f, dist};
                        if (tracer.isOccluded(shadowRay)) {
                            continue;
                        }

                        // Calculate radiance
                        glm::vec3 radiance;
                        if (pLight->type == LightType::Point) {
                            radiance = pLight->color * (pLight->intensity / distSq);
                        } else {
                            radiance = pLight->color * pLight->intensity;
                        }

                        // BRDF evaluation
                        const glm::vec3 brdf = mat->shade(hit.origin, hit.normal, viewDir, lightDir);

                        // Accumulate
                        color += brdf * radiance * cosAngle;
                    }

                    // Clamp after all lights
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