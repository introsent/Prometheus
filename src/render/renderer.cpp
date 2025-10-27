//
// Created by minaj on 10/27/2025.
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
    const auto amountOfPixels = static_cast<uint32_t>(m_width * m_height);
    constexpr uint32_t packetSize = 16;

    std::vector<uint32_t> packetStarts;
    packetStarts.reserve((amountOfPixels + packetSize - 1) / packetSize);
    for (uint32_t start = 0; start < amountOfPixels; start += packetSize)
        packetStarts.push_back(start);

#if defined(PARALLEL_EXECUTION)
    std::for_each(std::execution::par, packetStarts.begin(), packetStarts.end(),
                  [&](uint32_t start) {
#else
    for (uint32_t start : packetStarts) {
#endif
        std::vector<Ray> rays;
        rays.reserve(packetSize);
        std::vector<uint32_t> pixelIndices;
        pixelIndices.reserve(packetSize);

        for (uint32_t i = 0; i < packetSize; ++i) {
            uint32_t pixelIndex = start + i;
            if (pixelIndex >= amountOfPixels) break;
            const uint32_t px = pixelIndex % m_width;
            const uint32_t py = pixelIndex / m_width;

            const float u = static_cast<float>(px + 0.5f) / static_cast<float>(m_width);
            const float v = static_cast<float>(py + 0.5f) / static_cast<float>(m_height);

            Ray ray = camera.generateRay(u, v);
            rays.push_back(ray);
            pixelIndices.push_back(pixelIndex);
        }

        std::vector<HitResult> hits;
        tracer.intersectPacket16(rays, hits);

        for (size_t lane = 0; lane < hits.size(); ++lane) {
            const HitResult& hit = hits[lane];
            glm::vec3 color = {0.f, 0.f, 0.f};

            if (hit.didHit) {
                const unsigned char matId = scene.getGeometryMaterial(hit.geomID);
                const Material* mat = scene.getMaterial(matId);

                glm::vec3 viewDir = camera.getPosition() - hit.origin;
                viewDir = glm::normalize(viewDir);

                const auto& lights = scene.getLights();
                for (const auto& pLight : lights) {
                    glm::vec3 lightDirection = pLight->origin - hit.origin;
                    float distHitToLight = glm::length(lightDirection);
                    lightDirection = glm::normalize(lightDirection);

                    // Shadow check
                    const Ray shadowRay{hit.origin, lightDirection, 0.0001f, distHitToLight};
                    if (tracer.isOccluded(shadowRay)) {
                        continue;
                    }

                    // Calculate observed area (Lambert's cosine law)
                    float cosAngle = std::max(glm::dot(hit.normal, lightDirection), 0.f);
                    if (cosAngle <= 0.f) continue;

                    // Calculate radiance
                    glm::vec3 radiance;
                    if (pLight->type == LightType::Point) {
                        float distSq = distHitToLight * distHitToLight;
                        radiance = pLight->color * (pLight->intensity / distSq);
                    } else { // Directional
                        radiance = pLight->color * pLight->intensity;
                    }

                    // BRDF shading
                    glm::vec3 brdf = mat->shade(hit.origin, hit.normal, viewDir, lightDirection);

                    // Combined: BRDF * Radiance * cos(angle)
                    color += brdf * radiance * cosAngle;
                }

                // Clamp to valid range
                color = glm::clamp(color, 0.f, 1.f);
            }

            const auto r = static_cast<Uint8>(color.r * 255.f);
            const auto g = static_cast<Uint8>(color.g * 255.f);
            const auto b = static_cast<Uint8>(color.b * 255.f);

            uint32_t pixelIndex = pixelIndices[lane];
            m_pixels[pixelIndex] = (0xFFu << 24) | (r << 16) | (g << 8) | b;
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