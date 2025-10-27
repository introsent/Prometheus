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
    const RayTracer tracer(scene.getScene());
    const auto amountOfPixels = static_cast<uint32_t>(m_width * m_height);

#if defined(PARALLEL_EXECUTION)
    // Parallel execution path
    std::vector<uint32_t> pixelIndices;
    pixelIndices.reserve(amountOfPixels);

    for (uint32_t index = 0; index < amountOfPixels; ++index) {
        pixelIndices.emplace_back(index);
    }

    std::for_each(std::execution::par, pixelIndices.begin(), pixelIndices.end(),
                  [&](uint32_t i) {
        renderPixel(camera, scene, tracer, i);
    });
#else
    // Sequential execution path
    for (uint32_t pixelIndex = 0; pixelIndex < amountOfPixels; ++pixelIndex) {
        renderPixel(camera, scene, tracer, pixelIndex);
    }
#endif

    SDL_UpdateTexture(m_texture, nullptr, m_pixels.data(), m_width * sizeof(Uint32));
}

void Renderer::renderPixel(const Camera& camera, const SceneManager& scene,
                           const RayTracer& tracer, uint32_t pixelIndex) {
    const uint32_t px = pixelIndex % m_width;
    const uint32_t py = pixelIndex / m_width;

    const float u = static_cast<float>(px + 0.5f) / static_cast<float>(m_width);
    const float v = static_cast<float>(py + 0.5f) / static_cast<float>(m_height);

    Ray ray = camera.generateRay(u, v);
    const glm::vec3 color = traceRay(ray, tracer, scene);

    // Convert to ARGB8888
    const Uint8 r = static_cast<Uint8>(std::clamp(color.r * 255.f, 0.f, 255.f));
    const Uint8 g = static_cast<Uint8>(std::clamp(color.g * 255.f, 0.f, 255.f));
    const Uint8 b = static_cast<Uint8>(std::clamp(color.b * 255.f, 0.f, 255.f));

    m_pixels[py * m_width + px] = (0xFF << 24) | (r << 16) | (g << 8) | b;
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

glm::vec3 Renderer::traceRay(const Ray& ray, const RayTracer& tracer, const SceneManager& scene) {
    if (const HitResult hit = tracer.intersect(ray); hit.didHit) {

        const unsigned char matId = scene.getGeometryMaterial(hit.geomID);
        const Material* mat = scene.getMaterial(matId);
        glm::vec3 matColor = mat->getColor();

        const auto& lights = scene.getLights();


        for (const auto& pLight : lights) {
            // Configure light
            glm::vec3 lightDirection = pLight->origin - hit.origin;
            float distHitToLight =  glm::length(lightDirection);
            lightDirection = normalize(lightDirection);

            const Ray lightRay {hit.origin, lightDirection, 0.0001f,  distHitToLight};
            if (const HitResult lightHit = tracer.intersect(lightRay); lightHit.didHit) {
                return matColor * 0.5f;
            }
        }
        return matColor;
    }

    // Background color (black)
    return {0.f, 0.f, 0.f};
}

SDL_Color Renderer::toSDLColor(const glm::vec3& color) {
    return SDL_Color{
        static_cast<Uint8>(std::clamp(color.r * 255.f, 0.f, 255.f)),
        static_cast<Uint8>(std::clamp(color.g * 255.f, 0.f, 255.f)),
        static_cast<Uint8>(std::clamp(color.b * 255.f, 0.f, 255.f)),
        255
    };
}