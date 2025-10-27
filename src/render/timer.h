//
// Created by minaj on 10/27/2025.
//

#ifndef TIMER_H
#define TIMER_H
#include <cstdint>
#include <vector>

class Timer final {
public:
    Timer();
    virtual ~Timer() = default;

    Timer(const Timer&) = delete;
    Timer(Timer&&) noexcept = delete;
    Timer& operator=(const Timer&) = delete;
    Timer& operator=(Timer&&) noexcept = delete;

    void startBenchmark(int numFrames = 10);

    void reset();
    void start();
    void update();
    void stop();

    [[nodiscard]] uint32_t getFPS() const { return m_FPS; };
    [[nodiscard]] float getdFPS() const { return m_dFPS; };
    [[nodiscard]] float getElapsed() const { return m_elapsedTime; };
    [[nodiscard]] float getTotal() const { return m_totalTime; };
    [[nodiscard]] bool isRunning() const { return !m_isStopped; };

private:
    uint64_t m_baseTime = 0;
    uint64_t m_pausedTime = 0;
    uint64_t m_stopTime = 0;
    uint64_t m_previousTime = 0;
    uint64_t m_currentTime = 0;

    uint32_t m_FPS = 0;
    float m_dFPS = 0.0f;
    uint32_t m_FPSCount = 0;

    float m_totalTime = 0.0f;
    float m_elapsedTime = 0.0f;
    float m_secondsPerCount = 0.0f;
    float m_elapsedUpperBound = 0.03f;
    float m_FPSTimer = 0.0f;

    bool m_isStopped = true;
    bool m_forceElapsedUpperBound = false;

    bool m_benchmarkActive = false;
    float m_benchmarkHigh{ 0.f };
    float m_benchmarkLow{ 0.f };
    float m_benchmarkAvg{ 0.f };
    int m_benchmarkFrames{ 0 };
    int m_benchmarkCurrFrame{ 0 };
    std::vector<float> m_benchmarks{};
};

#endif //TIMER_H
