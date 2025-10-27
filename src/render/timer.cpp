//
// Created by minaj on 10/27/2025.
//

#include "timer.h"

#include <iostream>
#include <numeric>
#include <fstream>

#include "SDL3/SDL.h"

Timer::Timer()
{
	const uint64_t countsPerSecond = SDL_GetPerformanceFrequency();
	m_secondsPerCount = 1.0f / static_cast<float>(countsPerSecond);
}

void Timer::reset()
{
	const uint64_t currentTime = SDL_GetPerformanceCounter();

	m_baseTime = currentTime;
	m_previousTime = currentTime;
	m_stopTime = 0;
	m_FPSTimer = 0.0f;
	m_FPSCount = 0;
	m_isStopped = false;
}

void Timer::start()
{
	const uint64_t startTime = SDL_GetPerformanceCounter();

	if (m_isStopped)
	{
		m_pausedTime += (startTime - m_stopTime);

		m_previousTime = startTime;
		m_stopTime = 0;
		m_isStopped = false;
	}
}

void Timer::startBenchmark(int numFrames)
{
	if (m_benchmarkActive)
	{
		std::cout << "(Benchmark already running)";
		return;
	}

	m_benchmarkActive = true;

	m_benchmarkAvg = 0.f;
	m_benchmarkHigh = FLT_MIN;
	m_benchmarkLow = FLT_MAX;

	m_benchmarkFrames = numFrames;
	m_benchmarkCurrFrame = 0;

	m_benchmarks.clear();
	m_benchmarks.resize(m_benchmarkFrames);

	std::cout << "**BENCHMARK STARTED**\n";
}

void Timer::update()
{
	if (m_isStopped)
	{
		m_FPS = 0;
		m_elapsedTime = 0.0f;
		m_totalTime = static_cast<float>(((m_stopTime - m_pausedTime) - m_baseTime) * m_baseTime);
		return;
	}

	const uint64_t currentTime = SDL_GetPerformanceCounter();
	m_currentTime = currentTime;

	m_elapsedTime = static_cast<float>((m_currentTime - m_previousTime) * m_secondsPerCount);
	m_previousTime = m_currentTime;

	if (m_elapsedTime < 0.0f)
		m_elapsedTime = 0.0f;

	if (m_forceElapsedUpperBound && m_elapsedTime > m_elapsedUpperBound)
	{
		m_elapsedTime = m_elapsedUpperBound;
	}

	m_totalTime = static_cast<float>(((m_currentTime - m_pausedTime) - m_baseTime) * m_secondsPerCount);

	//FPS LOGIC
	m_FPSTimer += m_elapsedTime;
	++m_FPSCount;
	if (m_FPSTimer >= 1.0f)
	{
		m_dFPS = m_FPSCount / m_FPSTimer;
		m_FPS = m_FPSCount;
		m_FPSCount = 0;
		m_FPSTimer = 0.0f;

		if (m_benchmarkActive)
		{
			m_benchmarks[m_benchmarkCurrFrame] = m_dFPS;

			m_benchmarkLow = std::min(m_benchmarkLow, m_dFPS);
			m_benchmarkHigh = std::max(m_benchmarkHigh, m_dFPS);

			++m_benchmarkCurrFrame;
			if (m_benchmarkCurrFrame >= m_benchmarkFrames)
			{
				m_benchmarkActive = false;
				m_benchmarkAvg = std::accumulate(m_benchmarks.begin(), m_benchmarks.end(), 0.f) / static_cast<float>(m_benchmarkFrames);

				//print
				std::cout << "**BENCHMARK FINISHED**\n";
				std::cout << ">> HIGH = " << m_benchmarkHigh << std::endl;
				std::cout << ">> LOW = " << m_benchmarkLow << std::endl;
				std::cout << ">> AVG = " << m_benchmarkAvg << std::endl;

				//file save
				std::ofstream fileStream("benchmark.txt");
				fileStream << "FRAMES = " << m_benchmarkCurrFrame << std::endl;
				fileStream << "HIGH = " << m_benchmarkHigh << std::endl;
				fileStream << "LOW = " << m_benchmarkLow << std::endl;
				fileStream << "AVG = " << m_benchmarkAvg << std::endl;
				fileStream.close();
			}
		}
	}
}

void Timer::stop()
{
	if (!m_isStopped)
	{
		const uint64_t currentTime = SDL_GetPerformanceCounter();

		m_stopTime = currentTime;
		m_isStopped = true;
	}
}

