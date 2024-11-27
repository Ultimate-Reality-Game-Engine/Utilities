

#include <GameTimer.h>

namespace UltReality::Utilities
{
	GameTimer* GameTimer::GetInstance() noexcept
	{
		// Acquire a concurrent read lock (shared) for read access
		std::shared_lock<std::shared_mutex> readLock(s_sharedMutex);

		if (!GameTimer::s_instance)
		{
			// Release the read lock to prevent a deadlock as both read and write use the same mutex
			readLock.unlock();

			// Acquire an exclusive lock for write access
			std::unique_lock<std::shared_mutex> writeLock(s_sharedMutex);

			if (!GameTimer::s_instance)
			{
				// Create the instance
				GameTimer::s_instance = new GameTimer();
			}
		}

		// Return the instance pointer
		return GameTimer::s_instance;
	}

	void GameTimer::DestroyInstance() noexcept
	{
		// Acquire exclusive lock
		std::unique_lock<std::shared_mutex> writeLock(s_sharedMutex);

		if (s_instance)
		{
			delete s_instance;
			s_instance = nullptr;
		}
	}

	double GameTimer::GetTotalTime() const noexcept
	{
		if (!m_stopped)
		{
			return std::chrono::duration<double>(m_currTime - m_baseTime).count();
		}
		else
		{
			return std::chrono::duration<double>(m_stopTime - m_baseTime).count();
		}
	}

	double GameTimer::GetGameTime() const noexcept
	{
		using namespace std::chrono;

		duration<double> currDuration = duration_cast<duration<double>>(m_currTime.time_since_epoch());
		duration<double> pausedDuration = duration_cast<duration<double>>(m_pausedTime.time_since_epoch());
		duration<double> baseDuration = duration_cast<duration<double>>(m_baseTime.time_since_epoch());

		if (!m_stopped)
		{
			return (currDuration - pausedDuration - baseDuration).count();
		}
		else
		{
			duration<double> stopDuration = duration_cast<duration<double>>(m_stopTime.time_since_epoch());
			return (stopDuration - pausedDuration - baseDuration).count();
		}
	}

	double GameTimer::GetDeltaTime() const noexcept
	{
		return m_deltaTime;
	}

	void GameTimer::Initialize() noexcept
	{
		auto now = std::chrono::high_resolution_clock::now();
		m_baseTime = now;
		m_prevTime = now;
		m_stopTime = {};
		m_stopped = false;
	}

	void GameTimer::Start() noexcept
	{
		// Accumulate the time elapsed between stop and start pairs if we are not paused
		if (m_stopped)
		{
			auto startTime = std::chrono::high_resolution_clock::now();
			m_pausedTime += startTime - m_stopTime;
			m_prevTime = startTime;
			m_stopTime = {};
			m_stopped = false;
		}
	}

	void GameTimer::Stop() noexcept
	{
		// Pause the timer
		if (!m_stopped)
		{
			m_stopTime = std::chrono::high_resolution_clock::now();
			m_stopped = true;
		}
	}

	void GameTimer::Tick() noexcept
	{
		if (m_stopped)
		{
			m_deltaTime = 0.0;
			return;
		}

		m_currTime = std::chrono::high_resolution_clock::now();
		m_deltaTime = std::chrono::duration<double>(m_currTime - m_prevTime).count();
		m_prevTime = m_currTime;
	}
}