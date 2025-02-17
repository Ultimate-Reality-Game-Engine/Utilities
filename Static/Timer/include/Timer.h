#ifndef ULTREALITY_UTILITIES_GAME_TIMER_H
#define ULTREALITY_UTILITIES_GAME_TIMER_H

#include <stdint.h>
#include <chrono>

namespace UltReality::Utilities
{
	class Timer
	{
	private:
		double m_deltaTime = 0.0;

		std::chrono::time_point<std::chrono::high_resolution_clock> m_baseTime;		// Start time
		std::chrono::time_point<std::chrono::high_resolution_clock> m_pausedTime;	// Amount of time timer was paused
		std::chrono::time_point<std::chrono::high_resolution_clock> m_stopTime;		// Time at which the timer was stopped/paused
		std::chrono::time_point<std::chrono::high_resolution_clock> m_prevTime;		// Time at last tick
		std::chrono::time_point<std::chrono::high_resolution_clock> m_currTime;		// Time at this tick

		bool m_stopped = false;

	public:
		/// <summary>
		/// Create Timer instance
		/// </summary>
		Timer() noexcept = default;

		/// <summary>
		/// Destroy the Timer instance
		/// </summary>
		~Timer() noexcept = default;

		/// <summary>
		/// Gets the total time this Timer instance has been running
		/// </summary>
		/// <returns>Total running time</returns>
		double GetTotalTime() const noexcept;

		/// <summary>
		/// Gets the total time this Timer instance has been running and active (minus paused time)
		/// </summary>
		/// <returns>Active time since initialization</returns>
		double GetGameTime() const noexcept;

		/// <summary>
		/// Gets the amount of time that has elapsed since the last time step
		/// </summary>
		/// <returns>Elapsed time</returns>
		double GetDeltaTime() const noexcept;

		/// <summary>
		/// Call this method to set the timer's initial state and prepare it for update
		/// </summary>
		void Initialize() noexcept; // Call before message loop

		/// <summary>
		/// Start the timer up (un-pause)
		/// </summary>
		void Start() noexcept; // Call when unpaused

		/// <summary>
		/// Stop the timer (pause)
		/// </summary>
		void Stop() noexcept; // Call when paused

		/// <summary>
		/// Update timer state to measure elapsed time since last call
		/// </summary>
		void Tick() noexcept; // Call every frame
	};
}

#endif // !ULTREALITY_UTILITIES_GAME_TIMER_H
