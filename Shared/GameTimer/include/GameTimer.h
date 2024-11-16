#ifndef ULTREALITY_UTILITIES_GAME_TIMER_H
#define ULTREALITY_UTILITIES_GAME_TIMER_H

#include <stdint.h>
#include <chrono>
#include <shared_mutex>

#include <LibraryExport.h>

namespace UltReality::Utilities
{
	class LIBRARY_ABI GameTimer
	{
	private:
		double m_deltaTime = 0.0;

		std::chrono::time_point<std::chrono::high_resolution_clock> m_baseTime;		// Start time
		std::chrono::time_point<std::chrono::high_resolution_clock> m_pausedTime;	// Amount of time timer was paused
		std::chrono::time_point<std::chrono::high_resolution_clock> m_stopTime;		// Time at which the timer was stopped/paused
		std::chrono::time_point<std::chrono::high_resolution_clock> m_prevTime;		// Time at last tick
		std::chrono::time_point<std::chrono::high_resolution_clock> m_currTime;		// Time at this tick

		bool m_stopped = false;

		inline static std::shared_mutex s_sharedMutex;
		inline static GameTimer* s_instance = nullptr;

		/// <summary>
		/// Used internally by GetInstance
		/// </summary>
		GameTimer() noexcept;

	public:
		/// <summary>
		/// Deletes the singleton, frees memory
		/// </summary>
		LIBRARY_CALL ~GameTimer() noexcept;

		/// <summary>
		/// Gets a pointer to the singleton's instance. Creates the instance if it doesn't exist
		/// </summary>
		/// <returns>Pointer to instance</returns>
		static GameTimer* LIBRARY_CALL GetInstance() noexcept;

		/// <summary>
		/// Gets the total time the application has been running (including paused time)
		/// </summary>
		/// <returns>Total running time</returns>
		double LIBRARY_CALL GetTotalTime() const noexcept;

		/// <summary>
		/// Gets the total time the application has been running and active (minus paused time)
		/// </summary>
		/// <returns>Active time since initialization</returns>
		double LIBRARY_CALL GetGameTime() const noexcept;

		/// <summary>
		/// Gets the amount of time that has elapsed since the last time step
		/// </summary>
		/// <returns>Elapsed time</returns>
		double LIBRARY_CALL GetDeltaTime() const noexcept;

		/// <summary>
		/// Call this method to set the timer's initial state and prepare it for update
		/// </summary>
		void LIBRARY_CALL Initialize() noexcept; // Call before message loop

		/// <summary>
		/// Start the timer up/un-pause
		/// </summary>
		void LIBRARY_CALL Start() noexcept; // Call when unpaused

		/// <summary>
		/// Stop the timer/pause
		/// </summary>
		void LIBRARY_CALL Stop() noexcept; // Call when paused

		/// <summary>
		/// Update timer state to measure elapsed time since last call
		/// </summary>
		void LIBRARY_CALL Tick() noexcept; // Call every frame

		GameTimer(const GameTimer&) = delete;
		GameTimer& operator=(const GameTimer&) = delete;
		GameTimer(GameTimer&&) = delete;
		GameTimer& operator=(GameTimer&&) = delete;
	};
}

#endif // !ULTREALITY_UTILITIES_GAME_TIMER_H
