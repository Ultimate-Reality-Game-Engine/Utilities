#ifndef ULTREALITY_MATH_RANDOM_H
#define ULTREALITY_MATH_RANDOM_H

#include <stdint.h>

namespace UltReality::Math
{
	/// <summary>
	/// This is the core class that encapsulates a Pseudo Random Number Generator state and produces Pseudo Random Numbers using XORshift
	/// </summary>
	class RandomEngine
	{
	private:
		uint32_t m_state;

	public:
		/// <summary>
		/// Constructor that seeds with a default value
		/// </summary>
		/// <param name="seed">Value to seed engine with, default of 2463534242U</param>
		explicit RandomEngine(uint32_t seed = 2463534242U);

		/// <summary>
		/// Seed the engine with a manually provided value
		/// </summary>
		/// <param name="seed">Value to seed engine with</param>
		void Seed(uint32_t seed);

		/// <summary>
		/// Seed the engine with the current system time
		/// </summary>
		void SeedWithTime();

		/// <summary>
		/// XORshift32 algorithm
		/// </summary>
		/// <returns>Pseudo Random 32 bit unsinged integer</returns>
		uint32_t xorshift32();

		/// <summary>
		/// Return a random unsinged integer using the internal XORshift32
		/// </summary>
		uint32_t UInt();

		/// <summary>
		/// Return a random unsinged integer in the range of [min, max] using the internal XORshift
		/// </summary>
		/// <param name="min">Lower end of range, inclusive</param>
		/// <param name="max">Upper end of range, inclusive</param>
		uint32_t UInt(uint32_t min, uint32_t max);

		/// <summary>
		/// Return a random integer using the internal XORshift
		/// </summary>
		int32_t Int();

		/// <summary>
		/// Return a random integer in the range of [min, max] using the internal XORshift
		/// </summary>
		/// <param name="min">Lower end of range, inclusive</param>
		/// <param name="max">Upper end of range, inclusive</param>
		int32_t Int(int32_t min, int32_t max);

		/// <summary>
		/// Return a random float in the range of [0, 1) using the internal XORshift
		/// </summary>
		float Float();

		/// <summary>
		/// Return a random float in the range of [min, max) using the internal XORshift
		/// </summary>
		/// <param name="min">Lower end of range, inclusive</param>
		/// <param name="max">Upper end of range, exclusive</param>
		float Float(float min, float max);

		/// <summary>
		/// Return a random double in the range of [0, 1) using the internal XORshift
		/// </summary>
		double Double();

		/// <summary>
		/// Return a random double in the range of [min, max) using the internal XORshift
		/// </summary>
		/// <param name="min">Lower end of range, inclusive</param>
		/// <param name="max">Upper end of range, exclusive</param>
		double Double(double min, double max);
	};

	/// <summary>
	/// A static class which internally contains a thread_local instances of a <see cref="RandomEngine"/> which is used to generate Pseudo Random Numbers
	/// </summary>
	class Random
	{
	private:
		inline static thread_local RandomEngine s_globalEngine{ 2463534242U };

	public:
		/// <summary>
		/// Seed the global engine with a manually provided value
		/// </summary>
		/// <param name="seed">Value to seed engine with</param>
		static void Seed(uint32_t seed);

		/// <summary>
		/// Seed the global engine with the current system time
		/// </summary>
		static void SeedWithTime();

		/// <summary>
		/// Return a random unsinged integer using the global RandomEngine
		/// </summary>
		static uint32_t UInt();

		/// <summary>
		/// Return a random unsinged integer in the range of [min, max] using the global RandomEngine
		/// </summary>
		/// <param name="min">Lower end of range, inclusive</param>
		/// <param name="max">Upper end of range, inclusive</param>
		static uint32_t UInt(uint32_t min, uint32_t max);

		/// <summary>
		/// Return a random integer using the global RandomEngine
		/// </summary>
		static int32_t Int();

		/// <summary>
		/// Return a random integer in the range of [min, max] using the global RandomEngine
		/// </summary>
		/// <param name="min">Lower end of range, inclusive</param>
		/// <param name="max">Upper end of range, inclusive</param>
		static int32_t Int(int32_t min, int32_t max);

		/// <summary>
		/// Return a random float in the range [0, 1) using the global RandomEngine
		/// </summary>
		static float Float();

		/// <summary>
		/// Return a random float in the range of [min, max) using the global RandomEngine
		/// </summary>
		/// <param name="min">Lower end of range, inclusive</param>
		/// <param name="max">Upper end of range, exclusive</param>
		static float Float(float min, float max);

		/// <summary>
		/// Return a random double using the global RandomEngine
		/// </summary>
		static double Double();

		/// <summary>
		/// Return a random double in the range of [min, max) using the global RandomEngine
		/// </summary>
		/// <param name="min">Lower end of range, inclusive</param>
		/// <param name="max">Upper end of range, exclusive</param>
		static double Double(double min, double max);
	};
}

#include <Random.inl>

#endif // !ULTREALITY_MATH_RANDOM_H
