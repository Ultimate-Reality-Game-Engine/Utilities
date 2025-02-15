#include <gtest/gtest.h>

#include <Random.h>

#include <thread>
#include <atomic>

using namespace UltReality::Math;

TEST(RandomEngineTests, UIntInRange)
{
	RandomEngine engine(12345);

	// Test many values to verify they are within the full range
	for (int i = 0; i < 1000; i++)
	{
		uint32_t value = engine.UInt();

		EXPECT_GE(value, 0U);
		EXPECT_LE(value, UINT32_MAX);
	}
}

TEST(RandomEngineTests, UIntWithRange)
{
	RandomEngine engine(12345);
	const uint32_t min = 10;
	const uint32_t max = 20;

	for (int i = 0; i < 1000; i++)
	{
		uint32_t value = engine.UInt(min, max);

		// Inclusive range check
		EXPECT_GE(value, min);
		EXPECT_LE(value, max);
	}
}

TEST(RandomEngineTests, IntInRange)
{
	RandomEngine engine(12345);

	// Test many values to verify they are within the full range
	for (int i = 0; i < 1000; i++)
	{
		int32_t value = engine.Int();

		EXPECT_GE(value, INT32_MIN);
		EXPECT_LE(value, INT32_MAX);
	}
}

TEST(RandomEngineTests, IntWithRange)
{
	RandomEngine engine(12345);
	const int32_t min = -20;
	const int32_t max = 20;

	for (int i = 0; i < 1000; i++)
	{
		int32_t value = engine.Int(min, max);

		// Inclusive range check
		EXPECT_GE(value, min);
		EXPECT_LE(value, max);
	}
}

TEST(RandomEngineTests, FloatInRange) {
	RandomEngine engine(12345);
	
	for (int i = 0; i < 1000; i++)
	{
		float value = engine.Float();
		EXPECT_GE(value, 0.0f);
		// Since the float version is [min, max) the value should be strictly less than max.
		EXPECT_LT(value, 1.0f);
	}
}

TEST(RandomEngineTests, FloatWithRange) {
	RandomEngine engine(12345);
	const float min = 1.0f;
	const float max = 2.0f;
	for (int i = 0; i < 1000; i++)
	{
		float value = engine.Float(min, max);
		EXPECT_GE(value, min);
		// Since the float version is [min, max) the value should be strictly less than max.
		EXPECT_LT(value, max);
	}
}

TEST(RandomEngineTests, DoubleInRange) {
	RandomEngine engine(12345);

	for (int i = 0; i < 1000; i++)
	{
		double value = engine.Double();
		EXPECT_GE(value, 0.0);
		// Since the double version is [min, max) the value should be strictly less than max.
		EXPECT_LT(value, 1.0);
	}
}

TEST(RandomEngineTests, DoubleWithRange) {
	RandomEngine engine(12345);
	const double min = 1.0;
	const double max = 2.0;
	for (int i = 0; i < 1000; i++)
	{
		double value = engine.Double(min, max);
		EXPECT_GE(value, min);
		// Since the float version is [min, max) the value should be strictly less than max.
		EXPECT_LT(value, max);
	}
}

/*
	Note: Testing thread - local behavior is inherently a bit non-deterministic. 
	In this case, by seeding with different known values in two threads and comparing outputs, 
	we can check that each thread’s RNG is independent.
*/
TEST(RandomStaticTests, GlobalEngineThreadLocal)
{
	// We'll launch two threads that call static Random methods.
	// Each thread seeds its own thread_local engine with a known value and
	// generates a value. We expect that different seeds produce different results.

	std::atomic<uint32_t> value1{ 0 };
	std::atomic<uint32_t> value2{ 0 };

	auto threadFunc = [](uint32_t seed, std::atomic<uint32_t>& outVal) {
		Random::Seed(seed);  // This affects the thread_local global engine.
		outVal.store(Random::UInt());
		};

	std::thread t1(threadFunc, 11111, std::ref(value1));
	std::thread t2(threadFunc, 22222, std::ref(value2));

	t1.join();
	t2.join();

	// Since the seeds differ, the results should likely be different.
	EXPECT_NE(value1.load(), value2.load());
}