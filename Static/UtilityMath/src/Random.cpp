#include <Random.h>

#include <ctime>

namespace UltReality::Math
{
	RandomEngine::RandomEngine(uint32_t seed) : m_state(seed)
	{}

	void RandomEngine::Seed(uint32_t seed)
	{
		m_state = (seed != 0U) ? seed : 2463534242U;
	}

	void RandomEngine::SeedWithTime()
	{
		Seed(static_cast<uint32_t>(std::time(nullptr)));
	}

	void Random::Seed(uint32_t seed)
	{
		s_globalEngine.Seed(seed);
	}

	void Random::SeedWithTime()
	{
		s_globalEngine.SeedWithTime();
	}
}