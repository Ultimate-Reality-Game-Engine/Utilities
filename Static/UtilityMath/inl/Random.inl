#ifndef ULTREALITY_MATH_RANDOM_INL
#define ULTREALITY_MATH_RANDOM_INL

#if defined(__GNUC__) or defined(__clang__)
#define FORCE_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline
#endif

namespace UltReality::Math
{
	FORCE_INLINE uint32_t RandomEngine::xorshift32()
	{
		uint32_t x = m_state;
		x ^= x << 13;
		x ^= x >> 17;
		x ^= x << 5;
		m_state = x;
		return x;
	}

	FORCE_INLINE uint32_t RandomEngine::RandomUInt()
	{
		return xorshift32();
	}

	FORCE_INLINE uint32_t RandomEngine::RandomUInt(uint32_t min, uint32_t max)
	{
		uint32_t range = max - min + 1;
		return min + (xorshift32() % range);
	}

	FORCE_INLINE int32_t RandomEngine::RandomInt()
	{
		return static_cast<int32_t>(xorshift32());
	}

	FORCE_INLINE int32_t RandomEngine::RandomInt(int32_t min, int32_t max)
	{
		uint32_t range = static_cast<uint32_t>(max - min + 1);
		return min + static_cast<int32_t>(xorshift32() % range);
	}

	FORCE_INLINE float RandomEngine::RandomFloat()
	{
		return static_cast<float>(xorshift32()) / 4294967296.0f;
	}

	FORCE_INLINE float RandomEngine::RandomFloat(float min, float max)
	{
		return min + RandomFloat() * (max - min);
	}

	FORCE_INLINE double RandomEngine::RandomDouble()
	{
		return static_cast<double>(xorshift32()) / 4294967296.0;
	}

	FORCE_INLINE double RandomEngine::RandomDouble(double min, double max)
	{
		return min + RandomDouble() * (max - min);
	}

	FORCE_INLINE uint32_t Random::UInt()
	{
		return s_globalEngine.RandomUInt();
	}

	FORCE_INLINE uint32_t Random::UInt(uint32_t min, uint32_t max)
	{
		return s_globalEngine.RandomUInt(min, max);
	}

	FORCE_INLINE int32_t Random::Int()
	{
		return s_globalEngine.RandomInt();
	}

	FORCE_INLINE int32_t Random::Int(int32_t min, int32_t max)
	{
		return s_globalEngine.RandomInt(min, max);
	}

	FORCE_INLINE float Random::Float()
	{
		return s_globalEngine.RandomFloat();
	}

	FORCE_INLINE float Random::Float(float min, float max)
	{
		return s_globalEngine.RandomFloat(min, max);
	}

	FORCE_INLINE double Random::Double()
	{
		return s_globalEngine.RandomDouble();
	}

	FORCE_INLINE double Random::Double(double min, double max)
	{
		return s_globalEngine.RandomDouble(min, max);
	}
}
#endif // !ULTREALITY_MATH_RANDOM_INL
