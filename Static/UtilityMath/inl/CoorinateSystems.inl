#ifndef ULTREALITY_MATH_COORDINATE_SYSTEMS_INL
#define ULTREALITY_MATH_COORDINATE_SYSTEMS_INL

#if defined(__GNUC__) or defined(__clang__)
#define FORCE_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline
#endif

namespace UltReality::Math
{
	FORCE_INLINE VECTOR SphericalToCartesian(float radius, float theta, float phi)
	{
		return Vector::Set(
			radius * ScalarSine(phi) * ScalarCos(theta),	// x
			radius * ScalarCos(phi),						// y
			radius * ScalarSine(phi) * ScalarSine(theta),	// z
			1.0f											// w
		);
	}
}

#endif // !ULTREALITY_MATH_COORDINATE_SYSTEMS_INL
