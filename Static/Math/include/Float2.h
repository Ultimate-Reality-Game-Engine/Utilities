#ifndef ULTREALITY_MATH_FLOAT2_H
#define ULTREALITY_MATH_FLOAT2_H

#include <SSE2VectorConfig.h>

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4068 4201 4365 4324 4820)
	 // C4068: ignore unknown pragmas
	 // C4201: nonstandard extension used : nameless struct/union
	 // C4365: Off by default noise
	 // C4324/4820: padding warnings
#endif

#ifdef _PREFAST_
#pragma prefast(push)
#pragma prefast(disable : 25000, "A_VECTOR is 16 bytes")
#endif

namespace UltReality::Math
{
	/// <summary>
	/// 2D vector of 32-bit floating point components
	/// </summary>
	struct Float2
	{
		float x;
		float y;

		Float2() = default;

		Float2(const Float2&) = default;
		Float2& operator=(const Float2&) = default;

		Float2(Float2&&) = default;
		Float2& operator=(Float2&&) = default;

		constexpr Float2(float _x, float _y) noexcept : x(_x), y(_y) {}
		explicit Float2(_In_reads_(2) const float* pArray) noexcept : x(pArray[0]), y(pArray[1]) {}

		explicit Float2(_In_ A_VECTOR v) noexcept;
		Float2& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

#if(__cplusplus >= 202002L)
		bool operator==(const Float2&) const = default;
		auto operator<=>(const Float2&) const = default;
#endif

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	/// <summary>
	/// 2D vector of 32-bit floating point components aligned on a 16 byte boundary
	/// </summary>
	ALIGNED_STRUCT(16) AFloat2 : public Float2
	{
		using Float2::Float2;

		explicit AFloat2(_In_ A_VECTOR v) noexcept;
		AFloat2& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	namespace VEC
	{
		VECTOR VEC_CALLCONV LoadFloat2(_In_ const Float2* pSource) noexcept;
		VECTOR VEC_CALLCONV LoadAFloat2(_In_ const AFloat2* pSource) noexcept;

		void VEC_CALLCONV StoreFloat2(_Out_ Float2* pDestination, _In_ A_VECTOR v) noexcept;
		void VEC_CALLCONV StoreAFloat2(_Out_ AFloat2* pDestination, _In_ A_VECTOR v) noexcept;
	}
}

#include <Float2.inl>

#ifdef _PREFAST_
#pragma prefast(pop)
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // !ULTREALITY_MATH_FLOAT2_H
