#ifndef ULTREALITY_MATH_FLOAT4_H
#define ULTREALITY_MATH_FLOAT4_H

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
	/// 4D vector of 32-bit floating point components
	/// </summary>
	struct Float4
	{
		float x;
		float y;
		float z;
		float w;

		Float4() = default;

		Float4(const Float4&) = default;
		Float4& operator=(const Float4&) = default;

		Float4(Float4&&) = default;
		Float4& operator=(Float4&&) = default;

		constexpr Float4(float _x, float _y, float _z, float _w) : x(_x), y(_y), z(_z), w(_w) {}
		explicit Float4(_In_reads_(4) const float* pArray) noexcept : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}

		explicit Float4(_In_ A_VECTOR v) noexcept;
		Float4& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

#if (__cplusplus >= 202002L)
		bool operator==(const Float4&) const = default;
		auto operator<=>(const Float4&) const = default;
#endif

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	/// <summary>
	/// 4D vector of 32-bit floating point components aligned on a 16 byte boundary
	/// </summary>
	ALIGNED_STRUCT(16) AFloat4 : public Float4
	{
		using Float4::Float4;

		explicit AFloat4(_In_ A_VECTOR v) noexcept;
		AFloat4& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	namespace VEC
	{
		VECTOR VEC_CALLCONV LoadFloat4(_In_ const Float4* pSource) noexcept;
		VECTOR VEC_CALLCONV LoadAFloat4(_In_ const AFloat4* pSource) noexcept;

		void VEC_CALLCONV StoreFloat4(_Out_ Float4* pDestination, _In_ A_VECTOR v) noexcept;
		void VEC_CALLCONV StoreAFloat4(_Out_ AFloat4* pDestination, _In_ A_VECTOR v) noexcept;
	}
}

#include <Float4.inl>

#ifdef _PREFAST_
#pragma prefast(pop)
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // !ULTREALITY_MATH_FLOAT4_H
