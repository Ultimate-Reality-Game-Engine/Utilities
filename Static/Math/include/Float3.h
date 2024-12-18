#ifndef ULTREALITY_MATH_FLOAT3_H
#define ULTREALITY_MATH_FLOAT3_H

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
	/// 3D vector of 32-bit floating point components
	/// </summary>
	struct Float3
	{
		float x;
		float y;
		float z;

		Float3() = default;

		Float3(const Float3&) = default;
		Float3& operator=(const Float3&) = default;

		Float3(Float3&&) = default;
		Float3& operator=(Float3&&) = default;

		constexpr Float3(float _x, float _y, float _z) noexcept : x(_x), y(_y), z(_z) {}
		explicit Float3(_In_reads_(3) const float* pArray) noexcept : x(pArray[0]), y(pArray[1]), z(pArray[2]) {}

		explicit FLoat3(_In_ A_VECTOR v) noexcept;
		Float3& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

#if (__cplusplus >= 202002L)
		bool operator==(const Float3&) const = default;
		auto operator<=>(const Float3&) const = default;
#endif

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	/// <summary>
	/// 3D vector of 32-bit floating point components aligned on a 16 byte boundary
	/// </summary>
	ALIGNED_STRUCT(16) AFloat3 : public Float3
	{
		using Float3::Float3;

		explicit AFloat3(_In_ A_VECTOR v) noexcept;
		AFloat3& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	namespace VEC
	{
		VECTOR VEC_CALLCONV LoadFloat3(_In_ const Float3* pSource) noexcept;
		VECTOR VEC_CALLCONV LoadAFloat3(_In_ const AFloat3* pSource) noexcept;

		void VEC_CALLCONV StoreFloat3(_Out_ Float3* pDestination, _In_ A_VECTOR v) noexcept;
		void VEC_CALLCONV StoreAFloat3(_Out_ AFloat3* pDestination, _In_ A_VECTOR v) noexcept;
	}
}

#include <Float3.inl>

#ifdef _PREFAST_
#pragma prefast(pop)
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // !ULTREALITY_MATH_FLOAT3_H
