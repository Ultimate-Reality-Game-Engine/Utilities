#ifndef ULTREALITY_MATH_UINT3_H
#define ULTREALITY_MATH_UINT3_H

#include <SIMDVectorConfig.h>

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
	/// 3D vector of 32-bit unsigned integer components
	/// </summary>
	struct UInt3
	{
		uint32_t x;
		uint32_t y;
		uint32_t z;

		UInt3() = default;

		UInt3(const UInt3&) = default;
		UInt3& operator=(const UInt3&) = default;

		UInt3(UInt3&&) = default;
		UInt3& operator=(UInt3&&) = default;

		constexpr UInt3(uint32_t _x, uint32_t _y, uint32_t _z) noexcept : x(_x), y(_y), z(_z) {}
		explicit UInt3(_In_reads_(3) const uint32_t* pArray) noexcept : x(pArray[0]), y(pArray[1]), z(pArray[2]) {}

		explicit UInt3(_In_ A_VECTOR v) noexcept;
		UInt3& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

#if (__cplusplus >= 202002L)
		bool operator==(const UInt3&) const = default;
		auto operator<=>(const UInt3&) const = default;
#endif

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	/// <summary>
	/// 3D vector of 32-bit unsigned integer components aligned on a 16 byte boundary
	/// </summary>
	ALIGNED_STRUCT(16) AUInt3 : public UInt3
	{
		using UInt3::UInt3;

		explicit AUInt3(_In_ A_VECTOR v) noexcept;
		AUInt3& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	namespace Vector
	{
		VECTOR VEC_CALLCONV LoadUInt3(_In_ const UInt3* pSource) noexcept;
		VECTOR VEC_CALLCONV LoadAUInt3(_In_ const AUInt3* pSource) noexcept;

		void VEC_CALLCONV StoreUInt3(_Out_ UInt3* pDestination, _In_ A_VECTOR v) noexcept;
		void VEC_CALLCONV StoreAUInt3(_Out_ AUInt3* pDestination, _In_ A_VECTOR v) noexcept;
	}
}

#include <UInt3.inl>

#ifdef _PREFAST_
#pragma prefast(pop)
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // !ULTREALITY_MATH_UINT3_H
