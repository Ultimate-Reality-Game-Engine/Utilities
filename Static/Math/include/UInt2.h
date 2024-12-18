#ifndef ULTREALITY_MATH_UINT2_H
#define ULTREALITY_MATH_UINT2_H

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
	/// 2D vector of 32-bit unsigned integer components
	/// </summary>
	struct UInt2
	{
		uint32_t x;
		uint32_t y;

		UInt2() = default;
		
		UInt2(const UInt2&) = default;
		UInt2& operator=(const UInt2&) = default;

		UInt2(UInt2&&) = default;
		UInt2& operator=(UInt2&&) = default;

		constexpr UInt2(uint32_t _x, uint32_t _y) noexcept : x(_x), y(_y) {}
		explicit UInt2(_In_reads_(2) const uint32_t* pArray) noexcept : x(pArray[0]), y(pArray[1]) {}

		explicit UInt2(_In_ A_VECTOR v) noexcept;
		UInt2& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

#if (__cplusplus >= 202002L)
		bool operator==(const UInt2&) const = default;
		auto operator<=>(const UInt2&) const = default;
#endif

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	/// <summary>
	/// 2D vector of 32-bit unsigned integer components aligned on a 16 byte boundary
	/// </summary>
	ALIGNED_STRUCT(16) AUInt2 : public UInt2
	{
		using UInt2::UInt2;

		explicit AUInt2(_In_ A_VECTOR v) noexcept;
		AUInt2& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	namespace VEC
	{
		VECTOR VEC_CALLCONV LoadUInt2(_In_ const UInt2* pSource) noexcept;
		VECTOR VEC_CALLCONV LoadAUInt2(_In_ const AUInt2* pSource) noexcept;

		void VEC_CALLCONV StoreUInt2(_Out_ UInt2* pDestination, _In_ A_VECTOR v) noexcept;
		void VEC_CALLCONV StoreAUInt2(_Out_ AUInt2* pDestination, _In_ A_VECTOR v) noexcept;
	}
}

#include <UInt2.inl>

#ifdef _PREFAST_
#pragma prefast(pop)
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // !ULTREALITY_MATH_UINT2_H
