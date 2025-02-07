#ifndef ULTREALITY_MATH_UINT4_H
#define ULTREALITY_MATH_UINT4_H

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
	/// 4D vector of 32-bit unsigned integer components
	/// </summary>
	struct UInt4
	{
		uint32_t x;
		uint32_t y;
		uint32_t z;
		uint32_t w;

		UInt4() = default;

		UInt4(const UInt4&) = default;
		UInt4& operator=(const UInt4&) = default;

		UInt4(UInt4&&) = default;
		UInt4& operator=(UInt4&&) = default;

		constexpr UInt4(uint32_t _x, uint32_t _y, uint32_t _z, uint32_t _w) noexcept : x(_x), y(_y), z(_z), w(_w) {}
		explicit UInt4(_In_reads_(4) const uint32_t* pArray) noexcept : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}

		explicit UInt4(_In_ A_VECTOR v) noexcept;
		UInt4& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

#if (__cplusplus >= 202002L)
		bool operator==(const UInt4&) const = default;
		auto operator<=>(const UInt4&) const = default;
#endif

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	/// <summary>
	/// 4D vector of 32-bit unsigned integer components aligned on a 16 byte boundary
	/// </summary>
	ALIGNED_STRUCT(16) AUInt4 : public UInt4
	{
		using UInt4::UInt4;

		explicit AUInt4(_In_ A_VECTOR v) noexcept;
		AUInt4& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	namespace VEC
	{
		VECTOR VEC_CALLCONV LoadUInt4(_In_ const UInt4* pSource) noexcept;
		VECTOR VEC_CALLCONV LoadAUInt4(_In_ const AUInt4* pSource) noexcept;

		void VEC_CALLCONV StoreUInt4(_Out_ UInt4* pDestination, _In_ A_VECTOR v) noexcept;
		void VEC_CALLCONV StoreAUInt4(_Out_ AUInt4* pDestination, _In_ A_VECTOR v) noexcept;
	}
}
#include <UInt4.inl>

#ifdef _PREFAST_
#pragma prefast(pop)
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // !ULTREALITY_MATH_UINT4_H
