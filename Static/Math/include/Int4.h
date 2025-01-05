#ifndef ULTREALITY_MATH_INT4_H
#define ULTREALITY_MATH_INT4_H

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
	/// 4D vector of 32-bit integer components
	/// </summary>
	struct Int4
	{
		int32_t x;
		int32_t y;
		int32_t z;
		int32_t w;

		Int4() = default;

		Int4(const Int4&) = default;
		Int4& operator=(const Int4&) = default;

		Int4(Int4&&) = default;
		Int4& operator=(Int4&&) = default;

		constexpr Int4(int32_t _x, int32_t _y, int32_t _z, int32_t _w) : x(_x), y(_y), z(_z), w(_w) {}
		explicit Int4(_In_reads_(4) const int32_t* pArray) : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}

		explicit Int4(_In_ A_VECTOR v) noexcept;
		Int4& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

#if (__cplusplus >= 202002L)
		bool operator==(const Int4&) const = default;
		auto operator<=>(const Int4&) const = default;
#endif

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	/// <summary>
	/// 4D vector of 32-bit integer components aligned on a 16 byte boundary
	/// </summary>
	ALIGNED_STRUCT(16) AInt4 : public Int4
	{
		using Int4::Int4;

		explicit AInt4(_In_ A_VECTOR v) noexcept;
		AInt4& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	namespace VEC
	{
		VECTOR VEC_CALLCONV LoadInt4(_In_ const Int4* pSource) noexcept;
		VECTOR VEC_CALLCONV LoadAInt4(_In_ const AInt4* pSource) noexcept;

		void VEC_CALLCONV StoreInt4(_Out_ Int4* pDestination, _In_ A_VECTOR v) noexcept;
		void VEC_CALLCONV StoreAInt4(_Out_ AInt4* pDestination, _In_ A_VECTOR v) noexcept;
	}

	namespace VEC4
	{
		bool VEC_CALLCONV EqualInt(A_VECTOR V1, A_VECTOR V2) noexcept;
		uint32_t VEC_CALLCONV EqualIntR(A_VECTOR V1, A_VECTOR V2) noexcept;
		bool VEC_CALLCONV NotEqualInt(A_VECTOR V1, A_VECTOR V2) noexcept;
	}
}

#include <Int4.inl>

#ifdef _PREFAST_
#pragma prefast(pop)
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // !ULTREALITY_MATH_INT4_H
