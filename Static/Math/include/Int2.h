#ifndef ULTREALITY_MATH_INT2_H
#define ULTREALITY_MATH_INT2_H

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
	/// 2D vector of 32-bit signed integer components
	/// </summary>
	struct Int2
	{
		int32_t x;
		int32_t y;

		Int2() = default;

		Int2(const Int2&) = default;
		Int2& operator=(const Int2&) = default;

		Int2(Int2&&) = default;
		Int2& operator=(Int2&&) = default;

		constexpr Int2(int32_t _x, int32_t _y) noexcept : x(_x), y(_y) {}
		explicit Int2(_In_reads_(2) const int32_t* pArray) noexcept : x(pArray[0]), y(pArray[1]) {}

		explicit Int2(_In_ A_VECTOR v) noexcept;
		Int2& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

#if (__cplusplus >= 202002L)
		bool operator==(const Int2&) const = default;
		auto operator<=>(const Int2&) const = default;
#endif
		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	/// <summary>
	/// 2D vector of 32-bit signed integer components aligned on a 16 byte boundary
	/// </summary>
	ALIGNED_STRUCT(16) AInt2 : public Int2
	{	
		using Int2::Int2;

		explicit AInt2(_In_ A_VECTOR v) noexcept;
		AInt2& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	namespace VEC
	{
		VECTOR VEC_CALLCONV LoadInt2(_In_ const Int2* pSource) noexcept;
		VECTOR VEC_CALLCONV LoadAInt2(_In_ const AInt2* pSource) noexcept;

		void VEC_CALLCONV StoreInt2(_Out_ Int2* pDestination, _In_ A_VECTOR v) noexcept;
		void VEC_CALLCONV StoreAInt2(_Out_ AInt2* pDestination, _In_ A_VECTOR v) noexcept;
	}

	namespace VEC2
	{
		bool VEC_CALLCONV EqualInt(A_VECTOR V1, A_VECTOR V2) noexcept;
		uint32_t VEC_CALLCONV EqualIntR(A_VECTOR V1, A_VECTOR V2) noexcept;
		bool VEC_CALLCONV NotEqualInt(A_VECTOR V1, A_VECTOR V2) noexcept;
	}
}

#include <Int2.inl>

#ifdef _PREFAST_
#pragma prefast(pop)
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // !ULTREALITY_MATH_INT2_H
