#ifndef ULTREALITY_MATH_INT3_H
#define ULTREALITY_MATH_INT3_H

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
	/// 3D vector of 32-bit integer components
	/// </summary>
	struct Int3
	{
		int32_t x;
		int32_t y;
		int32_t z;

		Int3() = default;

		Int3(const Int3&) = default;
		Int3& operator=(const Int3&) = default;

		Int3(Int3&&) = default;
		Int3& operator=(Int3&&) = default;

		constexpr Int3(int32_t _x, int32_t _y, int32_t _z) noexcept : x(_x), y(_y), z(_z) {}
		explicit Int3(_In_reads_(3) const int32_t* pArray) noexcept : x(pArray[0]), y(pArray[1]), z(pArray[2]) {}

		explicit Int3(_In_ A_VECTOR v) noexcept;
		Int3& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

#if (__cplusplus >= 202002L)
		bool operator==(const Int3&) const = default;
		auto operator<=>(const Int3&) const = default;
#endif

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	/// <summary>
	/// 3D vector of 32-bit integer components aligned on a 16 byte boundary
	/// </summary>
	ALIGNED_STRUCT(16) AInt3 : public Int3
	{
		using Int3::Int3;

		explicit AInt3(_In_ A_VECTOR v) noexcept;
		AInt3& VEC_CALLCONV operator=(_In_ A_VECTOR v) noexcept;

		VECTOR VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_VECTOR v) noexcept;
	};

	namespace Vector
	{
		VECTOR VEC_CALLCONV LoadInt3(_In_ const Int3* pSource) noexcept;
		VECTOR VEC_CALLCONV LoadAInt3(_In_ const AInt3* pSource) noexcept;

		void VEC_CALLCONV StoreInt3(_Out_ Int3* pDestination, _In_ A_VECTOR v) noexcept;
		void VEC_CALLCONV StoreAInt3(_Out_ AInt3* pDestination, _In_ A_VECTOR v) noexcept;
	}

	namespace Vector3
	{
		bool VEC_CALLCONV EqualInt(A_VECTOR V1, A_VECTOR V2) noexcept;
		uint32_t VEC_CALLCONV EqualIntR(A_VECTOR V1, A_VECTOR V2) noexcept;
		bool VEC_CALLCONV NotEqualInt(A_VECTOR V1, A_VECTOR V2) noexcept;
	}
}

#include <Int3.inl>

#ifdef _PREFAST_
#pragma prefast(pop)
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // !ULTREALITY_MATH_INT3_H
