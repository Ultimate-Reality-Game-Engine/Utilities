#ifndef ULTREALITY_MATH_FLOAT4X3_H
#define ULTREALITY_MATH_FLOAT4X3_H

#include <MATRIX.h>

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

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wgnu-anonymous-struct"
#pragma clang diagnostic ignored "-Wnested-anon-types"
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif

namespace UltReality::Math
{
	/// <summary>
	/// 4x3 matrix of 32-bit floating point components
	/// </summary>
	struct Float4x3
	{
		union
		{
			struct
			{
				float _00, _01, _02;
				float _10, _11, _12;
				float _20, _21, _22;
				float _30, _31, _32;
			};
			float m[4][3];
		};

		Float4x3() = default;

		Float4x3(const Float4x3&) = default;
		Float4x3& operator=(const Float4x3&) = default;

		Float4x3(Float4x3&&) = default;
		Float4x3& operator=(Float4x3&&) = default;

		constexpr Float4x3(
			float m00, float m01, float m02, 
			float m10, float m11, float m12, 
			float m20, float m21, float m22, 
			float m30, float m31, float m32
		) noexcept : 
			_00(m00), _01(m01), _02(m02), 
			_10(m10), _11(m11), _12(m12), 
			_20(m20), _21(m21), _22(m22), 
			_30(m30), _31(m31), _32(m32)
		{}

		explicit Float4x3(_In_reads_(12) const float* pArray) noexcept;

		explicit Float4x3(_In_ A_MATRIX m) noexcept;
		Float4x3& VEC_CALLCONV operator=(_In_ A_MATRIX m) noexcept;

		float operator() (size_t row, size_t column) const noexcept { return m[row][column]; }
		float& operator() (size_t row, size_t column) noexcept { return m[row][column]; }

		const float* operator[](size_t row) const noexcept { return m[row]; }
		float* operator[](size_t row) noexcept { return m[row]; }

#if (__cplusplus >= 202002L)
		bool operator==(const Float4x3&) const = default;
		auto operator<=>(const Float4x3&) const = default;
#endif

		MATRIX VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_MATRIX m) noexcept;
	};

	/// <summary>
	/// 4x3 matrix of 32-bit floating point components aligned on a 16 byte boundary
	/// </summary>
	ALIGNED_STRUCT(16) AFloat4x3 : public Float4x3
	{
		using Float4x3::Float4x3;

		explicit AFloat4x3(_In_ A_MATRIX m) noexcept;
		AFloat4x3& VEC_CALLCONV operator=(_In_ A_MATRIX m) noexcept;

		MATRIX VEC_CALLCONV Load() noexcept;
		void VEC_CALLCONV Store(_In_ A_MATRIX m) noexcept;
	};

	namespace MAT
	{
		MATRIX VEC_CALLCONV LoadFloat4x3(_In_ const Float4x3* pSource) noexcept;
		MATRIX VEC_CALLCONV LoadAFloat4x3(_In_ const AFloat4x3* pSource) noexcept;

		void VEC_CALLCONV StoreFloat4x3(_Out_ Float4x3* pDestination, _In_ A_MATRIX m) noexcept;
		void VEC_CALLCONV StoreAFloat4x3(_Out_ AFloat4x3* pDestination, _In_ A_MATRIX m) noexcept;
	}
}

#include <Float4x3.inl>

#ifdef __clang__
#pragma clang diagnostic pop
#endif
#ifdef _PREFAST_
#pragma prefast(pop)
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // !ULTREALITY_MATH_FLOAT4X3_H
