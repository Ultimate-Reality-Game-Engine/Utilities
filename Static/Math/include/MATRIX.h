#ifndef ULTREALITY_MATH_SSE2_MATRIX_H
#define ULTREALITY_MATH_SSE2_MATRIX_H

#include <SSE2VectorConfig.h>

namespace UltReality::Math
{
	struct MATRIX;

	// Define alias to be used for 1st matrix type argument. Passed in registers for x86, and vector call convention; by reference otherwise
#if (defined(_M_IX86) || _VECTORCALL_ || __i386__) && !defined(_NO_INTRINSICS_)
	typedef const MATRIX A_MATRIX;
#else
	typedef const MATRIX& A_MATRIX;
#endif

	// Define alias to be used for 2nd+ matrix type arguments. Passed by reference
	typedef const MATRIX& B_MATRIX;

#if defined(_SSE2_INTRINSICS_) && !defined(_NO_INTRINSICS_)
	ALIGNED_STRUCT(16) MATRIX
#else
	struct MATRIX
#endif
	{
#if defined(_NO_INTRINSICS_)
		union
		{
			VECTOR r[4];
			struct
			{
				float _00, _01, _02, _03;
				float _10, _11, _12, _13;
				float _20, _21, _22, _23;
				float _30, _31, _32, _33;
			};
			float m[4][4];
		};
#else
		VECTOR r[4];
#endif

		MATRIX() = default;

		MATRIX(const MATRIX&) = default;

#if defined(_MSC_VER) && (_MSC_FULL_VER < 191426431)
		MATRIX& operator=(const MATRIX& m) noexcept
		{
			r[0] = m.r[0];
			r[1] = m.r[1];
			r[2] = m.r[2];
			r[3] = m.r[3];

			return *this;
		}
#else
		MATRIX& operator=(const MATRIX&) = default;

		MATRIX(MATRIX&&) = default;
		MATRIX& operator=(MATRIX&&) = default;
#endif

		constexpr MATRIX(A_VECTOR r0, A_VECTOR r1, A_VECTOR r2, B_VECTOR r3) noexcept : r{ r0, r1, r2, r3 } {}

		MATRIX(
			float m00, float m01, float m02, float m03,
			float m10, float m11, float m12, float m13,
			float m20, float m21, float m22, float m23,
			float m30, float m31, float m32, float m33
		) noexcept;

		explicit MATRIX(_In_reads_(16) const float* pArray) noexcept;

#if defined(_NO_INTRINSICS_)
		float operator() (size_t row, size_t column) const noexcept { return m[row][column]; }
		float& operator() (size_t row, size_t column) noexcept { return m[row][column]; }

		const float* operator[]size_t row) const noexcept { return m[row]; }
		float* operator[](size_t row) noexcept { return m[row]; }
#endif

		MATRIX operator+() const noexcept { return *this; }
		MATRIX operator-() const noexcept;

		MATRIX& VEC_CALLCONV operator+=(A_MATRIX m) noexcept;
		MATRIX& VEC_CALLCONV operator-=(A_MATRIX m) noexcept;
		MATRIX& VEC_CALLCONV operator*=(A_MATRIX m) noexcept;
		MATRIX& operator*=(float s) noexcept;
		MATRIX& operator/=(float s) noexcept;

		MATRIX VEC_CALLCONV operator+(A_MATRIX m) const noexcept;
		MATRIX VEC_CALLCONV operator-(A_MATRIX m) const noexcept;
		MATRIX VEC_CALLCONV operator*(A_MATRIX m) const noexcept;
		MATRIX operator*(float s) const noexcept;
		MATRIX operator/(float s) const noexcept;

		friend MATRIX VEC_CALLCONV operator*(float s, A_MATRIX m) noexcept;
	};

	namespace MAT
	{
		
	}
}

#include <MATRIX.inl>

#endif // !ULTREALITY_MATH_SSE2_MATRIX_H
