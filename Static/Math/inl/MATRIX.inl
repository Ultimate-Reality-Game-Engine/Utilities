#ifndef ULTREALITY_MATH_SSE2_MATRIX_INL
#define ULTREALITY_MATH_SSE2_MATRIX_INL

#if defined(__GNUC__) or defined(__clang__)
#define FORCE_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline
#endif

namespace UltReality::Math
{
	FORCE_INLINE MATRIX::MATRIX(
		float m00, float m01, float m02, float m03,
		float m10, float m11, float m12, float m13,
		float m20, float m21, float m22, float m23,
		float m30, float m31, float m32, float m33
	) noexcept
	{
		r[0] = VECTOR::VectorSet(m00, m01, m02, m03);
		r[1] = VECTOR::VectorSet(m10, m11, m12, m13);
		r[2] = VECTOR::VectorSet(m20, m21, m22, m23);
		r[3] = VECTOR::VectorSet(m30, m31, m32, m33);
	}

	struct Float4;

	_Use_decl_annotations_
		FORCE_INLINE MATRIX::MATRIX(const float* pArray) noexcept
	{
#if defined(DEBUG) || defined(_DEBUG)
		assert(pArray != nullptr);
#endif
		r[0] = VECTOR::LoadFloat4(reinterpret_cast<const Float4*>(pArray));
		r[1] = VECTOR::LoadFloat4(reinterpret_cast<const Float4*>(pArray + 4));
		r[2] = VECTOR::LoadFloat4(reinterpret_cast<const Float4*>(pArray + 8));
		r[3] = VECTOR::LoadFloat4(reinterpret_cast<const FLoat4*>(pArray + 12));
	}

	FORCE_INLINE MATRIX MATRIX::operator-() const noexcept
	{
		MATRIX m;
		m.r[0] = VECTOR::VectorNegate(r[0]);
		m.r[1] = VECTOR::VectorNegate(r[1]);
		m.r[2] = VECTOR::VectorNegate(r[2]);
		m.r[3] = VECTOR::VectorNegate(r[3]);

		return m;
	}

	FORCE_INLINE MATRIX& VEC_CALLCONV MATRIX::operator+=(A_MATRIX m) noexcept
	{
		r[0] = VECTOR::VectorAdd(r[0], m.r[0]);
		r[1] = VECTOR::VectorAdd(r[1], m.r[1]);
		r[2] = VECTOR::VectorAdd(r[2], m.r[2]);
		r[3] = VECTOR::VectorAdd(r[3], m.r[3]);

		return *this;
	}

	FORCE_INLINE MATRIX& VEC_CALLCONV MATRIX::operator-=(A_MATRIX m) noexcept
	{
		r[0] = VECTOR::VectorSubtract(r[0], m.r[0]);
		r[1] = VECTOR::VectorSubtract(r[1], m.r[1]);
		r[2] = VECTOR::VectorSubtract(r[2], m.r[2]);
		r[3] = VECTOR::VectorSubtract(r[3], m.r[3]);

		return *this;
	}

	FORCE_INLINE MATRIX& VEC_CALLCONV MATRIX::operator*=(A_MATRIX m) noexcept
	{
		*this = MATRIX::MatrixMultiply(*this, m);

		return *this;
	}

	FORCE_INLINE MATRIX& VEC_CALLCONV MATRIX::operator*=(float s) noexcept
	{
		r[0] = VECTOR::VectorScale(r[0], s);
		r[1] = VECTOR::VectorScale(r[1], s);
		r[2] = VECTOR::VectorScale(r[2], s);
		r[3] = VECTOR::VectorScale(r[3], s);

		return *this;
	}

	FORCE_INLINE MATRIX& VEC_CALLCONV MATRIX::operator/=(float s) noexcept
	{
#if defined(_NO_INTRINSICS_)
		VECTOR vS = VECTOR::VectorReplicate(s);
		r[0] = VECTOR::VectorDivide(r[0], vS);
		r[1] = VECTOR::VectorDivide(r[1], vS);
		r[2] = VECTOR::VectorDivide(r[2], vS);
		r[3] = VECTOR::VectorDivide(r[3], vS);

#elif defined(_SSE2_INTRINSICS_)
		__m128 vS = _mm_set_ps1(s);
		r[0] = _mm_div_ps(r[0], vS);
		r[1] = _mm_div_ps(r[1], vS);
		r[2] = _mm_div_ps(r[2], vS);
		r[3] = _mm_div_ps(r[3], vS);
#endif

		return *this;
	}

	FORCE_INLINE MATRIX VEC_CALLCONV MATRIX::operator+(A_MATRIX m) const noexcept
	{
		MATRIX R;
		R.r[0] = VECTOR::VectorAdd(r[0], m.r[0]);
		R.r[1] = VECTOR::VectorAdd(r[1], m.r[1]);
		R.r[2] = VECTOR::VectorAdd(r[2], m.r[2]);
		R.r[3] = VECTOR::VectorAdd(r[3], m.r[3]);

		return R;
	}

	FORCE_INLINE MATRIX VEC_CALLCONV MATRIX::operator-(A_MATRIX m) const noexcept
	{
		MATRIX R;
		R.r[0] = VECTOR::VectorSubtract(r[0], m.r[0]);
		R.r[1] = VECTOR::VectorSubtract(r[1], m.r[1]);
		R.r[2] = VECTOR::VectorSubtract(r[2], m.r[2]);
		R.r[3] = VECTOR::VectorSubtract(r[3], m.r[3]);

		return R;
	}

	FORCE_INLINE MATRIX VEC_CALLCONV MATRIX::operator*(A_MATRIX m) const noexcept
	{
		return MATRIX::MatrixMultiply(*this, M);
	}

	FORCE_INLINE MATRIX MATRIX::operator*(float s) const noexcept
	{
		MATRIX R;
		R.r[0] = VECTOR::VectorScale(r[0], s);
		R.r[1] = VECTOR::VectorScale(r[1], s);
		R.r[2] = VECTOR::VectorScale(r[2], s);
		R.r[3] = VECTOR::VectorScale(r[3], s);

		return R;
	}

	FORCE_INLINE MATRIX MATRIX::operator/(float s) const noexcept
	{
#if defined(_NO_INTRINSICS_)
		VECTOR vS = VECTOR::VectorReplicate(s);
		MATRIX R;
		R.r[0] = VECTOR::VectorDivide(r[0], vS);
		R.r[1] = VECTOR::VectorDivide(r[1], vS);
		R.r[2] = VECTOR::VectorDivide(r[2], vS);
		R.r[3] = VECTOR::VectorDivide(r[3], vS);

#elif defined(_SSE2_INTRINSICS_)
		__m128 vS = _mm_set_ps1(s);
		MATRIX R;
		R.r[0] = _mm_div_ps(r[0], vS);
		R.r[1] = _mm_div_ps(r[1], vS);
		R.r[2] = _mm_div_ps(r[2], vS);
		R.r[3] = _mm_div_ps(r[3], vS);
#endif

		return R;
	}

	FORCE_INLINE MATRIX VEC_CALLCONV operator*(float s, A_MATRIX m) noexcept
	{
		MATRIX R;
		R.r[0] = VECTOR::VectorScale(m.r[0], s);
		R.r[1] = VECTOR::VectorScale(m.r[1], s);
		R.r[2] = VECTOR::VectorScale(m.r[2], s);
		R.r[3] = VECTOR::VectorScale(m.r[3], s);

		return R;
	}
}

#endif // !ULTREALITY_MATH_SSE2_MATRIX_INL
