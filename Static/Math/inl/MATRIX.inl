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
		return MATRIX::MatrixMultiply(*this, m);
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

#if !defined(__NO_INTRINSICS_) && defined(_MSC_VER) && !defined(__INTEL_COMPILER)
#pragma float_control(push)
#pragma float_control(precise, on)
#endif

	namespace MAT
	{
		FORCE_INLINE bool VEC_CALLCONV IsNaN(A_MATRIX m) noexcept
		{
#if defined(_NO_INTRINSICS_)
			size_t i = 16;
			const uint32_t* pWork = reinterpret_cast<const uint32_t*>(&m.m[0][0]);
			do
			{
				// Fetch value into integer unit
				uint32_t uTest = pWork[0];
				// Remove sign
				uTest &= 0x7FFFFFFFU;
				// NaN is 0x7F800001 through 0x7FFFFFFF inclusive
				uTest -= 0x7F800001U;
				if (uTest < 0x007FFFFFU)
				{
					break; // NaN detected
				}
				++pWork; // Advance to next entry
			} while (--i);

			return i != 0; // i == 0 if nothing matched

#elif defined(_SSE2_INTRINSICS_)
			// Load in registers
			VECTOR vX = m.r[0];
			VECTOR vY = m.r[1];
			VECTOR vZ = m.r[2];
			VECTOR vW = m.r[3];
			// Test themselves to check for NaN
			vX = _mm_cmpneq_ps(vX, vX);
			vY = _mm_cmpneq_ps(vY, vY);
			vZ = _mm_cmpneq_ps(vZ, vZ);
			vW = _mm_cmpneq_ps(vW, vW);
			// Or all the results
			vX = _mm_or_ps(vX, vZ);
			vY = _mm_or_ps(vY, vW);
			vX = _mm_or_ps(vX, vY);
			
			// If any tested true, return true
			return (_mm_movemask_ps(vX) != 0);
#endif
		}

#if !defined(__NO_INTRINSICS_) && defined(_MSC_VER) && !defined(__INTEL_COMPILER)
#pragma float_control(pop)
#endif

		FORCE_INLINE bool VEC_CALLCONV IsInfinite(A_MATRIX m) noexcept
		{
#if defined(_NO_INTRINSICS_)
			size_t i = 16;
			const uint32_t* pWork = reinterpret_cast<const uint32_t*>(&m.m[0][0]);
			do {
				// Fetch value into integer unit
				uint32_t uTest = pWork[0];
				// Remove sign
				uTest &= 0x7FFFFFFFU;
				// INF is 0x7F800000
				if (uTest == 0x7F800000U)
				{
					break;      // INF found
				}
				++pWork;        // Next entry
			} while (--i);
			return (i != 0);      // i == 0 if nothing matched

#elif defined(_SSE2_INTRINSICS_)
			// Mask off the sign bits
			VECTOR vTemp1 = _mm_and_ps(m.r[0], g_AbsMask);
			VECTOR vTemp2 = _mm_and_ps(m.r[1], g_AbsMask);
			VECTOR vTemp3 = _mm_and_ps(m.r[2], g_AbsMask);
			VECTOR vTemp4 = _mm_and_ps(m.r[3], g_AbsMask);
			// Compare to infinity
			vTemp1 = _mm_cmpeq_ps(vTemp1, g_Infinity);
			vTemp2 = _mm_cmpeq_ps(vTemp2, g_Infinity);
			vTemp3 = _mm_cmpeq_ps(vTemp3, g_Infinity);
			vTemp4 = _mm_cmpeq_ps(vTemp4, g_Infinity);
			// Or the answers together
			vTemp1 = _mm_or_ps(vTemp1, vTemp2);
			vTemp3 = _mm_or_ps(vTemp3, vTemp4);
			vTemp1 = _mm_or_ps(vTemp1, vTemp3);
			
			// If any are infinity, the signs are true.
			return (_mm_movemask_ps(vTemp1) != 0);
#endif
		}

		FORCE_INLINE bool VEC_CALLCONV IsIdentity(A_MATRIX m) noexcept
		{
#if defined(_NO_INTRINSICS_)
			// Use the integer pipeline to reduce branching to a minimum
			const uint32_t* pWork = reinterpret_cast<const uint32_t*>(&m.m[0][0]);
			// Convert 1.0f to zero and or them together
			uint32_t uOne = pWork[0] ^ 0x3F800000U;
			// Or all the 0.0f entries together
			uint32_t uZero = pWork[1];
			uZero |= pWork[2];
			uZero |= pWork[3];
			// 2nd row
			uZero |= pWork[4];
			uOne |= pWork[5] ^ 0x3F800000U;
			uZero |= pWork[6];
			uZero |= pWork[7];
			// 3rd row
			uZero |= pWork[8];
			uZero |= pWork[9];
			uOne |= pWork[10] ^ 0x3F800000U;
			uZero |= pWork[11];
			// 4th row
			uZero |= pWork[12];
			uZero |= pWork[13];
			uZero |= pWork[14];
			uOne |= pWork[15] ^ 0x3F800000U;
			// If all zero entries are zero, the uZero==0
			uZero &= 0x7FFFFFFF;    // Allow -0.0f
			// If all 1.0f entries are 1.0f, then uOne==0
			uOne |= uZero;
			
			return (uOne == 0);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp1 = _mm_cmpeq_ps(m.r[0], g_IdentityR0);
			VECTOR vTemp2 = _mm_cmpeq_ps(m.r[1], g_IdentityR1);
			VECTOR vTemp3 = _mm_cmpeq_ps(m.r[2], g_IdentityR2);
			VECTOR vTemp4 = _mm_cmpeq_ps(m.r[3], g_IdentityR3);
			vTemp1 = _mm_and_ps(vTemp1, vTemp2);
			vTemp3 = _mm_and_ps(vTemp3, vTemp4);
			vTemp1 = _mm_and_ps(vTemp1, vTemp3);
			
			return (_mm_movemask_ps(vTemp1) == 0x0f);
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV Multiply(A_MATRIX M1, B_MATRIX M2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			MATRIX mResult;
			// Cache the invariants in registers
			float x = M1.m[0][0];
			float y = M1.m[0][1];
			float z = M1.m[0][2];
			float w = M1.m[0][3];
			// Perform the operation on the first row
			mResult.m[0][0] = (M2.m[0][0] * x) + (M2.m[1][0] * y) + (M2.m[2][0] * z) + (M2.m[3][0] * w);
			mResult.m[0][1] = (M2.m[0][1] * x) + (M2.m[1][1] * y) + (M2.m[2][1] * z) + (M2.m[3][1] * w);
			mResult.m[0][2] = (M2.m[0][2] * x) + (M2.m[1][2] * y) + (M2.m[2][2] * z) + (M2.m[3][2] * w);
			mResult.m[0][3] = (M2.m[0][3] * x) + (M2.m[1][3] * y) + (M2.m[2][3] * z) + (M2.m[3][3] * w);
			// Repeat for all the other rows
			x = M1.m[1][0];
			y = M1.m[1][1];
			z = M1.m[1][2];
			w = M1.m[1][3];
			mResult.m[1][0] = (M2.m[0][0] * x) + (M2.m[1][0] * y) + (M2.m[2][0] * z) + (M2.m[3][0] * w);
			mResult.m[1][1] = (M2.m[0][1] * x) + (M2.m[1][1] * y) + (M2.m[2][1] * z) + (M2.m[3][1] * w);
			mResult.m[1][2] = (M2.m[0][2] * x) + (M2.m[1][2] * y) + (M2.m[2][2] * z) + (M2.m[3][2] * w);
			mResult.m[1][3] = (M2.m[0][3] * x) + (M2.m[1][3] * y) + (M2.m[2][3] * z) + (M2.m[3][3] * w);
			x = M1.m[2][0];
			y = M1.m[2][1];
			z = M1.m[2][2];
			w = M1.m[2][3];
			mResult.m[2][0] = (M2.m[0][0] * x) + (M2.m[1][0] * y) + (M2.m[2][0] * z) + (M2.m[3][0] * w);
			mResult.m[2][1] = (M2.m[0][1] * x) + (M2.m[1][1] * y) + (M2.m[2][1] * z) + (M2.m[3][1] * w);
			mResult.m[2][2] = (M2.m[0][2] * x) + (M2.m[1][2] * y) + (M2.m[2][2] * z) + (M2.m[3][2] * w);
			mResult.m[2][3] = (M2.m[0][3] * x) + (M2.m[1][3] * y) + (M2.m[2][3] * z) + (M2.m[3][3] * w);
			x = M1.m[3][0];
			y = M1.m[3][1];
			z = M1.m[3][2];
			w = M1.m[3][3];
			mResult.m[3][0] = (M2.m[0][0] * x) + (M2.m[1][0] * y) + (M2.m[2][0] * z) + (M2.m[3][0] * w);
			mResult.m[3][1] = (M2.m[0][1] * x) + (M2.m[1][1] * y) + (M2.m[2][1] * z) + (M2.m[3][1] * w);
			mResult.m[3][2] = (M2.m[0][2] * x) + (M2.m[1][2] * y) + (M2.m[2][2] * z) + (M2.m[3][2] * w);
			mResult.m[3][3] = (M2.m[0][3] * x) + (M2.m[1][3] * y) + (M2.m[2][3] * z) + (M2.m[3][3] * w);
			
			return mResult;

#elif defined(_AVX2_INTRINSICS_)
			__m256 t0 = _mm256_castps128_ps256(M1.r[0]);
			t0 = _mm256_insertf128_ps(t0, M1.r[1], 1);
			__m256 t1 = _mm256_castps128_ps256(M1.r[2]);
			t1 = _mm256_insertf128_ps(t1, M1.r[3], 1);

			__m256 u0 = _mm256_castps128_ps256(M2.r[0]);
			u0 = _mm256_insertf128_ps(u0, M2.r[1], 1);
			__m256 u1 = _mm256_castps128_ps256(M2.r[2]);
			u1 = _mm256_insertf128_ps(u1, M2.r[3], 1);

			__m256 a0 = _mm256_shuffle_ps(t0, t0, _MM_SHUFFLE(0, 0, 0, 0));
			__m256 a1 = _mm256_shuffle_ps(t1, t1, _MM_SHUFFLE(0, 0, 0, 0));
			__m256 b0 = _mm256_permute2f128_ps(u0, u0, 0x00);
			__m256 c0 = _mm256_mul_ps(a0, b0);
			__m256 c1 = _mm256_mul_ps(a1, b0);

			a0 = _mm256_shuffle_ps(t0, t0, _MM_SHUFFLE(1, 1, 1, 1));
			a1 = _mm256_shuffle_ps(t1, t1, _MM_SHUFFLE(1, 1, 1, 1));
			b0 = _mm256_permute2f128_ps(u0, u0, 0x11);
			__m256 c2 = _mm256_fmadd_ps(a0, b0, c0);
			__m256 c3 = _mm256_fmadd_ps(a1, b0, c1);

			a0 = _mm256_shuffle_ps(t0, t0, _MM_SHUFFLE(2, 2, 2, 2));
			a1 = _mm256_shuffle_ps(t1, t1, _MM_SHUFFLE(2, 2, 2, 2));
			__m256 b1 = _mm256_permute2f128_ps(u1, u1, 0x00);
			__m256 c4 = _mm256_mul_ps(a0, b1);
			__m256 c5 = _mm256_mul_ps(a1, b1);

			a0 = _mm256_shuffle_ps(t0, t0, _MM_SHUFFLE(3, 3, 3, 3));
			a1 = _mm256_shuffle_ps(t1, t1, _MM_SHUFFLE(3, 3, 3, 3));
			b1 = _mm256_permute2f128_ps(u1, u1, 0x11);
			__m256 c6 = _mm256_fmadd_ps(a0, b1, c4);
			__m256 c7 = _mm256_fmadd_ps(a1, b1, c5);

			t0 = _mm256_add_ps(c2, c6);
			t1 = _mm256_add_ps(c3, c7);

			MATRIX mResult;
			mResult.r[0] = _mm256_castps256_ps128(t0);
			mResult.r[1] = _mm256_extractf128_ps(t0, 1);
			mResult.r[2] = _mm256_castps256_ps128(t1);
			mResult.r[3] = _mm256_extractf128_ps(t1, 1);
			return mResult;

#elif defined(_SSE2_INTRINSICS_)
    		MATRIX mResult;
    		// Splat the component X,Y,Z then W
#if defined(_AVX_INTRINSICS_)
			VECTOR vX = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[0]) + 0);
			VECTOR vY = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[0]) + 1);
			VECTOR vZ = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[0]) + 2);
			VECTOR vW = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[0]) + 3);
#else
			// Use vW to hold the original row
			VECTOR vW = M1.r[0];
			VECTOR vX = PERMUTE_PS(vW, _MM_SHUFFLE(0, 0, 0, 0));
			VECTOR vY = PERMUTE_PS(vW, _MM_SHUFFLE(1, 1, 1, 1));
			VECTOR vZ = PERMUTE_PS(vW, _MM_SHUFFLE(2, 2, 2, 2));
			vW = PERMUTE_PS(vW, _MM_SHUFFLE(3, 3, 3, 3));
#endif
			// Perform the operation on the first row
			vX = _mm_mul_ps(vX, M2.r[0]);
			vY = _mm_mul_ps(vY, M2.r[1]);
			vZ = _mm_mul_ps(vZ, M2.r[2]);
			vW = _mm_mul_ps(vW, M2.r[3]);
			// Perform a binary add to reduce cumulative errors
			vX = _mm_add_ps(vX, vZ);
			vY = _mm_add_ps(vY, vW);
			vX = _mm_add_ps(vX, vY);
			mResult.r[0] = vX;
			// Repeat for the other 3 rows
#if defined(_AVX_INTRINSICS_)
			vX = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[1]) + 0);
			vY = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[1]) + 1);
			vZ = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[1]) + 2);
			vW = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[1]) + 3);
#else
			vW = M1.r[1];
			vX = PERMUTE_PS(vW, _MM_SHUFFLE(0, 0, 0, 0));
			vY = PERMUTE_PS(vW, _MM_SHUFFLE(1, 1, 1, 1));
			vZ = PERMUTE_PS(vW, _MM_SHUFFLE(2, 2, 2, 2));
			vW = PERMUTE_PS(vW, _MM_SHUFFLE(3, 3, 3, 3));
#endif
    		vX = _mm_mul_ps(vX, M2.r[0]);
			vY = _mm_mul_ps(vY, M2.r[1]);
			vZ = _mm_mul_ps(vZ, M2.r[2]);
			vW = _mm_mul_ps(vW, M2.r[3]);
			vX = _mm_add_ps(vX, vZ);
			vY = _mm_add_ps(vY, vW);
			vX = _mm_add_ps(vX, vY);
			mResult.r[1] = vX;
#if defined(_AVX_INTRINSICS_)
			vX = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[2]) + 0);
			vY = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[2]) + 1);
			vZ = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[2]) + 2);
			vW = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[2]) + 3);
#else
			vW = M1.r[2];
			vX = PERMUTE_PS(vW, _MM_SHUFFLE(0, 0, 0, 0));
			vY = PERMUTE_PS(vW, _MM_SHUFFLE(1, 1, 1, 1));
			vZ = PERMUTE_PS(vW, _MM_SHUFFLE(2, 2, 2, 2));
			vW = PERMUTE_PS(vW, _MM_SHUFFLE(3, 3, 3, 3));
#endif
			vX = _mm_mul_ps(vX, M2.r[0]);
			vY = _mm_mul_ps(vY, M2.r[1]);
			vZ = _mm_mul_ps(vZ, M2.r[2]);
			vW = _mm_mul_ps(vW, M2.r[3]);
			vX = _mm_add_ps(vX, vZ);
			vY = _mm_add_ps(vY, vW);
			vX = _mm_add_ps(vX, vY);
			mResult.r[2] = vX;
#if defined(_AVX_INTRINSICS_)
			vX = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[3]) + 0);
			vY = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[3]) + 1);
			vZ = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[3]) + 2);
			vW = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[3]) + 3);
#else
			vW = M1.r[3];
			vX = PERMUTE_PS(vW, _MM_SHUFFLE(0, 0, 0, 0));
			vY = PERMUTE_PS(vW, _MM_SHUFFLE(1, 1, 1, 1));
			vZ = PERMUTE_PS(vW, _MM_SHUFFLE(2, 2, 2, 2));
			vW = PERMUTE_PS(vW, _MM_SHUFFLE(3, 3, 3, 3));
#endif
			vX = _mm_mul_ps(vX, M2.r[0]);
			vY = _mm_mul_ps(vY, M2.r[1]);
			vZ = _mm_mul_ps(vZ, M2.r[2]);
			vW = _mm_mul_ps(vW, M2.r[3]);
			vX = _mm_add_ps(vX, vZ);
			vY = _mm_add_ps(vY, vW);
			vX = _mm_add_ps(vX, vY);
			mResult.r[3] = vX;
			
			return mResult;
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV MultiplyTranspose(A_MATRIX M1, B_MATRIX M2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			MATRIX mResult;
			// Cache the invariants in registers
			float x = M2.m[0][0];
			float y = M2.m[1][0];
			float z = M2.m[2][0];
			float w = M2.m[3][0];
			// Perform the operation on the first row
			mResult.m[0][0] = (M1.m[0][0] * x) + (M1.m[0][1] * y) + (M1.m[0][2] * z) + (M1.m[0][3] * w);
			mResult.m[0][1] = (M1.m[1][0] * x) + (M1.m[1][1] * y) + (M1.m[1][2] * z) + (M1.m[1][3] * w);
			mResult.m[0][2] = (M1.m[2][0] * x) + (M1.m[2][1] * y) + (M1.m[2][2] * z) + (M1.m[2][3] * w);
			mResult.m[0][3] = (M1.m[3][0] * x) + (M1.m[3][1] * y) + (M1.m[3][2] * z) + (M1.m[3][3] * w);
			// Repeat for all the other rows
			x = M2.m[0][1];
			y = M2.m[1][1];
			z = M2.m[2][1];
			w = M2.m[3][1];
			mResult.m[1][0] = (M1.m[0][0] * x) + (M1.m[0][1] * y) + (M1.m[0][2] * z) + (M1.m[0][3] * w);
			mResult.m[1][1] = (M1.m[1][0] * x) + (M1.m[1][1] * y) + (M1.m[1][2] * z) + (M1.m[1][3] * w);
			mResult.m[1][2] = (M1.m[2][0] * x) + (M1.m[2][1] * y) + (M1.m[2][2] * z) + (M1.m[2][3] * w);
			mResult.m[1][3] = (M1.m[3][0] * x) + (M1.m[3][1] * y) + (M1.m[3][2] * z) + (M1.m[3][3] * w);
			x = M2.m[0][2];
			y = M2.m[1][2];
			z = M2.m[2][2];
			w = M2.m[3][2];
			mResult.m[2][0] = (M1.m[0][0] * x) + (M1.m[0][1] * y) + (M1.m[0][2] * z) + (M1.m[0][3] * w);
			mResult.m[2][1] = (M1.m[1][0] * x) + (M1.m[1][1] * y) + (M1.m[1][2] * z) + (M1.m[1][3] * w);
			mResult.m[2][2] = (M1.m[2][0] * x) + (M1.m[2][1] * y) + (M1.m[2][2] * z) + (M1.m[2][3] * w);
			mResult.m[2][3] = (M1.m[3][0] * x) + (M1.m[3][1] * y) + (M1.m[3][2] * z) + (M1.m[3][3] * w);
			x = M2.m[0][3];
			y = M2.m[1][3];
			z = M2.m[2][3];
			w = M2.m[3][3];
			mResult.m[3][0] = (M1.m[0][0] * x) + (M1.m[0][1] * y) + (M1.m[0][2] * z) + (M1.m[0][3] * w);
			mResult.m[3][1] = (M1.m[1][0] * x) + (M1.m[1][1] * y) + (M1.m[1][2] * z) + (M1.m[1][3] * w);
			mResult.m[3][2] = (M1.m[2][0] * x) + (M1.m[2][1] * y) + (M1.m[2][2] * z) + (M1.m[2][3] * w);
			mResult.m[3][3] = (M1.m[3][0] * x) + (M1.m[3][1] * y) + (M1.m[3][2] * z) + (M1.m[3][3] * w);
			
			return mResult;

#elif defined(_AVX2_INTRINSICS_)
			__m256 t0 = _mm256_castps128_ps256(M1.r[0]);
			t0 = _mm256_insertf128_ps(t0, M1.r[1], 1);
			__m256 t1 = _mm256_castps128_ps256(M1.r[2]);
			t1 = _mm256_insertf128_ps(t1, M1.r[3], 1);

			__m256 u0 = _mm256_castps128_ps256(M2.r[0]);
			u0 = _mm256_insertf128_ps(u0, M2.r[1], 1);
			__m256 u1 = _mm256_castps128_ps256(M2.r[2]);
			u1 = _mm256_insertf128_ps(u1, M2.r[3], 1);

			__m256 a0 = _mm256_shuffle_ps(t0, t0, _MM_SHUFFLE(0, 0, 0, 0));
			__m256 a1 = _mm256_shuffle_ps(t1, t1, _MM_SHUFFLE(0, 0, 0, 0));
			__m256 b0 = _mm256_permute2f128_ps(u0, u0, 0x00);
			__m256 c0 = _mm256_mul_ps(a0, b0);
			__m256 c1 = _mm256_mul_ps(a1, b0);

			a0 = _mm256_shuffle_ps(t0, t0, _MM_SHUFFLE(1, 1, 1, 1));
			a1 = _mm256_shuffle_ps(t1, t1, _MM_SHUFFLE(1, 1, 1, 1));
			b0 = _mm256_permute2f128_ps(u0, u0, 0x11);
			__m256 c2 = _mm256_fmadd_ps(a0, b0, c0);
			__m256 c3 = _mm256_fmadd_ps(a1, b0, c1);

			a0 = _mm256_shuffle_ps(t0, t0, _MM_SHUFFLE(2, 2, 2, 2));
			a1 = _mm256_shuffle_ps(t1, t1, _MM_SHUFFLE(2, 2, 2, 2));
			__m256 b1 = _mm256_permute2f128_ps(u1, u1, 0x00);
			__m256 c4 = _mm256_mul_ps(a0, b1);
			__m256 c5 = _mm256_mul_ps(a1, b1);

			a0 = _mm256_shuffle_ps(t0, t0, _MM_SHUFFLE(3, 3, 3, 3));
			a1 = _mm256_shuffle_ps(t1, t1, _MM_SHUFFLE(3, 3, 3, 3));
			b1 = _mm256_permute2f128_ps(u1, u1, 0x11);
			__m256 c6 = _mm256_fmadd_ps(a0, b1, c4);
			__m256 c7 = _mm256_fmadd_ps(a1, b1, c5);

			t0 = _mm256_add_ps(c2, c6);
			t1 = _mm256_add_ps(c3, c7);

			// Transpose result
			__m256 vTemp = _mm256_unpacklo_ps(t0, t1);
			__m256 vTemp2 = _mm256_unpackhi_ps(t0, t1);
			__m256 vTemp3 = _mm256_permute2f128_ps(vTemp, vTemp2, 0x20);
			__m256 vTemp4 = _mm256_permute2f128_ps(vTemp, vTemp2, 0x31);
			vTemp = _mm256_unpacklo_ps(vTemp3, vTemp4);
			vTemp2 = _mm256_unpackhi_ps(vTemp3, vTemp4);
			t0 = _mm256_permute2f128_ps(vTemp, vTemp2, 0x20);
			t1 = _mm256_permute2f128_ps(vTemp, vTemp2, 0x31);

			MATRIX mResult;
			mResult.r[0] = _mm256_castps256_ps128(t0);
			mResult.r[1] = _mm256_extractf128_ps(t0, 1);
			mResult.r[2] = _mm256_castps256_ps128(t1);
			mResult.r[3] = _mm256_extractf128_ps(t1, 1);
			
			return mResult;

#elif defined(_SSE2_INTRINSICS_)
    		// Splat the component X,Y,Z then W
#if defined(_AVX_INTRINSICS_)
			VECTOR vX = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[0]) + 0);
			VECTOR vY = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[0]) + 1);
			VECTOR vZ = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[0]) + 2);
			VECTOR vW = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[0]) + 3);
#else
			// Use vW to hold the original row
			VECTOR vW = M1.r[0];
			VECTOR vX = PERMUTE_PS(vW, _MM_SHUFFLE(0, 0, 0, 0));
			VECTOR vY = PERMUTE_PS(vW, _MM_SHUFFLE(1, 1, 1, 1));
			VECTOR vZ = PERMUTE_PS(vW, _MM_SHUFFLE(2, 2, 2, 2));
			vW = PERMUTE_PS(vW, _MM_SHUFFLE(3, 3, 3, 3));
#endif
			// Perform the operation on the first row
			vX = _mm_mul_ps(vX, M2.r[0]);
			vY = _mm_mul_ps(vY, M2.r[1]);
			vZ = _mm_mul_ps(vZ, M2.r[2]);
			vW = _mm_mul_ps(vW, M2.r[3]);
			// Perform a binary add to reduce cumulative errors
			vX = _mm_add_ps(vX, vZ);
			vY = _mm_add_ps(vY, vW);
			vX = _mm_add_ps(vX, vY);
			VECTOR r0 = vX;
			// Repeat for the other 3 rows
#if defined(_AVX_INTRINSICS_)
			vX = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[1]) + 0);
			vY = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[1]) + 1);
			vZ = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[1]) + 2);
			vW = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[1]) + 3);
#else
			vW = M1.r[1];
			vX = PERMUTE_PS(vW, _MM_SHUFFLE(0, 0, 0, 0));
			vY = PERMUTE_PS(vW, _MM_SHUFFLE(1, 1, 1, 1));
			vZ = PERMUTE_PS(vW, _MM_SHUFFLE(2, 2, 2, 2));
			vW = PERMUTE_PS(vW, _MM_SHUFFLE(3, 3, 3, 3));
#endif
			vX = _mm_mul_ps(vX, M2.r[0]);
			vY = _mm_mul_ps(vY, M2.r[1]);
			vZ = _mm_mul_ps(vZ, M2.r[2]);
			vW = _mm_mul_ps(vW, M2.r[3]);
			vX = _mm_add_ps(vX, vZ);
			vY = _mm_add_ps(vY, vW);
			vX = _mm_add_ps(vX, vY);
			VECTOR r1 = vX;
#if defined(_AVX_INTRINSICS_)
			vX = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[2]) + 0);
			vY = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[2]) + 1);
			vZ = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[2]) + 2);
			vW = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[2]) + 3);
#else
			vW = M1.r[2];
			vX = PERMUTE_PS(vW, _MM_SHUFFLE(0, 0, 0, 0));
			vY = PERMUTE_PS(vW, _MM_SHUFFLE(1, 1, 1, 1));
			vZ = PERMUTE_PS(vW, _MM_SHUFFLE(2, 2, 2, 2));
			vW = PERMUTE_PS(vW, _MM_SHUFFLE(3, 3, 3, 3));
#endif
			vX = _mm_mul_ps(vX, M2.r[0]);
			vY = _mm_mul_ps(vY, M2.r[1]);
			vZ = _mm_mul_ps(vZ, M2.r[2]);
			vW = _mm_mul_ps(vW, M2.r[3]);
			vX = _mm_add_ps(vX, vZ);
			vY = _mm_add_ps(vY, vW);
			vX = _mm_add_ps(vX, vY);
			VECTOR r2 = vX;
#if defined(_AVX_INTRINSICS_)
			vX = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[3]) + 0);
			vY = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[3]) + 1);
			vZ = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[3]) + 2);
			vW = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[3]) + 3);
#else
			vW = M1.r[3];
			vX = PERMUTE_PS(vW, _MM_SHUFFLE(0, 0, 0, 0));
			vY = PERMUTE_PS(vW, _MM_SHUFFLE(1, 1, 1, 1));
			vZ = PERMUTE_PS(vW, _MM_SHUFFLE(2, 2, 2, 2));
			vW = PERMUTE_PS(vW, _MM_SHUFFLE(3, 3, 3, 3));
#endif
			vX = _mm_mul_ps(vX, M2.r[0]);
			vY = _mm_mul_ps(vY, M2.r[1]);
			vZ = _mm_mul_ps(vZ, M2.r[2]);
			vW = _mm_mul_ps(vW, M2.r[3]);
			vX = _mm_add_ps(vX, vZ);
			vY = _mm_add_ps(vY, vW);
			vX = _mm_add_ps(vX, vY);
			VECTOR r3 = vX;

			// Transpose result
			// x.x,x.y,y.x,y.y
			VECTOR vTemp1 = _mm_shuffle_ps(r0, r1, _MM_SHUFFLE(1, 0, 1, 0));
			// x.z,x.w,y.z,y.w
			VECTOR vTemp3 = _mm_shuffle_ps(r0, r1, _MM_SHUFFLE(3, 2, 3, 2));
			// z.x,z.y,w.x,w.y
			VECTOR vTemp2 = _mm_shuffle_ps(r2, r3, _MM_SHUFFLE(1, 0, 1, 0));
			// z.z,z.w,w.z,w.w
			VECTOR vTemp4 = _mm_shuffle_ps(r2, r3, _MM_SHUFFLE(3, 2, 3, 2));

			MATRIX mResult;
			// x.x,y.x,z.x,w.x
			mResult.r[0] = _mm_shuffle_ps(vTemp1, vTemp2, _MM_SHUFFLE(2, 0, 2, 0));
			// x.y,y.y,z.y,w.y
			mResult.r[1] = _mm_shuffle_ps(vTemp1, vTemp2, _MM_SHUFFLE(3, 1, 3, 1));
			// x.z,y.z,z.z,w.z
			mResult.r[2] = _mm_shuffle_ps(vTemp3, vTemp4, _MM_SHUFFLE(2, 0, 2, 0));
			// x.w,y.w,z.w,w.w
			mResult.r[3] = _mm_shuffle_ps(vTemp3, vTemp4, _MM_SHUFFLE(3, 1, 3, 1));
			
			return mResult;
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV Transpose(A_MATRIX m) noexcept
		{
#if defined(_NO_INTRINSICS_)
			// Original matrix:
			//
			//     m00m01m02m03
			//     m10m11m12m13
			//     m20m21m22m23
			//     m30m31m32m33

			MATRIX P;
			P.r[0] = VEC::MergeXY(m.r[0], m.r[2]); // m00m20m01m21
			P.r[1] = VEC::MergeXY(m.r[1], m.r[3]); // m10m30m11m31
			P.r[2] = VEC::MergeZW(m.r[0], m.r[2]); // m02m22m03m23
			P.r[3] = VEC::MergeZW(m.r[1], m.r[3]); // m12m32m13m33

			MATRIX MT;
			MT.r[0] = VectorMergeXY(P.r[0], P.r[1]); // m00m10m20m30
			MT.r[1] = VectorMergeZW(P.r[0], P.r[1]); // m01m11m21m31
			MT.r[2] = VectorMergeXY(P.r[2], P.r[3]); // m02m12m22m32
			MT.r[3] = VectorMergeZW(P.r[2], P.r[3]); // m03m13m23m33
			
			return MT;

#elif defined(_AVX2_INTRINSICS_)
			__m256 t0 = _mm256_castps128_ps256(m.r[0]);
			t0 = _mm256_insertf128_ps(t0, m.r[1], 1);
			__m256 t1 = _mm256_castps128_ps256(m.r[2]);
			t1 = _mm256_insertf128_ps(t1, m.r[3], 1);

			__m256 vTemp = _mm256_unpacklo_ps(t0, t1);
			__m256 vTemp2 = _mm256_unpackhi_ps(t0, t1);
			__m256 vTemp3 = _mm256_permute2f128_ps(vTemp, vTemp2, 0x20);
			__m256 vTemp4 = _mm256_permute2f128_ps(vTemp, vTemp2, 0x31);
			vTemp = _mm256_unpacklo_ps(vTemp3, vTemp4);
			vTemp2 = _mm256_unpackhi_ps(vTemp3, vTemp4);
			t0 = _mm256_permute2f128_ps(vTemp, vTemp2, 0x20);
			t1 = _mm256_permute2f128_ps(vTemp, vTemp2, 0x31);

			MATRIX mResult;
			mResult.r[0] = _mm256_castps256_ps128(t0);
			mResult.r[1] = _mm256_extractf128_ps(t0, 1);
			mResult.r[2] = _mm256_castps256_ps128(t1);
			mResult.r[3] = _mm256_extractf128_ps(t1, 1);
			
			return mResult;

#elif defined(_SSE2_INTRINSICS_)
			// x.x,x.y,y.x,y.y
			VECTOR vTemp1 = _mm_shuffle_ps(m.r[0], m.r[1], _MM_SHUFFLE(1, 0, 1, 0));
			// x.z,x.w,y.z,y.w
			VECTOR vTemp3 = _mm_shuffle_ps(m.r[0], m.r[1], _MM_SHUFFLE(3, 2, 3, 2));
			// z.x,z.y,w.x,w.y
			VECTOR vTemp2 = _mm_shuffle_ps(m.r[2], m.r[3], _MM_SHUFFLE(1, 0, 1, 0));
			// z.z,z.w,w.z,w.w
			VECTOR vTemp4 = _mm_shuffle_ps(m.r[2], m.r[3], _MM_SHUFFLE(3, 2, 3, 2));

			MATRIX mResult;
			// x.x,y.x,z.x,w.x
			mResult.r[0] = _mm_shuffle_ps(vTemp1, vTemp2, _MM_SHUFFLE(2, 0, 2, 0));
			// x.y,y.y,z.y,w.y
			mResult.r[1] = _mm_shuffle_ps(vTemp1, vTemp2, _MM_SHUFFLE(3, 1, 3, 1));
			// x.z,y.z,z.z,w.z
			mResult.r[2] = _mm_shuffle_ps(vTemp3, vTemp4, _MM_SHUFFLE(2, 0, 2, 0));
			// x.w,y.w,z.w,w.w
			mResult.r[3] = _mm_shuffle_ps(vTemp3, vTemp4, _MM_SHUFFLE(3, 1, 3, 1));
			
			return mResult;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE MATRIX VEC_CALLCONV Inverse(VECTOR* pDeterminant, A_MATRIX m) noexcept
		{
#if defined(_NO_INTRINSICS_)
			MATRIX MT = Transpose(m);

			VECTOR V0[4], V1[4];
			V0[0] = VEC::Swizzle<SWIZZLE_X, SWIZZLE_X, SWIZZLE_Y, SWIZZLE_Y>(MT.r[2]);
			V1[0] = VEC::Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_Z, SWIZZLE_W>(MT.r[3]);
			V0[1] = VEC::Swizzle<SWIZZLE_X, SWIZZLE_X, SWIZZLE_Y, SWIZZLE_Y>(MT.r[0]);
			V1[1] = VEC::Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_Z, SWIZZLE_W>(MT.r[1]);
			V0[2] = VEC::Permute<PERMUTE_0X, PERMUTE_0Z, PERMUTE_1X, PERMUTE_1Z>(MT.r[2], MT.r[0]);
			V1[2] = VEC::Permute<PERMUTE_0Y, PERMUTE_0W, PERMUTE_1Y, PERMUTE_1W>(MT.r[3], MT.r[1]);

			VECTOR D0 = VEC::Multiply(V0[0], V1[0]);
			VECTOR D1 = VEC::Multiply(V0[1], V1[1]);
			VECTOR D2 = VEC::Multiply(V0[2], V1[2]);

			V0[0] = VEC::Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_Z, SWIZZLE_W>(MT.r[2]);
			V1[0] = VEC::Swizzle<SWIZZLE_X, SWIZZLE_X, SWIZZLE_Y, SWIZZLE_Y>(MT.r[3]);
			V0[1] = VEC::Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_Z, SWIZZLE_W>(MT.r[0]);
			V1[1] = VEC::Swizzle<SWIZZLE_X, SWIZZLE_X, SWIZZLE_Y, SWIZZLE_Y>(MT.r[1]);
			V0[2] = VEC::Permute<PERMUTE_0Y, PERMUTE_0W, PERMUTE_1Y, PERMUTE_1W>(MT.r[2], MT.r[0]);
			V1[2] = VEC::Permute<PERMUTE_0X, PERMUTE_0Z, PERMUTE_1X, PERMUTE_1Z>(MT.r[3], MT.r[1]);

			D0 = VEC::NegativeMultiplySubtract(V0[0], V1[0], D0);
			D1 = VEC::NegativeMultiplySubtract(V0[1], V1[1], D1);
			D2 = VEC::NegativeMultiplySubtract(V0[2], V1[2], D2);

			V0[0] = VEC::Swizzle<SWIZZLE_Y, SWIZZLE_Z, SWIZZLE_X, SWIZZLE_Y>(MT.r[1]);
			V1[0] = VEC::Permute<PERMUTE_1Y, PERMUTE_0Y, PERMUTE_0W, PERMUTE_0X>(D0, D2);
			V0[1] = VEC::Swizzle<SWIZZLE_Z, SWIZZLE_X, SWIZZLE_Y, SWIZZLE_X>(MT.r[0]);
			V1[1] = VEC::Permute<PERMUTE_0W, PERMUTE_1Y, PERMUTE_0Y, PERMUTE_0Z>(D0, D2);
			V0[2] = VEC::Swizzle<SWIZZLE_Y, SWIZZLE_Z, SWIZZLE_X, SWIZZLE_Y>(MT.r[3]);
			V1[2] = VEC::Permute<PERMUTE_1W, PERMUTE_0Y, PERMUTE_0W, PERMUTE_0X>(D1, D2);
			V0[3] = VEC::Swizzle<SWIZZLE_Z, SWIZZLE_X, SWIZZLE_Y, SWIZZLE_X>(MT.r[2]);
			V1[3] = VEC::Permute<PERMUTE_0W, PERMUTE_1W, PERMUTE_0Y, PERMUTE_0Z>(D1, D2);

			VECTOR C0 = VEC::Multiply(V0[0], V1[0]);
			VECTOR C2 = VEC::Multiply(V0[1], V1[1]);
			VECTOR C4 = VEC::Multiply(V0[2], V1[2]);
			VECTOR C6 = VEC::Multiply(V0[3], V1[3]);

			V0[0] = VEC::Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_Y, SWIZZLE_Z>(MT.r[1]);
			V1[0] = VEC::Permute<PERMUTE_0W, PERMUTE_0X, PERMUTE_0Y, PERMUTE_1X>(D0, D2);
			V0[1] = VEC::Swizzle<SWIZZLE_W, SWIZZLE_Z, SWIZZLE_W, SWIZZLE_Y>(MT.r[0]);
			V1[1] = VEC::Permute<PERMUTE_0Z, PERMUTE_0Y, PERMUTE_1X, PERMUTE_0X>(D0, D2);
			V0[2] = VEC::Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_Y, SWIZZLE_Z>(MT.r[3]);
			V1[2] = VEC::Permute<PERMUTE_0W, PERMUTE_0X, PERMUTE_0Y, PERMUTE_1Z>(D1, D2);
			V0[3] = VEC::Swizzle<SWIZZLE_W, SWIZZLE_Z, SWIZZLE_W, SWIZZLE_Y>(MT.r[2]);
			V1[3] = VEC::Permute<PERMUTE_0Z, PERMUTE_0Y, PERMUTE_1Z, PERMUTE_0X>(D1, D2);

			C0 = VEC::NegativeMultiplySubtract(V0[0], V1[0], C0);
			C2 = VEC::NegativeMultiplySubtract(V0[1], V1[1], C2);
			C4 = VEC::NegativeMultiplySubtract(V0[2], V1[2], C4);
			C6 = VEC::NegativeMultiplySubtract(V0[3], V1[3], C6);

			V0[0] = VEC::Swizzle<SWIZZLE_W, SWIZZLE_X, SWIZZLE_W, SWIZZLE_X>(MT.r[1]);
			V1[0] = VEC::Permute<PERMUTE_0Z, PERMUTE_1Y, PERMUTE_1X, PERMUTE_0Z>(D0, D2);
			V0[1] = VEC::Swizzle<SWIZZLE_Y, SWIZZLE_W, SWIZZLE_X, SWIZZLE_Z>(MT.r[0]);
			V1[1] = VEC::Permute<PERMUTE_1Y, PERMUTE_0X, PERMUTE_0W, PERMUTE_1X>(D0, D2);
			V0[2] = VEC::Swizzle<SWIZZLE_W, SWIZZLE_X, SWIZZLE_W, SWIZZLE_X>(MT.r[3]);
			V1[2] = VEC::Permute<PERMUTE_0Z, PERMUTE_1W, PERMUTE_1Z, PERMUTE_0Z>(D1, D2);
			V0[3] = VEC::Swizzle<SWIZZLE_Y, SWIZZLE_W, SWIZZLE_X, SWIZZLE_Z>(MT.r[2]);
			V1[3] = VEC::Permute<PERMUTE_1W, PERMUTE_0X, PERMUTE_0W, PERMUTE_1Z>(D1, D2);

			VECTOR C1 = VEC::NegativeMultiplySubtract(V0[0], V1[0], C0);
			C0 = VEC::MultiplyAdd(V0[0], V1[0], C0);
			VECTOR C3 = VEC::MultiplyAdd(V0[1], V1[1], C2);
			C2 = VEC::NegativeMultiplySubtract(V0[1], V1[1], C2);
			VECTOR C5 = VEC::NegativeMultiplySubtract(V0[2], V1[2], C4);
			C4 = VEC::MultiplyAdd(V0[2], V1[2], C4);
			VECTOR C7 = VEC::MultiplyAdd(V0[3], V1[3], C6);
			C6 = VEC::NegativeMultiplySubtract(V0[3], V1[3], C6);

			MATRIX R;
			R.r[0] = VEC::Select(C0, C1, g_Select0101.v);
			R.r[1] = VEC::Select(C2, C3, g_Select0101.v);
			R.r[2] = VEC::Select(C4, C5, g_Select0101.v);
			R.r[3] = VEC::Select(C6, C7, g_Select0101.v);

			VECTOR determinant = VEC4::Dot(R.r[0], MT.r[0]);

			if (pDeterminant != nullptr)
				*pDeterminant = determinant;

			VECTOR reciprocal = VEC::Reciprocal(determinant);

			MATRIX Result;
			Result.r[0] = VEC::Multiply(R.r[0], reciprocal);
			Result.r[1] = VEC::Multiply(R.r[1], reciprocal);
			Result.r[2] = VEC::Multiply(R.r[2], reciprocal);
			Result.r[3] = VEC::Multiply(R.r[3], reciprocal);
			return Result;

#elif defined(_SSE2_INTRINSICS_)
			// Transpose matrix
			VECTOR vTemp1 = _mm_shuffle_ps(m.r[0], m.r[1], _MM_SHUFFLE(1, 0, 1, 0));
			VECTOR vTemp3 = _mm_shuffle_ps(m.r[0], m.r[1], _MM_SHUFFLE(3, 2, 3, 2));
			VECTOR vTemp2 = _mm_shuffle_ps(m.r[2], m.r[3], _MM_SHUFFLE(1, 0, 1, 0));
			VECTOR vTemp4 = _mm_shuffle_ps(m.r[2], m.r[3], _MM_SHUFFLE(3, 2, 3, 2));

			MATRIX MT;
			MT.r[0] = _mm_shuffle_ps(vTemp1, vTemp2, _MM_SHUFFLE(2, 0, 2, 0));
			MT.r[1] = _mm_shuffle_ps(vTemp1, vTemp2, _MM_SHUFFLE(3, 1, 3, 1));
			MT.r[2] = _mm_shuffle_ps(vTemp3, vTemp4, _MM_SHUFFLE(2, 0, 2, 0));
			MT.r[3] = _mm_shuffle_ps(vTemp3, vTemp4, _MM_SHUFFLE(3, 1, 3, 1));

			VECTOR V00 = PERMUTE_PS(MT.r[2], _MM_SHUFFLE(1, 1, 0, 0));
			VECTOR V10 = PERMUTE_PS(MT.r[3], _MM_SHUFFLE(3, 2, 3, 2));
			VECTOR V01 = PERMUTE_PS(MT.r[0], _MM_SHUFFLE(1, 1, 0, 0));
			VECTOR V11 = PERMUTE_PS(MT.r[1], _MM_SHUFFLE(3, 2, 3, 2));
			VECTOR V02 = _mm_shuffle_ps(MT.r[2], MT.r[0], _MM_SHUFFLE(2, 0, 2, 0));
			VECTOR V12 = _mm_shuffle_ps(MT.r[3], MT.r[1], _MM_SHUFFLE(3, 1, 3, 1));

			VECTOR D0 = _mm_mul_ps(V00, V10);
			VECTOR D1 = _mm_mul_ps(V01, V11);
			VECTOR D2 = _mm_mul_ps(V02, V12);

			V00 = PERMUTE_PS(MT.r[2], _MM_SHUFFLE(3, 2, 3, 2));
			V10 = PERMUTE_PS(MT.r[3], _MM_SHUFFLE(1, 1, 0, 0));
			V01 = PERMUTE_PS(MT.r[0], _MM_SHUFFLE(3, 2, 3, 2));
			V11 = PERMUTE_PS(MT.r[1], _MM_SHUFFLE(1, 1, 0, 0));
			V02 = _mm_shuffle_ps(MT.r[2], MT.r[0], _MM_SHUFFLE(3, 1, 3, 1));
			V12 = _mm_shuffle_ps(MT.r[3], MT.r[1], _MM_SHUFFLE(2, 0, 2, 0));

			D0 = FNMADD_PS(V00, V10, D0);
			D1 = FNMADD_PS(V01, V11, D1);
			D2 = FNMADD_PS(V02, V12, D2);
			// V11 = D0Y,D0W,D2Y,D2Y
			V11 = _mm_shuffle_ps(D0, D2, _MM_SHUFFLE(1, 1, 3, 1));
			V00 = PERMUTE_PS(MT.r[1], _MM_SHUFFLE(1, 0, 2, 1));
			V10 = _mm_shuffle_ps(V11, D0, _MM_SHUFFLE(0, 3, 0, 2));
			V01 = PERMUTE_PS(MT.r[0], _MM_SHUFFLE(0, 1, 0, 2));
			V11 = _mm_shuffle_ps(V11, D0, _MM_SHUFFLE(2, 1, 2, 1));
			// V13 = D1Y,D1W,D2W,D2W
			VECTOR V13 = _mm_shuffle_ps(D1, D2, _MM_SHUFFLE(3, 3, 3, 1));
			V02 = PERMUTE_PS(MT.r[3], _MM_SHUFFLE(1, 0, 2, 1));
			V12 = _mm_shuffle_ps(V13, D1, _MM_SHUFFLE(0, 3, 0, 2));
			VECTOR V03 = PERMUTE_PS(MT.r[2], _MM_SHUFFLE(0, 1, 0, 2));
			V13 = _mm_shuffle_ps(V13, D1, _MM_SHUFFLE(2, 1, 2, 1));

			VECTOR C0 = _mm_mul_ps(V00, V10);
			VECTOR C2 = _mm_mul_ps(V01, V11);
			VECTOR C4 = _mm_mul_ps(V02, V12);
			VECTOR C6 = _mm_mul_ps(V03, V13);

			// V11 = D0X,D0Y,D2X,D2X
			V11 = _mm_shuffle_ps(D0, D2, _MM_SHUFFLE(0, 0, 1, 0));
			V00 = PERMUTE_PS(MT.r[1], _MM_SHUFFLE(2, 1, 3, 2));
			V10 = _mm_shuffle_ps(D0, V11, _MM_SHUFFLE(2, 1, 0, 3));
			V01 = PERMUTE_PS(MT.r[0], _MM_SHUFFLE(1, 3, 2, 3));
			V11 = _mm_shuffle_ps(D0, V11, _MM_SHUFFLE(0, 2, 1, 2));
			// V13 = D1X,D1Y,D2Z,D2Z
			V13 = _mm_shuffle_ps(D1, D2, _MM_SHUFFLE(2, 2, 1, 0));
			V02 = PERMUTE_PS(MT.r[3], _MM_SHUFFLE(2, 1, 3, 2));
			V12 = _mm_shuffle_ps(D1, V13, _MM_SHUFFLE(2, 1, 0, 3));
			V03 = PERMUTE_PS(MT.r[2], _MM_SHUFFLE(1, 3, 2, 3));
			V13 = _mm_shuffle_ps(D1, V13, _MM_SHUFFLE(0, 2, 1, 2));

			C0 = FNMADD_PS(V00, V10, C0);
			C2 = FNMADD_PS(V01, V11, C2);
			C4 = FNMADD_PS(V02, V12, C4);
			C6 = FNMADD_PS(V03, V13, C6);

			V00 = PERMUTE_PS(MT.r[1], _MM_SHUFFLE(0, 3, 0, 3));
			// V10 = D0Z,D0Z,D2X,D2Y
			V10 = _mm_shuffle_ps(D0, D2, _MM_SHUFFLE(1, 0, 2, 2));
			V10 = PERMUTE_PS(V10, _MM_SHUFFLE(0, 2, 3, 0));
			V01 = PERMUTE_PS(MT.r[0], _MM_SHUFFLE(2, 0, 3, 1));
			// V11 = D0X,D0W,D2X,D2Y
			V11 = _mm_shuffle_ps(D0, D2, _MM_SHUFFLE(1, 0, 3, 0));
			V11 = PERMUTE_PS(V11, _MM_SHUFFLE(2, 1, 0, 3));
			V02 = PERMUTE_PS(MT.r[3], _MM_SHUFFLE(0, 3, 0, 3));
			// V12 = D1Z,D1Z,D2Z,D2W
			V12 = _mm_shuffle_ps(D1, D2, _MM_SHUFFLE(3, 2, 2, 2));
			V12 = PERMUTE_PS(V12, _MM_SHUFFLE(0, 2, 3, 0));
			V03 = PERMUTE_PS(MT.r[2], _MM_SHUFFLE(2, 0, 3, 1));
			// V13 = D1X,D1W,D2Z,D2W
			V13 = _mm_shuffle_ps(D1, D2, _MM_SHUFFLE(3, 2, 3, 0));
			V13 = PERMUTE_PS(V13, _MM_SHUFFLE(2, 1, 0, 3));

			V00 = _mm_mul_ps(V00, V10);
			V01 = _mm_mul_ps(V01, V11);
			V02 = _mm_mul_ps(V02, V12);
			V03 = _mm_mul_ps(V03, V13);
			VECTOR C1 = _mm_sub_ps(C0, V00);
			C0 = _mm_add_ps(C0, V00);
			VECTOR C3 = _mm_add_ps(C2, V01);
			C2 = _mm_sub_ps(C2, V01);
			VECTOR C5 = _mm_sub_ps(C4, V02);
			C4 = _mm_add_ps(C4, V02);
			VECTOR C7 = _mm_add_ps(C6, V03);
			C6 = _mm_sub_ps(C6, V03);

			C0 = _mm_shuffle_ps(C0, C1, _MM_SHUFFLE(3, 1, 2, 0));
			C2 = _mm_shuffle_ps(C2, C3, _MM_SHUFFLE(3, 1, 2, 0));
			C4 = _mm_shuffle_ps(C4, C5, _MM_SHUFFLE(3, 1, 2, 0));
			C6 = _mm_shuffle_ps(C6, C7, _MM_SHUFFLE(3, 1, 2, 0));
			C0 = PERMUTE_PS(C0, _MM_SHUFFLE(3, 1, 2, 0));
			C2 = PERMUTE_PS(C2, _MM_SHUFFLE(3, 1, 2, 0));
			C4 = PERMUTE_PS(C4, _MM_SHUFFLE(3, 1, 2, 0));
			C6 = PERMUTE_PS(C6, _MM_SHUFFLE(3, 1, 2, 0));
			// Get the determinant
			VECTOR vTemp = VEC4::Dot(C0, MT.r[0]);
			if (pDeterminant != nullptr)
				*pDeterminant = vTemp;
			vTemp = _mm_div_ps(g_One, vTemp);
			MATRIX mResult;
			mResult.r[0] = _mm_mul_ps(C0, vTemp);
			mResult.r[1] = _mm_mul_ps(C2, vTemp);
			mResult.r[2] = _mm_mul_ps(C4, vTemp);
			mResult.r[3] = _mm_mul_ps(C6, vTemp);
			
			return mResult;
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV VectorTensorProduct(A_VECTOR V1, A_VECTOR V2) noexcept
		{
			MATRIX mResult;
			mResult.r[0] = VEC::Multiply(VEC::Swizzle<0, 0, 0, 0>(V1), V2);
			mResult.r[1] = VEC::Multiply(VEC::Swizzle<1, 1, 1, 1>(V1), V2);
			mResult.r[2] = VEC::Multiply(VEC::Swizzle<2, 2, 2, 2>(V1), V2);
			mResult.r[3] = VEC::Multiply(VEC::Swizzle<3, 3, 3, 3>(V1), V2);
			
			return mResult;
		}

		FORCE_INLINE MATRIX VEC_CALLCONV Determinant(A_MATRIX m) noexcept
		{
			static const VECTOR_F32 sign = { { { 1.0f, -1.0f, 1.0f, -1.0f } } };

			VECTOR V0 = VEC::Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_X, SWIZZLE_X>(m.r[2]);
			VECTOR V1 = VEC::Swizzle<SWIZZLE_Z, SWIZZLE_Z, SWIZZLE_Y, SWIZZLE_Y>(m.r[3]);
			VECTOR V2 = VEC::Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_X, SWIZZLE_X>(m.r[2]);
			VECTOR V3 = VEC::Swizzle<SWIZZLE_W, SWIZZLE_W, SWIZZLE_W, SWIZZLE_Z>(m.r[3]);
			VECTOR V4 = VEC::Swizzle<SWIZZLE_Z, SWIZZLE_Z, SWIZZLE_Y, SWIZZLE_Y>(m.r[2]);
			VECTOR V5 = VEC::Swizzle<SWIZZLE_W, SWIZZLE_W, SWIZZLE_W, SWIZZLE_Z>(m.r[3]);

			VECTOR P0 = VEC::Multiply(V0, V1);
			VECTOR P1 = VEC::Multiply(V2, V3);
			VECTOR P2 = VEC::Multiply(V4, V5);

			V0 = VEC::Swizzle<SWIZZLE_Z, SWIZZLE_Z, SWIZZLE_Y, SWIZZLE_Y>(m.r[2]);
			V1 = VEC::Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_X, SWIZZLE_X>(m.r[3]);
			V2 = VEC::Swizzle<SWIZZLE_W, SWIZZLE_W, SWIZZLE_W, SWIZZLE_Z>(m.r[2]);
			V3 = VEC::Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_X, SWIZZLE_X>(m.r[3]);
			V4 = VEC::Swizzle<SWIZZLE_W, SWIZZLE_W, SWIZZLE_W, SWIZZLE_Z>(m.r[2]);
			V5 = VEC::Swizzle<SWIZZLE_Z, SWIZZLE_Z, SWIZZLE_Y, SWIZZLE_Y>(m.r[3]);

			P0 = VEC::NegativeMultiplySubtract(V0, V1, P0);
			P1 = VEC::NegativeMultiplySubtract(V2, V3, P1);
			P2 = VEC::NegativeMultiplySubtract(V4, V5, P2);

			V0 = VEC::Swizzle<SWIZZLE_W, SWIZZLE_W, SWIZZLE_W, SWIZZLE_Z>(m.r[1]);
			V1 = VEC::Swizzle<SWIZZLE_Z, SWIZZLE_Z, SWIZZLE_Y, SWIZZLE_Y>(m.r[1]);
			V2 = VEC::Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_X, SWIZZLE_X>(m.r[1]);

			VECTOR S = VEC::Multiply(m.r[0], sign.v);
			VECTOR R = VEC::Multiply(V0, P0);
			R = VEC::NegativeMultiplySubtract(V1, P1, R);
			R = VEC::MultiplyAdd(V2, P2, R);

			return VEC4::Dot(S, R);
		}

#define _3RANKDECOMPOSE(a, b, c, x, y, z)      \
    if((x) < (y))                   \
    {                               \
        if((y) < (z))               \
        {                           \
            (a) = 2;                \
            (b) = 1;                \
            (c) = 0;                \
        }                           \
        else                        \
        {                           \
            (a) = 1;                \
                                    \
            if((x) < (z))           \
            {                       \
                (b) = 2;            \
                (c) = 0;            \
            }                       \
            else                    \
            {                       \
                (b) = 0;            \
                (c) = 2;            \
            }                       \
        }                           \
    }                               \
    else                            \
    {                               \
        if((x) < (z))               \
        {                           \
            (a) = 2;                \
            (b) = 0;                \
            (c) = 1;                \
        }                           \
        else                        \
        {                           \
            (a) = 0;                \
                                    \
            if((y) < (z))           \
            {                       \
                (b) = 2;            \
                (c) = 1;            \
            }                       \
            else                    \
            {                       \
                (b) = 1;            \
                (c) = 2;            \
            }                       \
        }                           \
    }

#define _3_DECOMP_EPSILON 0.0001f

		_Use_decl_annotations_
		FORCE_INLINE bool VEC_CALLCONV Decompose(VECTOR* outScale, VECTOR* outRotationQuaternion, VECTOR* outTranslation, A_MATRIX m) noexcept
		{
			static const VECTOR* pvCanonicalBasis[3] = {
				&g_IdentityR0.v,
				&g_IdentityR1.v,
				&g_IdentityR2.v
			};

#if defined(DEBUG) || defined(_DEBUG)
			assert(outScale != nullptr);
			assert(outRotQuat != nullptr);
			assert(outTrans != nullptr);
#endif

			// Get the translation
			outTranslation[0] = m.r[3];

			VECTOR* ppvBasis[3];
			MATRIX matTemp;
			ppvBasis[0] = &matTemp.r[0];
			ppvBasis[1] = &matTemp.r[1];
			ppvBasis[2] = &matTemp.r[2];

			matTemp.r[0] = m.r[0];
			matTemp.r[1] = m.r[1];
			matTemp.r[2] = m.r[2];
			matTemp.r[3] = g_IdentityR3.v;

			auto pfScales = reinterpret_cast<float*>(outScale);

			size_t a, b, c;
			VEC::GetXPtr(&pfScales[0], VEC3::Length(ppvBasis[0][0]));
			VEC::GetXPtr(&pfScales[1], VEC3::Length(ppvBasis[1][0]));
			VEC::GetXPtr(&pfScales[2], VEC3::Length(ppvBasis[2][0]));
			pfScales[3] = 0.f;

			_3RANKDECOMPOSE(a, b, c, pfScales[0], pfScales[1], pfScales[2])

				if (pfScales[a] < _3_DECOMP_EPSILON)
				{
					ppvBasis[a][0] = pvCanonicalBasis[a][0];
				}
			ppvBasis[a][0] = VEC3::Normalize(ppvBasis[a][0]);

			if (pfScales[b] < _3_DECOMP_EPSILON)
			{
				size_t aa, bb, cc;
				float fAbsX, fAbsY, fAbsZ;

				fAbsX = fabsf(VEC::GetX(ppvBasis[a][0]));
				fAbsY = fabsf(VEC::GetY(ppvBasis[a][0]));
				fAbsZ = fabsf(VEC::GetZ(ppvBasis[a][0]));

				_3RANKDECOMPOSE(aa, bb, cc, fAbsX, fAbsY, fAbsZ)

					ppvBasis[b][0] = VEC3::Cross(ppvBasis[a][0], pvCanonicalBasis[cc][0]);
			}

			ppvBasis[b][0] = VEC3::Normalize(ppvBasis[b][0]);

			if (pfScales[c] < _3_DECOMP_EPSILON)
			{
				ppvBasis[c][0] = VEC3::Cross(ppvBasis[a][0], ppvBasis[b][0]);
			}

			ppvBasis[c][0] = VEC3::Normalize(ppvBasis[c][0]);

			float fDet = VEC::GetX(Determinant(matTemp));

			// use Kramer's rule to check for handedness of coordinate system
			if (fDet < 0.0f)
			{
				// switch coordinate system by negating the scale and inverting the basis vector on the x-axis
				pfScales[a] = -pfScales[a];
				ppvBasis[a][0] = VEC::Negate(ppvBasis[a][0]);

				fDet = -fDet;
			}

			fDet -= 1.0f;
			fDet *= fDet;

			if (_3_DECOMP_EPSILON < fDet)
			{
				// Non-SRT matrix encountered
				return false;
			}

			// generate the quaternion from the matrix
			outRotQuat[0] = Quat::RotationMatrix(matTemp);
			return true;
		}

#undef _3_DECOMP_EPSILON
#undef _3RANKDECOMPOSE

		FORCE_INLINE MATRIX VEC_CALLCONV Identity() noexcept
		{
			MATRIX m;
			m.r[0] = g_IdentityR0.v;
			m.r[1] = g_IdentityR1.v;
			m.r[2] = g_IdentityR2.v;
			m.r[3] = g_IdentityR3.v;

			return m;
		}

		FORCE_INLINE MATRIX VEC_CALLCONV Set(
			float m00, float m01, float m02, float m03, 
			float m10, float m11, float m12, float m13, 
			float m20, float m21, float m22, float m23, 
			float m30, float m31, float m32, float m33
		) noexcept
		{
			MATRIX M;
#if defined(_NO_INTRINSICS_)
			M.m[0][0] = m00; M.m[0][1] = m01; M.m[0][2] = m02; M.m[0][3] = m03;
			M.m[1][0] = m10; M.m[1][1] = m11; M.m[1][2] = m12; M.m[1][3] = m13;
			M.m[2][0] = m20; M.m[2][1] = m21; M.m[2][2] = m22; M.m[2][3] = m23;
			M.m[3][0] = m30; M.m[3][1] = m31; M.m[3][2] = m32; M.m[3][3] = m33;
#else
			M.r[0] = VEC::Set(m00, m01, m02, m03);
			M.r[1] = VEC::Set(m10, m11, m12, m13);
			M.r[2] = VEC::Set(m20, m21, m22, m23);
			M.r[3] = VEC::Set(m30, m31, m32, m33);
#endif
    		return M;
		}

		FORCE_INLINE MATRIX VEC_CALLCONV Translation(float offsetX, float offsetY, float offsetZ) noexcept
		{
#if defined(_XM_NO_INTRINSICS_)
			MATRIX M;
			M.m[0][0] = 1.0f;
			M.m[0][1] = 0.0f;
			M.m[0][2] = 0.0f;
			M.m[0][3] = 0.0f;

			M.m[1][0] = 0.0f;
			M.m[1][1] = 1.0f;
			M.m[1][2] = 0.0f;
			M.m[1][3] = 0.0f;

			M.m[2][0] = 0.0f;
			M.m[2][1] = 0.0f;
			M.m[2][2] = 1.0f;
			M.m[2][3] = 0.0f;

			M.m[3][0] = offsetX;
			M.m[3][1] = offsetY;
			M.m[3][2] = offsetZ;
			M.m[3][3] = 1.0f;
			
			return M;

#elif defined(_SSE2_INTRINSICS_)
			MATRIX M;
			M.r[0] = g_IdentityR0.v;
			M.r[1] = g_IdentityR1.v;
			M.r[2] = g_IdentityR2.v;
			M.r[3] = VEC::Set(offsetX, offsetY, offsetZ, 1.0f);
			
			return M;
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV TranslationFromVector(A_VECTOR offset) noexcept
		{
#if defined(_NO_INTRINSICS_)
			MATRIX M;
			M.m[0][0] = 1.0f;
			M.m[0][1] = 0.0f;
			M.m[0][2] = 0.0f;
			M.m[0][3] = 0.0f;

			M.m[1][0] = 0.0f;
			M.m[1][1] = 1.0f;
			M.m[1][2] = 0.0f;
			M.m[1][3] = 0.0f;

			M.m[2][0] = 0.0f;
			M.m[2][1] = 0.0f;
			M.m[2][2] = 1.0f;
			M.m[2][3] = 0.0f;

			M.m[3][0] = offset.vector4_f32[0];
			M.m[3][1] = offset.vector4_f32[1];
			M.m[3][2] = offset.vector4_f32[2];
			M.m[3][3] = 1.0f;
			
			return M;

#elif defined(_SSE2_INTRINSICS_)
			MATRIX M;
			M.r[0] = g_IdentityR0.v;
			M.r[1] = g_IdentityR1.v;
			M.r[2] = g_IdentityR2.v;
			M.r[3] = VEC::Select(g_IdentityR3.v, offset, g_Select1110.v);
			
			return M;
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV Scaling(float scaleX, float scaleY, float scaleZ) noexcept
		{
#if defined(_NO_INTRINSICS_)
			MATRIX M;
			M.m[0][0] = scaleX;
			M.m[0][1] = 0.0f;
			M.m[0][2] = 0.0f;
			M.m[0][3] = 0.0f;

			M.m[1][0] = 0.0f;
			M.m[1][1] = scaleY;
			M.m[1][2] = 0.0f;
			M.m[1][3] = 0.0f;

			M.m[2][0] = 0.0f;
			M.m[2][1] = 0.0f;
			M.m[2][2] = scaleZ;
			M.m[2][3] = 0.0f;

			M.m[3][0] = 0.0f;
			M.m[3][1] = 0.0f;
			M.m[3][2] = 0.0f;
			M.m[3][3] = 1.0f;
			
			return M;

#elif defined(_SSE2_INTRINSICS_)
			MATRIX M;
			M.r[0] = _mm_set_ps(0, 0, 0, scaleX);
			M.r[1] = _mm_set_ps(0, 0, scaleY, 0);
			M.r[2] = _mm_set_ps(0, scaleZ, 0, 0);
			M.r[3] = g_IdentityR3.v;
			
			return M;
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV ScalingFromVector(A_VECTOR scale) noexcept
		{
#if defined(_NO_INTRINSICS_)
			MATRIX M;
			M.m[0][0] = scale.vector4_f32[0];
			M.m[0][1] = 0.0f;
			M.m[0][2] = 0.0f;
			M.m[0][3] = 0.0f;

			M.m[1][0] = 0.0f;
			M.m[1][1] = scale.vector4_f32[1];
			M.m[1][2] = 0.0f;
			M.m[1][3] = 0.0f;

			M.m[2][0] = 0.0f;
			M.m[2][1] = 0.0f;
			M.m[2][2] = scale.vector4_f32[2];
			M.m[2][3] = 0.0f;

			M.m[3][0] = 0.0f;
			M.m[3][1] = 0.0f;
			M.m[3][2] = 0.0f;
			M.m[3][3] = 1.0f;
			
			return M;

#elif defined(_SSE2_INTRINSICS_)
			MATRIX M;
			M.r[0] = _mm_and_ps(scale, g_MaskX);
			M.r[1] = _mm_and_ps(scale, g_MaskY);
			M.r[2] = _mm_and_ps(scale, g_MaskZ);
			M.r[3] = g_IdentityR3.v;
			
			return M;
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV RotationX(float angle) noexcept
		{
#if defined(_NO_INTRINSICS_)
			float    fSinAngle;
			float    fCosAngle;
			VEC::ScalarSinCos(&fSinAngle, &fCosAngle, angle);

			MATRIX M;
			M.m[0][0] = 1.0f;
			M.m[0][1] = 0.0f;
			M.m[0][2] = 0.0f;
			M.m[0][3] = 0.0f;

			M.m[1][0] = 0.0f;
			M.m[1][1] = fCosAngle;
			M.m[1][2] = fSinAngle;
			M.m[1][3] = 0.0f;

			M.m[2][0] = 0.0f;
			M.m[2][1] = -fSinAngle;
			M.m[2][2] = fCosAngle;
			M.m[2][3] = 0.0f;

			M.m[3][0] = 0.0f;
			M.m[3][1] = 0.0f;
			M.m[3][2] = 0.0f;
			M.m[3][3] = 1.0f;
			
			return M;

#elif defined(_SSE2_INTRINSICS_)
			float    SinAngle;
			float    CosAngle;
			VEC::ScalarSinCos(&SinAngle, &CosAngle, angle);

			VECTOR vSin = _mm_set_ss(SinAngle);
			VECTOR vCos = _mm_set_ss(CosAngle);
			// x = 0,y = cos,z = sin, w = 0
			vCos = _mm_shuffle_ps(vCos, vSin, _MM_SHUFFLE(3, 0, 0, 3));
			MATRIX M;
			M.r[0] = g_IdentityR0;
			M.r[1] = vCos;
			// x = 0,y = sin,z = cos, w = 0
			vCos = PERMUTE_PS(vCos, _MM_SHUFFLE(3, 1, 2, 0));
			// x = 0,y = -sin,z = cos, w = 0
			vCos = _mm_mul_ps(vCos, g_NegateY);
			M.r[2] = vCos;
			M.r[3] = g_IdentityR3;
			
			return M;
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV RotationY(float angle) noexcept
		{
#if defined(_NO_INTRINSICS_)
			float    fSinAngle;
			float    fCosAngle;
			VEC::ScalarSinCos(&fSinAngle, &fCosAngle, angle);

			MATRIX M;
			M.m[0][0] = fCosAngle;
			M.m[0][1] = 0.0f;
			M.m[0][2] = -fSinAngle;
			M.m[0][3] = 0.0f;

			M.m[1][0] = 0.0f;
			M.m[1][1] = 1.0f;
			M.m[1][2] = 0.0f;
			M.m[1][3] = 0.0f;

			M.m[2][0] = fSinAngle;
			M.m[2][1] = 0.0f;
			M.m[2][2] = fCosAngle;
			M.m[2][3] = 0.0f;

			M.m[3][0] = 0.0f;
			M.m[3][1] = 0.0f;
			M.m[3][2] = 0.0f;
			M.m[3][3] = 1.0f;
			
			return M;

#elif defined(_SSE2_INTRINSICS_)
			float    SinAngle;
			float    CosAngle;
			VEC::ScalarSinCos(&SinAngle, &CosAngle, angle);

			VECTOR vSin = _mm_set_ss(SinAngle);
			VECTOR vCos = _mm_set_ss(CosAngle);
			// x = sin,y = 0,z = cos, w = 0
			vSin = _mm_shuffle_ps(vSin, vCos, _MM_SHUFFLE(3, 0, 3, 0));
			MATRIX M;
			M.r[2] = vSin;
			M.r[1] = g_IdentityR1;
			// x = cos,y = 0,z = sin, w = 0
			vSin = PERMUTE_PS(vSin, _MM_SHUFFLE(3, 0, 1, 2));
			// x = cos,y = 0,z = -sin, w = 0
			vSin = _mm_mul_ps(vSin, g_NegateZ);
			M.r[0] = vSin;
			M.r[3] = g_IdentityR3;
			
			return M;
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV RotationZ(angle) noexcept
		{
#if defined(_NO_INTRINSICS_)
			float    fSinAngle;
			float    fCosAngle;
			VEC::ScalarSinCos(&fSinAngle, &fCosAngle, angle);

			MATRIX M;
			M.m[0][0] = fCosAngle;
			M.m[0][1] = fSinAngle;
			M.m[0][2] = 0.0f;
			M.m[0][3] = 0.0f;

			M.m[1][0] = -fSinAngle;
			M.m[1][1] = fCosAngle;
			M.m[1][2] = 0.0f;
			M.m[1][3] = 0.0f;

			M.m[2][0] = 0.0f;
			M.m[2][1] = 0.0f;
			M.m[2][2] = 1.0f;
			M.m[2][3] = 0.0f;

			M.m[3][0] = 0.0f;
			M.m[3][1] = 0.0f;
			M.m[3][2] = 0.0f;
			M.m[3][3] = 1.0f;
			
			return M;

#elif defined(_SSE2_INTRINSICS_)
			float    SinAngle;
			float    CosAngle;
			VEC::ScalarSinCos(&SinAngle, &CosAngle, angle);

			VECTOR vSin = _mm_set_ss(SinAngle);
			VECTOR vCos = _mm_set_ss(CosAngle);
			// x = cos,y = sin,z = 0, w = 0
			vCos = _mm_unpacklo_ps(vCos, vSin);
			MATRIX M;
			M.r[0] = vCos;
			// x = sin,y = cos,z = 0, w = 0
			vCos = PERMUTE_PS(vCos, _MM_SHUFFLE(3, 2, 0, 1));
			// x = cos,y = -sin,z = 0, w = 0
			vCos = _mm_mul_ps(vCos, g_NegateX);
			M.r[1] = vCos;
			M.r[2] = g_IdentityR2;
			M.r[3] = g_IdentityR3;
			
			return M;
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV RotationPitchYawRoll(float pitch, float yaw, float roll) noexcept
		{
#if defined(_NO_INTRINSICS_)
			float cp = cosf(pitch);
			float sp = sinf(pitch);

			float cy = cosf(yaw);
			float sy = sinf(yaw);

			float cr = cosf(roll);
			float sr = sinf(roll);

			MATRIX M;
			M.m[0][0] = cr * cy + sr * sp * sy;
			M.m[0][1] = sr * cp;
			M.m[0][2] = sr * sp * cy - cr * sy;
			M.m[0][3] = 0.0f;

			M.m[1][0] = cr * sp * sy - sr * cy;
			M.m[1][1] = cr * cp;
			M.m[1][2] = sr * sy + cr * sp * cy;
			M.m[1][3] = 0.0f;

			M.m[2][0] = cp * sy;
			M.m[2][1] = -sp;
			M.m[2][2] = cp * cy;
			M.m[2][3] = 0.0f;

			M.m[3][0] = 0.0f;
			M.m[3][1] = 0.0f;
			M.m[3][2] = 0.0f;
			M.m[3][3] = 1.0f;
			
			return M;

#else
			VECTOR angles = VEC::Set(pitch, yaw, roll, 0.0f);
    		
			return RotationRollPitchYawFromVector(angles);
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV RotationPitchYawRollFromVector(A_VECTOR angles) noexcept
		{
#if defined(_NO_INTRINSICS_)
			float cp = cosf(angles.vector4_f32[0]);
			float sp = sinf(angles.vector4_f32[0]);

			float cy = cosf(angles.vector4_f32[1]);
			float sy = sinf(angles.vector4_f32[1]);

			float cr = cosf(angles.vector4_f32[2]);
			float sr = sinf(angles.vector4_f32[2]);

			MATRIX M;
			M.m[0][0] = cr * cy + sr * sp * sy;
			M.m[0][1] = sr * cp;
			M.m[0][2] = sr * sp * cy - cr * sy;
			M.m[0][3] = 0.0f;

			M.m[1][0] = cr * sp * sy - sr * cy;
			M.m[1][1] = cr * cp;
			M.m[1][2] = sr * sy + cr * sp * cy;
			M.m[1][3] = 0.0f;

			M.m[2][0] = cp * sy;
			M.m[2][1] = -sp;
			M.m[2][2] = cp * cy;
			M.m[2][3] = 0.0f;

			M.m[3][0] = 0.0f;
			M.m[3][1] = 0.0f;
			M.m[3][2] = 0.0f;
			M.m[3][3] = 1.0f;
			
			return M;

#else 
			static const VECTOR_F32  sign = { { { 1.0f, -1.0f, -1.0f, 1.0f } } };

			VECTOR SinAngles, CosAngles;
			VEC::SinCos(&SinAngles, &CosAngles, angles);

			VECTOR P0 = VEC::Permute<PERMUTE_1X, PERMUTE_0Z, PERMUTE_1Z, PERMUTE_1X>(SinAngles, CosAngles);
			VECTOR Y0 = VEC::Permute<PERMUTE_0Y, PERMUTE_1X, PERMUTE_1X, PERMUTE_1Y>(SinAngles, CosAngles);
			VECTOR P1 = VEC::Permute<PERMUTE_1Z, PERMUTE_0Z, PERMUTE_1Z, PERMUTE_0Z>(SinAngles, CosAngles);
			VECTOR Y1 = VEC::Permute<PERMUTE_1Y, PERMUTE_1Y, PERMUTE_0Y, PERMUTE_0Y>(SinAngles, CosAngles);
			VECTOR P2 = VEC::Permute<PERMUTE_0Z, PERMUTE_1Z, PERMUTE_0Z, PERMUTE_1Z>(SinAngles, CosAngles);
			VECTOR P3 = VEC::Permute<PERMUTE_0Y, PERMUTE_0Y, PERMUTE_1Y, PERMUTE_1Y>(SinAngles, CosAngles);
			VECTOR Y2 = VEC::SplatX(SinAngles);
			VECTOR NS = VEC::Negate(SinAngles);

			VECTOR Q0 = VEC::(P0, Y0);
			VECTOR Q1 = VEC::Multiply(P1, sign.v);
			Q1 = VEC::Multiply(Q1, Y1);
			VECTOR Q2 = VEC::Multiply(P2, Y2);
			Q2 = VEC::MultiplyAdd(Q2, P3, Q1);

			VECTOR V0 = VEC::Permute<PERMUTE_1X, PERMUTE_0Y, PERMUTE_1Z, PERMUTE_0W>(Q0, Q2);
			VECTOR V1 = VEC::Permute<PERMUTE_1Y, PERMUTE_0Z, PERMUTE_1W, PERMUTE_0W>(Q0, Q2);
			VECTOR V2 = VEC::Permute<PERMUTE_0X, PERMUTE_1X, PERMUTE_0W, PERMUTE_0W>(Q0, NS);

			MATRIX M;
			M.r[0] = VEC::Select(g_Zero, V0, g_Select1110.v);
			M.r[1] = VEC::Select(g_Zero, V1, g_Select1110.v);
			M.r[2] = VEC::Select(g_Zero, V2, g_Select1110.v);
			M.r[3] = g_IdentityR3;
			
			return M;
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV RotationNormal(A_VECTOR normalAxis, float angle) noexcept
		{
#if defined(_NO_INTRINSICS_)
			float    fSinAngle;
			float    fCosAngle;
			VEC::ScalarSinCos(&fSinAngle, &fCosAngle, angle);

			VECTOR A = VEC::Set(fSinAngle, fCosAngle, 1.0f - fCosAngle, 0.0f);

			VECTOR C2 = VEC::SplatZ(A);
			VECTOR C1 = VEC::SplatY(A);
			VECTOR C0 = VEC::SplatX(A);

			VECTOR N0 = VEC::Swizzle<SWIZZLE_Y, SWIZZLE_Z, SWIZZLE_X, SWIZZLE_W>(normalAxis);
			VECTOR N1 = VEC::Swizzle<SWIZZLE_Z, SWIZZLE_X, SWIZZLE_Y, SWIZZLE_W>(normalAxis);

			VECTOR V0 = VEC::Multiply(C2, N0);
			V0 = VEC::Multiply(V0, N1);

			VECTOR R0 = VEC::Multiply(C2, normalAxis);
			R0 = VE::MultiplyAdd(R0, normalAxis, C1);

			VECTOR R1 = VEC::MultiplyAdd(C0, normalAxis, V0);
			VECTOR R2 = VEC::NegativeMultiplySubtract(C0, normalAxis, V0);

			V0 = VEC::Select(A, R0, g_Select1110.v);
			VECTOR V1 = VEC::Permute<PERMUTE_0Z, PERMUTE_1Y, PERMUTE_1Z, PERMUTE_0X>(R1, R2);
			VECTOR V2 = VEC::Permute<PERMUTE_0Y, PERMUTE_1X, PERMUTE_0Y, PERMUTE_1X>(R1, R2);

			MATRIX M;
			M.r[0] = VEC::Permute<PERMUTE_0X, PERMUTE_1X, PERMUTE_1Y, PERMUTE_0W>(V0, V1);
			M.r[1] = VEC::Permute<PERMUTE_1Z, PERMUTE_0Y, PERMUTE_1W, PERMUTE_0W>(V0, V1);
			M.r[2] = VEC::Permute<PERMUTE_1X, PERMUTE_1Y, PERMUTE_0Z, PERMUTE_0W>(V0, V2);
			M.r[3] = g_IdentityR3.v;
			
			return M;

#elif defined(_SSE2_INTRINSICS_)
			float    fSinAngle;
			float    fCosAngle;
			VEC::ScalarSinCos(&fSinAngle, &fCosAngle, angle);

			VECTOR C2 = _mm_set_ps1(1.0f - fCosAngle);
			VECTOR C1 = _mm_set_ps1(fCosAngle);
			VECTOR C0 = _mm_set_ps1(fSinAngle);

			VECTOR N0 = PERMUTE_PS(normalAxis, _MM_SHUFFLE(3, 0, 2, 1));
			VECTOR N1 = PERMUTE_PS(normalAxis, _MM_SHUFFLE(3, 1, 0, 2));

			VECTOR V0 = _mm_mul_ps(C2, N0);
			V0 = _mm_mul_ps(V0, N1);

			VECTOR R0 = _mm_mul_ps(C2, normalAxis);
			R0 = _mm_mul_ps(R0, normalAxis);
			R0 = _mm_add_ps(R0, C1);

			VECTOR R1 = _mm_mul_ps(C0, normalAxis);
			R1 = _mm_add_ps(R1, V0);
			VECTOR R2 = _mm_mul_ps(C0, normalAxis);
			R2 = _mm_sub_ps(V0, R2);

			V0 = _mm_and_ps(R0, g_Mask3);
			VECTOR V1 = _mm_shuffle_ps(R1, R2, _MM_SHUFFLE(2, 1, 2, 0));
			V1 = PERMUTE_PS(V1, _MM_SHUFFLE(0, 3, 2, 1));
			VECTOR V2 = _mm_shuffle_ps(R1, R2, _MM_SHUFFLE(0, 0, 1, 1));
			V2 = PERMUTE_PS(V2, _MM_SHUFFLE(2, 0, 2, 0));

			R2 = _mm_shuffle_ps(V0, V1, _MM_SHUFFLE(1, 0, 3, 0));
			R2 = PERMUTE_PS(R2, _MM_SHUFFLE(1, 3, 2, 0));

			MATRIX M;
			M.r[0] = R2;

			R2 = _mm_shuffle_ps(V0, V1, _MM_SHUFFLE(3, 2, 3, 1));
			R2 = PERMUTE_PS(R2, _MM_SHUFFLE(1, 3, 0, 2));
			M.r[1] = R2;

			V2 = _mm_shuffle_ps(V2, V0, _MM_SHUFFLE(3, 2, 1, 0));
			M.r[2] = V2;
			M.r[3] = g_IdentityR3.v;
			
			return M;
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV RotationAxis(A_VECTOR axis, float angle) noexcept
		{
#if defined(_NO_INTRINSICS_)
			assert(!VEC3::Equal(axis, VEC::Zero()));
			assert(!VEC3::IsInfinite(axis));
#endif

			VECTOR normal = VEC3::Normalize(axis);
			
			return RotationNormal(normal, angle);
		}

		FORCE_INLINE MATRIX VEC_CALLCONV RotationQuaternion(A_VECTOR quaternion) noexcept
		{
#if defined(_NO_INTRINSICS_)
			float qx = quaternion.vector4_f32[0];
			float qxx = qx * qx;

			float qy = quaternion.vector4_f32[1];
			float qyy = qy * qy;

			float qz = quaternion.vector4_f32[2];
			float qzz = qz * qz;

			float qw = quaternion.vector4_f32[3];

			MATRIX M;
			M.m[0][0] = 1.f - 2.f * qyy - 2.f * qzz;
			M.m[0][1] = 2.f * qx * qy + 2.f * qz * qw;
			M.m[0][2] = 2.f * qx * qz - 2.f * qy * qw;
			M.m[0][3] = 0.f;

			M.m[1][0] = 2.f * qx * qy - 2.f * qz * qw;
			M.m[1][1] = 1.f - 2.f * qxx - 2.f * qzz;
			M.m[1][2] = 2.f * qy * qz + 2.f * qx * qw;
			M.m[1][3] = 0.f;

			M.m[2][0] = 2.f * qx * qz + 2.f * qy * qw;
			M.m[2][1] = 2.f * qy * qz - 2.f * qx * qw;
			M.m[2][2] = 1.f - 2.f * qxx - 2.f * qyy;
			M.m[2][3] = 0.f;

			M.m[3][0] = 0.f;
			M.m[3][1] = 0.f;
			M.m[3][2] = 0.f;
			M.m[3][3] = 1.0f;
			
			return M;

#elif defined(_SSE2_INTRINSICS_)
			static const VECTOR_F32  Constant1110 = { { { 1.0f, 1.0f, 1.0f, 0.0f } } };

			VECTOR Q0 = _mm_add_ps(quaternion, quaternion);
			VECTOR Q1 = _mm_mul_ps(quaternion, Q0);

			VECTOR V0 = PERMUTE_PS(Q1, _MM_SHUFFLE(3, 0, 0, 1));
			V0 = _mm_and_ps(V0, g_Mask3);
			VECTOR V1 = PERMUTE_PS(Q1, _MM_SHUFFLE(3, 1, 2, 2));
			V1 = _mm_and_ps(V1, g_Mask3);
			VECTOR R0 = _mm_sub_ps(Constant1110, V0);
			R0 = _mm_sub_ps(R0, V1);

			V0 = PERMUTE_PS(quaternion, _MM_SHUFFLE(3, 1, 0, 0));
			V1 = PERMUTE_PS(Q0, _MM_SHUFFLE(3, 2, 1, 2));
			V0 = _mm_mul_ps(V0, V1);

			V1 = PERMUTE_PS(quaternion, _MM_SHUFFLE(3, 3, 3, 3));
			VECTOR V2 = PERMUTE_PS(Q0, _MM_SHUFFLE(3, 0, 2, 1));
			V1 = _mm_mul_ps(V1, V2);

			VECTOR R1 = _mm_add_ps(V0, V1);
			VECTOR R2 = _mm_sub_ps(V0, V1);

			V0 = _mm_shuffle_ps(R1, R2, _MM_SHUFFLE(1, 0, 2, 1));
			V0 = PERMUTE_PS(V0, _MM_SHUFFLE(1, 3, 2, 0));
			V1 = _mm_shuffle_ps(R1, R2, _MM_SHUFFLE(2, 2, 0, 0));
			V1 = PERMUTE_PS(V1, _MM_SHUFFLE(2, 0, 2, 0));

			Q1 = _mm_shuffle_ps(R0, V0, _MM_SHUFFLE(1, 0, 3, 0));
			Q1 = PERMUTE_PS(Q1, _MM_SHUFFLE(1, 3, 2, 0));

			MATRIX M;
			M.r[0] = Q1;

			Q1 = _mm_shuffle_ps(R0, V0, _MM_SHUFFLE(3, 2, 3, 1));
			Q1 = PERMUTE_PS(Q1, _MM_SHUFFLE(1, 3, 0, 2));
			M.r[1] = Q1;

			Q1 = _mm_shuffle_ps(V1, R0, _MM_SHUFFLE(3, 2, 1, 0));
			M.r[2] = Q1;
			M.r[3] = g_IdentityR3;
			
			return M;
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV Transformation2D(
			A_VECTOR scalingOrigin, 
			float scalingOrientation, 
			A_VECTOR scaling, 
			A_VECTOR rotationOrigin, 
			float rotation, 
			B_VECTOR translation
		) noexcept
		{
			// M = Inverse(MScalingOrigin) * Transpose(MScalingOrientation) * MScaling * MScalingOrientation *
    		//         MScalingOrigin * Inverse(MRotationOrigin) * MRotation * MRotationOrigin * MTranslation;

			VECTOR VScalingOrigin = VEC::Select(g_Select1100.v, scalingOrigin, g_Select1100.v);
			VECTOR NegScalingOrigin = VEC::Negate(VScalingOrigin);

			MATRIX MScalingOriginI = TranslationFromVector(NegScalingOrigin);
			MATRIX MScalingOrientation = RotationZ(scalingOrientation);
			MATRIX MScalingOrientationT = Transpose(scalingOrientation);
			VECTOR VScaling = VEC::Select(g_One.v, scaling, g_Select1100.v);
			MATRIX MScaling = ScalingFromVector(VScaling);
			VECTOR VRotationOrigin = VEC::Select(g_Select1100.v, rotationOrigin, g_Select1100.v);
			MATRIX MRotation = RotationZ(rotation);
			VECTOR VTranslation = VEC::Select(g_Select1100.v, translation, g_Select1100.v);

			MATRIX M = Multiply(MScalingOriginI, MScalingOrientationT);
			M = Multiply(M, MScaling);
			M = Multiply(M, MScalingOrientation);
			M.r[3] = VEC::Add(M.r[3], VScalingOrigin);
			M.r[3] = VEC::Subtract(M.r[3], VRotationOrigin);
			M = Multiply(M, MRotation);
			M.r[3] = VEC::Add(M.r[3], VRotationOrigin);
			M.r[3] = VEC::Add(M.r[3], VTranslation);

			return M;
		}

		FORCE_INLINE MATRIX VEC_CALLCONV Transformation(
			A_VECTOR scalingOrigin, 
			A_VECTOR scalingOrientationQuaternion, 
			A_VECTOR scaling, 
			B_VECTOR rotationOrigin, 
			C_VECTOR rotationQuaternion, 
			C_VECTOR translation
		) noexcept
		{
			// M = Inverse(MScalingOrigin) * Transpose(MScalingOrientation) * MScaling * MScalingOrientation *
			//         MScalingOrigin * Inverse(MRotationOrigin) * MRotation * MRotationOrigin * MTranslation;

			VECTOR VScalingOrigin = VEC::Select(g_Select1110.v, scalingOrigin, g_Select1110.v);
			VECTOR NegScalingOrigin = VEC::Negate(scalingOrigin);

			MATRIX MScalingOriginI = TranslationFromVector(NegScalingOrigin);
			MATRIX MScalingOrientation = RotationQuaternion(scalingOrientationQuaternion);
			MATRIX MScalingOrientationT = Transpose(MScalingOrientation);
			MATRIX MScaling = ScalingFromVector(scaling);
			VECTOR VRotationOrigin = VEC::Select(g_Select1110.v, rotationOrigin, g_Select1110.v);
			MATRIX MRotation = RotationQuaternion(rotationQuaternion);
			VECTOR VTranslation = VEC::Select(g_Select1110.v, translation, g_Select1110.v);

			MATRIX M;
			M = Multiply(MScalingOriginI, MScalingOrientationT);
			M = Multiply(M, MScaling);
			M = Multiply(M, MScalingOrientation);
			M.r[3] = VEC::Add(M.r[3], VScalingOrigin);
			M.r[3] = VEC::Subtract(M.r[3], VRotationOrigin);
			M = Multiply(M, MRotation);
			M.r[3] = VEC::Add(M.r[3], VRotationOrigin);
			M.r[3] = VEC::Add(M.r[3], VTranslation);
			
			return M;
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV AffineTransformation2D(
			A_VECTOR scaling, 
			A_VECTOR rotationOrigin, 
			float rotation, 
			A_VECTOR translation
		) noexcept
		{
			// M = MScaling * Inverse(MRotationOrigin) * MRotation * MRotationOrigin * MTranslation;

			VECTOR VScaling = VEC::Select(g_One.v, scaling, g_Select1100.v);
			MATRIX MScaling = ScalingFromVector(VScaling);
			VECTOR VRotationOrigin = VEC::Select(g_Select1100.v, rotationOrigin, g_Select1100.v);
			MATRIX MRotation = RotationZ(rotation);
			VECTOR VTranslation = VEC::Select(g_Select1100.v, translation, g_Select1100.v);

			MATRIX M;
			M = MScaling;
			M.r[3] = VEC::Subtract(M.r[3], VRotationOrigin);
			M = Multiply(M, MRotation);
			M.r[3] = VEC::Add(M.r[3], VRotationOrigin);
			M.r[3] = VEC::Add(M.r[3], VTranslation);
			
			return M;
		}

		FORCE_INLINE MATRIX VEC_CALLCONV AffineTransformation(
			A_VECTOR scaling, 
			A_VECTOR rotationOrigin, 
			A_VECTOR rotationQuaternion, 
			B_VECTOR translation
		) noexcept
		{
			// M = MScaling * Inverse(MRotationOrigin) * MRotation * MRotationOrigin * MTranslation;

			MATRIX MScaling = ScalingFromVector(scaling);
			VECTOR VRotationOrigin = VEC::Select(g_Select1110.v, rotationOrigin, g_Select1110.v);
			MATRIX MRotation = RotationQuaternion(rotationQuaternion);
			VECTOR VTranslation = VEC::Select(g_Select1110.v, translation, g_Select1110.v);

			MATRIX M;
			M = MScaling;
			M.r[3] = VEC::Subtract(M.r[3], VRotationOrigin);
			M = Multiply(M, MRotation);
			M.r[3] = VEC::Add(M.r[3], VRotationOrigin);
			M.r[3] = VEC::Add(M.r[3], VTranslation);
			
			return M;
		}

		FORCE_INLINE MATRIX VEC_CALLCONV Reflect(A_VECTOR reflectionPlane) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(!VEC3::Equal(reflectionPlane, VEC::Zero()));
			assert(!VEC::PlaneIsInfinite(reflectionPlane));
#endif

			static const VECTOR_F32 NegativeTwo = { { { -2.0f, -2.0f, -2.0f, 0.0f } } };

			VECTOR P = VEC::PlaneNormalize(reflectionPlane);
			XMVECTOR S = VEC::Multiply(P, NegativeTwo);

			VECTOR A = VEC::SplatX(P);
			VECTOR B = VEC::SplatY(P);
			VECTOR C = VEC::SplatZ(P);
			VECTOR D = VEC::SplatW(P);

			MATRIX M;
			M.r[0] = VEC::MultiplyAdd(A, S, g_IdentityR0.v);
			M.r[1] = VEC::MultiplyAdd(B, S, g_IdentityR1.v);
			M.r[2] = VEC::MultiplyAdd(C, S, g_IdentityR2.v);
			M.r[3] = VEC::MultiplyAdd(D, S, g_IdentityR3.v);
			
			return M;
		}

		FORCE_INLINE MATRIX VEC_CALLCONV Shadow(A_VECTOR shadowPlane, A_VECTOR lightPosition) noexcept
		{
			static const VECTOR_U32 Select0001 = { { { SELECT_0, SELECT_0, SELECT_0, SELECT_1 } } };

#if defined(DEBUG) || defined(_DEBUG)
			assert(!VEC3::Equal(shadowPlane, VEC::Zero()));
			assert(!VEC::PlaneIsInfinite(shadowPlane));
#endif

			VECTOR P = VEC::PlaneNormalize(shadowPlane);
			VECTOR dot = VEC::PlaneDot(P, lightPosition);
			P = VEC::Negate(P);
			VECTOR D = VEC::SplatW(P);
			VECTOR C = VEC::SplatZ(P);
			VECTOR B = VEC::SplatY(P);
			VECTOR A = VEC::SplatX(P);
			dot = VEC::Select(Select0001.v, dot, Select0001.v);

			MATRIX M;
			M.r[3] = VEC::MultiplyAdd(D, lightPosition, dot);
			dot = VEC::RotateLeft(dot, 1);
			M.r[2] = VEC::MultiplyAdd(C, lightPosition, dot);
			Dot = VEC::RotateLeft(dot, 1);
			M.r[1] = VEC::MultiplyAdd(B, lightPosition, dot);
			dot = VEC::RotateLeft(dot, 1);
			M.r[0] = VEC::MultiplyAdd(A, lightPosition, dot);
			
			return M;
		}

		FORCE_INLINE MATRIX VEC_CALLCONV LootAtLH(
			A_VECTOR eyePosition, 
			A_VECTOR focusPosition, 
			A_VECTOR upDirection
		) noexcept
		{
			VECTOR eyeDirection = VEC::Subtract(focusPosition, eyePosition);

			return LookToLH(eyePosition, eyeDirection, upDirection);
		}

		FORCE_INLINE MATRIX VEC_CALLCONV LookAtRH(
			A_VECTOR eyePosition, 
			A_VECTOR focusPosition, 
			A_VECTOR upDirection
		) noexcept
		{
			VECTOR negEyeDirection = VEC::Subtract(eyePosition, focusPosition);

			return LookToLH(eyePosition, negEyeDirection, upDirection);
		}

		FORCE_INLINE MATRIX VEC_CALLCONV LookToLH(
			A_VECTOR eyePosition, 
			A_VECTOR eyeDirection, 
			A_VECTOR upDirection
		) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(!VEC3::Equal(eyeDirection, VEC::Zero()));
			assert(!VEC3::IsInfinite(eyeDirection));
			assert(!VEC3::Equal(upDirection, VEC::Zero()));
			assert(!VE3::IsInfinite(upDirection));
#endif

			VECTOR R2 = VEC3::Normalize(eyeDirection);

			VECTOR R0 = VEC3::Cross(upDirection, R2);
			R0 = VEC3::Normalize(R0);

			VECTOR R1 = VEC3::Cross(R2, R0);

			VECTOR negEyePosition = VEC::Negate(eyePosition);

			VECTOR D0 = VEC3::Dot(R0, negEyePosition);
			VECTOR D1 = VEC3::Dot(R1, negEyePosition);
			VECTOR D2 = VEC3::Dot(R2, negEyePosition);

			MATRIX M;
			M.r[0] = VEC::Select(D0, R0, g_Select1110.v);
			M.r[1] = VEC::Select(D1, R1, g_Select1110.v);
			M.r[2] = VEC::Select(D2, R2, g_Select1110.v);
			M.r[3] = g_IdentityR3.v;

			M = Transpose(M);

			return M;
		}

		FORCE_INLINE MATRIX VEC_CALLCONV LookToRH(
			A_VECTOR eyePosition, 
			A_VECTOR eyeDirection, 
			A_VECTOR upDirection
		) noexcept
		{
			VECTOR negEyeDirection = VEC::Negate(eyeDirection);

			return LookToLH(eyePosition, negEyeDirection, upDirection);
		}

#ifdef _PREFAST_
#pragma prefast(push)
#pragma prefast(disable:28931, "PREfast noise: Esp:1266")
#endif

		FORCE_INLINE MATRIX VEC_CALLCONV PerspectiveLH(
			float viewWidth, 
			float viewHeight, 
			float nearZ, 
			float farZ
		) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(nearZ > 0.f && farZ > 0.f);
			assert(!ScalarNearEqual(viewWidth, 0.0f, 0.00001f));
			assert(!ScalarNearEqual(viewHeight, 0.0f, 0.00001f));
			assert(!ScalarNearEqual(farZ, nearZ, 0.00001f));
#endif

#if defined(_NO_INTRINSICS_)
			float twoNearZ = nearZ + nearZ;
			float fRange = farZ / (farZ - nearZ);

			MATRIX M;
			M.m[0][0] = twoNearZ / viewWidth;
			M.m[0][1] = 0.0f;
			M.m[0][2] = 0.0f;
			M.m[0][3] = 0.0f;

			M.m[1][0] = 0.0f;
			M.m[1][1] = twoNearZ / viewHeight;
			M.m[1][2] = 0.0f;
			M.m[1][3] = 0.0f;

			M.m[2][0] = 0.0f;
			M.m[2][1] = 0.0f;
			M.m[2][2] = fRange;
			M.m[2][3] = 1.0f;

			M.m[3][0] = 0.0f;
			M.m[3][1] = 0.0f;
			M.m[3][2] = -fRange * nearZ;
			M.m[3][3] = 0.0f;
			
			return M;

#elif defined(_SSE2_INTRINSICS_)
			MATRIX M;
			float twoNearZ = nearZ + nearZ;
			float fRange = farZ / (farZ - nearZ);
			// Note: This is recorded on the stack
			VECTOR rMem = {
				twoNearZ / viewWidth,
				twoNearZ / viewHeight,
				fRange,
				-fRange * nearZ
			};
			// Copy from memory to SSE register
			VECTOR vValues = rMem;
			VECTOR vTemp = _mm_setzero_ps();
			// Copy x only
			vTemp = _mm_move_ss(vTemp, vValues);
			// TwoNearZ / ViewWidth,0,0,0
			M.r[0] = vTemp;
			// 0,TwoNearZ / ViewHeight,0,0
			vTemp = vValues;
			vTemp = _mm_and_ps(vTemp, g_MaskY);
			M.r[1] = vTemp;
			// x=fRange,y=-fRange * NearZ,0,1.0f
			vValues = _mm_shuffle_ps(vValues, g_IdentityR3, _MM_SHUFFLE(3, 2, 3, 2));
			// 0,0,fRange,1.0f
			vTemp = _mm_setzero_ps();
			vTemp = _mm_shuffle_ps(vTemp, vValues, _MM_SHUFFLE(3, 0, 0, 0));
			M.r[2] = vTemp;
			// 0,0,-fRange * NearZ,0
			vTemp = _mm_shuffle_ps(vTemp, vValues, _MM_SHUFFLE(2, 1, 0, 0));
			M.r[3] = vTemp;
			
			return M;
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV PerspectiveRH(
			float viewWidth, 
			float viewHeight, 
			float nearZ, 
			float farZ
		) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(nearZ > 0.f && farZ > 0.f);
			assert(!ScalarNearEqual(viewWidth, 0.0f, 0.00001f));
			assert(!ScalarNearEqual(viewHeight, 0.0f, 0.00001f));
			assert(!ScalarNearEqual(farZ, nearZ, 0.00001f));
#endif

#if defined(_NO_INTRINSICS_)

			float twoNearZ = nearZ + nearZ;
			float fRange = farZ / (nearZ - farZ);

			MATRIX M;
			M.m[0][0] = twoNearZ / viewWidth;
			M.m[0][1] = 0.0f;
			M.m[0][2] = 0.0f;
			M.m[0][3] = 0.0f;

			M.m[1][0] = 0.0f;
			M.m[1][1] = twoNearZ / viewHeight;
			M.m[1][2] = 0.0f;
			M.m[1][3] = 0.0f;

			M.m[2][0] = 0.0f;
			M.m[2][1] = 0.0f;
			M.m[2][2] = fRange;
			M.m[2][3] = -1.0f;

			M.m[3][0] = 0.0f;
			M.m[3][1] = 0.0f;
			M.m[3][2] = fRange * nearZ;
			M.m[3][3] = 0.0f;
			
			return M;

#elif defined(_SSE2_INTRINSICS_)
			MATRIX M;
			float twoNearZ = nearZ + nearZ;
			float fRange = farZ / (nearZ - farZ);
			// Note: This is recorded on the stack
			VECTOR rMem = {
				twoNearZ / viewWidth,
				twoNearZ / viewHeight,
				fRange,
				fRange * nearZ
			};
			// Copy from memory to SSE register
			VECTOR vValues = rMem;
			VECTOR vTemp = _mm_setzero_ps();
			// Copy x only
			vTemp = _mm_move_ss(vTemp, vValues);
			// TwoNearZ / ViewWidth,0,0,0
			M.r[0] = vTemp;
			// 0,TwoNearZ / ViewHeight,0,0
			vTemp = vValues;
			vTemp = _mm_and_ps(vTemp, g_MaskY);
			M.r[1] = vTemp;
			// x=fRange,y=-fRange * NearZ,0,-1.0f
			vValues = _mm_shuffle_ps(vValues, g_NegIdentityR3, _MM_SHUFFLE(3, 2, 3, 2));
			// 0,0,fRange,-1.0f
			vTemp = _mm_setzero_ps();
			vTemp = _mm_shuffle_ps(vTemp, vValues, _MM_SHUFFLE(3, 0, 0, 0));
			M.r[2] = vTemp;
			// 0,0,-fRange * NearZ,0
			vTemp = _mm_shuffle_ps(vTemp, vValues, _MM_SHUFFLE(2, 1, 0, 0));
			M.r[3] = vTemp;
			
			return M;
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV PerspectiveFovLH(
			float fovAngleY, 
			float aspectRatio, 
			float nearZ, 
			float farZ
		) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(nearZ > 0.f && farZ > 0.f);
			assert(!ScalarNearEqual(fovAngleY, 0.0f, 0.00001f * 2.0f));
			assert(!ScalarNearEqual(aspectRatio, 0.0f, 0.00001f));
			assert(!ScalarNearEqual(farZ, nearZ, 0.00001f));
#endif

#if defined(_NO_INTRINSICS_)

			float    SinFov;
			float    CosFov;
			ScalarSinCos(&SinFov, &CosFov, 0.5f * fovAngleY);

			float Height = CosFov / SinFov;
			float Width = Height / aspectRatio;
			float fRange = farZ / (farZ - nearZ);

			MATRIX M;
			M.m[0][0] = Width;
			M.m[0][1] = 0.0f;
			M.m[0][2] = 0.0f;
			M.m[0][3] = 0.0f;

			M.m[1][0] = 0.0f;
			M.m[1][1] = Height;
			M.m[1][2] = 0.0f;
			M.m[1][3] = 0.0f;

			M.m[2][0] = 0.0f;
			M.m[2][1] = 0.0f;
			M.m[2][2] = fRange;
			M.m[2][3] = 1.0f;

			M.m[3][0] = 0.0f;
			M.m[3][1] = 0.0f;
			M.m[3][2] = -fRange * nearZ;
			M.m[3][3] = 0.0f;
			
			return M;

#elif defined(_SSE2_INTRINSICS_)
			float    SinFov;
			float    CosFov;
			ScalarSinCos(&SinFov, &CosFov, 0.5f * fovAngleY);

			float fRange = farZ / (farZ - nearZ);
			// Note: This is recorded on the stack
			float Height = CosFov / SinFov;
			VECTOR rMem = {
				Height / aspectRatio,
				Height,
				fRange,
				-fRange * nearZ
			};
			// Copy from memory to SSE register
			VECTOR vValues = rMem;
			VECTOR vTemp = _mm_setzero_ps();
			// Copy x only
			vTemp = _mm_move_ss(vTemp, vValues);
			// Height / AspectRatio,0,0,0
			MATRIX M;
			M.r[0] = vTemp;
			// 0,Height,0,0
			vTemp = vValues;
			vTemp = _mm_and_ps(vTemp, g_MaskY);
			M.r[1] = vTemp;
			// x=fRange,y=-fRange * NearZ,0,1.0f
			vTemp = _mm_setzero_ps();
			vValues = _mm_shuffle_ps(vValues, g_IdentityR3, _MM_SHUFFLE(3, 2, 3, 2));
			// 0,0,fRange,1.0f
			vTemp = _mm_shuffle_ps(vTemp, vValues, _MM_SHUFFLE(3, 0, 0, 0));
			M.r[2] = vTemp;
			// 0,0,-fRange * NearZ,0.0f
			vTemp = _mm_shuffle_ps(vTemp, vValues, _MM_SHUFFLE(2, 1, 0, 0));
			M.r[3] = vTemp;
			
			return M;
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV PerspectiveFovRH(
			float fovAngleY, 
			float aspectRatio, 
			float nearZ, 
			float farZ
		) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(nearZ > 0.f && farZ > 0.f);
			assert(!ScalarNearEqual(fovAngleY, 0.0f, 0.00001f * 2.0f));
			assert(!XMScalarNearEqual(aspectRatio, 0.0f, 0.00001f));
			assert(!XMScalarNearEqual(farZ, nearZ, 0.00001f));
#endif

#if defined(_NO_INTRINSICS_)
			float    SinFov;
			float    CosFov;
			ScalarSinCos(&SinFov, &CosFov, 0.5f * fovAngleY);

			float Height = CosFov / SinFov;
			float Width = Height / aspectRatio;
			float fRange = farZ / (nearZ - farZ);

			MATRIX M;
			M.m[0][0] = Width;
			M.m[0][1] = 0.0f;
			M.m[0][2] = 0.0f;
			M.m[0][3] = 0.0f;

			M.m[1][0] = 0.0f;
			M.m[1][1] = Height;
			M.m[1][2] = 0.0f;
			M.m[1][3] = 0.0f;

			M.m[2][0] = 0.0f;
			M.m[2][1] = 0.0f;
			M.m[2][2] = fRange;
			M.m[2][3] = -1.0f;

			M.m[3][0] = 0.0f;
			M.m[3][1] = 0.0f;
			M.m[3][2] = fRange * nearZ;
			M.m[3][3] = 0.0f;
			
			return M;

#elif defined(_SSE2_INTRINSICS_)
			float    SinFov;
			float    CosFov;
			ScalarSinCos(&SinFov, &CosFov, 0.5f * fovAngleY);
			float fRange = farZ / (nearZ - farZ);
			// Note: This is recorded on the stack
			float Height = CosFov / SinFov;
			VECTOR rMem = {
				Height / aspectRatio,
				Height,
				fRange,
				fRange * nearZ
			};
			// Copy from memory to SSE register
			VECTOR vValues = rMem;
			VECTOR vTemp = _mm_setzero_ps();
			// Copy x only
			vTemp = _mm_move_ss(vTemp, vValues);
			// Height / AspectRatio,0,0,0
			MATRIX M;
			M.r[0] = vTemp;
			// 0,Height,0,0
			vTemp = vValues;
			vTemp = _mm_and_ps(vTemp, g_MaskY);
			M.r[1] = vTemp;
			// x=fRange,y=-fRange * NearZ,0,-1.0f
			vTemp = _mm_setzero_ps();
			vValues = _mm_shuffle_ps(vValues, g_NegIdentityR3, _MM_SHUFFLE(3, 2, 3, 2));
			// 0,0,fRange,-1.0f
			vTemp = _mm_shuffle_ps(vTemp, vValues, _MM_SHUFFLE(3, 0, 0, 0));
			M.r[2] = vTemp;
			// 0,0,fRange * NearZ,0.0f
			vTemp = _mm_shuffle_ps(vTemp, vValues, _MM_SHUFFLE(2, 1, 0, 0));
			M.r[3] = vTemp;
			
			return M;
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV PerspectiveOffCenterLH(
			float viewLeft, 
			float viewRight, 
			float viewBottom, 
			float viewTop, 
			float nearZ, 
			float farZ
		) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(nearZ > 0.f && farZ > 0.f);
			assert(!ScalarNearEqual(viewRight, viewLeft, 0.00001f));
			assert(!ScalarNearEqual(viewTop, viewBottom, 0.00001f));
			assert(!ScalarNearEqual(farZ, nearZ, 0.00001f));
#endif

#if defined(_NO_INTRINSICS_)
			float TwoNearZ = nearZ + nearZ;
			float ReciprocalWidth = 1.0f / (viewRight - viewLeft);
			float ReciprocalHeight = 1.0f / (viewTop - viewBottom);
			float fRange = farZ / (farZ - nearZ);

			MATRIX M;
			M.m[0][0] = TwoNearZ * ReciprocalWidth;
			M.m[0][1] = 0.0f;
			M.m[0][2] = 0.0f;
			M.m[0][3] = 0.0f;

			M.m[1][0] = 0.0f;
			M.m[1][1] = TwoNearZ * ReciprocalHeight;
			M.m[1][2] = 0.0f;
			M.m[1][3] = 0.0f;

			M.m[2][0] = -(viewLeft + viewRight) * ReciprocalWidth;
			M.m[2][1] = -(viewTop + viewBottom) * ReciprocalHeight;
			M.m[2][2] = fRange;
			M.m[2][3] = 1.0f;

			M.m[3][0] = 0.0f;
			M.m[3][1] = 0.0f;
			M.m[3][2] = -fRange * nearZ;
			M.m[3][3] = 0.0f;
			
			return M;

#elif defined(_SSE2_INTRINSICS_)
			MATRIX M;
			float TwoNearZ = nearZ + nearZ;
			float ReciprocalWidth = 1.0f / (viewRight - viewLeft);
			float ReciprocalHeight = 1.0f / (viewTop - viewBottom);
			float fRange = farZ / (farZ - nearZ);
			// Note: This is recorded on the stack
			VECTOR rMem = {
				TwoNearZ * ReciprocalWidth,
				TwoNearZ * ReciprocalHeight,
				-fRange * nearZ,
				0
			};
			// Copy from memory to SSE register
			VECTOR vValues = rMem;
			VECTOR vTemp = _mm_setzero_ps();
			// Copy x only
			vTemp = _mm_move_ss(vTemp, vValues);
			// TwoNearZ*ReciprocalWidth,0,0,0
			M.r[0] = vTemp;
			// 0,TwoNearZ*ReciprocalHeight,0,0
			vTemp = vValues;
			vTemp = _mm_and_ps(vTemp, g_MaskY);
			M.r[1] = vTemp;
			// 0,0,fRange,1.0f
			M.r[2] = VEC::Set(-(viewLeft + viewRight) * ReciprocalWidth,
				-(viewTop + viewBottom) * ReciprocalHeight,
				fRange,
				1.0f);
			// 0,0,-fRange * NearZ,0.0f
			vValues = _mm_and_ps(vValues, g_MaskZ);
			M.r[3] = vValues;
			
			return M;
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV PerspectiveOffCenterRH(
			float viewLeft, 
			float viewRight, 
			float viewBottom, 
			float viewTop, 
			float nearZ, 
			float farZ
		) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(nearZ > 0.f && farZ > 0.f);
			assert(!ScalarNearEqual(viewRight, viewLeft, 0.00001f));
			assert(!ScalarNearEqual(viewTop, viewBottom, 0.00001f));
			assert(!ScalarNearEqual(farZ, nearZ, 0.00001f));
#endif

#if defined(_NO_INTRINSICS_)
			float TwoNearZ = nearZ + nearZ;
			float ReciprocalWidth = 1.0f / (viewRight - viewLeft);
			float ReciprocalHeight = 1.0f / (viewTop - viewBottom);
			float fRange = farZ / (nearZ - farZ);

			MATRIX M;
			M.m[0][0] = TwoNearZ * ReciprocalWidth;
			M.m[0][1] = 0.0f;
			M.m[0][2] = 0.0f;
			M.m[0][3] = 0.0f;

			M.m[1][0] = 0.0f;
			M.m[1][1] = TwoNearZ * ReciprocalHeight;
			M.m[1][2] = 0.0f;
			M.m[1][3] = 0.0f;

			M.m[2][0] = (viewLeft + viewRight) * ReciprocalWidth;
			M.m[2][1] = (viewTop + viewBottom) * ReciprocalHeight;
			M.m[2][2] = fRange;
			M.m[2][3] = -1.0f;

			M.m[3][0] = 0.0f;
			M.m[3][1] = 0.0f;
			M.m[3][2] = fRange * nearZ;
			M.m[3][3] = 0.0f;
			
			return M;

#elif defined(_SSE2_INTRINSICS_)
			MATRIX M;
			float TwoNearZ = nearZ + nearZ;
			float ReciprocalWidth = 1.0f / (viewRight - viewLeft);
			float ReciprocalHeight = 1.0f / (viewTop - viewBottom);
			float fRange = farZ / (nearZ - farZ);
			// Note: This is recorded on the stack
			VECTOR rMem = {
				TwoNearZ * ReciprocalWidth,
				TwoNearZ * ReciprocalHeight,
				fRange * nearZ,
				0
			};
			// Copy from memory to SSE register
			VECTOR vValues = rMem;
			VECTOR vTemp = _mm_setzero_ps();
			// Copy x only
			vTemp = _mm_move_ss(vTemp, vValues);
			// TwoNearZ*ReciprocalWidth,0,0,0
			M.r[0] = vTemp;
			// 0,TwoNearZ*ReciprocalHeight,0,0
			vTemp = vValues;
			vTemp = _mm_and_ps(vTemp, g_MaskY);
			M.r[1] = vTemp;
			// 0,0,fRange,1.0f
			M.r[2] = VEC::Set((viewLeft + viewRight) * ReciprocalWidth,
				(viewTop + viewBottom) * ReciprocalHeight,
				fRange,
				-1.0f);
			// 0,0,-fRange * NearZ,0.0f
			vValues = _mm_and_ps(vValues, g_MaskZ);
			M.r[3] = vValues;
			
			return M;
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV OrthographicLH(
			float viewWidth, 
			float viewHeight, 
			float nearZ, 
			float farZ
		) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(!ScalarNearEqual(viewWidth, 0.0f, 0.00001f));
			assert(!ScalarNearEqual(viewHeight, 0.0f, 0.00001f));
			assert(!ScalarNearEqual(farZ, nearZ, 0.00001f));
#endif

#if defined(_NO_INTRINSICS_)
			float fRange = 1.0f / (farZ - nearZ);

			MATRIX M;
			M.m[0][0] = 2.0f / viewWidth;
			M.m[0][1] = 0.0f;
			M.m[0][2] = 0.0f;
			M.m[0][3] = 0.0f;

			M.m[1][0] = 0.0f;
			M.m[1][1] = 2.0f / viewHeight;
			M.m[1][2] = 0.0f;
			M.m[1][3] = 0.0f;

			M.m[2][0] = 0.0f;
			M.m[2][1] = 0.0f;
			M.m[2][2] = fRange;
			M.m[2][3] = 0.0f;

			M.m[3][0] = 0.0f;
			M.m[3][1] = 0.0f;
			M.m[3][2] = -fRange * nearZ;
			M.m[3][3] = 1.0f;
			
			return M;

#elif defined(_SSE2_INTRINSICS_)
			MATRIX M;
			float fRange = 1.0f / (farZ - nearZ);
			// Note: This is recorded on the stack
			VECTOR rMem = {
				2.0f / viewWidth,
				2.0f / viewHeight,
				fRange,
				-fRange * nearZ
			};
			// Copy from memory to SSE register
			VECTOR vValues = rMem;
			VECTOR vTemp = _mm_setzero_ps();
			// Copy x only
			vTemp = _mm_move_ss(vTemp, vValues);
			// 2.0f / ViewWidth,0,0,0
			M.r[0] = vTemp;
			// 0,2.0f / ViewHeight,0,0
			vTemp = vValues;
			vTemp = _mm_and_ps(vTemp, g_MaskY);
			M.r[1] = vTemp;
			// x=fRange,y=-fRange * NearZ,0,1.0f
			vTemp = _mm_setzero_ps();
			vValues = _mm_shuffle_ps(vValues, g_IdentityR3, _MM_SHUFFLE(3, 2, 3, 2));
			// 0,0,fRange,0.0f
			vTemp = _mm_shuffle_ps(vTemp, vValues, _MM_SHUFFLE(2, 0, 0, 0));
			M.r[2] = vTemp;
			// 0,0,-fRange * NearZ,1.0f
			vTemp = _mm_shuffle_ps(vTemp, vValues, _MM_SHUFFLE(3, 1, 0, 0));
			M.r[3] = vTemp;
			
			return M;
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV OrthographicRH(
			float viewWidth, 
			float viewHeight, 
			float nearZ, 
			float farZ
		) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(!ScalarNearEqual(viewWidth, 0.0f, 0.00001f));
			assert(!ScalarNearEqual(viewHeight, 0.0f, 0.00001f));
			assert(!ScalarNearEqual(farZ, nearZ, 0.00001f));
#endif

#if defined(_NO_INTRINSICS_)
			float fRange = 1.0f / (nearZ - farZ);

			MATRIX M;
			M.m[0][0] = 2.0f / viewWidth;
			M.m[0][1] = 0.0f;
			M.m[0][2] = 0.0f;
			M.m[0][3] = 0.0f;

			M.m[1][0] = 0.0f;
			M.m[1][1] = 2.0f / viewHeight;
			M.m[1][2] = 0.0f;
			M.m[1][3] = 0.0f;

			M.m[2][0] = 0.0f;
			M.m[2][1] = 0.0f;
			M.m[2][2] = fRange;
			M.m[2][3] = 0.0f;

			M.m[3][0] = 0.0f;
			M.m[3][1] = 0.0f;
			M.m[3][2] = fRange * nearZ;
			M.m[3][3] = 1.0f;
			
			return M;

#elif defined(_SSE2_INTRINSICS_)
			MATRIX M;
			float fRange = 1.0f / (nearZ - farZ);
			// Note: This is recorded on the stack
			VECTOR rMem = {
				2.0f / viewWidth,
				2.0f / viewHeight,
				fRange,
				fRange * nearZ
			};
			// Copy from memory to SSE register
			VECTOR vValues = rMem;
			VECTOR vTemp = _mm_setzero_ps();
			// Copy x only
			vTemp = _mm_move_ss(vTemp, vValues);
			// 2.0f / ViewWidth,0,0,0
			M.r[0] = vTemp;
			// 0,2.0f / ViewHeight,0,0
			vTemp = vValues;
			vTemp = _mm_and_ps(vTemp, g_MaskY);
			M.r[1] = vTemp;
			// x=fRange,y=fRange * NearZ,0,1.0f
			vTemp = _mm_setzero_ps();
			vValues = _mm_shuffle_ps(vValues, g_IdentityR3, _MM_SHUFFLE(3, 2, 3, 2));
			// 0,0,fRange,0.0f
			vTemp = _mm_shuffle_ps(vTemp, vValues, _MM_SHUFFLE(2, 0, 0, 0));
			M.r[2] = vTemp;
			// 0,0,fRange * NearZ,1.0f
			vTemp = _mm_shuffle_ps(vTemp, vValues, _MM_SHUFFLE(3, 1, 0, 0));
			M.r[3] = vTemp;
			
			return M;
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV OrthographicOffCenterLH(
			float viewLeft, 
			float viewRight, 
			float viewBottom, 
			float viewTop, 
			float nearZ, 
			float farZ
		) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(!ScalarNearEqual(viewRight, viewLeft, 0.00001f));
			assert(!ScalarNearEqual(viewTop, viewBottom, 0.00001f));
			assert(!ScalarNearEqual(farZ, nearZ, 0.00001f));
#endif

#if defined(_NO_INTRINSICS_)
			float ReciprocalWidth = 1.0f / (viewRight - viewLeft);
			float ReciprocalHeight = 1.0f / (viewTop - viewBottom);
			float fRange = 1.0f / (farZ - nearZ);

			MATRIX M;
			M.m[0][0] = ReciprocalWidth + ReciprocalWidth;
			M.m[0][1] = 0.0f;
			M.m[0][2] = 0.0f;
			M.m[0][3] = 0.0f;

			M.m[1][0] = 0.0f;
			M.m[1][1] = ReciprocalHeight + ReciprocalHeight;
			M.m[1][2] = 0.0f;
			M.m[1][3] = 0.0f;

			M.m[2][0] = 0.0f;
			M.m[2][1] = 0.0f;
			M.m[2][2] = fRange;
			M.m[2][3] = 0.0f;

			M.m[3][0] = -(viewLeft + viewRight) * ReciprocalWidth;
			M.m[3][1] = -(viewTop + viewBottom) * ReciprocalHeight;
			M.m[3][2] = -fRange * nearZ;
			M.m[3][3] = 1.0f;
			
			return M;

#elif defined(_SSE2_INTRINSICS_)
			MATRIX M;
			float fReciprocalWidth = 1.0f / (viewRight - viewLeft);
			float fReciprocalHeight = 1.0f / (viewTop - viewBottom);
			float fRange = 1.0f / (farZ - nearZ);
			// Note: This is recorded on the stack
			VECTOR rMem = {
				fReciprocalWidth,
				fReciprocalHeight,
				fRange,
				1.0f
			};
			VECTOR rMem2 = {
				-(viewLeft + viewRight),
				-(viewTop + viewBottom),
				-nearZ,
				1.0f
			};
			// Copy from memory to SSE register
			VECTOR vValues = rMem;
			VECTOR vTemp = _mm_setzero_ps();
			// Copy x only
			vTemp = _mm_move_ss(vTemp, vValues);
			// fReciprocalWidth*2,0,0,0
			vTemp = _mm_add_ss(vTemp, vTemp);
			M.r[0] = vTemp;
			// 0,fReciprocalHeight*2,0,0
			vTemp = vValues;
			vTemp = _mm_and_ps(vTemp, g_MaskY);
			vTemp = _mm_add_ps(vTemp, vTemp);
			M.r[1] = vTemp;
			// 0,0,fRange,0.0f
			vTemp = vValues;
			vTemp = _mm_and_ps(vTemp, g_MaskZ);
			M.r[2] = vTemp;
			// -(ViewLeft + ViewRight)*fReciprocalWidth,-(ViewTop + ViewBottom)*fReciprocalHeight,fRange*-NearZ,1.0f
			vValues = _mm_mul_ps(vValues, rMem2);
			M.r[3] = vValues;
			
			return M;
#endif
		}

		FORCE_INLINE MATRIX VEC_CALLCONV OrthographicOffCenterRH(
			float viewLeft, 
			float viewRight, 
			float viewBottom, 
			float viewTop, 
			float nearZ, 
			float farZ
		) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(!ScalarNearEqual(viewRight, viewLeft, 0.00001f));
			assert(!ScalarNearEqual(viewTop, viewBottom, 0.00001f));
			assert(!ScalarNearEqual(farZ, nearZ, 0.00001f));
#endif

#if defined(_NO_INTRINSICS_)
			float ReciprocalWidth = 1.0f / (viewRight - viewLeft);
			float ReciprocalHeight = 1.0f / (viewTop - viewBottom);
			float fRange = 1.0f / (nearZ - farZ);

			MATRIX M;
			M.m[0][0] = ReciprocalWidth + ReciprocalWidth;
			M.m[0][1] = 0.0f;
			M.m[0][2] = 0.0f;
			M.m[0][3] = 0.0f;

			M.m[1][0] = 0.0f;
			M.m[1][1] = ReciprocalHeight + ReciprocalHeight;
			M.m[1][2] = 0.0f;
			M.m[1][3] = 0.0f;

			M.m[2][0] = 0.0f;
			M.m[2][1] = 0.0f;
			M.m[2][2] = fRange;
			M.m[2][3] = 0.0f;

			M.r[3] = VEC::Set(-(viewLeft + viewRight) * ReciprocalWidth,
				-(viewTop + viewBottom) * ReciprocalHeight,
				fRange * nearZ,
				1.0f);
			
			return M;
#elif defined(_SSE2_INTRINSICS_)
			MATRIX M;
			float fReciprocalWidth = 1.0f / (viewRight - viewLeft);
			float fReciprocalHeight = 1.0f / (viewTop - viewBottom);
			float fRange = 1.0f / (nearZ - farZ);
			// Note: This is recorded on the stack
			VECTOR rMem = {
				fReciprocalWidth,
				fReciprocalHeight,
				fRange,
				1.0f
			};
			VECTOR rMem2 = {
				-(viewLeft + viewRight),
				-(viewTop + viewBottom),
				nearZ,
				1.0f
			};
			// Copy from memory to SSE register
			VECTOR vValues = rMem;
			VECTOR vTemp = _mm_setzero_ps();
			// Copy x only
			vTemp = _mm_move_ss(vTemp, vValues);
			// fReciprocalWidth*2,0,0,0
			vTemp = _mm_add_ss(vTemp, vTemp);
			M.r[0] = vTemp;
			// 0,fReciprocalHeight*2,0,0
			vTemp = vValues;
			vTemp = _mm_and_ps(vTemp, g_MaskY);
			vTemp = _mm_add_ps(vTemp, vTemp);
			M.r[1] = vTemp;
			// 0,0,fRange,0.0f
			vTemp = vValues;
			vTemp = _mm_and_ps(vTemp, g_MaskZ);
			M.r[2] = vTemp;
			// -(ViewLeft + ViewRight)*fReciprocalWidth,-(ViewTop + ViewBottom)*fReciprocalHeight,fRange*-NearZ,1.0f
			vValues = _mm_mul_ps(vValues, rMem2);
			M.r[3] = vValues;
			
			return M;
#endif
		}

#ifdef _PREFAST_
#pragma prefast(pop)
#endif

		
	}
}

#endif // !ULTREALITY_MATH_SSE2_MATRIX_INL
