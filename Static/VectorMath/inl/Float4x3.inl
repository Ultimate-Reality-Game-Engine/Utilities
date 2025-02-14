#ifndef ULTREALITY_MATH_FLOAT4X3_INL
#define ULTREALITY_MATH_FLOAT4X3_INL

#if defined(__GNUC__) or defined(__clang__)
#define FORCE_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline
#endif

namespace UltReality::Math
{
	_Use_decl_annotations_
	FORCE_INLINE Float4x3::Float4x3(const float* pArray) noexcept
	{
#if defined(DEBUG) || defined(_DEBUG)
		assert(pArray != nullptr);
#endif

		m[0][0] = pArray[0];
		m[0][1] = pArray[1];
		m[0][2] = pArray[2];

		m[1][0] = pArray[3];
		m[1][1] = pArray[4];
		m[1][2] = pArray[5];

		m[2][0] = pArray[6];
		m[2][1] = pArray[7];
		m[2][2] = pArray[8];

		m[3][0] = pArray[9];
		m[3][1] = pArray[10];
		m[3][2] = pArray[11];
	}

	_Use_decl_annotations_
	FORCE_INLINE Float4x3::Float4x3(A_MATRIX m) noexcept
	{
		Matrix::StoreFloat4x3(this, m);
	}

	_Use_decl_annotations_
	FORCE_INLINE Float4x3& VEC_CALLCONV Float4x3::operator=(A_MATRIX m) noexcept
	{
		Matrix::StoreFloat4x3(this, m);

		return *this;
	}

	FORCE_INLINE MATRIX VEC_CALLCONV Float4x3::Load() noexcept
	{
		return Matrix::LoadFloat4x3(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV Float4x3::Store(A_MATRIX m) noexcept
	{
		Matrix::StoreFloat4x3(this, m);
	}

	_Use_decl_annotations_
	FORCE_INLINE AFloat4x3::AFloat4x3(A_MATRIX m) noexcept
	{
		Matrix::StoreAFloat4x3(this, m);
	}

	_Use_decl_annotations_
	FORCE_INLINE AFloat4x3& VEC_CALLCONV AFloat4x3::operator=(A_MATRIX m) noexcept
	{
		Matrix::StoreAFloat4x3(this, m);

		return *this;
	}

	FORCE_INLINE MATRIX VEC_CALLCONV AFloat4x3::Load() noexcept
	{
		return Matrix::LoadAFloat4x3(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV AFloat4x3::Store(A_MATRIX m) noexcept
	{
		Matrix::StoreAFloat4x3(this, m);
	}

	namespace Matrix
	{
		_Use_decl_annotations_
		FORCE_INLINE MATRIX VEC_CALLCONV LoadFloat4x3(const Float4x3* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			MATRIX M;

			M.r[0].vector4_f32[0] = pSource->m[0][0];
			M.r[0].vector4_f32[1] = pSource->m[0][1];
			M.r[0].vector4_f32[2] = pSource->m[0][2];
			M.r[0].vector4_f32[3] = 0.0f;

			M.r[1].vector4_f32[0] = pSource->m[1][0];
			M.r[1].vector4_f32[1] = pSource->m[1][1];
			M.r[1].vector4_f32[2] = pSource->m[1][2];
			M.r[1].vector4_f32[3] = 0.0f;

			M.r[2].vector4_f32[0] = pSource->m[2][0];
			M.r[2].vector4_f32[1] = pSource->m[2][1];
			M.r[2].vector4_f32[2] = pSource->m[2][2];
			M.r[2].vector4_f32[3] = 0.0f;

			M.r[3].vector4_f32[0] = pSource->m[3][0];
			M.r[3].vector4_f32[1] = pSource->m[3][1];
			M.r[3].vector4_f32[2] = pSource->m[3][2];
			M.r[3].vector4_f32[3] = 1.0f;

			return M;
#elif defined(_SSE2_INTRINSICS_)
			// Use unaligned loads to load the 12 floats
			// V1 = x1, y1, z1, x2
			__m128 V1 = _mm_loadu_ps(&pSource->m[0][0]);

			// V2 = y2, z2, x3, y3
			__m128 V2 = _mm_loadu_ps(&pSource->m[1][1]);

			// V4 = z3, x4, y4, z4
			__m128 V4 = _mm_loadu_ps(&pSource->m[2][2]);

			// V3 = x3, y3, z3, z3
			__m128 V3 = _mm_shuffle_ps(V2, V4, _MM_SHUFFLE(0, 0, 3, 2));

			// V2 = y2, z2, x2, x2
			V2 = _mm_shuffle_ps(V2, V1, _MM_SHUFFLE(3, 3, 1, 0));
			// V2 = x2, y2, z2, z2
			V2 = PERMUTE_PS(V2, _MM_SHUFFLE(1, 1, 0, 2));

			// V1 = x1, y1, z1, 0
			V1 = _mm_and_ps(V1, g_Mask3);

			// V2 = x2, y2, z2, 0
			V2 = _mm_and_ps(V2, g_Mask3);

			// V3 = x3, y3, z3, 0
			V3 = _mm_and_ps(V3, g_Mask3);

			// V4i = x4, y4, z4, 0
			__m128i V4i = _mm_srli_si128(_mm_castps_si128(V4), 32 / 8);

			// V4i = x4, y4, z4, 1.0f
			V4i = _mm_or_si128(V4i, g_IdentityR3);

			MATRIX M(V1, V2, V3, _mm_castsi128_ps(V4i));

			return M;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE MATRIX VEC_CALLCONV LoadAFloat4x3(const AFloat4x3* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
			assert((reinterpret_cast<uintptr_t>(pSource) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			MATRIX M;

			M.r[0].vector4_f32[0] = pSource->m[0][0];
			M.r[0].vector4_f32[1] = pSource->m[0][1];
			M.r[0].vector4_f32[2] = pSource->m[0][2];
			M.r[0].vector4_f32[3] = 0.0f;

			M.r[1].vector4_f32[0] = pSource->m[1][0];
			M.r[1].vector4_f32[1] = pSource->m[1][1];
			M.r[1].vector4_f32[2] = pSource->m[1][2];
			M.r[1].vector4_f32[3] = 0.0f;

			M.r[2].vector4_f32[0] = pSource->m[2][0];
			M.r[2].vector4_f32[1] = pSource->m[2][1];
			M.r[2].vector4_f32[2] = pSource->m[2][2];
			M.r[2].vector4_f32[3] = 0.0f;

			M.r[3].vector4_f32[0] = pSource->m[3][0];
			M.r[3].vector4_f32[1] = pSource->m[3][1];
			M.r[3].vector4_f32[2] = pSource->m[3][2];
			M.r[3].vector4_f32[3] = 1.0f;

			return M;
#elif defined(_SSE2_INTRINSICS_)
			// Use aligned loads to load the 12 floats
			// V1 = x1, y1, z1, x2
			__m128 V1 = _mm_load_ps(&pSource->m[0][0]);

			// V2 = y2, z2, x3, y3
			__m128 V2 = _mm_load_ps(&pSource->m[1][1]);

			// V4 = z3, x4, y4, z4
			__m128 V4 = _mm_load_ps(&pSource->m[2][2]);

			// V3 = x3, y3, z3, z3
			__m128 V3 = _mm_shuffle_ps(V2, V4, _MM_SHUFFLE(0, 0, 3, 2));

			// V2 = y2, z2, x2, x2
			V2 = _mm_shuffle_ps(V2, V1, _MM_SHUFFLE(3, 3, 1, 0));
			// V2 = x2, y2, z2, z2
			V2 = PERMUTE_PS(V2, _MM_SHUFFLE(1, 1, 0, 2));

			// V1 = x1, y1, z1, 0
			V1 = _mm_and_ps(V1, g_Mask3);

			// V2 = x2, y2, z2, 0
			V2 = _mm_and_ps(V2, g_Mask3);

			// V3 = x3, y3, z3, 0
			V3 = _mm_and_ps(V3, g_Mask3);

			// V4i = x4, y4, z4, 0
			__m128i V4i = _mm_srli_si128(_mm_castps_si128(V4), 32 / 8);

			// V4i = x4, y4, z4, 1.0f
			V4i = _mm_or_si128(V4i, g_IdentityR3);

			MATRIX M(V1, V2, V3, _mm_castsi128_ps(V4i));

			return M;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreFloat4x3(Float4x3* pDestination, A_MATRIX m) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->m[0][0] = m.r[0].vector4_f32[0];
			pDestination->m[0][1] = m.r[0].vector4_f32[1];
			pDestination->m[0][2] = m.r[0].vector4_f32[2];
									
			pDestination->m[1][0] = m.r[1].vector4_f32[0];
			pDestination->m[1][1] = m.r[1].vector4_f32[1];
			pDestination->m[1][2] = m.r[1].vector4_f32[2];
									
			pDestination->m[2][0] = m.r[2].vector4_f32[0];
			pDestination->m[2][1] = m.r[2].vector4_f32[1];
			pDestination->m[2][2] = m.r[2].vector4_f32[2];
									
			pDestination->m[3][0] = m.r[3].vector4_f32[0];
			pDestination->m[3][1] = m.r[3].vector4_f32[1];
			pDestination->m[3][2] = m.r[3].vector4_f32[2];

#elif defined(_SSE2_INTRINSICS_)
			// x1, y1, z1, w1
			VECTOR V1 = m.r[0];
			// x2, y2, z2, w2
			VECTOR V2 = m.r[1];
			// x3, y3, z3, w3
			VECTOR V3 = m.r[2];
			// x4, y4, z4, w4
			VECTOR V4 = m.r[3];
			
			// y2, z2, x3, y3
			VECTOR vTemp2x = _mm_shuffle_ps(V2, V3, _MM_SHUFFLE(1, 0, 2, 1));

			// y2, z2, x3, y3
			V2 = _mm_shuffle_ps(V2, V1, _MM_SHUFFLE(2, 2, 0, 0));
			// x1, y1, z1, x2
			V1 = _mm_shuffle_ps(V1, V2, _MM_SHUFFLE(0, 2, 1, 0));
			// z3, z3, x4, x4
			V3 = _mm_shuffle_ps(V3, V4, _MM_SHUFFLE(0, 0, 2, 2));
			// z3, x4, y4, z4 (final)
			V3 = _mm_shuffle_ps(V3, V4, _MM_SHUFFLE(2, 1, 2, 0));

			_mm_storeu_ps(&pDestination->m[0][0], V1);
			_mm_storeu_ps(&pDestination->m[1][1], vTemp2x);
			_mm_storeu_ps(&pDestination->m[2][2], V3);

			// x1, y1, z1
			// x2, y2, z2
			// x3, y3, z3
			// x4, y4, z4
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreAFloat4x3(AFloat4x3* pDestination, A_MATRIX m) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
			assert((reinterpret_cast<uintptr_t>(pDestination) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->m[0][0] = m.r[0].vector4_f32[0];
			pDestination->m[0][1] = m.r[0].vector4_f32[1];
			pDestination->m[0][2] = m.r[0].vector4_f32[2];

			pDestination->m[1][0] = m.r[1].vector4_f32[0];
			pDestination->m[1][1] = m.r[1].vector4_f32[1];
			pDestination->m[1][2] = m.r[1].vector4_f32[2];

			pDestination->m[2][0] = m.r[2].vector4_f32[0];
			pDestination->m[2][1] = m.r[2].vector4_f32[1];
			pDestination->m[2][2] = m.r[2].vector4_f32[2];

			pDestination->m[3][0] = m.r[3].vector4_f32[0];
			pDestination->m[3][1] = m.r[3].vector4_f32[1];
			pDestination->m[3][2] = m.r[3].vector4_f32[2];

#elif defined(_SSE2_INTRINSICS_)
			// x1, y1, z1, w1
			VECTOR V1 = m.r[0];
			// x2, y2, z2, w2
			VECTOR V2 = m.r[1];
			// x3, y3, z3, w3
			VECTOR V3 = m.r[2];
			// x4, y4, z4, w4
			VECTOR V4 = m.r[3];
			
			// z1, z1, x2, y2
			VECTOR vTemp2x = _mm_shuffle_ps(V1, V2, _MM_SHUFFLE(1, 0, 2, 2));

			// y2, z2, x3, y3 (final)
			V2 = _mm_shuffle_ps(V2, V3, _MM_SHUFFLE(1, 0, 2, 1));
			// x1, y1, z1, x2 (final)
			V1 = _mm_shuffle_ps(V1, V2, _MM_SHUFFLE(2, 0, 1, 0));
			// z3, z3, x4, x4
			V3 = _mm_shuffle_ps(V3, V4, _MM_SHUFFLE(0, 0, 2, 2));
			// z3, x4, y4, z4 (final)
			V3 = _mm_shuffle_ps(V3, V4, _MM_SHUFFLE(2, 1, 2, 0));

			// x1, y1, z1, x2
			// y2, z2, x3, y3
			// z3, x4, y4, z4

			_mm_store_ps(&pDestination->m[0][0], V1);
			_mm_store_ps(&pDestination->m[1][1], V2);
			_mm_store_ps(&pDestination->m[2][2], V3);

			// x1, y1, z1
			// x2, y2, z2
			// x3, y3, z3
			// x4, y4, z4
#endif
		}
	}
}

#endif // !ULTREALITY_MATH_FLOAT4X3_INL
