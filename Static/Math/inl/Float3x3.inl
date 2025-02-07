#ifndef ULTREALITY_MATH_FLOAT3X3_INL
#define ULTREALITY_MATH_FLOAT3X3_INL

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
	FORCE_INLINE Float3x3::Float3x3(const float* pArray) noexcept
	{
#if defined(DEBUG) || defined(_DEBUG)
		assert(pArray != nullptr);
#endif

		for (size_t row = 0; row < 3; row++)
		{
			for (size_t column = 0; column < 3; column++)
			{
				m[row][column] = pArray[row * 3 + column];
			}
		}
	}

	_Use_decl_annotations_
	FORCE_INLINE Float3x3::Float3x3(A_MATRIX m) noexcept
	{
		MAT::StoreFloat3x3(this, m);
	}

	_Use_decl_annotations_
	FORCE_INLINE Float3x3& VEC_CALLCONV Float3x3::operator=(A_MATRIX m) noexcept
	{
		MAT::StoreFloat3x3(this, m);

		return *this;
	}

	FORCE_INLINE MATRIX VEC_CALLCONV Float3x3::Load() noexcept
	{
		return MAT::LoadFloat3x3(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV Float3x3::Store(A_MATRIX m) noexcept
	{
		MAT::StoreFloat3x3(this, m);
	}

	_Use_decl_annotations_
	FORCE_INLINE AFloat3x3::AFloat3x3(A_MATRIX m) noexcept
	{
		MAT::StoreAFloat3x3(this, m);
	}

	_Use_decl_annotations_
	FORCE_INLINE AFloat3x3& VEC_CALLCONV AFloat3x3::operator=(A_MATRIX m) noexcept
	{
		MAT::StoreAFloat3x3(this, m);

		return *this;
	}

	FORCE_INLINE MATRIX VEC_CALLCONV AFloat3x3::Load() noexcept
	{
		return MAT::LoadAFloat3x3(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV AFloat3x3::Store(A_MATRIX m) noexcept
	{
		MAT::StoreAFloat3x3(this, m);
	}

	namespace MAT
	{
		_Use_decl_annotations_
		FORCE_INLINE MATRIX VEC_CALLCONV LoadFloat3x3(const Float3x3* pSource) noexcept
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

			M.r[2].vector4_f32[0] = pSource->[2][0];
			M.r[2].vector4_f32[1] = pSource->[2][1];
			M.r[2].vector4_f32[2] = pSource->[2][2];
			M.r[2].vector4_f32[3] = 0.0f;

			M.r[3].vector4_f32[0] = 0.0f;
			M.r[3].vector4_f32[1] = 0.0f;
			M.r[3].vector4_f32[2] = 0.0f;
			M.r[3].vector4_f32[3] = 1.0f;

			return M;
#elif defined(_SSE2_INTRINSICS_)
			__m128 Z = _mm_setzero_ps();

			__m128 V1 = _mm_loadu_ps(&pSource->m[0][0]);
			__m128 V2 = _mm_loadu_ps(&pSource->m[1][1]);
			__m128 V3 = _mm_loadu_ps(&pSource->m[2][2]);

			__m128 T1 = _mm_unpackhi_ps(V1, Z);
			__m128 T2 = _mm_unpacklo_ps(V2, Z);
			__m128 T3 = _mm_shuffle_ps(V3, T2, _MM_SHUFFLE(0, 1, 0, 0));
			__m128 T4 = _mm_movehl_ps(T2, T3);
			__m128 T5 = _mm_movehl_ps(Z, T1);

			MATRIX M;
			M.r[0] = _mm_movelh_ps(V1, T1);
			M.r[1] = _mm_add_ps(T4, T5);
			M.r[2] = _mm_shuffle_ps(V2, V3, _MM_SHUFFLE(1, 0, 3, 2));
			M.r[3] = g_IdentityR3;

			return M;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE MATRIX VEC_CALLCONV LoadAFloat3x3(const AFloat3x3* pSource) noexcept
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

			M.r[2].vector4_f32[0] = pSource->[2][0];
			M.r[2].vector4_f32[1] = pSource->[2][1];
			M.r[2].vector4_f32[2] = pSource->[2][2];
			M.r[2].vector4_f32[3] = 0.0f;

			M.r[3].vector4_f32[0] = 0.0f;
			M.r[3].vector4_f32[1] = 0.0f;
			M.r[3].vector4_f32[2] = 0.0f;
			M.r[3].vector4_f32[3] = 1.0f;

			return M;
#elif defined(_SSE2_INTRINSICS_)
			__m128 row0 = _mm_load_ps(&pSource->m[0][0]);
			__m128 row1 = _mm_load_ps(&pSource->m[1][0]);
			__m128 row2 = _mm_load_ps(&pSource->m[2][0]);

			// mask the rows to clear the w-components
			row0 = _mm_and_ps(row0, g_Mask3);
			row1 = _mm_and_ps(row1, g_Mask3);
			row2 = _mm_and_ps(row2, g_Mask3);

			MATRIX M;
			M.r[0] = row0;
			M.r[1] = row1;
			M.r[2] = row2;
			M.r[3] = g_IdentityR3;

			return M;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreFloat3x3(Float3x3* pDestination, A_MATRIX m) noexcept
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

#elif defined(_SSE2_INTRINSICS_)
			VECTOR V1 = m.r[0];
			VECTOR V2 = m.r[1];
			VECTOR V3 = m.r[2];
			VECTOR vWork = _mm_shuffle_ps(V1, V2, _MM_SHUFFLE(0, 0, 2, 2));
			
			V1 = _mm_shuffle_ps(V1, vWork, _MM_SHUFFLE(2, 0, 1, 0));
			_mm_storeu_ps(&pDestination->m[0][0], V1);
			
			V2 = _mm_shuffle_ps(V2, V3, _MM_SHUFFLE(1, 0, 2, 1));
			_mm_storeu_ps(&pDestination->m[1][1], V2);

			V3 = PERMUTE_PS(V3, _MM_SHUFFLE(2, 2, 2, 2));
			_mm_store_ss(&pDestination->m[2][2], V3);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreAFloat3x3(AFloat3x3* pDestination, A_MATRIX m) noexcept
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

#elif defined(_SSE2_INTRINSICS_)
			VECTOR V1 = m.r[0];
			VECTOR V2 = m.r[1];
			VECTOR V3 = m.r[2];
			VECTOR vWork = _mm_shuffle_ps(V1, V2, _MM_SHUFFLE(0, 0, 2, 2));

			V1 = _mm_shuffle_ps(V1, vWork, _MM_SHUFFLE(2, 0, 1, 0));
			_mm_store_ps(&pDestination->m[0][0], V1);

			V2 = _mm_shuffle_ps(V2, V3, _MM_SHUFFLE(1, 0, 2, 1));
			_mm_store_ps(&pDestination->m[1][1], V2);

			V3 = _mm_movehl_ps(V3, V3);
			_mm_store_ss(&pDestination->m[2][2], V3);
#endif
		}
	}
}

#endif // !ULTREALITY_MATH_FLOAT3X3_INL
