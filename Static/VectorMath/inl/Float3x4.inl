#ifndef ULTREALITY_MATH_FLOAT3X4_INL
#define ULTREALITY_MATH_FLOAT3X4_INL

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
	FORCE_INLINE Float3x4::Float3x4(const float* pArray) noexcept
	{
#if defined(DEBUG) || defined(_DEBUG)
		assert(pArray != nullptr);
#endif

		m[0][0] = pArray[0];
		m[0][1] = pArray[1];
		m[0][2] = pArray[2];
		m[0][3] = pArray[3];

		m[1][0] = pArray[4];
		m[1][1] = pArray[5];
		m[1][2] = pArray[6];
		m[1][3] = pArray[7];

		m[2][0] = pArray[8];
		m[2][1] = pArray[9];
		m[2][2] = pArray[10];
		m[2][3] = pArray[11];
	}

	_Use_decl_annotations_
	FORCE_INLINE Float3x4::Float3x4(A_MATRIX m) noexcept
	{
		Matrix::StoreFloat3x4(this, m);
	}

	_Use_decl_annotations_
	FORCE_INLINE Float3x4& VEC_CALLCONV Float3x4::operator=(A_MATRIX m) noexcept
	{
		Matrix::StoreFloat3x4(this, m);

		return *this;
	}

	FORCE_INLINE MATRIX VEC_CALLCONV Float3x4::Load() noexcept
	{
		return Matrix::LoadFloat3x4(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV Float3x4::Store(A_MATRIX m) noexcept
	{
		Matrix::StoreFloat3x4(this, m);
	}

	_Use_decl_annotations_
	FORCE_INLINE AFloat3x4::AFloat3x4(A_MATRIX m) noexcept
	{
		Matrix::StoreAFloat3x4(this, m);
	}

	_Use_decl_annotations_
	FORCE_INLINE AFloat3x4& VEC_CALLCONV AFloat3x4::operator=(A_MATRIX m) noexcept
	{
		Matrix::StoreAFloat3x4(this, m);

		return *this;
	}

	FORCE_INLINE MATRIX VEC_CALLCONV AFloat3x4::Load() noexcept
	{
		return Matrix::LoadAFloat3x4(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV AFloat3x4::Store(A_MATRIX m) noexcept
	{
		Matrix::StoreAFloat3x4(this, m);
	}

	namespace Matrix
	{
		_Use_decl_annotations_
		FORCE_INLINE MATRIX VEC_CALLCONV LoadFloat3x4(const Float3x4* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			MATRIX M;
			M.r[0].vector4_f32[0] = pSource->m[0][0];
			M.r[0].vector4_f32[1] = pSource->m[1][0];
			M.r[0].vector4_f32[2] = pSource->m[2][0];
			M.r[0].vector4_f32[3] = 0.0f;

			M.r[1].vector4_f32[0] = pSource->m[0][1];
			M.r[1].vector4_f32[1] = pSource->m[1][1];
			M.r[1].vector4_f32[2] = pSource->m[2][1];
			M.r[1].vector4_f32[3] = 0.0f;

			M.r[2].vector4_f32[0] = pSource->m[0][2];
			M.r[2].vector4_f32[1] = pSource->m[1][2];
			M.r[2].vector4_f32[2] = pSource->m[2][2];
			M.r[2].vector4_f32[3] = 0.0f;

			M.r[3].vector4_f32[0] = pSource->m[0][3];
			M.r[3].vector4_f32[1] = pSource->m[1][3];
			M.r[3].vector4_f32[2] = pSource->m[2][3];
			M.r[3].vector4_f32[3] = 1.0f;

			return M;
#elif defined(_SSE2_INTRINSICS_)
			MATRIX M;
			M.r[0] = _mm_loadu_ps(&pSource->_00);
			M.r[1] = _mm_loadu_ps(&pSource->_10);
			M.r[2] = _mm_loadu_ps(&pSource->_20);
			M.r[3] = g_IdentityR3;

			// x.x, x.y, y.x, y.y
			__m128 V1 = _mm_shuffle_ps(M.r[0], M.r[1], _MM_SHUFFLE(1, 0, 1, 0));

			// x.z, x.w, y.z, y.w
			__m128 V3 = _mm_shuffle_ps(M.r[0], M.r[1], _MM_SHUFFLE(3, 2, 3, 2));

			// z.x, z.y, w.x, w.y
			__m128 V2 = _mm_shuffle_ps(M.r[2], M.r[3], _MM_SHUFFLE(1, 0, 1, 0));

			// z.z, z.w, w.z, w.w
			__m128 V4 = _mm_shuffle_ps(M.r[2], M.r[3], _MM_SHUFFLE(3, 2, 3, 2));

			MATRIX mResult;
			// x.x, y.x, z.x, w.x
			mResult.r[0] = _mm_shuffle_ps(V1, V2, _MM_SHUFFLE(2, 0, 2, 0));

			// x.y, y.y, z.y, w.y
			mResult.r[1] = _mm_shuffle_ps(V1, V2, _MM_SHUFFLE(3, 1, 3, 1));

			// x.z, y.z, z.z, w.z
			mResult.r[2] = _mm_shuffle_ps(V3, V4, _MM_SHUFFLE(2, 0, 2, 0));

			// x.w, y.w, z.w, w.w
			mResult.r[3] = _mm_shuffle_ps(V3, V4, _MM_SHUFFLE(3, 1, 3, 1));

			return mResult;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE MATRIX VEC_CALLCONV LoadAFloat3x4(const AFloat3x4* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
			assert((reinterpret_cast<uintptr_t>(pSource) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			MATRIX M;
			M.r[0].vector4_f32[0] = pSource->m[0][0];
			M.r[0].vector4_f32[1] = pSource->m[1][0];
			M.r[0].vector4_f32[2] = pSource->m[2][0];
			M.r[0].vector4_f32[3] = 0.0f;

			M.r[1].vector4_f32[0] = pSource->m[0][1];
			M.r[1].vector4_f32[1] = pSource->m[1][1];
			M.r[1].vector4_f32[2] = pSource->m[2][1];
			M.r[1].vector4_f32[3] = 0.0f;

			M.r[2].vector4_f32[0] = pSource->m[0][2];
			M.r[2].vector4_f32[1] = pSource->m[1][2];
			M.r[2].vector4_f32[2] = pSource->m[2][2];
			M.r[2].vector4_f32[3] = 0.0f;

			M.r[3].vector4_f32[0] = pSource->m[0][3];
			M.r[3].vector4_f32[1] = pSource->m[1][3];
			M.r[3].vector4_f32[2] = pSource->m[2][3];
			M.r[3].vector4_f32[3] = 1.0f;

			return M;
#elif defined(_SSE2_INTRINSICS_)
			MATRIX M;
			M.r[0] = _mm_load_ps(&pSource->_00);
			M.r[1] = _mm_load_ps(&pSource->_10);
			M.r[2] = _mm_load_ps(&pSource->_20);
			M.r[3] = g_IdentityR3;

			// x.x, x.y, y.x, y.y
			__m128 V1 = _mm_shuffle_ps(M.r[0], M.r[1], _MM_SHUFFLE(1, 0, 1, 0));

			// x.z, x.w, y.z, y.w
			__m128 V3 = _mm_shuffle_ps(M.r[0], M.r[1], _MM_SHUFFLE(3, 2, 3, 2));

			// z.x, z.y, w.x, w.y
			__m128 V2 = _mm_shuffle_ps(M.r[2], M.r[3], _MM_SHUFFLE(1, 0, 1, 0));

			// z.z, z.w, w.z, w.w
			__m128 V4 = _mm_shuffle_ps(M.r[2], M.r[3], _MM_SHUFFLE(3, 2, 3, 2));

			MATRIX mResult;
			// x.x, y.x, z.x, w.x
			mResult.r[0] = _mm_shuffle_ps(V1, V2, _MM_SHUFFLE(2, 0, 2, 0));

			// x.y, y.y, z.y, w.y
			mResult.r[1] = _mm_shuffle_ps(V1, V2, _MM_SHUFFLE(3, 1, 3, 1));

			// x.z, y.z, z.z, w.z
			mResult.r[2] = _mm_shuffle_ps(V3, V4, _MM_SHUFFLE(2, 0, 2, 0));

			// x.w, y.w, z.w, w.w
			mResult.r[3] = _mm_shuffle_ps(V3, V4, _MM_SHUFFLE(3, 1, 3, 1));

			return mResult;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreFloat3x4(Float3x4* pDestination, A_MATRIX m) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->m[0][0] = m.r[0].vector4_f32[0];
			pDestination->m[0][1] = m.r[1].vector4_f32[0];
			pDestination->m[0][2] = m.r[2].vector4_f32[0];
			pDestination->m[0][3] = m.r[3].vector4_f32[0];
									
			pDestination->m[1][0] = m.r[0].vector4_f32[1];
			pDestination->m[1][1] = m.r[1].vector4_f32[1];
			pDestination->m[1][2] = m.r[2].vector4_f32[1];
			pDestination->m[1][3] = m.r[3].vector4_f32[1];
									
			pDestination->m[2][0] = m.r[0].vector4_f32[2];
			pDestination->m[2][1] = m.r[1].vector4_f32[2];
			pDestination->m[2][2] = m.r[2].vector4_f32[2];
			pDestination->m[2][3] = m.r[3].vector4_f32[2];

#elif defined(_SSE2_INTRINSICS_)
			// x.x, x.y, y.x, y.y
			VECTOR V1 = _mm_shuffle_ps(m.r[0], m.r[1], _MM_SHUFFLE(1, 0, 1, 0));
			// x.z, x.w, y.z, y.w
			VECTOR V3 = _mm_shuffle_ps(m.r[0], m.r[1], _MM_SHUFFLE(3, 2, 3, 2));
			// z.x, z.y, w.x, w.y
			VECTOR V2 = _mm_shuffle_ps(m.r[2], m.r[3], _MM_SHUFFLE(1, 0, 1, 0));
			// z.z, z.w, w.z, w.w
			VECTOR V4 = _mm_shuffle_ps(m.r[2], m.r[3], _MM_SHUFFLE(3, 2, 3, 2));

			// x.x, y.x, z.x, w.x
			VECTOR RO = _mm_shuffle_ps(V1, V2, _MM_SHUFFLE(2, 0, 2, 0));
			// x.y, y.y, z.y, w.y
			VECTOR R1 = _mm_shuffle_ps(V1, V2, _MM_SHUFFLE(3, 1, 3, 1));
			// x.z, y.z, z.z, w.z
			VECTOR R2 = _mm_shuffle_ps(V3, V4, _MM_SHUFFLE(2, 0, 2, 0));

			_mm_storeu_ps(&pDestination->m[0][0], R0);
			_mm_storeu_ps(&pDestination->m[1][0], R1);
			_mm_storeu_ps(&pDestination->m[2][0], R2);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreAFloat3x4(AFloat3x4* pDestination, A_MATRIX m) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
			assert((reinterpret_cast<uintptr_t>(pDestination) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->m[0][0] = m.r[0].vector4_f32[0];
			pDestination->m[0][1] = m.r[1].vector4_f32[0];
			pDestination->m[0][2] = m.r[2].vector4_f32[0];
			pDestination->m[0][3] = m.r[3].vector4_f32[0];

			pDestination->m[1][0] = m.r[0].vector4_f32[1];
			pDestination->m[1][1] = m.r[1].vector4_f32[1];
			pDestination->m[1][2] = m.r[2].vector4_f32[1];
			pDestination->m[1][3] = m.r[3].vector4_f32[1];

			pDestination->m[2][0] = m.r[0].vector4_f32[2];
			pDestination->m[2][1] = m.r[1].vector4_f32[2];
			pDestination->m[2][2] = m.r[2].vector4_f32[2];
			pDestination->m[2][3] = m.r[3].vector4_f32[2];

#elif defined(_SSE2_INTRINSICS_)
			// x.x, x.y, y.x, y.y
			VECTOR V1 = _mm_shuffle_ps(m.r[0], m.r[1], _MM_SHUFFLE(1, 0, 1, 0));
			// x.z, x.w, y.z, y.w
			VECTOR V3 = _mm_shuffle_ps(m.r[0], m.r[1], _MM_SHUFFLE(3, 2, 3, 2));
			// z.x, z.y, w.x, w.y
			VECTOR V2 = _mm_shuffle_ps(m.r[2], m.r[3], _MM_SHUFFLE(1, 0, 1, 0));
			// z.z, z.w, w.z, w.w
			VECTOR V4 = _mm_shuffle_ps(m.r[2], m.r[3], _MM_SHUFFLE(3, 2, 3, 2));

			// x.x, y.x, z.x, w.x
			VECTOR RO = _mm_shuffle_ps(V1, V2, _MM_SHUFFLE(2, 0, 2, 0));
			// x.y, y.y, z.y, w.y
			VECTOR R1 = _mm_shuffle_ps(V1, V2, _MM_SHUFFLE(3, 1, 3, 1));
			// x.z, y.z, z.z, w.z
			VECTOR R2 = _mm_shuffle_ps(V3, V4, _MM_SHUFFLE(2, 0, 2, 0));

			_mm_store_ps(&pDestination->m[0][0], R0);
			_mm_store_ps(&pDestination->m[1][0], R1);
			_mm_store_ps(&pDestination->m[2][0], R2);
#endif
		}
	}
}

#endif // !ULTREALITY_MATH_FLOAT3X4_INL
