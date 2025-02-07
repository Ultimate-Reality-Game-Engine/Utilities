#ifndef ULTREALITY_MATH_FLOAT4X4_INL
#define ULTREALITY_MATH_FLOAT4X4_INL

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
	FORCE_INLINE Float4x4::Float4x4(const float* pArray) noexcept
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

		m[3][0] = pArray[12];
		m[3][1] = pArray[13];
		m[3][2] = pArray[14];
		m[3][3] = pArray[15];
	}

	_Use_decl_annotations_
	FORCE_INLINE Float4x4::Float4x4(A_MATRIX m) noexcept
	{
		MAT::StoreFloat4x4(this, m);
	}

	_Use_decl_annotations_
	FORCE_INLINE Float4x4& VEC_CALLCONV Float4x4::operator=(A_MATRIX m) noexcept
	{
		MAT::StoreFloat4x4(this, m);

		return *this;
	}

	FORCE_INLINE MATRIX VEC_CALLCONV Float4x4::Load() noexcept
	{
		return MAT::LoadFloat4x4(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV Float4x4::Store(A_MATRIX m) noexcept
	{
		MAT::StoreFloat4x4(this, m);
	}

	_Use_decl_annotations_
	FORCE_INLINE AFloat4x4::AFloat4x4(A_MATRIX m) noexcept
	{
		MAT::StoreAFloat4x4(this, m);
	}

	_Use_decl_annotations_
	FORCE_INLINE AFloat4x4& VEC_CALLCONV AFloat4x4::operator=(A_MATRIX m) noexcept
	{
		MAT::StoreAFloat4x4(this, m);

		return *this;
	}

	FORCE_INLINE MATRIX VEC_CALLCONV AFloat4x4::Load() noexcept
	{
		return MAT::LoadAFloat4x4(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV AFloat4x4::Store(A_MATRIX m) noexcept
	{
		MAT::StoreAFloat4x4(this, m);
	}

	namespace MAT
	{
		_Use_decl_annotations_
			FORCE_INLINE MATRIX VEC_CALLCONV LoadFloat4x4(const Float4x4* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			MATRIX M;

			M.r[0].vector4_f32[0] = pSource->m[0][0];
			M.r[0].vector4_f32[1] = pSource->m[0][1];
			M.r[0].vector4_f32[2] = pSource->m[0][2];
			M.r[0].vector4_f32[3] = pSource->m[0][3];

			M.r[1].vector4_f32[0] = pSource->m[1][0];
			M.r[1].vector4_f32[1] = pSource->m[1][1];
			M.r[1].vector4_f32[2] = pSource->m[1][2];
			M.r[1].vector4_f32[3] = pSource->m[1][3];

			M.r[2].vector4_f32[0] = pSource->m[2][0];
			M.r[2].vector4_f32[1] = pSource->m[2][1];
			M.r[2].vector4_f32[2] = pSource->m[2][2];
			M.r[2].vector4_f32[3] = pSource->m[2][3];

			M.r[3].vector4_f32[0] = pSource->m[3][0];
			M.r[3].vector4_f32[1] = pSource->m[3][1];
			M.r[3].vector4_f32[2] = pSource->m[3][2];
			M.r[3].vector4_f32[3] = pSource->m[3][3];

			return M;
#elif defined(_SSE2_INTRINSICS_)
			MATRIX M;

			M.r[0] = _mm_loadu_ps(&pSource->_00);
			M.r[1] = _mm_loadu_ps(&pSource->_10);
			M.r[2] = _mm_loadu_ps(&pSource->_20);
			M.r[3] = _mm_loadu_ps(&pSource->_30);

			return M;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE MATRIX VEC_CALLCONV LoadAFloat4x4(const AFloat4x4* pSource) noexcept
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
			M.r[0].vector4_f32[3] = pSource->m[0][3];

			M.r[1].vector4_f32[0] = pSource->m[1][0];
			M.r[1].vector4_f32[1] = pSource->m[1][1];
			M.r[1].vector4_f32[2] = pSource->m[1][2];
			M.r[1].vector4_f32[3] = pSource->m[1][3];

			M.r[2].vector4_f32[0] = pSource->m[2][0];
			M.r[2].vector4_f32[1] = pSource->m[2][1];
			M.r[2].vector4_f32[2] = pSource->m[2][2];
			M.r[2].vector4_f32[3] = pSource->m[2][3];

			M.r[3].vector4_f32[0] = pSource->m[3][0];
			M.r[3].vector4_f32[1] = pSource->m[3][1];
			M.r[3].vector4_f32[2] = pSource->m[3][2];
			M.r[3].vector4_f32[3] = pSource->m[3][3];

			return M;
#elif defined(_SSE2_INTRINSICS_)
			MATRIX M;

			M.r[0] = _mm_load_ps(&pSource->_00);
			M.r[1] = _mm_load_ps(&pSource->_10);
			M.r[2] = _mm_load_ps(&pSource->_20);
			M.r[3] = _mm_load_ps(&pSource->_30);

			return M;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreFloat4x4(Float4x4* pDestination, A_MATRIX m) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->m[0][0] = m.r[0].vector4_f32[0];
			pDestination->m[0][1] = m.r[0].vector4_f32[1];
			pDestination->m[0][2] = m.r[0].vector4_f32[2];
			pDestination->m[0][3] = m.r[0].vector4_f32[3];
									
			pDestination->m[1][0] = m.r[1].vector4_f32[0];
			pDestination->m[1][1] = m.r[1].vector4_f32[1];
			pDestination->m[1][2] = m.r[1].vector4_f32[2];
			pDestination->m[1][3] = m.r[1].vector4_f32[3];
									
			pDestination->m[2][0] = m.r[2].vector4_f32[0];
			pDestination->m[2][1] = m.r[2].vector4_f32[1];
			pDestination->m[2][2] = m.r[2].vector4_f32[2];
			pDestination->m[2][3] = m.r[2].vector4_f32[3];
									
			pDestination->m[3][0] = m.r[3].vector4_f32[0];
			pDestination->m[3][1] = m.r[3].vector4_f32[1];
			pDestination->m[3][2] = m.r[3].vector4_f32[2];
			pDestination->m[3][3] = m.r[3].vector4_f32[3];

#elif defined(_SSE2_INTRINSICS_)
			_mm_storeu_ps(&pDestination->_00, m.r[0]);
			_mm_storeu_ps(&pDestination->_10, m.r[1]);
			_mm_storeu_ps(&pDestination->_20, m.r[2]);
			_mm_storeu_ps(&pDestination->_30, m.r[3]);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreAFloat4x4(AFloat4x4* pDestination, A_MATRIX m) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
			assert((reinterpret_cast<uintptr_t>(pDestination) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->m[0][0] = m.r[0].vector4_f32[0];
			pDestination->m[0][1] = m.r[0].vector4_f32[1];
			pDestination->m[0][2] = m.r[0].vector4_f32[2];
			pDestination->m[0][3] = m.r[0].vector4_f32[3];

			pDestination->m[1][0] = m.r[1].vector4_f32[0];
			pDestination->m[1][1] = m.r[1].vector4_f32[1];
			pDestination->m[1][2] = m.r[1].vector4_f32[2];
			pDestination->m[1][3] = m.r[1].vector4_f32[3];

			pDestination->m[2][0] = m.r[2].vector4_f32[0];
			pDestination->m[2][1] = m.r[2].vector4_f32[1];
			pDestination->m[2][2] = m.r[2].vector4_f32[2];
			pDestination->m[2][3] = m.r[2].vector4_f32[3];

			pDestination->m[3][0] = m.r[3].vector4_f32[0];
			pDestination->m[3][1] = m.r[3].vector4_f32[1];
			pDestination->m[3][2] = m.r[3].vector4_f32[2];
			pDestination->m[3][3] = m.r[3].vector4_f32[3];

#elif defined(_SSE2_INTRINSICS_)
			_mm_store_ps(&pDestination->_00, m.r[0]);
			_mm_store_ps(&pDestination->_10, m.r[1]);
			_mm_store_ps(&pDestination->_20, m.r[2]);
			_mm_store_ps(&pDestination->_30, m.r[3]);
#endif
		}
	}
}

#endif // !ULTREALITY_MATH_FLOAT4X4_INL
