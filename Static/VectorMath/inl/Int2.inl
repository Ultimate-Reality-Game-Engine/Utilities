#ifndef ULTREALITY_MATH_INT2_INL
#define ULTREALITY_MATH_INT2_INL

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
	FORCE_INLINE Int2::Int2(A_VECTOR v) noexcept
	{
		Vector::StoreInt2(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE Int2& VEC_CALLCONV Int2::operator=(A_VECTOR v) noexcept
	{
		Vector::StoreInt2(this, v);

		return *this;
	}

	FORCE_INLINE VECTOR VEC_CALLCONV Int2::Load() noexcept
	{
		return Vector::LoadInt2(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV Int2::Store(A_VECTOR v) noexcept
	{
		Vector::StoreInt2(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE AInt2::AInt2(A_VECTOR v) noexcept
	{
		Vector::StoreAInt2(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE AInt2& VEC_CALLCONV AInt2::operator=(A_VECTOR v) noexcept
	{
		Vector::StoreAInt2(this, v);

		return *this;
	}

	FORCE_INLINE VECTOR VEC_CALLCONV AInt2::Load() noexcept
	{
		return Vector::LoadAInt2(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV AInt2::Store(A_VECTOR v) noexcept
	{
		Vector::StoreAInt2(this, v);
	}

	namespace Vector
	{
		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadInt2(const Int2* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR V;
			V.vector4_f32[0] = static_cast<float>(pSource->x);
			V.vector4_f32[1] = static_cast<float>(pSource->y);
			V.vector4_f32[2] = 0.0f;
			V.vector4_f32[3] = 0.0f;

			return V;
#elif defined(_SSE2_INTRINSICS_)
			__m128 V = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(pSource)));
			return _mm_cvtepi32_ps(_mm_castps_si128(V));
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadAInt2(const AInt2* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
			assert((reinterpret_cast<uintptr_t>(pSource) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR V;
			V.vector4_f32[0] = static_cast<float>(pSource->x);
			V.vector4_f32[1] = static_cast<float>(pSource->y);
			V.vector4_f32[2] = 0.0f;
			V.vector4_f32[3] = 0.0f;

			return V;
#elif defined(_SSE2_INTRINSICS_)
			// Use aligned load for 16-byte aligned AInt2
			__m128i V = _mm_load_si128(reinterpret_cast<const __m128i*>(pSource));
			// Convert 32-bit integers to 32-bit floats
			return _mm_cvtepi32_ps(V);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreInt2(Int2* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->x = static_cast<int32_t>(v.vector4_f32[0]);
			pDestination->y = static_cast<int32_t>(v.vector4_f32[1]);

#elif defined(_SSE2_INTRINSICS_)
			// Detect positive overflow
			VECTOR  vOverflow = _mm_cmpgt_ps(v, g_MaxInt);

			// Convert from float to int
			__m128i vResulti = _mm_cvttps_epi32(v);

			// Set positive overflow elements to 0x7FFFFFFF
			VECTOR vResult = _mm_and_ps(vOverflow, g_AbsMask);
			vOverflow = _mm_andnot_ps(vOverflow, _mm_castsi128_ps(vResulti));
			vOverflow = _mm_or_ps(vOverflow, vResult);

			// Write x, y to Int2
			_mm_store_sd(reinterpret_cast<double*>(pDestination), _mm_castps_pd(vOverflow));
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreAInt2(AInt2* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
			assert((reinterpret_cast<uintptr_t>(pDestination) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->x = static_cast<int32_t>(v.vector4_f32[0]);
			pDestination->y = static_cast<int32_t>(v.vector4_f32[1]);

#elif defined(_SSE2_INTRINSICS_)
			// Detect positive overflow
			VECTOR  vOverflow = _mm_cmpgt_ps(v, g_MaxInt);

			// Convert from float to int
			__m128i vResulti = _mm_cvttps_epi32(v);

			// Set positive overflow elements to 0x7FFFFFFF
			VECTOR vResult = _mm_and_ps(vOverflow, g_AbsMask);
			vOverflow = _mm_andnot_ps(vOverflow, _mm_castsi128_ps(vResulti));
			vOverflow = _mm_or_ps(vOverflow, vResult);

			// Write x, y to Int2
			_mm_store_sd(reinterpret_cast<double*>(pDestination), _mm_castps_pd(vOverflow));
#endif
		}
	}

	namespace Vector2
	{
		FORCE_INLINE bool VEC_CALLCONV EqualInt(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return (((V1.vector4_u32[0] == V2.vector4_u32[0]) && (V1.vector4_u32[1] == V2.vector4_u32[1])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			__m128i vTempi = _mm_cmpeq_epi32(_mm_castps_si128(V1), _mm_castps_si128(V2));

			return (((_mm_movemask_ps(_mm_castsi128_ps(vTempi)) & 3) == 3) != 0);
#endif
		}

		FORCE_INLINE uint32_t VEC_CALLCONV EqualIntR(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			uint32_t CR = 0;
			if ((V1.vector4_u32[0] == V2.vector4_u32[0]) &&
				(V1.vector4_u32[1] == V2.vector4_u32[1]))
			{
				CR = CRMASK_CR6TRUE;
			}
			else if ((V1.vector4_u32[0] != V2.vector4_u32[0]) &&
				(V1.vector4_u32[1] != V2.vector4_u32[1]))
			{
				CR = CRMASK_CR6FALSE;
			}

			return CR;

#elif defined(_SSE2_INTRINSICS_)
			__m128i vTempi = _mm_cmpeq_epi32(_mm_castps_si128(V1), _mm_castps_si128(V2));
			int iTest = _mm_movemask_ps(_mm_castsi128_ps(vTempi)) & 3;
			uint32_t CR = 0;
			if (iTest == 3)
			{
				CR = CRMASK_CR6TRUE;
			}
			else if (!iTest)
			{
				CR = CRMASK_CR6FALSE;
			}

			return CR;
#endif
		}

		FORCE_INLINE bool VEC_CALLCONV NotEqualInt(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return (((V1.vector4_u32[0] != V2.vector4_u32[0]) || (V1.vector4_u32[1] != V2.vector4_u32[1])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			__m128i vTempi = _mm_cmpeq_epi32(_mm_castps_si128(V1), _mm_castps_si128(V2));

			return (((_mm_movemask_ps(_mm_castsi128_ps(vTempi)) & 3) != 3) != 0);
#endif
		}
	}
}

#endif // !ULTREALITY_MATH_INT2_INL
