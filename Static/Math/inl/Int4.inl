#ifndef ULTREALITY_MATH_INT4_INL
#define ULTREALITY_MATH_INT4_INL

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
	FORCE_INLINE Int4::Int4(A_VECTOR v) noexcept
	{
		VEC::StoreInt4(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE Int4& VEC_CALLCONV Int4::operator=(A_VECTOR v) noexcept
	{
		VEC::StoreInt4(this, v);

		return *this;
	}

	FORCE_INLINE VECTOR VEC_CALLCONV Int4::Load() noexcept
	{
		return VEC::LoadInt4(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV Int4::Store(A_VECTOR v) noexcept
	{
		VEC::StoreInt4(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE AInt4::AInt4(A_VECTOR v) noexcept
	{
		VEC::StoreAInt4(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE AInt4& VEC_CALLCONV AInt4::operator=(A_VECTOR v) noexcept
	{
		VEC::StoreAInt4(this, v);

		return *this;
	}

	FORCE_INLINE VECTOR VEC_CALLCONV AInt4::Load() noexcept
	{
		return VEC::LoadAInt4(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV AInt4::Store(A_VECTOR v) noexcept
	{
		VEC::StoreAInt4(this, v);
	}

	namespace VEC
	{
		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadInt4(const Int4* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR V;
			V.vector4_f32[0] = static_cast<float>(pSource->x);
			V.vector4_f32[1] = static_cast<float>(pSource->y);
			V.vector4_f32[2] = static_cast<float>(pSource->z);
			V.vector4_f32[3] = static_cast<float>(pSource->w);

			return V;
#elif defined(_SSE2_INTRINSICS_)
			__m128i V = _mm_loadu_si128(reinterpret_cast<const __m128i*>(pSource));
			return _mm_cvtepi32_ps(V);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadAInt4(const AInt4* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
			assert((reinterpret_cast<uintptr_t>(pSource) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR V;
			V.vector4_f32[0] = static_cast<float>(pSource->x);
			V.vector4_f32[1] = static_cast<float>(pSource->y);
			V.vector4_f32[2] = static_cast<float>(pSource->z);
			V.vector4_f32[3] = static_cast<float>(pSource->w);

			return V;
#elif defined(_SSE2_INTRINSICS_)
			__m128i V = _mm_load_si128(reinterpret_cast<const __m128i*>(pSource));
			return _mm_cvtepi32_ps(V);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreInt4(Int4* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->x = static_cast<int32_t>(v.vector4_f32[0]);
			pDestination->y = static_cast<int32_t>(v.vector4_f32[1]);
			pDestination->z = static_cast<int32_t>(v.vector4_f32[2]);
			pDestination->w = static_cast<int32_t>(v.vector4_f32[3]);

#elif defined(_SSE2_INTRINSICS_)
			// Detect positive overflow in elements
			VECTOR vOverflow = _mm_cmpgt_ps(v, g_MaxInt);

			// Convert float to int
			__m128i vResulti = _mm_cvttps_epi32(v);

			// Any elements with overflow, set to 0x7FFFFFFF
			VECTOR vResult = _mm_and_ps(vOverflow, g_AbsMask);
			vOverflow = _mm_andnot_ps(vOverflow, _mm_castsi128_ps(vResulti));
			vOverflow = _mm_or_ps(vOverflow, vResult);

			// Write x, y, z, w
			_mm_storeu_si128(reinterpret_cast<__m128i*>(pDestination), _mm_castps_si128(vOverflow));
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreAInt4(AInt4* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
			assert((reinterpret_cast<uintptr_t>(pDestination) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->x = static_cast<int32_t>(v.vector4_f32[0]);
			pDestination->y = static_cast<int32_t>(v.vector4_f32[1]);
			pDestination->z = static_cast<int32_t>(v.vector4_f32[2]);
			pDestination->w = static_cast<int32_t>(v.vector4_f32[3]);

#elif defined(_SSE2_INTRINSICS_)
			// Detect positive overflow in elements
			VECTOR vOverflow = _mm_cmpgt_ps(v, g_MaxInt);

			// Convert float to int
			__m128i vResulti = _mm_cvttps_epi32(v);

			// Any elements with overflow, set to 0x7FFFFFFF
			VECTOR vResult = _mm_and_ps(vOverflow, g_AbsMask);
			vOverflow = _mm_andnot_ps(vOverflow, _mm_castsi128_ps(vResulti));
			vOverflow = _mm_or_ps(vOverflow, vResult);

			// Write x, y, z, w
			_mm_store_si128(reinterpret_cast<__m128i*>(pDestination), _mm_castps_si128(vOverflow));
#endif
		}
	}

	namespace VEC4
	{
		FORCE_INLINE bool VEC_CALLCONV EqualInt(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
    		return (((V1.vector4_u32[0] == V2.vector4_u32[0]) && (V1.vector4_u32[1] == V2.vector4_u32[1]) && (V1.vector4_u32[2] == V2.vector4_u32[2]) && (V1.vector4_u32[3] == V2.vector4_u32[3])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			__m128i vTemp = _mm_cmpeq_epi32(_mm_castps_si128(V1), _mm_castps_si128(V2));

			return ((_mm_movemask_ps(_mm_castsi128_ps(vTemp)) == 0xf) != 0);
#endif
		}

		FORCE_INLINE uint32_t VEC_CALLCONV EqualIntR(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			uint32_t CR = 0;
			if (V1.vector4_u32[0] == V2.vector4_u32[0] &&
				V1.vector4_u32[1] == V2.vector4_u32[1] &&
				V1.vector4_u32[2] == V2.vector4_u32[2] &&
				V1.vector4_u32[3] == V2.vector4_u32[3])
			{
				CR = CRMASK_CR6TRUE;
			}
			else if (V1.vector4_u32[0] != V2.vector4_u32[0] &&
				V1.vector4_u32[1] != V2.vector4_u32[1] &&
				V1.vector4_u32[2] != V2.vector4_u32[2] &&
				V1.vector4_u32[3] != V2.vector4_u32[3])
			{
				CR = CRMASK_CR6FALSE;
			}
			
			return CR;

#elif defined(_SSE2_INTRINSICS_)
			__m128i vTemp = _mm_cmpeq_epi32(_mm_castps_si128(V1), _mm_castps_si128(V2));
			int iTest = _mm_movemask_ps(_mm_castsi128_ps(vTemp));
			uint32_t CR = 0;
			if (iTest == 0xf)
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
    		return (((V1.vector4_u32[0] != V2.vector4_u32[0]) || (V1.vector4_u32[1] != V2.vector4_u32[1]) || (V1.vector4_u32[2] != V2.vector4_u32[2]) || (V1.vector4_u32[3] != V2.vector4_u32[3])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			__m128i vTemp = _mm_cmpeq_epi32(_mm_castps_si128(V1), _mm_castps_si128(V2));
    		
			return ((_mm_movemask_ps(_mm_castsi128_ps(vTemp)) != 0xF) != 0);
#endif
		}
	}
}

#endif // !ULTREALITY_MATH_INT4_INL
