#ifndef ULTREALITY_MATH_UINT4_INL
#define ULTREALITY_MATH_UINT4_INL

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
	FORCE_INLINE UInt4(A_VECTOR v) noexcept
	{
		Vector::StoreUInt4(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE UInt4& VEC_CALLCONV UInt4::operator=(A_VECTOR v) noexcept
	{
		Vector::StoreUInt4(this, v);

		return *this;
	}

	FORCE_INLINE VECTOR VEC_CALLCONV UInt4::Load() noexcept
	{
		return Vector::LoadUInt4(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV UInt4::Store(A_VECTOR v) noexcept
	{
		Vector::StoreUInt4(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE AUInt4::AUInt4(A_VECTOR v) noexcept
	{
		Vector::StoreAUInt4(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE AUInt4& VEC_CALLCONV AUInt4::operator=(A_VECTOR v) noexcept
	{
		Vector::StoreAUInt4(this, v);

		return *this;
	}

	FORCE_INLINE VECTOR VEC_CALLCONV AUInt4::Load() noexcept
	{
		return Vector::LoadAUInt4(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV AUInt4::Store(A_VECTOR v) noexcept
	{
		Vector::StoreAUInt4(this, v);

		return *this;
	}

	namespace Vector
	{
		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadUInt4(const UInt4* pSource) noexcept
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

			// For values higher than 0x7FFFFFFF a fix is needed
			// Determine which elements need the fix
			VECTOR vMask = _mm_and_ps(_mm_castsi128_ps(V), g_NegativeZero);

			// Force all elements positive
			VECTOR vResult = _mm_xor_ps(_mm_castsi128_ps(V), vMask);

			// Convert to floats
			vResult = _mm_cvtepi32_ps(_mm_castps_si128(vResult));

			// Convert 0x80000000 -> 0xFFFFFFFF
			__m128i iMask = _mm_srai_epi32(_mm_castps_si128(vMask), 31);

			// For the elements that are too big apply the fix
			vMask = _mm_and_ps(_mm_castsi128_ps(iMask), g_FixUnsigned);
			vResult = _mm_and_ps(vResult, vMask);
			return vResult;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadAUInt4(const AUInt4* pSource) noexcept
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

			// For values higher than 0x7FFFFFFF a fix is needed
			// Determine which elements need the fix
			VECTOR vMask = _mm_and_ps(_mm_castsi128_ps(V), g_NegativeZero);

			// Force all elements positive
			VECTOR vResult = _mm_xor_ps(_mm_castsi128_ps(V), vMask);

			// Convert to floats
			vResult = _mm_cvtepi32_ps(_mm_castps_si128(vResult));

			// Convert 0x80000000 -> 0xFFFFFFFF
			__m128i iMask = _mm_srai_epi32(_mm_castps_si128(vMask), 31);

			// For the elements that are too big apply the fix
			vMask = _mm_and_ps(_mm_castsi128_ps(iMask), g_FixUnsigned);
			vResult = _mm_and_ps(vResult, vMask);
			return vResult;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreUInt4(UInt4* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->x = v.vector4_u32[0];
			pDestination->y = v.vector4_u32[1];
			pDestination->z = v.vector4_u32[2];
			pDestination->w = v.vector4_u32[3];

#elif defined(_SSE2_INTRINSICS_)
			// Clamp to >= 0
			VECTOR vResult = _mm_max_ps(v, g_Zero);

			// Any element that is too big, set to 0xFFFFFFFFU
			VECTOR vOverflow = _mm_cmpgt_ps(vResult, g_MaxUInt);
			
			// Find elements too large for signed integer
			VECTOR vMask = _mm_cmpge_ps(vResult, g_UnsignedFix);

			// Zero for nubbers lower than 0x80000000, 32768.0f*65536.0f otherwise
			vValue = _mm_and_ps(vValue, vMask);

			// Perform fix only on elements that are too large
			vResult = _mm_sub_ps(vResult, vValue);
			__m128i vResulti = _mm_cvttps_epi32(vResult);

			// Convert from singed to unsigned for elements greater than 0x80000000
			vMask = _mm_and_ps(vMask, g_NegativeZero);
			vResult = _mm_xor_ps(_mm_castsi128_ps(vResulti), vMask);

			// On elements that are too large, set to 0xFFFFFFFF
			vResult = _mm_or_ps(vResult, vOverflow);
			_mm_storeu_si128(reinterpret_cast<__m128i*>(pDestination), _mm_castps_si128(vResult));
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreAUInt4(AUInt4* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
			assert((reinterpret_cast<uintptr_t>(pDestination) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->x = v.vector4_u32[0];
			pDestination->y = v.vector4_u32[1];
			pDestination->z = v.vector4_u32[2];
			pDestination->w = v.vector4_u32[3];

#elif defined(_SSE2_INTRINSICS_)
			// Clamp to >= 0
			VECTOR vResult = _mm_max_ps(v, g_Zero);

			// Any element that is too big, set to 0xFFFFFFFFU
			VECTOR vOverflow = _mm_cmpgt_ps(vResult, g_MaxUInt);

			// Find elements too large for signed integer
			VECTOR vMask = _mm_cmpge_ps(vResult, g_UnsignedFix);

			// Zero for nubbers lower than 0x80000000, 32768.0f*65536.0f otherwise
			vValue = _mm_and_ps(vValue, vMask);

			// Perform fix only on elements that are too large
			vResult = _mm_sub_ps(vResult, vValue);
			__m128i vResulti = _mm_cvttps_epi32(vResult);

			// Convert from singed to unsigned for elements greater than 0x80000000
			vMask = _mm_and_ps(vMask, g_NegativeZero);
			vResult = _mm_xor_ps(_mm_castsi128_ps(vResulti), vMask);

			// On elements that are too large, set to 0xFFFFFFFF
			vResult = _mm_or_ps(vResult, vOverflow);
			_mm_store_si128(reinterpret_cast<__m128i*>(pDestination), _mm_castps_si128(vResult));
#endif
		}
	}
}

#endif // !ULTREALITY_MATH_UINT4_INL
