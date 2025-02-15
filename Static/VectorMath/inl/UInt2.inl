#ifndef ULTREALITY_MATH_UINT2_INL
#define ULTREALITY_MATH_UINT2_INL

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
	FORCE_INLINE UInt2::UInt2(A_VECTOR v) noexcept
	{
		Vector::StoreUInt2(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE UInt2& VEC_CALLCONV UInt2::operator=(A_VECTOR v) noexcept
	{
		Vector::StoreUInt2(this, v);

		return *this;
	}

	FORCE_INLINE VECTOR VEC_CALLCONV UInt2::Load() noexcept
	{
		return Vector::LoadUInt2(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV UInt2::Store(A_VECTOR v) noexcept
	{
		Vector::StoreUInt2(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE AUInt2::AUInt2(A_VECTOR v) noexcept
	{
		Vector::StoreAUInt2(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE AUInt2& VEC_CALLCONV AUInt2::operator=(A_VECTOR v) noexcept
	{
		Vector::StoreAUInt2(this, v);

		return *this;
	}

	FORCE_INLINE VECTOR VEC_CALLCONV AUInt2::Load() noexcept
	{
		return Vector::LoadAUInt2(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV AUInt2::Store(A_VECTOR v) noexcept
	{
		Vector::StoreAUInt2(this, v);
	}

	namespace Vector
	{
		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadUInt2(const UInt2* pSource) noexcept
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
			// For values higher than 0x7FFFFFFF, a fix is needed
			VECTOR vMask = _mm_and_ps(V, g_NegativeZero);
			// Force all values positive
			VECTOR vResult = _mm_xor_ps(V, vMask);
			// Convert to floats
			vResult = _mm_cvtepi32_ps(_mm_castps_si128(vResult));
			// Convert 0x80000000 -> 0xFFFFFFFF
			__m128i iMask = _mm_srai_epi32(_mm_castps_si128(vMask), 31);
			// For only the ones that are too big, add the fix
			vMask = _mm_and_ps(_mm_castsi128_ps(iMask), g_FixUnsigned);
			vResult = _mm_and_ps(vResult, vMask);
			return vResult;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadAUInt2(const AUInt2* pSource) noexcept
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
			// Load aligned 16-byte data
			__m128 V = _mm_castsi128_ps(_mm_load_si128(reinterpret_cast<const __m128i*>(pSource)));

			// For values higher than 0x7FFFFFFF, a fix is needed
			VECTOR vMask = _mm_and_ps(V, g_NegativeZero);
			// Force all values positive
			VECTOR vResult = _mm_xor_ps(V, vMask);
			// Convert to floats
			vResult = _mm_cvtepi32_ps(_mm_castps_si128(vResult));
			// Convert 0x80000000 -> 0xFFFFFFFF
			__m128i iMask = _mm_srai_epi32(_mm_castps_si128(vMask), 31);
			// For only the ones that are too big, add the fix
			vMask = _mm_and_ps(_mm_castsi128_ps(iMask), g_FixUnsigned);
			vResult = _mm_and_ps(vResult, vMask);
			return vResult;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreUInt2(UInt2* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->x = v.vector4_u32[0];
			pDestination->y = v.vector4_u32[1];

#elif defined(_SSE2_INTRINSICS_)
			// Clamp elements to >= 0
			VECTOR vResult = _mm_max_ps(v, g_Zero);

			// Any numbers that are too large, set to 0x7FFFFFFFU
			VECTOR vOverflow = _mm_cmpgt_ps(vResult, g_MaxUInt);
			VECTOR vValue = g_UnsignedFix;

			// Detect elements too large for signed int
			VECTOR vMask = _mm_cmpge_ps(vResult, vValue);

			// Zero numbers lower than 0x80000000, 32768.0f*65536.0f otherwise
			vValue = _mm_and_ps(vValue, vMask);

			// Perform fix for the numbers that are too large
			vResult = _mm_sub_ps(vValue, vMask);
			__m128i vResulti = _mm_cvttps_epi32(vResult);

			// Convert from signed to unsigned if element greater than 0x80000000
			vMask = _mm_and_ps(vMask, g_NegativeZero);
			vResult = _mm_xor_ps(_mm_castsi128_ps(vResulti), vMask);

			// On elements that are too large, set to 0xFFFFFFFF
			vResult = _mm_or_ps(vResult, vOverflow);

			// Write x and y to UInt2
			_mm_store_sd(reinterpret_cast<double*>(pDestination), _mm_castps_pd(vResult));
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreAUInt2(AUInt2* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
			assert((reinterpret_cast<uintptr_t>(pDestination) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->x = v.vector4_u32[0];
			pDestination->y = v.vector4_u32[1];

#elif defined(_SSE2_INTRINSICS_)
			// Clamp elements to >= 0
			VECTOR vResult = _mm_max_ps(v, g_Zero);

			// Any numbers that are too large, set to 0x7FFFFFFFU
			VECTOR vOverflow = _mm_cmpgt_ps(vResult, g_MaxUInt);
			VECTOR vValue = g_UnsignedFix;

			// Detect elements too large for signed int
			VECTOR vMask = _mm_cmpge_ps(vResult, vValue);

			// Zero numbers lower than 0x80000000, 32768.0f*65536.0f otherwise
			vValue = _mm_and_ps(vValue, vMask);

			// Perform fix for the numbers that are too large
			vResult = _mm_sub_ps(vValue, vMask);
			__m128i vResulti = _mm_cvttps_epi32(vResult);

			// Convert from signed to unsigned if element greater than 0x80000000
			vMask = _mm_and_ps(vMask, g_NegativeZero);
			vResult = _mm_xor_ps(_mm_castsi128_ps(vResulti), vMask);

			// On elements that are too large, set to 0xFFFFFFFF
			vResult = _mm_or_ps(vResult, vOverflow);

			// Write x and y to UInt2
			_mm_store_sd(reinterpret_cast<double*>(pDestination), _mm_castps_pd(vResult));
#endif
		}
	}
}

#endif // !ULTREALITY_MATH_UINT2_INL
