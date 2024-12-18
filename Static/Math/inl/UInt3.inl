#ifndef ULTREALITY_MATH_UINT3_INL
#define ULTREALITY_MATH UINT3_INL

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
	FORCE_INLINE UInt3::UInt3(A_VECTOR v) noexcept
	{
		VEC::StoreUInt3(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE UInt3& VEC_CALLCONV UInt3::operator=(A_VECTOR v) noexcept
	{
		VEC::StoreUInt3(this, v);

		return *this;
	}

	FORCE_INLINE VECTOR VEC_CALLCONV UInt3::Load() noexcept
	{
		return VEC::LoadUInt3(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV UInt3::Store(A_VECTOR v) noexcept
	{
		VEC::StoreUInt3(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE AUInt3::AUInt3(A_VECTOR v) noexcept
	{
		VEC::StoreAUInt3(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE AUInt3& VEC_CALLCONV AUInt3::operator=(A_VECTOR v) noexcept
	{
		VEC::StoreAUInt3(this, v);
	}

	FORCE_INLINE VECTOR VEC_CALLCONV AUInt3::Load() noexcept
	{
		return VEC::LoadAUInt3(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV AUInt3::Store(A_VECTOR v) noexcept
	{
		VEC::StoreAUInt3(this, v);
	}

	namespace VEC
	{
		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadUInt3(const UInt3* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR V;
			V.vector4_f32[0] = static_cast<float>(pSource->x);
			V.vector4_f32[1] = static_cast<float>(pSource->y);
			V.vector4_f32[2] = static_cast<float>(pSource->z);
			V.vector4_f32[3] = 0.0f;

			return V;
#elif defined(_SSE2_INTRINSICS_)
			__m128 xy = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(pSource)));
			__m128 z = _mm_load_ss(reinterpret_cast<const float*>(&pSource->z));
			__m128 V = _mm_movelh_ps(xy, z);

			// For values greater tha 0x7FFFFFFF, a fix is needed
			// Determine which elements need the fix
			VECTOR vMask = _mm_and_ps(V, g_NegativeZero);

			// Force all values positive
			VECTOR vResult = _mm_xor_ps(V, vMask);

			// Convert to floats
			vResult = _mm_cvtepi32_ps(_mm_castps_si128(vResult));

			// Convert 0x80000000 -> 0xFFFFFFFF
			__m128i iMask = _mm_srai_epi32(_mm_castps_si128(vMask), 31);

			// For the elements that are too big, apply the fix
			vMask = _mm_and_ps(_mm_castsi128_ps(iMask), g_FixUnsigned);

			vResult = _mm_add_ps(vResult, vMask);
			
			return vResult;
#endif
		}

		_Use_dec_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadAUInt3(const AUInt3* pSource) noexcept
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
			V.vector4_f32[3] = 0.0f;

			return V;
#elif defined(_SSE2_INTRINSICS_)
			// Load all 16 bytes (4 integers) into a __m128i
			__m128i alignedInts = _mm_load_si128(reinterpret_cast<const __m128i*>(pSource));

			// Convert to floats
			__m128 V = _mm_cvtepi32_ps(alignedInts);

			// Mask off the w-component (keep x, y, z; zero w)
			V = _mm_and_ps(V, g_Mask3);

			// Handle unsigned integer conversion for large values (similar fix to the unaligned version)
			__m128 vMask = _mm_and_ps(V, g_NegativeZero); // Identify large unsigned values
			__m128 vResult = _mm_xor_ps(V, vMask);        // Force values positive

			// Convert large values (0x80000000 -> 0xFFFFFFFF)
			__m128i iMask = _mm_srai_epi32(_mm_castps_si128(vMask), 31); // Shift to create signed masks
			__m128 vFix = _mm_and_ps(_mm_castsi128_ps(iMask), g_FixUnsigned);

			// Apply fix for large values
			vResult = _mm_add_ps(vResult, vFix);

			return vResult;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreUInt3(UInt3* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->x = v.vector4_u32[0];
			pDestination->y = v.vector4_u32[1];
			pDestination->z = v.vector4_u32[2];

#elif defined(_SSE2_INTRINSICS_)
			// Clamp elements to >= 0
			VECTOR vResult = _mm_max_ps(v, g_Zero);

			// Elements that are too but, set to 0xFFFFFFFFU
			VECTOR vOverflow = _mm_cmpgt_ps(vResult, g_MaxUInt);
			
			// Check for elements too large for signed integer
			VECTOR vMask = _mm_cmpge_ps(vResult, g_UnsignedFix);

			// Zero for numbers lower than 0x80000000, 32768.0f*65536.0f otherwise
			vValue = _mm_and_ps(vValue, vMask);

			// Perform fix on numbers too large
			vResult = _mm_sub_ps(vResult, vValue);
			__m128i vResulti = _mm_cvttps_epi32(vResult);

			// Convert from signed to unsigned only of element is greater than 0x80000000
			vMask = _mm_and_ps(vMask, g_NegativeZero);
			vResult = _mm_xor_ps(_mm_castsi128_ps(vResulti), vMask);

			// On elements too large, set to 0xFFFFFFFF
			vResult = _mm_or_ps(vResult, vOverflow);

			// Write x, y, z
			_mm_store_sd(reinterpret_cast<double*>(pDestination), _mm_castps_pd(vResult));
			__m128 z = PERMUTE_PS(vResult, _MM_SHUFFLE(2, 2, 2, 2));
			_mm_store_ss(reinterpret_cast<float*>(&pDestination->z), z);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreAUInt3(AUInt3* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
			assert((reinterpret_cast<uintptr_t>(pDestination) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->x = v.vector4_u32[0];
			pDestination->y = v.vector4_u32[1];
			pDestination->z = v.vector4_u32[2];

#elif defined(_SSE2_INTRINSICS_)
			// Clamp elements to >= 0
			VECTOR vResult = _mm_max_ps(v, g_Zero);

			// Elements that are too but, set to 0xFFFFFFFFU
			VECTOR vOverflow = _mm_cmpgt_ps(vResult, g_MaxUInt);

			// Check for elements too large for signed integer
			VECTOR vMask = _mm_cmpge_ps(vResult, g_UnsignedFix);

			// Zero for numbers lower than 0x80000000, 32768.0f*65536.0f otherwise
			vValue = _mm_and_ps(vValue, vMask);

			// Perform fix on numbers too large
			vResult = _mm_sub_ps(vResult, vValue);
			__m128i vResulti = _mm_cvttps_epi32(vResult);

			// Convert from signed to unsigned only of element is greater than 0x80000000
			vMask = _mm_and_ps(vMask, g_NegativeZero);
			vResult = _mm_xor_ps(_mm_castsi128_ps(vResulti), vMask);

			// On elements too large, set to 0xFFFFFFFF
			vResult = _mm_or_ps(vResult, vOverflow);

			// Write x, y, z
			_mm_store_sd(reinterpret_cast<double*>(pDestination), _mm_castps_pd(vResult));
			__m128 z = _mm_movehl_ps(v, v);
			_mm_store_ss(reinterpret_cast<float*>(&pDestination->z), z);
#endif
		}
	}
}

#endif // !ULTREALITY_MATH_UINT3_INL
