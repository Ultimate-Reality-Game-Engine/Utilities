#ifndef ULTREALITY_MATH_INT3_INL
#define ULTREALITY_MATH_INT3_INL

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
	FORCE_INLINE Int3::Int3(A_VECTOR v) noexcept
	{
		VEC::StoreInt3(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE Int3& VEC_CALLCONV Int3::operator=(A_VECTOR v) noexcept
	{
		VEC::StoreInt3(this, v);

		return *this;
	}

	FORCE_INLINE VECTOR VEC_CALLCONV Int3::Load() noexcept
	{
		return VEC::LoadInt3(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV Int3::Store(A_VECTOR v) noexcept
	{
		VEC::StoreInt3(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE AInt3::AInt3(A_VECTOR v) noexcept
	{
		VEC::StoreAInt3(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE AInt3& VEC_CALLCONV AInt3::operator=(A_VECTOR v) noexcept
	{
		VEC::StoreAInt3(this, v);

		return *this;
	}

	FORCE_INLINE VECTOR VEC_CALLCONV AInt3::Load() noexcept
	{
		return VEC::LoadAInt3(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV AInt3::Store(A_VECTOR v) noexcept
	{
		VEC::StoreAInt3(this, v);
	}

	namespace VEC
	{
		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadInt3(const Int3* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR V;
			V.vector4_f32[0] = static_cast<float>(pSource->x);
			V.vector4_f32[1] = static_cast<float>(pSource->y);
			V.vector4_f32[2] = static_cast<float>(pSource->z);
			V.vector4_f32[0] = 0.0f;

			return V;
#elif defined(_SSE2_INTRINSICS_)
			__m128 xy = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(pSource)));
			__m128 x = _mm_load_ss(reinterpret_cast<const float*>(&pSource->z));
			__m128 V = _mm_movelh_ps(xy, z);
			
			return _mm_cvtepi32_ps(_mm_castps_si128(V));
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadAInt3(const AInt3* pSource) noexcept
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
			V.vector4_f32[0] = 0.0f;

			return V;
#elif defined(_SSE2_INTRINSICS_)
			// Load all 16 bytes directly as aligned memory
			__m128i alignedInts = _mm_load_si128(reinterpret_cast<const __m128i*>(pSource));

			// Convert to floats
			__m128 V = _mm_cvtepi32_ps(alignedInts);

			// Mask off the w-component
			return _mm_and_ps(V, g_Mask3);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreInt3(Int3* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->x = static_cast<int32_t>(v.vector4_f32[0]);
			pDestination->y = static_cast<int32_t>(v.vector4_f32[1]);
			pDestination->z = static_cast<int32_t>(v.vector4_f32[2]);

#elif defined(_SSE2_INTRINSICS_)
			// Detect positive overflow
			VECTOR vOverflow = _mm_cmpgt_ps(v, g_MaxInt);

			// Convert float to int
			__m128i vResulti = _mm_cvttps_epi32(v);

			// Set overflow elements to 0x7FFFFFFF
			VECTOR vResult = _mm_and_ps(vOverflow, g_AbsMask);
			vOverflow = _mm_andnot_ps(vOverflow, _mm_castsi128_ps(vResulti));
			vOverflow = _mm_or_ps(vOverflow, vResult);

			// Write x, y, z
			_mm_store_sd(reinterpret_cast<double*>(pDestination), _mm_castps_pd(vOverflow));
			__m128 z = PERMUTE_PS(vOverflow, _MM_SHUFFLE(2, 2, 2, 2));
			_mm_store_ss(reinterpret_cast<float*>(&pDestination->z), z);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreAInt3(AInt3* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
			assert((reinterpret_cast<uintptr_t>(pDestination) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->x = static_cast<int32_t>(v.vector4_f32[0]);
			pDestination->y = static_cast<int32_t>(v.vector4_f32[1]);
			pDestination->z = static_cast<int32_t>(v.vector4_f32[2]);

#elif defined(_SSE2_INTRINSICS_)
			// Detect positive overflow
			VECTOR vOverflow = _mm_cmpgt_ps(v, g_MaxInt);

			// Convert float to int
			__m128i vResulti = _mm_cvttps_epi32(v);

			// Set overflow elements to 0x7FFFFFFF
			VECTOR vResult = _mm_and_ps(vOverflow, g_AbsMask);
			vOverflow = _mm_andnot_ps(vOverflow, _mm_castsi128_ps(vResulti));
			vOverflow = _mm_or_ps(vOverflow, vResult);

			// Write x, y, z
			_mm_store_sd(reinterpret_cast<double*>(pDestination), _mm_castps_pd(vOverflow));
			__m128 z = _mm_movehl_ps(v, v);
			_mm_store_ss(reinterpret_cast<float*>(&pDestination->z), z);
#endif
		}
	}
}

#endif // !ULTREALITY_MATH_INT3_INL
