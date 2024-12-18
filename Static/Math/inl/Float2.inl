#ifndef ULTREALITY_MATH_FLOAT2_INL
#define ULTREALITY_MATH_FLOAT2_INL

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
	FORCE_INLINE Float2::Float2(A_VECTOR v) noexcept
	{
		VEC::StoreFloat2(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE Float2& VEC_CALLCONV Float2::operator=(A_VECTOR v) noexcept
	{
		VEC::StoreFloat2(this, v);

		return *this;
	}

	FORCE_INLINE VECTOR VEC_CALLCONV Float2::Load() noexcept
	{
		return VEC::LoadFloat2(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV Float2::Store(A_VECTOR v) noexcept
	{
		VEC::StoreFloat2(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE AFloat2::AFloat2(A_VECTOR v) noexcept
	{
		VEC::StoreAFloat2(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE AFloat2& VEC_CALLCONV AFloat2::operator=(A_VECTOR v) noexcept
	{
		VEC::StoreAFloat2(this, v);

		return *this;
	}

	FORCE_INLINE VECTOR VEC_CALLCONV AFloat2::Load() noexcept
	{
		return VEC::LoadAFloat2(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV AFloat2::Store(A_VECTOR v) noexcept
	{
		VEC::StoreAFloat2(this, v);
	}

	namespace VEC
	{
		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadFloat2(const Float2* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR V;
			V.vector4_f32[0] = pSource->x;
			V.vector4_f32[1] = pSource->y;
			V.vector4_f32[2] = 0.0f;
			V.vector4_f32[3] = 0.0f;

			return V;
#elif defined(_SSE2_INTRINSICS_)
			return _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(pSource)));
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadAFloat2(const AFloat2* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
			assert((reinterpret_cast<uintptr_t>(pSource) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR V;
			V.vector4_f32[0] = pSource->x;
			V.vector4_f32[1] = pSource->y;
			V.vector4_f32[2] = 0.0f;
			V.vector4_f32[3] = 0.0f;

			return V;
#elif defined(_SSE2_INTRINSICS_)
			// Use aligned load for 16-byte aligned AFLoat2
			return _mm_load_ps(reinterpret_cast<const float*>(pSource));
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreFloat2(Float2* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->x = v.vector4_f32[0];
			pDestination->y = v.vector4_f32[1];

#elif defined(_SSE2_INTRINSICS_)
			_mm_store_sd(reinterpret_cast<double*>(pDestination), _mm_castps_pd(v));
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreAFloat2(AFloat2* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
			assert((reinterpret_cast<uintptr_t>(pDestination) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->x = v.vector4_f32[0];
			pDestination->y = v.vector4_f32[1];

#elif defined(_SSE2_INTRINSICS_)
			_mm_store_sd(reinterpret_cast<double*>(pDestination), _mm_castps_pd(v));
#endif
		}
	}
}

#endif // !ULTREALITY_MATH_FLOAT2_INL
