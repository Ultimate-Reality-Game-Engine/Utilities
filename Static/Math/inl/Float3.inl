#ifndef ULTREALITY_MATH_FLOAT3_INL
#define ULTREALITY_MATH_FLOAT3_INL

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
	FORCE_INLINE Float3::Float3(A_VECTOR v) noexcept
	{
		VEC::StoreFloat3(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE Float3& VEC_CALLCONV Float3::operator=(A_VECTOR v) noexcept
	{
		VEC::StoreFloat3(this, v);

		return *this;
	}

	FORCE_INLINE VECTOR VEC_CALLCONV Float3::Load() noexcept
	{
		return VEC::LoadFloat3(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV Float3::Store(A_VECTOR v) noexcept
	{
		VEC::StoreFloat3(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE AFloat3::AFloat3(A_VECTOR v) noexcept
	{
		VEC::StoreAFloat3(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE AFloat3& VEC_CALLCONV AFloat3::operator=(A_VECTOR v) noexcept
	{
		VEC::StoreAFloat3(this, v);

		return *this;
	}

	FORCE_INLINE VECTOR VEC_CALLCONV AFLoat3::Load() noexcept
	{
		return VEC::LoadAFloat3(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV AFloat3::Store(A_VECTOR v) noexcept
	{
		VEC::StoreAFloat3(this, v);
	}

	namespace VEC
	{
		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadFloat3(const Float3* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR V;
			V.vector4_f32[0] = pSource->x;
			V.vector4_f32[1] = pSource->y;
			V.vector4_f32[2] = pSource->z;
			V.vector4_f32[3] = 0.0f;

			return V;
#elif defined(_SSE2_INTRINSICS_)
			__m128 xy = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(pSource)));
			__m128 z = _mm_load_ss(&pSource->z);
			return _mm_movelh_ps(xy, z);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadAFloat3(const AFloat3* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
			assert((reinterpret_cast<uintptr_t>(pSource) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR V;
			V.vector4_f32[0] = pSource->x;
			V.vector4_f32[1] = pSource->y;
			V.vector4_f32[2] = pSource->z;
			V.vector4_f32[3] = 0.0f;

			return V;
#elif defined(_SSE2_INTRINSICS_)
			// Load all 16 bytes into __m128
			__m128 V = _mm_load_ps(&pSource->x);
			return _mm_and_ps(V, g_Mask3);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreFloat3(Float3* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->x = v.vector4_f32[0];
			pDestination->y = v.vector4_f32[1];
			pDestination->z = v.vector4_f32[2];

#elif defined(_SSE2_INTRINSICS_)
			_mm_store_sd(reinterpret_cast<double*>(pDestination)m _mm_castps_pd(v));
			__m128 z = PERMUTE_PS(v, _MM_SHUFFLE(2, 2, 2, 2));
			_mm_store_ss(&pDestination->z, z);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreAFloat3(AFloat3* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
			assert((reinterpret_cast<uintptr_t>(pDestination) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->x = v.vector4_f32[0];
			pDestination->y = v.vector4_f32[1];
			pDestination->z = v.vector4_f32[2];

#elif defined(_SSE2_INTRINSICS_)
			_mm_store_sd(reinterpret_cast<double*>(pDestination)m _mm_castps_pd(v));
			__m128 z = _mm_movehl_ps(v, v);
			_mm_store_ss(&pDestination->z, z);
#endif
		}
	}
}

#endif // !ULTREALITY_MATH_FLOAT3_INL
