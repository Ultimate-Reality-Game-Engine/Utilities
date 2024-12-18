#ifndef ULTREALITY_MATH_FLOAT4_INL
#define ULTREALITY_MATH_FLOAT4_INL

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
	FORCE_INLINE Float4::Float4(A_VECTOR v) noexcept
	{
		VEC::StoreFloat4(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE Float4& VEC_CALLCONV Float4::operator=(A_VECTOR v) noexcept
	{
		VEC::StoreFloat4(this, v);

		return *this;
	}

	FORCE_INLINE VECTOR VEC_CALLCONV Float4::Load() noexcept
	{
		return VEC::LoadFloat4(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV Float4::Store(A_VECTOR v) noexcept
	{
		VEC::StoreFloat4(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE AFloat4::AFloat4(A_VECTOR v) noexcept
	{
		VEC::StoreAFloat4(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE AFloat4& VEC_CALLCONV AFloat4::operator=(A_VECTOR v) noexcept
	{
		VEC::StoreAFloat4(this, v);

		return *this;
	}

	FORCE_INLINE VECTOR VEC_CALLCONV AFloat4::Load() noexcept
	{
		return VEC::LoadAFloat4(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV AFloat4::Store(A_VECTOR v) noexcept
	{
		VEC::StoreAFloat4(this, v);
	}

	namespace VEC
	{
		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadFloat4(const Float4* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR V;
			V.vector4_f32[0] = pSource->x;
			V.vector4_f32[1] = pSource->y;
			V.vector4_f32[2] = pSource->z;
			V.vector4_f32[3] = pSource->w;

			return V;
#elif defined(_SSE2_INTRINSICS_)
			return _mm_loadu_ps(&pSource->x);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadAFloat4(const AFloat4* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
			assert((reinterpret_cast<uintptr_t>(pSource) & 0xf) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR V;
			V.vector4_f32[0] = pSource->x;
			V.vector4_f32[1] = pSource->y;
			V.vector4_f32[2] = pSource->z;
			V.vector4_f32[3] = pSource->w;

			return V;
#elif defined(_SSE2_INTRINSICS_)
			return _mm_load_ps(&pSource->x);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreFloat4(Float4* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->x = v.vector4_f32[0];
			pDestination->y = v.vector4_f32[1];
			pDestination->z = v.vector4_f32[2];
			pDestination->w = v.vector4_f32[3];

#elif defined(_SSE2_INTRINSICS_)
			_mm_storeu_ps(&pDestination->x, v);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreAFloat4(AFloat4* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
			assert((reinterpret_cast<uintptr_t>(pDestination) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination->x = v.vector4_f32[0];
			pDestination->y = v.vector4_f32[1];
			pDestination->z = v.vector4_f32[2];
			pDestination->w = v.vector4_f32[3];

#elif defined(_SSE2_INTRINSICS_)
			_mm_store_ps(&pDestination->x, v);
#endif
		}
	}
}

#endif // !ULTREALITY_MATH_FLOAT4_INL
