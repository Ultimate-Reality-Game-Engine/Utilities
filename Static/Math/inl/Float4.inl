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

	namespace VEC4
	{
		FORCE_INLINE bool VEC_CALLCONV Equal(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return (((V1.vector4_f32[0] == V2.vector4_f32[0]) && (V1.vector4_f32[1] == V2.vector4_f32[1]) && (V1.vector4_f32[2] == V2.vector4_f32[2]) && (V1.vector4_f32[3] == V2.vector4_f32[3])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmpeq_ps(V1, V2);

			return ((_mm_movemask_ps(vTemp) == 0x0f) != 0);
#endif
		}

		FORCE_INLINE uint32_t VEC_CALLCONV EqualR(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			uint32_t CR = 0;

			if ((V1.vector4_f32[0] == V2.vector4_f32[0]) &&
				(V1.vector4_f32[1] == V2.vector4_f32[1]) &&
				(V1.vector4_f32[2] == V2.vector4_f32[2]) &&
				(V1.vector4_f32[3] == V2.vector4_f32[3]))
			{
				CR = CRMASK_CR6TRUE;
			}
			else if ((V1.vector4_f32[0] != V2.vector4_f32[0]) &&
				(V1.vector4_f32[1] != V2.vector4_f32[1]) &&
				(V1.vector4_f32[2] != V2.vector4_f32[2]) &&
				(V1.vector4_f32[3] != V2.vector4_f32[3]))
			{
				CR = CRMASK_CR6FALSE;
			}

			return CR;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmpeq_ps(V1, V2);
			int iTest = _mm_movemask_ps(vTemp);
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

		FORCE_INLINE bool VEC_CALLCONV NearEqual(A_VECTOR V1, A_VECTOR V2, A_VECTOR epsilon) noexcept
		{
#if defined(_NO_INTRINSICS_)
			float dx, dy, dz, dw;

			dx = fabsf(V1.vector4_f32[0] - V2.vector4_f32[0]);
			dy = fabsf(V1.vector4_f32[1] - V2.vector4_f32[1]);
			dz = fabsf(V1.vector4_f32[2] - V2.vector4_f32[2]);
			dw = fabsf(V1.vector4_f32[3] - V2.vector4_f32[3]);
			
			return (((dx <= epsilon.vector4_f32[0]) &&
				(dy <= epsilon.vector4_f32[1]) &&
				(dz <= epsilon.vector4_f32[2]) &&
				(dw <= epsilon.vector4_f32[3])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			// Get the difference
			VECTOR vDelta = _mm_sub_ps(V1, V2);
			// Get the absolute value of the difference
			VECTOR vTemp = _mm_setzero_ps();
			vTemp = _mm_sub_ps(vTemp, vDelta);
			vTemp = _mm_max_ps(vTemp, vDelta);
			vTemp = _mm_cmple_ps(vTemp, epsilon);

			return ((_mm_movemask_ps(vTemp) == 0xf) != 0);
#endif
		}

		FORCE_INLINE bool VEC_CALLCONV NotEqual(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
    		return (((V1.vector4_f32[0] != V2.vector4_f32[0]) || (V1.vector4_f32[1] != V2.vector4_f32[1]) || (V1.vector4_f32[2] != V2.vector4_f32[2]) || (V1.vector4_f32[3] != V2.vector4_f32[3])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmpneq_ps(V1, V2);

			return (_mm_movemask_ps(vTemp) != 0);
#endif
		}

		FORCE_INLINE bool VEC_CALLCONV Greater(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
    		return (((V1.vector4_f32[0] > V2.vector4_f32[0]) && (V1.vector4_f32[1] > V2.vector4_f32[1]) && (V1.vector4_f32[2] > V2.vector4_f32[2]) && (V1.vector4_f32[3] > V2.vector4_f32[3])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmpgt_ps(V1, V2);

			return ((_mm_movemask_ps(vTemp) == 0x0f) != 0);
#endif
		}

		FORCE_INLINE uint32_t VEC_CALLCONV GreaterR(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_XM_NO_INTRINSICS_)
			uint32_t CR = 0;
			if (V1.vector4_f32[0] > V2.vector4_f32[0] &&
				V1.vector4_f32[1] > V2.vector4_f32[1] &&
				V1.vector4_f32[2] > V2.vector4_f32[2] &&
				V1.vector4_f32[3] > V2.vector4_f32[3])
			{
				CR = CRMASK_CR6TRUE;
			}
			else if (V1.vector4_f32[0] <= V2.vector4_f32[0] &&
				V1.vector4_f32[1] <= V2.vector4_f32[1] &&
				V1.vector4_f32[2] <= V2.vector4_f32[2] &&
				V1.vector4_f32[3] <= V2.vector4_f32[3])
			{
				CR = CRMASK_CR6FALSE;
			}

			return CR;

#elif defined(_SSE2_INTRINSICS_)
			uint32_t CR = 0;
			VECTOR vTemp = _mm_cmpgt_ps(V1, V2);
			int iTest = _mm_movemask_ps(vTemp);
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

		FORCE_INLINE bool VEC_CALLCONV GreaterOrEqual(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
    		return (((V1.vector4_f32[0] >= V2.vector4_f32[0]) && (V1.vector4_f32[1] >= V2.vector4_f32[1]) && (V1.vector4_f32[2] >= V2.vector4_f32[2]) && (V1.vector4_f32[3] >= V2.vector4_f32[3])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmpge_ps(V1, V2);
    		
			return ((_mm_movemask_ps(vTemp) == 0x0f) != 0);
#endif
		}

		FORCE_INLINE uint32_t GreaterOrEqualR(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			uint32_t CR = 0;
			if ((V1.vector4_f32[0] >= V2.vector4_f32[0]) &&
				(V1.vector4_f32[1] >= V2.vector4_f32[1]) &&
				(V1.vector4_f32[2] >= V2.vector4_f32[2]) &&
				(V1.vector4_f32[3] >= V2.vector4_f32[3]))
			{
				CR = CRMASK_CR6TRUE;
			}
			else if ((V1.vector4_f32[0] < V2.vector4_f32[0]) &&
				(V1.vector4_f32[1] < V2.vector4_f32[1]) &&
				(V1.vector4_f32[2] < V2.vector4_f32[2]) &&
				(V1.vector4_f32[3] < V2.vector4_f32[3]))
			{
				CR = CRMASK_CR6FALSE;
			}

			return CR;

#elif defined(_SSE2_INTRINSICS_)
			uint32_t CR = 0;
			VECTOR vTemp = _mm_cmpge_ps(V1, V2);
			int iTest = _mm_movemask_ps(vTemp);
			if (iTest == 0x0f)
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

		FORCE_INLINE bool VEC_CALLCONV Less(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
    		return (((V1.vector4_f32[0] < V2.vector4_f32[0]) && (V1.vector4_f32[1] < V2.vector4_f32[1]) && (V1.vector4_f32[2] < V2.vector4_f32[2]) && (V1.vector4_f32[3] < V2.vector4_f32[3])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmplt_ps(V1, V2);
    		
			return ((_mm_movemask_ps(vTemp) == 0x0f) != 0);
#endif
		}

		FORCE_INLINE bool VEC_CALLCONV LessOrEqual(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
    		return (((V1.vector4_f32[0] <= V2.vector4_f32[0]) && (V1.vector4_f32[1] <= V2.vector4_f32[1]) && (V1.vector4_f32[2] <= V2.vector4_f32[2]) && (V1.vector4_f32[3] <= V2.vector4_f32[3])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmple_ps(V1, V2);
    		
			return ((_mm_movemask_ps(vTemp) == 0x0f) != 0);
#endif
		}

		FORCE_INLINE bool VEC_CALLCONV InBounds(A_VECTOR v, A_VECTOR bounds) noexcept
		{
#if defined(_XM_NO_INTRINSICS_)
			return (((v.vector4_f32[0] <= bounds.vector4_f32[0] && v.vector4_f32[0] >= -bounds.vector4_f32[0]) &&
				(v.vector4_f32[1] <= bounds.vector4_f32[1] && v.vector4_f32[1] >= -bounds.vector4_f32[1]) &&
				(v.vector4_f32[2] <= bounds.vector4_f32[2] && v.vector4_f32[2] >= -bounds.vector4_f32[2]) &&
				(v.vector4_f32[3] <= bounds.vector4_f32[3] && v.vector4_f32[3] >= -bounds.vector4_f32[3])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			// Test if less than or equal
			VECTOR vTemp1 = _mm_cmple_ps(v, bounds);
			// Negate the bounds
			VECTOR vTemp2 = _mm_mul_ps(bounds, g_NegativeOne);
			// Test if greater or equal (Reversed)
			vTemp2 = _mm_cmple_ps(vTemp2, v);
			// Blend answers
			vTemp1 = _mm_and_ps(vTemp1, vTemp2);
			// All in bounds?
			return ((_mm_movemask_ps(vTemp1) == 0x0f) != 0);
#endif
		}

#if !defined(_NO_INTRINSICS_) && defined(_MSC_VER) && !defined(__INTEL_COMPILER)
#pragma float_control(push)
#pragma float_control(precise, on)
#endif

		FORCE_INLINE bool VEC_CALLCONV IsNaN(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return (ISNAN(v.vector4_f32[0]) ||
				ISNAN(v.vector4_f32[1]) ||
				ISNAN(v.vector4_f32[2]) ||
				ISNAN(v.vector4_f32[3]));

#elif defined(_SSE2_INTRINSICS_)
#if defined(__clang__) && defined(__FINITE_MATH_ONLY__)
			ALIGNED(16) float tmp[4];
			_mm_store_ps(tmp, v);
			
			return isnan(tmp[0]) || isnan(tmp[1]) || isnan(tmp[2]) || isnan(tmp[3]);
#else
			// Test against itself. NaN is always not equal
			VECTOR vTempNan = _mm_cmpneq_ps(v, v);
			
			// If any are NaN, the mask is non-zero
			return (_mm_movemask_ps(vTempNan) != 0);
#endif
#endif
		}

#if !defined(_NO_INTRINSICS_) && defined(_MSC_VER) && !defined(__INTEL_COMPILER)
#pragma float_control(pop)
#endif

		FORCE_INLINE bool VEC_CALLCONV IsInfinite(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return (ISINF(V.vector4_f32[0]) ||
				ISINF(V.vector4_f32[1]) ||
				ISINF(V.vector4_f32[2]) ||
				ISINF(V.vector4_f32[3]));

#elif defined(_SSE2_INTRINSICS_)
			// Mask off the sign bit
			VECTOR vTemp = _mm_and_ps(v, g_AbsMask);
			// Compare to infinity
			vTemp = _mm_cmpeq_ps(vTemp, g_Infinity);
			
			// If any are infinity, the signs are true.
			return (_mm_movemask_ps(vTemp) != 0);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Dot(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 Result;
			Result.f[0] =
				Result.f[1] =
				Result.f[2] =
				Result.f[3] = V1.vector4_f32[0] * V2.vector4_f32[0] + V1.vector4_f32[1] * V2.vector4_f32[1] + V1.vector4_f32[2] * V2.vector4_f32[2] + V1.vector4_f32[3] * V2.vector4_f32[3];
			
			return Result.v;

#elif defined(_SSE4_INTRINSICS_)
			return _mm_dp_ps(V1, V2, 0xff);

#elif defined(_SSE3_INTRINSICS_)
			VECTOR vTemp = _mm_mul_ps(V1, V2);
			vTemp = _mm_hadd_ps(vTemp, vTemp);

			return _mm_hadd_ps(vTemp, vTemp);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp2 = V2;
			VECTOR vTemp = _mm_mul_ps(V1, vTemp2);
			vTemp2 = _mm_shuffle_ps(vTemp2, vTemp, _MM_SHUFFLE(1, 0, 0, 0)); // Copy X to the Z position and Y to the W position
			vTemp2 = _mm_add_ps(vTemp2, vTemp); // Add Z = X+Z; W = Y+W;
			vTemp = _mm_shuffle_ps(vTemp, vTemp2, _MM_SHUFFLE(0, 3, 0, 0)); // Copy W to the Z position
			vTemp = _mm_add_ps(vTemp, vTemp2); // Add Z and W together
			
			return PERMUTE_PS(vTemp, _MM_SHUFFLE(2, 2, 2, 2)); // Splat Z and return
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Cross(A_VECTOR V1, A_VECTOR V2, A_VECTOR V3) noexcept
		{
			// [ ((v2.z*v3.w-v2.w*v3.z)*v1.y)-((v2.y*v3.w-v2.w*v3.y)*v1.z)+((v2.y*v3.z-v2.z*v3.y)*v1.w),
			//   ((v2.w*v3.z-v2.z*v3.w)*v1.x)-((v2.w*v3.x-v2.x*v3.w)*v1.z)+((v2.z*v3.x-v2.x*v3.z)*v1.w),
			//   ((v2.y*v3.w-v2.w*v3.y)*v1.x)-((v2.x*v3.w-v2.w*v3.x)*v1.y)+((v2.x*v3.y-v2.y*v3.x)*v1.w),
			//   ((v2.z*v3.y-v2.y*v3.z)*v1.x)-((v2.z*v3.x-v2.x*v3.z)*v1.y)+((v2.y*v3.x-v2.x*v3.y)*v1.z) ]
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 Result = { { {
				(((V2.vector4_f32[2] * V3.vector4_f32[3]) - (V2.vector4_f32[3] * V3.vector4_f32[2])) * V1.vector4_f32[1]) - (((V2.vector4_f32[1] * V3.vector4_f32[3]) - (V2.vector4_f32[3] * V3.vector4_f32[1])) * V1.vector4_f32[2]) + (((V2.vector4_f32[1] * V3.vector4_f32[2]) - (V2.vector4_f32[2] * V3.vector4_f32[1])) * V1.vector4_f32[3]),
				(((V2.vector4_f32[3] * V3.vector4_f32[2]) - (V2.vector4_f32[2] * V3.vector4_f32[3])) * V1.vector4_f32[0]) - (((V2.vector4_f32[3] * V3.vector4_f32[0]) - (V2.vector4_f32[0] * V3.vector4_f32[3])) * V1.vector4_f32[2]) + (((V2.vector4_f32[2] * V3.vector4_f32[0]) - (V2.vector4_f32[0] * V3.vector4_f32[2])) * V1.vector4_f32[3]),
				(((V2.vector4_f32[1] * V3.vector4_f32[3]) - (V2.vector4_f32[3] * V3.vector4_f32[1])) * V1.vector4_f32[0]) - (((V2.vector4_f32[0] * V3.vector4_f32[3]) - (V2.vector4_f32[3] * V3.vector4_f32[0])) * V1.vector4_f32[1]) + (((V2.vector4_f32[0] * V3.vector4_f32[1]) - (V2.vector4_f32[1] * V3.vector4_f32[0])) * V1.vector4_f32[3]),
				(((V2.vector4_f32[2] * V3.vector4_f32[1]) - (V2.vector4_f32[1] * V3.vector4_f32[2])) * V1.vector4_f32[0]) - (((V2.vector4_f32[2] * V3.vector4_f32[0]) - (V2.vector4_f32[0] * V3.vector4_f32[2])) * V1.vector4_f32[1]) + (((V2.vector4_f32[1] * V3.vector4_f32[0]) - (V2.vector4_f32[0] * V3.vector4_f32[1])) * V1.vector4_f32[2]),
			} } };
			
			return Result.v;

#elif defined(_SSE2_INTRINSICS_)
			// V2zwyz * V3wzwy
			VECTOR vResult = PERMUTE_PS(V2, _MM_SHUFFLE(2, 1, 3, 2));
			VECTOR vTemp3 = PERMUTE_PS(V3, _MM_SHUFFLE(1, 3, 2, 3));
			vResult = _mm_mul_ps(vResult, vTemp3);
			// - V2wzwy * V3zwyz
			VECTOR vTemp2 = PERMUTE_PS(V2, _MM_SHUFFLE(1, 3, 2, 3));
			vTemp3 = PERMUTE_PS(vTemp3, _MM_SHUFFLE(1, 3, 0, 1));
			vResult = FNMADD_PS(vTemp2, vTemp3, vResult);
			// term1 * V1yxxx
			VECTOR vTemp1 = PERMUTE_PS(V1, _MM_SHUFFLE(0, 0, 0, 1));
			vResult = _mm_mul_ps(vResult, vTemp1);

			// V2ywxz * V3wxwx
			vTemp2 = PERMUTE_PS(V2, _MM_SHUFFLE(2, 0, 3, 1));
			vTemp3 = PERMUTE_PS(V3, _MM_SHUFFLE(0, 3, 0, 3));
			vTemp3 = _mm_mul_ps(vTemp3, vTemp2);
			// - V2wxwx * V3ywxz
			vTemp2 = PERMUTE_PS(vTemp2, _MM_SHUFFLE(2, 1, 2, 1));
			vTemp1 = PERMUTE_PS(V3, _MM_SHUFFLE(2, 0, 3, 1));
			vTemp3 = FNMADD_PS(vTemp2, vTemp1, vTemp3);
			// vResult - temp * V1zzyy
			vTemp1 = PERMUTE_PS(V1, _MM_SHUFFLE(1, 1, 2, 2));
			vResult = FNMADD_PS(vTemp1, vTemp3, vResult);

			// V2yzxy * V3zxyx
			vTemp2 = PERMUTE_PS(V2, _MM_SHUFFLE(1, 0, 2, 1));
			vTemp3 = PERMUTE_PS(V3, _MM_SHUFFLE(0, 1, 0, 2));
			vTemp3 = _mm_mul_ps(vTemp3, vTemp2);
			// - V2zxyx * V3yzxy
			vTemp2 = PERMUTE_PS(vTemp2, _MM_SHUFFLE(2, 0, 2, 1));
			vTemp1 = PERMUTE_PS(V3, _MM_SHUFFLE(1, 0, 2, 1));
			vTemp3 = FNMADD_PS(vTemp1, vTemp2, vTemp3);
			// vResult + term * V1wwwz
			vTemp1 = PERMUTE_PS(V1, _MM_SHUFFLE(2, 3, 3, 3));
			vResult = FMADD_PS(vTemp3, vTemp1, vResult);
			
			return vResult;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV LengthSq(A_VECTOR v) noexcept
		{
			return Dot(v, v);
		}

		FORCE_INLINE VECTOR VEC_CALLCONV ReciprocalLengthEst(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			XMVECTOR Result;

			Result = LengthSq(v);
			Result = VEC::ReciprocalSqrtEst(Result);

			return Result;

#elif defined(_SSE4_INTRINSICS_)
			VECTOR vTemp = _mm_dp_ps(v, v, 0xff);

			return _mm_rsqrt_ps(vTemp);

#elif defined(_SSE3_INTRINSICS_)
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
			vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
			
			return _mm_rsqrt_ps(vLengthSq);

#elif defined(_SSE2_INTRINSICS_)
			// Perform the dot product on x,y,z and w
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			// vTemp has z and w
			VECTOR vTemp = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(3, 2, 3, 2));
			// x+z, y+w
			vLengthSq = _mm_add_ps(vLengthSq, vTemp);
			// x+z,x+z,x+z,y+w
			vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(1, 0, 0, 0));
			// ??,??,y+w,y+w
			vTemp = _mm_shuffle_ps(vTemp, vLengthSq, _MM_SHUFFLE(3, 3, 0, 0));
			// ??,??,x+z+y+w,??
			vLengthSq = _mm_add_ps(vLengthSq, vTemp);
			// Splat the length
			vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(2, 2, 2, 2));
			
			// Get the reciprocal
			return _mm_rsqrt_ps(vLengthSq);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV ReciprocalLength(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR Result;

			Result = LengthSq(v);
			Result = VEC::ReciprocalSqrt(Result);

			return Result;

#elif defined(_SSE4_INTRINSICS_)
			VECTOR vTemp = _mm_dp_ps(v, v, 0xff);
			VECTOR vLengthSq = _mm_sqrt_ps(vTemp);

			return _mm_div_ps(g_One, vLengthSq);

#elif defined(_SSE3_INTRINSICS_)
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
			vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
			vLengthSq = _mm_sqrt_ps(vLengthSq);
			
			return _mm_div_ps(g_One, vLengthSq);

#elif defined(_SSE2_INTRINSICS_)
			// Perform the dot product on x,y,z and w
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			// vTemp has z and w
			VECTOR vTemp = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(3, 2, 3, 2));
			// x+z, y+w
			vLengthSq = _mm_add_ps(vLengthSq, vTemp);
			// x+z,x+z,x+z,y+w
			vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(1, 0, 0, 0));
			// ??,??,y+w,y+w
			vTemp = _mm_shuffle_ps(vTemp, vLengthSq, _MM_SHUFFLE(3, 3, 0, 0));
			// ??,??,x+z+y+w,??
			vLengthSq = _mm_add_ps(vLengthSq, vTemp);
			// Splat the length
			vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(2, 2, 2, 2));
			// Get the reciprocal
			vLengthSq = _mm_sqrt_ps(vLengthSq);
			
			// Accurate!
			return _mm_div_ps(g_One, vLengthSq);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV LengthEst(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR Result = LengthSq(v);
			
			return VEC::SqrtEst(Result);

#elif defined(_SSE4_INTRINSICS_)
			VECTOR vTemp = _mm_dp_ps(v, v, 0xff);

			return _mm_sqrt_ps(vTemp);

#elif defined(_SSE3_INTRINSICS_)
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
			vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);

			return _mm_sqrt_ps(vLengthSq);

#elif defined(_SSE2_INTRINSICS_)
			// Perform the dot product on x,y,z and w
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			// vTemp has z and w
			VECTOR vTemp = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(3, 2, 3, 2));
			// x+z, y+w
			vLengthSq = _mm_add_ps(vLengthSq, vTemp);
			// x+z,x+z,x+z,y+w
			vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(1, 0, 0, 0));
			// ??,??,y+w,y+w
			vTemp = _mm_shuffle_ps(vTemp, vLengthSq, _MM_SHUFFLE(3, 3, 0, 0));
			// ??,??,x+z+y+w,??
			vLengthSq = _mm_add_ps(vLengthSq, vTemp);
			// Splat the length
			vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(2, 2, 2, 2));
			
			// Get the length
			return _mm_sqrt_ps(vLengthSq);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Length(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR Result = LengthSq(v);

			return VEC::Sqrt(Result);

#elif defined(_SSE4_INTRINSICS_)
			VECTOR vTemp = _mm_dp_ps(v, v, 0xff);

			return _mm_sqrt_ps(vTemp);

#elif defined(_SSE3_INTRINSICS_)
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
			vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);

			return _mm_sqrt_ps(vLengthSq);

#elif defined(_SSE2_INTRINSICS_)
			// Perform the dot product on x,y,z and w
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			// vTemp has z and w
			VECTOR vTemp = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(3, 2, 3, 2));
			// x+z, y+w
			vLengthSq = _mm_add_ps(vLengthSq, vTemp);
			// x+z,x+z,x+z,y+w
			vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(1, 0, 0, 0));
			// ??,??,y+w,y+w
			vTemp = _mm_shuffle_ps(vTemp, vLengthSq, _MM_SHUFFLE(3, 3, 0, 0));
			// ??,??,x+z+y+w,??
			vLengthSq = _mm_add_ps(vLengthSq, vTemp);
			// Splat the length
			vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(2, 2, 2, 2));
			
			// Get the length
			return _mm_sqrt_ps(vLengthSq);
#endif
		}

		
	}
}

#endif // !ULTREALITY_MATH_FLOAT4_INL
