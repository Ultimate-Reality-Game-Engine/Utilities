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
		Vector::StoreFloat4(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE Float4& VEC_CALLCONV Float4::operator=(A_VECTOR v) noexcept
	{
		Vector::StoreFloat4(this, v);

		return *this;
	}

	FORCE_INLINE VECTOR VEC_CALLCONV Float4::Load() noexcept
	{
		return Vector::LoadFloat4(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV Float4::Store(A_VECTOR v) noexcept
	{
		Vector::StoreFloat4(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE AFloat4::AFloat4(A_VECTOR v) noexcept
	{
		Vector::StoreAFloat4(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE AFloat4& VEC_CALLCONV AFloat4::operator=(A_VECTOR v) noexcept
	{
		Vector::StoreAFloat4(this, v);

		return *this;
	}

	FORCE_INLINE VECTOR VEC_CALLCONV AFloat4::Load() noexcept
	{
		return Vector::LoadAFloat4(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV AFloat4::Store(A_VECTOR v) noexcept
	{
		Vector::StoreAFloat4(this, v);
	}

	namespace Vector
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

	namespace Vector4
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

		FORCE_INLINE uint32_t VEC_CALLCONV GreaterOrEqualR(A_VECTOR V1, A_VECTOR V2) noexcept
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
			Result = Vector::ReciprocalSqrtEst(Result);

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
			Result = Vector::ReciprocalSqrt(Result);

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
			
			return Vector::SqrtEst(Result);

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

			return Vector::Sqrt(Result);

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

		FORCE_INLINE VECTOR VEC_CALLCONV NormalizeEst(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR Result = ReciprocalLength(v);

			return Vector::Multiply(v, Result);

#elif defined(_SSE4_INTRINSICS_)
			VECTOR vTemp = _mm_dp_ps(v, v, 0xff);
			VECTOR vResult = _mm_rsqrt_ps(vTemp);

			return _mm_mul_ps(vResult, v);

#elif defined(_SSE3_INTRINSICS_)
			VECTOR vDot = _mm_mul_ps(v, v);
			vDot = _mm_hadd_ps(vDot, vDot);
			vDot = _mm_hadd_ps(vDot, vDot);
			vDot = _mm_rsqrt_ps(vDot);

			return _mm_mul_ps(vDot, v);

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
			VECTOR vResult = _mm_rsqrt_ps(vLengthSq);
			
			// Reciprocal mul to perform the normalization
			return _mm_mul_ps(vResult, v);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Normalize(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			float fLength;
			VECTOR vResult;

			vResult = Length(v);
			fLength = vResult.vector4_f32[0];

			// Prevent divide by zero
			if (fLength > 0)
			{
				fLength = 1.0f / fLength;
			}

			vResult.vector4_f32[0] = v.vector4_f32[0] * fLength;
			vResult.vector4_f32[1] = v.vector4_f32[1] * fLength;
			vResult.vector4_f32[2] = v.vector4_f32[2] * fLength;
			vResult.vector4_f32[3] = v.vector4_f32[3] * fLength;
			
			return vResult;

#elif defined(_SSE4_INTRINSICS_)
			VECTOR vLengthSq = _mm_dp_ps(v, v, 0xff);
			// Prepare for the division
			VECTOR vResult = _mm_sqrt_ps(vLengthSq);
			// Create zero with a single instruction
			VECTOR vZeroMask = _mm_setzero_ps();
			// Test for a divide by zero (Must be FP to detect -0.0)
			vZeroMask = _mm_cmpneq_ps(vZeroMask, vResult);
			// Failsafe on zero (Or epsilon) length planes
			// If the length is infinity, set the elements to zero
			vLengthSq = _mm_cmpneq_ps(vLengthSq, g_Infinity);
			// Divide to perform the normalization
			vResult = _mm_div_ps(v, vResult);
			// Any that are infinity, set to zero
			vResult = _mm_and_ps(vResult, vZeroMask);
			// Select qnan or result based on infinite length
			VECTOR vTemp1 = _mm_andnot_ps(vLengthSq, g_QNaN);
			VECTOR vTemp2 = _mm_and_ps(vResult, vLengthSq);
			
			return _mm_or_ps(vTemp1, vTemp2);

#elif defined(_SSE3_INTRINSICS_)
			// Perform the dot product on x,y,z and w
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
			vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
			// Prepare for the division
			VECTOR vResult = _mm_sqrt_ps(vLengthSq);
			// Create zero with a single instruction
			VECTOR vZeroMask = _mm_setzero_ps();
			// Test for a divide by zero (Must be FP to detect -0.0)
			vZeroMask = _mm_cmpneq_ps(vZeroMask, vResult);
			// Failsafe on zero (Or epsilon) length planes
			// If the length is infinity, set the elements to zero
			vLengthSq = _mm_cmpneq_ps(vLengthSq, g_Infinity);
			// Divide to perform the normalization
			vResult = _mm_div_ps(v, vResult);
			// Any that are infinity, set to zero
			vResult = _mm_and_ps(vResult, vZeroMask);
			// Select qnan or result based on infinite length
			VECTOR vTemp1 = _mm_andnot_ps(vLengthSq, g_QNaN);
			VECTOR vTemp2 = _mm_and_ps(vResult, vLengthSq);
			
			return _mm_or_ps(vTemp1, vTemp2);

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
			// Prepare for the division
			VECTOR vResult = _mm_sqrt_ps(vLengthSq);
			// Create zero with a single instruction
			VECTOR vZeroMask = _mm_setzero_ps();
			// Test for a divide by zero (Must be FP to detect -0.0)
			vZeroMask = _mm_cmpneq_ps(vZeroMask, vResult);
			// Failsafe on zero (Or epsilon) length planes
			// If the length is infinity, set the elements to zero
			vLengthSq = _mm_cmpneq_ps(vLengthSq, g_Infinity);
			// Divide to perform the normalization
			vResult = _mm_div_ps(v, vResult);
			// Any that are infinity, set to zero
			vResult = _mm_and_ps(vResult, vZeroMask);
			// Select qnan or result based on infinite length
			VECTOR vTemp1 = _mm_andnot_ps(vLengthSq, g_QNaN);
			VECTOR vTemp2 = _mm_and_ps(vResult, vLengthSq);
			
			return _mm_or_ps(vTemp1, vTemp2);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV ClampLength(A_VECTOR v, float lengthMin, float lengthMax) noexcept
		{
			VECTOR clampMax = Vector::Replicate(lengthMax);
			VECTOR clampMin = Vector::Replicate(lengthMin);

			return ClampLengthV(v, clampMin, clampMax);
		}

		FORCE_INLINE VECTOR VEC_CALLCONV ClampLengthV(A_VECTOR v, A_VECTOR lengthMin, A_VECTOR lengthMax) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert((Vector::GetY(lengthMin) == Vector::GetX(lengthMin)) && (Vector::GetZ(lengthMin) == Vector::GetX(lengthMin)) && (Vector::GetW(lengthMin) == Vector::GetX(lengthMin)));
			assert((Vector::GetY(lengthMax) == Vector::GetX(lengthMax)) && (Vector::GetZ(lengthMax) == Vector::GetX(lengthMax)) && (Vector::GetW(lengthMax) == Vector::GetX(lengthMax)));
			assert(GreaterOrEqual(lengthMin, Vector::Zero()));
			assert(GreaterOrEqual(lengthMax, Vector::Zero()));
			assert(GreaterOrEqual(lengthMax, lengthMin));
#endif

			VECTOR lengthSq = LengthSq(v);

			const VECTOR zero = Vector::Zero();

			VECTOR rcpLength = Vector::ReciprocalSqrt(lengthSq);

			VECTOR infiniteLength = Vector::EqualInt(lengthSq, g_Infinity.v);
			VECTOR zeroLength = Vector::Equal(lengthSq, zero);

			VECTOR normal = Vector::Multiply(v, rcpLength);

			VECTOR length = Vector::Multiply(lengthSq, rcpLength);

			VECTOR select = Vector::EqualInt(infiniteLength, zeroLength);
			length = Vector::Select(lengthSq, length, select);
			normal = Vector::Select(lengthSq, normal, select);

			VECTOR controlMax = Vector::Greater(length, lengthMax);
			VECTOR controlMin = Vector::Less(length, lengthMin);

			VECTOR clampLength = Vector::Select(length, lengthMax, controlMax);
			clampLength = Vector::Select(clampLength, lengthMin, controlMin);

			VECTOR Result = Vector::Multiply(normal, clampLength);

			// Preserve the original vector (with no precision loss) if the length falls within the given range
			VECTOR control = Vector::EqualInt(controlMax, controlMin);
			Result = Vector::Select(Result, v, control);

			return Result;
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Reflect(A_VECTOR incident, A_VECTOR normal) noexcept
		{
			// Result = incident - (2 - dot(incident, normal)) * normal
			VECTOR Result = Dot(incident, normal);
			Result = Vector::Add(Result, Result);

			return Vector::NegativeMultiplySubtract(Result, normal, incident);
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Refract(A_VECTOR incident, A_VECTOR normal, float refractionIndex) noexcept
		{
			VECTOR index = Vector::Replicate(refractionIndex);

			return RefractV(incident, normal, index);
		}

		FORCE_INLINE VECTOR VEC_CALLCONV RefractV(A_VECTOR incident, A_VECTOR normal, A_VECTOR refractionIndex) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR IDotN = Dot(incident, normal);
			
			VECTOR R = Vector::NegativeMultiplySubtract(IDotN, IDotN, g_One.v);
			R = Vector::Multiply(R, refractionIndex);
			R = Vector::NegativeMultiplySubtract(R, refractionIndex, g_One.v);
			
			const VECTOR zero = Vector::Zero();

			if (LessOrEqual(R, zero))
			{
				// Total internal reflection
				return zero;
			}
			else
			{
				R = Vector::Sqrt(R);
				R = Vector::MultiplyAdd(refractionIndex, IDotN, R);

				VECTOR Result = Vector::Multiply(refractionIndex, incident);
				Result = Vector::NegativeMultiplySubtract(normal, R, Result);

				return Result;
			}

#elif defined(_SSE2_INTRINSICS_)
			VECTOR IDotN = Dot(incident, normal);

			// R = 1.0f - RefractionIndex * RefractionIndex * (1.0f - IDotN * IDotN)
			VECTOR R = FNMADD_PS(IDotN, IDotN, g_One);
			VECTOR R2 = _mm_mul_ps(refractionIndex, refractionIndex);
			R = FNMADD_PS(R, R2, g_One);

			VECTOR vResult = _mm_cmple_ps(R, g_Zero);
			if (_mm_movemask_ps(vResult) == 0x0f)
			{
				// Total internal reflection
				vResult = g_Zero;
			}
			else
			{
				// R = RefractionIndex * IDotN + sqrt(R)
				R = _mm_sqrt_ps(R);
				R = FMADD_PS(refractionIndex, IDotN, R);
				// Result = RefractionIndex * Incident - Normal * R
				vResult = _mm_mul_ps(refractionIndex, incident);
				vResult = FNMADD_PS(R, normal, vResult);
			}

			return vResult;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Orthogonal(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 Result = { { {
				v.vector_f32[2], 
				v.vector_f32[3], 
				-v.vector_f32[0], 
				-v.vector_f32[1]
			} } };

			return Result.v;

#elif defined(_SSE2_INTRINSICS_)
			static const VECTOR_F32 flipZW = { { { 1.0f, 1.0f, -1.0f, -1.0f } } };
			VECTOR vResult = PERMUTE_PS(v, _MM_SHUFFLE(1, 0, 3, 2));

			return _mm_mul_ps(vResult, flipZW);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV AngleBetweenNormalsEst(A_VECTOR N1, A_VECTOR N2) noexcept
		{
			VECTOR Result = Dot(N1, N2);
			Result = Vector::Clamp(Result, g_NegativeOne.v, g_One.v);

			return Vector::ACosEst(Result);
		}

		FORCE_INLINE VECTOR VEC_CALLCONV AngleBetweenNormals(A_VECTOR N1, A_VECTOR N2) noexcept
		{
			VECTOR Result = Dot(N1, N2);
			Result = Vector::Clamp(Result, g_NegativeOne.v, g_One.v);

			return Vector::ACos(Result);
		}

		FORCE_INLINE VECTOR VEC_CALLCONV AngleBetweenVectors(A_VECTOR V1, A_VECTOR V2) noexcept
		{
			VECTOR L1 = ReciprocalLength(V1);
			VECTOR L2 = ReciprocalLength(V2);

			VECTOR dot = Dot(V1, V2);

			L1 = Vector::Multiply(L1, L2);

			VECTOR cosAngle = Vector::Multiply(dot, L1);
			cosAngle = Vector::Clamp(cosAngle, g_NegativeOne.v, g_One.v);

			return Vector::ACos(cosAngle);
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Transform(A_VECTOR v, A_MATRIX m) noexcept
		{
#if defined(_XM_NO_INTRINSICS_)
			float fX = (m.m[0][0] * v.vector4_f32[0]) + (m.m[1][0] * v.vector4_f32[1]) + (m.m[2][0] * v.vector4_f32[2]) + (m.m[3][0] * v.vector4_f32[3]);
			float fY = (m.m[0][1] * v.vector4_f32[0]) + (m.m[1][1] * v.vector4_f32[1]) + (m.m[2][1] * v.vector4_f32[2]) + (m.m[3][1] * v.vector4_f32[3]);
			float fZ = (m.m[0][2] * v.vector4_f32[0]) + (m.m[1][2] * v.vector4_f32[1]) + (m.m[2][2] * v.vector4_f32[2]) + (m.m[3][2] * v.vector4_f32[3]);
			float fW = (m.m[0][3] * v.vector4_f32[0]) + (m.m[1][3] * v.vector4_f32[1]) + (m.m[2][3] * v.vector4_f32[2]) + (m.m[3][3] * v.vector4_f32[3]);
			VECTOR_F32 vResult = { { { fX, fY, fZ, fW } } };
			
			return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vResult = PERMUTE_PS(v, _MM_SHUFFLE(3, 3, 3, 3)); // W
			vResult = _mm_mul_ps(vResult, m.r[3]);
			VECTOR vTemp = PERMUTE_PS(v, _MM_SHUFFLE(2, 2, 2, 2)); // Z
			vResult = FMADD_PS(vTemp, m.r[2], vResult);
			vTemp = PERMUTE_PS(v, _MM_SHUFFLE(1, 1, 1, 1)); // Y
			vResult = FMADD_PS(vTemp, m.r[1], vResult);
			vTemp = PERMUTE_PS(v, _MM_SHUFFLE(0, 0, 0, 0)); // X
			vResult = FMADD_PS(vTemp, m.r[0], vResult);
			
			return vResult;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE Float4* VEC_CALLCONV TransformStream(Float4* pOutputStream, 
			size_t outputStride, 
			const Float4* pInputStream, 
			size_t inputStride, 
			size_t vectorCount, 
			A_MATRIX m) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pOutputStream != nullptr);
			assert(pInputStream != nullptr);

			assert(inputStride >= sizeof(Float4));
			_Analysis_assume_(inputStride >= sizeof(Float4));

			assert(outputStride >= sizeof(Float4));
			_Analysis_assume_(outputStride >= sizeof(Float4));
#endif

#if defined(_NO_INTRINSICS_)
			const uint8_t* pInputVector = reinterpret_cast<const uint8_t*>(pInputStream);
			uint8_t* pOutputVector = reinterpret_cast<uint8_t*>(pOutputStream);

			const VECTOR row0 = m.r[0];
			const VECTOR row1 = m.r[1];
			const VECTOR row2 = m.r[2];
			const VECTOR row3 = m.r[3];

			for (size_t i = 0; i < vectorCount; i++)
			{
				VECTOR V = Vector::LoadFloat4(reinterpret_cast<const Float4*>(pInputVector));
				VECTOR W = Vector::SplatW(v);
				VECTOR Z = Vector::SplatZ(v);
				VECTOR Y = Vector::SplatY(v);
				VECTOR X = Vector::SplatX(v);

				VECTOR Result = Vector::Multiply(W, row3);
				Result = Vector::MultiplyAdd(Z, row2, Result);
				Result = Vector::MultiplyAdd(Y, row1, Result);
				Result = Vector::MultiplyAdd(X, row0, Result);

#ifdef _PREFAST_
#pragma prefast(push)
#pragma prefast(disable : 26015, "PREfast noise: Esp:1307" )
#endif

        		Vector::StoreFloat4(reinterpret_cast<Float4*>(pOutputVector), Result);

#ifdef _PREFAST_
#pragma prefast(pop)
#endif

				pInputVector += inputStride;
				pOutputVector += outputStride;
			}

			return pOutputStream;

#elif defined(_AVX2_INTRINSICS_)
			const uint8_t* pInputVector = reinterpret_cast<const uint8_t*>(pInputStream);
			uint8_t* pOutputVector = reinterpret_cast<uint8_t*>(pOutputStream);

			size_t i = 0;
			size_t two = vectorCount >> 1;
			if (two > 0)
			{
				__m256 row0 = _mm256_broadcast_ps(&m.r[0]);
				__m256 row1 = _mm256_broadcast_ps(&m.r[1]);
				__m256 row2 = _mm256_broadcast_ps(&m.r[2]);
				__m256 row3 = _mm256_broadcast_ps(&m.r[3]);

				if (inputStride == sizeof(Float4))
				{
					if (outputStride == sizeof(Float4))
					{
						if (!(reinterpret_cast<uintptr_t>(pOutputStream) & 0x1F))
						{
							// Packed input, aligned & packed output
							for (size_t j = 0; j < two; ++j)
							{
								__m256 VV = _mm256_loadu_ps(reinterpret_cast<const float*>(pInputVector));
								pInputVector += sizeof(Float4) * 2;

								__m256 vTempX = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(0, 0, 0, 0));
								__m256 vTempY = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(1, 1, 1, 1));
								__m256 vTempZ = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(2, 2, 2, 2));
								__m256 vTempW = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(3, 3, 3, 3));

								vTempX = _mm256_mul_ps(vTempX, row0);
								vTempY = _mm256_mul_ps(vTempY, row1);
								vTempZ = _mm256_fmadd_ps(vTempZ, row2, vTempX);
								vTempW = _mm256_fmadd_ps(vTempW, row3, vTempY);
								vTempX = _mm256_add_ps(vTempZ, vTempW);

								_256_STREAM_PS(reinterpret_cast<float*>(pOutputVector), vTempX);
								pOutputVector += sizeof(Float4) * 2;

								i += 2;
							}
						}
						else
						{
							// Packed input, packed output
							for (size_t j = 0; j < two; ++j)
							{
								__m256 VV = _mm256_loadu_ps(reinterpret_cast<const float*>(pInputVector));
								pInputVector += sizeof(Float4) * 2;

								__m256 vTempX = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(0, 0, 0, 0));
								__m256 vTempY = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(1, 1, 1, 1));
								__m256 vTempZ = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(2, 2, 2, 2));
								__m256 vTempW = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(3, 3, 3, 3));

								vTempX = _mm256_mul_ps(vTempX, row0);
								vTempY = _mm256_mul_ps(vTempY, row1);
								vTempZ = _mm256_fmadd_ps(vTempZ, row2, vTempX);
								vTempW = _mm256_fmadd_ps(vTempW, row3, vTempY);
								vTempX = _mm256_add_ps(vTempZ, vTempW);

								_mm256_storeu_ps(reinterpret_cast<float*>(pOutputVector), vTempX);
								pOutputVector += sizeof(Float4) * 2;

								i += 2;
							}
						}
					}
					else
					{
						// Packed input, unpacked output
						for (size_t j = 0; j < two; ++j)
						{
							__m256 VV = _mm256_loadu_ps(reinterpret_cast<const float*>(pInputVector));
							pInputVector += sizeof(Float4) * 2;

							__m256 vTempX = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(0, 0, 0, 0));
							__m256 vTempY = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(1, 1, 1, 1));
							__m256 vTempZ = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(2, 2, 2, 2));
							__m256 vTempW = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(3, 3, 3, 3));

							vTempX = _mm256_mul_ps(vTempX, row0);
							vTempY = _mm256_mul_ps(vTempY, row1);
							vTempZ = _mm256_fmadd_ps(vTempZ, row2, vTempX);
							vTempW = _mm256_fmadd_ps(vTempW, row3, vTempY);
							vTempX = _mm256_add_ps(vTempZ, vTempW);

							_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), _mm256_castps256_ps128(vTempX));
							pOutputVector += outputStride;

							_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), _mm256_extractf128_ps(vTempX, 1));
							pOutputVector += outputStride;
							i += 2;
						}
					}
				}
			}

			if (i < vectorCount)
			{
				const VECTOR row0 = m.r[0];
				const VECTOR row1 = m.r[1];
				const VECTOR row2 = m.r[2];
				const VECTOR row3 = m.r[3];

				for (; i < vectorCount; i++)
				{
					__m128 V = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
					pInputVector += inputStride;

					VECTOR vTempX = PERMUTE_PS(V, _MM_SHUFFLE(0, 0, 0, 0));
					VECTOR vTempY = PERMUTE_PS(V, _MM_SHUFFLE(1, 1, 1, 1));
					VECTOR vTempZ = PERMUTE_PS(V, _MM_SHUFFLE(2, 2, 2, 2));
					VECTOR vTempW = PERMUTE_PS(V, _MM_SHUFFLE(3, 3, 3, 3));

					vTempX = _mm_mul_ps(vTempX, row0);
					vTempY = _mm_mul_ps(vTempY, row1);
					vTempZ = FMADD_PS(vTempZ, row2, vTempX);
					vTempW = FMADD_PS(vTempW, row3, vTempY);
					vTempX = _mm_add_ps(vTempZ, vTempW);

					_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), vTempX);
					pOutputVector += outputStride;
				}
			}

			SFENCE();

			return pOutputStream;

#elif defined(_SSE2_INTRINSICS_)
			const uint8_t* pInputVector = reinterpret_cast<const uint8_t*>(pInputStream);
			uint8_t* pOutputVector = reinterpret_cast<uint8_t*>(pOutputStream);

			const VECTOR row0 = m.r[0];
			const VECTOR row1 = m.r[1];
			const VECTOR row2 = m.r[2];
			const VECTOR row3 = m.r[3];

			if (!(reinterpret_cast<uintptr_t>(pOutputStream) & 0xF) && !(outputStride & 0xF))
			{
				if (!(reinterpret_cast<uintptr_t>(pInputStream) & 0xF) && !(inputStride & 0xF))
				{
					// Aligned input, aligned output
					for (size_t i = 0; i < vectorCount; i++)
					{
						__m128 v = _mm_load_ps(reinterpret_cast<const float*>(pInputVector));
						pInputVector += inputStride;

						VECTOR vTempX = PERMUTE_PS(v, _MM_SHUFFLE(0, 0, 0, 0));
						VECTOR vTempY = PERMUTE_PS(v, _MM_SHUFFLE(1, 1, 1, 1));
						VECTOR vTempZ = PERMUTE_PS(v, _MM_SHUFFLE(2, 2, 2, 2));
						VECTOR vTempW = PERMUTE_PS(v, _MM_SHUFFLE(3, 3, 3, 3));

						vTempX = _mm_mul_ps(vTempX, row0);
						vTempY = _mm_mul_ps(vTempY, row1);
						vTempZ = FMADD_PS(vTempZ, row2, vTempX);
						vTempW = FMADD_PS(vTempW, row3, vTempY);
						vTempX = _mm_add_ps(vTempZ, vTempW);

						STREAM_PS(reinterpret_cast<float*>(pOutputVector), vTempX);
						pOutputVector += outputStride;
					}
				}
				else
				{
					// Unaligned input, aligned output
					for (size_t i = 0; i < vectorCount; i++)
					{
						__m128 v = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
						pInputVector += inputStride;

						VECTOR vTempX = PERMUTE_PS(v, _MM_SHUFFLE(0, 0, 0, 0));
						VECTOR vTempY = PERMUTE_PS(v, _MM_SHUFFLE(1, 1, 1, 1));
						VECTOR vTempZ = PERMUTE_PS(v, _MM_SHUFFLE(2, 2, 2, 2));
						VECTOR vTempW = PERMUTE_PS(v, _MM_SHUFFLE(3, 3, 3, 3));

						vTempX = _mm_mul_ps(vTempX, row0);
						vTempY = _mm_mul_ps(vTempY, row1);
						vTempZ = FMADD_PS(vTempZ, row2, vTempX);
						vTempW = FMADD_PS(vTempW, row3, vTempY);
						vTempX = _mm_add_ps(vTempZ, vTempW);

						STREAM_PS(reinterpret_cast<float*>(pOutputVector), vTempX);
						pOutputVector += outputStride;
					}
				}
			}
			else
			{
				if (!(reinterpret_cast<uintptr_t>(pInputStream) & 0xF) && !(inputStride & 0xF))
				{
					// Aligned input, unaligned output
					for (size_t i = 0; i < vectorCount; i++)
					{
						__m128 v = _mm_load_ps(reinterpret_cast<const float*>(pInputVector));
						pInputVector += inputStride;

						VECTOR vTempX = PERMUTE_PS(v, _MM_SHUFFLE(0, 0, 0, 0));
						VECTOR vTempY = PERMUTE_PS(v, _MM_SHUFFLE(1, 1, 1, 1));
						VECTOR vTempZ = PERMUTE_PS(v, _MM_SHUFFLE(2, 2, 2, 2));
						VECTOR vTempW = PERMUTE_PS(v, _MM_SHUFFLE(3, 3, 3, 3));

						vTempX = _mm_mul_ps(vTempX, row0);
						vTempY = _mm_mul_ps(vTempY, row1);
						vTempZ = FMADD_PS(vTempZ, row2, vTempX);
						vTempW = FMADD_PS(vTempW, row3, vTempY);
						vTempX = _mm_add_ps(vTempZ, vTempW);

						_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), vTempX);
						pOutputVector += outputStride;
					}
				}
				else
				{
					// Unaligned input, unaligned output
					for (size_t i = 0; i < vectorCount; i++)
					{
						__m128 v = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
						pInputVector += inputStride;

						VECTOR vTempX = PERMUTE_PS(v, _MM_SHUFFLE(0, 0, 0, 0));
						VECTOR vTempY = PERMUTE_PS(v, _MM_SHUFFLE(1, 1, 1, 1));
						VECTOR vTempZ = PERMUTE_PS(v, _MM_SHUFFLE(2, 2, 2, 2));
						VECTOR vTempW = PERMUTE_PS(v, _MM_SHUFFLE(3, 3, 3, 3));

						vTempX = _mm_mul_ps(vTempX, row0);
						vTempY = _mm_mul_ps(vTempY, row1);
						vTempZ = FMADD_PS(vTempZ, row2, vTempX);
						vTempW = FMADD_PS(vTempW, row3, vTempY);
						vTempX = _mm_add_ps(vTempZ, vTempW);

						_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), vTempX);
						pOutputVector += outputStride;
					}
				}
			}

			SFENCE();

			return pOutputStream;
#endif
		}
	}
}

#endif // !ULTREALITY_MATH_FLOAT4_INL
