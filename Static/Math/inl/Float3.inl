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

	namespace VEC3
	{
		FORCE_INLINE bool VEC_CALLCONV Equal(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return (((V1.vector4_f32[0] == V2.vector4_f32[0]) && (V1.vector4_f32[1] == V2.vector4_f32[1]) && (V1.vector4_f32[2] == V2.vector4_f32[2])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmpeq_ps(V1, V2);

			return (((_mm_movemask_ps(vTemp) & 7) == 7) != 0);
#endif
		}

		FORCE_INLINE uint32_t VEC_CALLCONV EqualR(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			uint32_t CR = 0;
			if ((V1.vector4_f32[0] == V2.vector4_f32[0]) &&
				(V1.vector4_f32[1] == V2.vector4_f32[1]) &&
				(V1.vector4_f32[2] == V2.vector4_f32[2]))
			{
				CR = CRMASK_CR6TRUE;
			}
			else if ((V1.vector4_f32[0] != V2.vector4_f32[0]) &&
				(V1.vector4_f32[1] != V2.vector4_f32[1]) &&
				(V1.vector4_f32[2] != V2.vector4_f32[2]))
			{
				CR = CRMASK_CR6FALSE;
			}

			return CR;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmpeq_ps(V1, V2);
			int iTest = _mm_movemask_ps(vTemp) & 7;
			uint32_t CR = 0;
			if (iTest == 7)
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
			float dx, dy, dz;

			dx = fabsf(V1.vector4_f32[0] - V2.vector4_f32[0]);
			dy = fabsf(V1.vector4_f32[1] - V2.vector4_f32[1]);
			dz = fabsf(V1.vector4_f32[2] - V2.vector4_f32[2]);
			
			return (((dx <= epsilon.vector4_f32[0]) &&
					 (dy <= epsilon.vector4_f32[1]) &&
					 (dz <= epsilon.vector4_f32[2])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			// Get the difference
			VECTOR vDelta = _mm_sub_ps(V1, V2);
			// Get the absolute value of the difference
			VECTOR vTemp = _mm_setzero_ps();
			vTemp = _mm_sub_ps(vTemp, vDelta);
			vTemp = _mm_max_ps(vTemp, vDelta);
			vTemp = _mm_cmple_ps(vTemp, epsilon);
			
			// w is don't care
			return (((_mm_movemask_ps(vTemp) & 7) == 0x7) != 0);
#endif
		}

		FORCE_INLINE bool VEC_CALLCONV NotEqual(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return (((V1.vector4_f32[0] != V2.vector4_f32[0]) || (V1.vector4_f32[1] != V2.vector4_f32[1]) || (V1.vector4_f32[2] != V2.vector4_f32[2])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmpeq_ps(V1, V2);
			
			return (((_mm_movemask_ps(vTemp) & 7) != 7) != 0);
#endif
		}

		FORCE_INLINE bool VEC_CALLCONV Greater(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return (((V1.vector4_f32[0] > V2.vector4_f32[0]) && (V1.vector4_f32[1] > V2.vector4_f32[1]) && (V1.vector4_f32[2] > V2.vector4_f32[2])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmpgt_ps(V1, V2);
    		
			return (((_mm_movemask_ps(vTemp) & 7) == 7) != 0);
#endif
		}

		FORCE_INLINE uint32_t VEC_CALLCONV GreaterR(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			uint32_t CR = 0;
			if ((V1.vector4_f32[0] > V2.vector4_f32[0]) &&
				(V1.vector4_f32[1] > V2.vector4_f32[1]) &&
				(V1.vector4_f32[2] > V2.vector4_f32[2]))
			{
				CR = CRMASK_CR6TRUE;
			}
			else if ((V1.vector4_f32[0] <= V2.vector4_f32[0]) &&
				(V1.vector4_f32[1] <= V2.vector4_f32[1]) &&
				(V1.vector4_f32[2] <= V2.vector4_f32[2]))
			{
				CR = CRMASK_CR6FALSE;
			}
			
			return CR;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmpgt_ps(V1, V2);
			uint32_t CR = 0;
			int iTest = _mm_movemask_ps(vTemp) & 7;
			if (iTest == 7)
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
			return (((V1.vector4_f32[0] >= V2.vector4_f32[0]) && (V1.vector4_f32[1] >= V2.vector4_f32[1]) && (V1.vector4_f32[2] >= V2.vector4_f32[2])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmpge_ps(V1, V2);
    		
			return (((_mm_movemask_ps(vTemp) & 7) == 7) != 0);
#endif
		}

		FORCE_INLINE uint32_t VEC_CALLCONV GreaterOrEqualR(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			uint32_t CR = 0;
			if ((V1.vector4_f32[0] >= V2.vector4_f32[0]) &&
				(V1.vector4_f32[1] >= V2.vector4_f32[1]) &&
				(V1.vector4_f32[2] >= V2.vector4_f32[2]))
			{
				CR = CRMASK_CR6TRUE;
			}
			else if ((V1.vector4_f32[0] < V2.vector4_f32[0]) &&
				(V1.vector4_f32[1] < V2.vector4_f32[1]) &&
				(V1.vector4_f32[2] < V2.vector4_f32[2]))
			{
				CR = CRMASK_CR6FALSE;
			}

			return CR;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmpge_ps(V1, V2);
			uint32_t CR = 0;
			int iTest = _mm_movemask_ps(vTemp) & 7;
			if (iTest == 7)
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
			return (((V1.vector4_f32[0] < V2.vector4_f32[0]) && (V1.vector4_f32[1] < V2.vector4_f32[1]) && (V1.vector4_f32[2] < V2.vector4_f32[2])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmplt_ps(V1, V2);

			return (((_mm_movemask_ps(vTemp) & 7) == 7) != 0);
#endif
		}

		FORCE_INLINE bool VEC_CALLCONV LessOrEqual(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return (((V1.vector4_f32[0] <= V2.vector4_f32[0]) && (V1.vector4_f32[1] <= V2.vector4_f32[1]) && (V1.vector4_f32[2] <= V2.vector4_f32[2])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmple_ps(V1, V2);

			return (((_mm_movemask_ps(vTemp) & 7) == 7) != 0);
#endif
		}

		FORCE_INLINE bool VEC_CALLCONV InBounds(A_VECTOR v, A_VECTOR bounds) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return (((v.vector4_f32[0] <= bounds.vector4_f32[0] && v.vector4_f32[0] >= -bounds.vector4_f32[0]) &&
					 (v.vector4_f32[1] <= bounds.vector4_f32[1] && v.vector4_f32[1] >= -bounds.vector4_f32[1]) &&
					 (v.vector4_f32[2] <= bounds.vector4_f32[2] && v.vector4_f32[2] >= -bounds.vector4_f32[2])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			// Test if less than or equal
			VECTOR vTemp1 = _mm_cmple_ps(v, bounds);
			// Negate the bounds
			VECTOR vTemp2 = _mm_mul_ps(bounds, g_NegativeOne);
			// Test if greater or equal (Reversed)
			vTemp2 = _mm_cmple_ps(vTemp2, v);
			// Blend answers
			vTemp1 = _mm_and_ps(vTemp1, vTemp2);
			
			// x,y and z in bounds? (w is don't care)
			return (((_mm_movemask_ps(vTemp1) & 0x7) == 0x7) != 0);
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
        			ISNAN(v.vector4_f32[2]));

#elif defined(_SSE2_INTRINSICS_)
#if defined(__clang__) && defined(__FINITE_MATH_ONLY__)
    		ALIGNED(16) float tmp[4];
    		_mm_store_ps(tmp, v);
    		
			return isnan(tmp[0]) || isnan(tmp[1]) || isnan(tmp[2]);
#else
    		// Test against itself. NaN is always not equal
   		 	VECTOR vTempNan = _mm_cmpneq_ps(v, v);
    		
			// If x or y or z are NaN, the mask is non-zero
    		return ((_mm_movemask_ps(vTempNan) & 7) != 0);
#endif
#endif
		}

#if !defined(_NO_INTRINSICS_) && defined(_MSC_VER) && !defined(__INTEL_COMPILER)
#pragma float_control(pop)
#endif

		FORCE_INLINE bool VEC_CALLCONV IsInfinite(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return (ISINF(v.vector4_f32[0]) ||
        			ISINF(v.vector4_f32[1]) ||
        			ISINF(v.vector4_f32[2]));

#elif defined(_SSE2_INTRINSICS_)
			// Mask off the sign bit
			__m128 vTemp = _mm_and_ps(v, g_AbsMask);
			// Compare to infinity
			vTemp = _mm_cmpeq_ps(vTemp, g_Infinity);
			
			// If x,y or z are infinity, the signs are true.
			return ((_mm_movemask_ps(vTemp) & 7) != 0);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Dot(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			float fValue = V1.vector4_f32[0] * V2.vector4_f32[0] + V1.vector4_f32[1] * V2.vector4_f32[1] + V1.vector4_f32[2] * V2.vector4_f32[2];
			VECTOR_F32 vResult;
			vResult.f[0] =
				vResult.f[1] =
				vResult.f[2] =
				vResult.f[3] = fValue;
			
			return vResult.v;

#elif defined(_SSE4_INTRINSICS_)
    		return _mm_dp_ps(V1, V2, 0x7f);
			
#elif defined(_SSE3_INTRINSICS_)
			VECTOR vTemp = _mm_mul_ps(V1, V2);
			vTemp = _mm_and_ps(vTemp, g_Mask3);
			vTemp = _mm_hadd_ps(vTemp, vTemp);
			
			return _mm_hadd_ps(vTemp, vTemp);

#elif defined(_SSE2_INTRINSICS_)
			// Perform the dot product
			VECTOR vDot = _mm_mul_ps(V1, V2);
			// x=Dot.vector4_f32[1], y=Dot.vector4_f32[2]
			VECTOR vTemp = PERMUTE_PS(vDot, _MM_SHUFFLE(2, 1, 2, 1));
			// Result.vector4_f32[0] = x+y
			vDot = _mm_add_ss(vDot, vTemp);
			// x=Dot.vector4_f32[2]
			vTemp = PERMUTE_PS(vTemp, _MM_SHUFFLE(1, 1, 1, 1));
			// Result.vector4_f32[0] = (x+y)+z
			vDot = _mm_add_ss(vDot, vTemp);
			
			// Splat x
			return PERMUTE_PS(vDot, _MM_SHUFFLE(0, 0, 0, 0));
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Cross(A_VECTOR V1, A_VECTOR V2) noexcept
		{
			// [V1.y*V2.z - V1.z*V2.y, V1.z*V2.x - V1.x*V2.z, V1.x*V2.y - V1.y*V2.x]
#if defined(_NO_INTRINSICS_)
			VECTORF_32 vResult = { { {
				(V1.vector4_f32[1] * V2.vector4_f32[2]) - (V1.vector4_f32[2] * V2.vector4_f32[1]),
				(V1.vector4_f32[2] * V2.vector4_f32[0]) - (V1.vector4_f32[0] * V2.vector4_f32[2]),
				(V1.vector4_f32[0] * V2.vector4_f32[1]) - (V1.vector4_f32[1] * V2.vector4_f32[0]),
				0.0f
			} } };

			return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
			// y1,z1,x1,w1
			VECTOR vTemp1 = PERMUTE_PS(V1, _MM_SHUFFLE(3, 0, 2, 1));
			// z2,x2,y2,w2
			VECTOR vTemp2 = PERMUTE_PS(V2, _MM_SHUFFLE(3, 1, 0, 2));
			// Perform the left operation
			VECTOR vResult = _mm_mul_ps(vTemp1, vTemp2);
			// z1,x1,y1,w1
			vTemp1 = PERMUTE_PS(vTemp1, _MM_SHUFFLE(3, 0, 2, 1));
			// y2,z2,x2,w2
			vTemp2 = PERMUTE_PS(vTemp2, _MM_SHUFFLE(3, 1, 0, 2));
			// Perform the right operation
			vResult = FNMADD_PS(vTemp1, vTemp2, vResult);
			// Set w to zero
			return _mm_and_ps(vResult, g_Mask3);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV LengthSq(A_VECTOR v) noexcept
		{
			return Dot(v, v);
		}

		FORCE_INLINE VECTOR VEC_CALLCONV ReciprocalLengthEst(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR Result;

			Result = LengthSq(v);
			Result = ReciprocalSqrtEst(Result);

			return Result;

#elif defined(_SSE4_INTRINSICS_)
			VECTOR vTemp = _mm_dp_ps(v, v, 0x7f);
			
			return _mm_rsqrt_ps(vTemp);

#elif defined(_SSE3_INTRINSICS_)
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			vLengthSq = _mm_and_ps(vLengthSq, g_Mask3);
			vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
			vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
			vLengthSq = _mm_rsqrt_ps(vLengthSq);
			
			return vLengthSq;

#elif defined(_SSE2_INTRINSICS_)
			// Perform the dot product on x,y and z
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			// vTemp has z and y
			VECTOR vTemp = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(1, 2, 1, 2));
			// x+z, y
			vLengthSq = _mm_add_ss(vLengthSq, vTemp);
			// y,y,y,y
			vTemp = PERMUTE_PS(vTemp, _MM_SHUFFLE(1, 1, 1, 1));
			// x+z+y,??,??,??
			vLengthSq = _mm_add_ss(vLengthSq, vTemp);
			// Splat the length squared
			vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(0, 0, 0, 0));
			// Get the reciprocal
			vLengthSq = _mm_rsqrt_ps(vLengthSq);
			
			return vLengthSq;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV ReciprocalLength(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR Result = LengthSq(v);
			
			return VEC::ReciprocalSqrt(Result);

#elif defined(_SSE4_INTRINSICS_)
			VECTOR vTemp = _mm_dp_ps(v, v, 0x7f);
			VECTOR vLengthSq = _mm_sqrt_ps(vTemp);
			
			return _mm_div_ps(g_One, vLengthSq);

#elif defined(_SSE3_INTRINSICS_)
			VECTOR vDot = _mm_mul_ps(v, v);
			vDot = _mm_and_ps(vDot, g_Mask3);
			vDot = _mm_hadd_ps(vDot, vDot);
			vDot = _mm_hadd_ps(vDot, vDot);
			vDot = _mm_sqrt_ps(vDot);
			vDot = _mm_div_ps(g_One, vDot);
			
			return vDot;

#elif defined(_SSE2_INTRINSICS_)
			// Perform the dot product
			VECTOR vDot = _mm_mul_ps(v, v);
			// x=Dot.y, y=Dot.z
			VECTOR vTemp = PERMUTE_PS(vDot, _MM_SHUFFLE(2, 1, 2, 1));
			// Result.x = x+y
			vDot = _mm_add_ss(vDot, vTemp);
			// x=Dot.z
			vTemp = PERMUTE_PS(vTemp, _MM_SHUFFLE(1, 1, 1, 1));
			// Result.x = (x+y)+z
			vDot = _mm_add_ss(vDot, vTemp);
			// Splat x
			vDot = PERMUTE_PS(vDot, _MM_SHUFFLE(0, 0, 0, 0));
			// Get the reciprocal
			vDot = _mm_sqrt_ps(vDot);
			// Get the reciprocal
			vDot = _mm_div_ps(g_One, vDot);
			
			return vDot;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV LengthEst(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR Result = LengthSq(v);

			return VEC::SqrtEst(Result);

#elif defined(_SSE4_INTRINSICS_)
			VECTOR vTemp = _mm_dp_ps(v, v, 0x7f);
			
			return _mm_sqrt_ps(vTemp);

#elif defined(_SSE3_INTRINSICS_)
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			vLengthSq = _mm_and_ps(vLengthSq, g_Mask3);
			vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
			vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
			vLengthSq = _mm_sqrt_ps(vLengthSq);
			
			return vLengthSq;

#elif defined(_SSE2_INTRINSICS_)
			// Perform the dot product on x,y and z
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			// vTemp has z and y
			VECTOR vTemp = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(1, 2, 1, 2));
			// x+z, y
			vLengthSq = _mm_add_ss(vLengthSq, vTemp);
			// y,y,y,y
			vTemp = PERMUTE_PS(vTemp, _MM_SHUFFLE(1, 1, 1, 1));
			// x+z+y,??,??,??
			vLengthSq = _mm_add_ss(vLengthSq, vTemp);
			// Splat the length squared
			vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(0, 0, 0, 0));
			// Get the length
			vLengthSq = _mm_sqrt_ps(vLengthSq);
			
			return vLengthSq;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Length(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR Result = LengthSq(v);

			return VEC::Sqrt(Result);

#elif defined(_SSE4_INTRINSICS_)
			VECTOR vTemp = _mm_dp_ps(v, v, 0x7f);
			
			return _mm_sqrt_ps(vTemp);

#elif defined(_SSE3_INTRINSICS_)
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			vLengthSq = _mm_and_ps(vLengthSq, g_Mask3);
			vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
			vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
			vLengthSq = _mm_sqrt_ps(vLengthSq);
			
			return vLengthSq;

#elif defined(_SSE2_INTRINSICS_)
			// Perform the dot product on x,y and z
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			// vTemp has z and y
			VECTOR vTemp = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(1, 2, 1, 2));
			// x+z, y
			vLengthSq = _mm_add_ss(vLengthSq, vTemp);
			// y,y,y,y
			vTemp = PERMUTE_PS(vTemp, _MM_SHUFFLE(1, 1, 1, 1));
			// x+z+y,??,??,??
			vLengthSq = _mm_add_ss(vLengthSq, vTemp);
			// Splat the length squared
			vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(0, 0, 0, 0));
			// Get the length
			vLengthSq = _mm_sqrt_ps(vLengthSq);
			
			return vLengthSq;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV NormalizeEst(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR Result = ReciprocalLength(v);

			return VEC::Multiply(v, Result);

#elif defined(_SSE4_INTRINSICS_)
			VECTOR vTemp = _mm_dp_ps(v, v, 0x7f);
			VECTOR vResult = _mm_rsqrt_ps(vTemp);
			
			return _mm_mul_ps(vResult, v);

#elif defined(_SSE3_INTRINSICS_)
			VECTOR vDot = _mm_mul_ps(v, v);
			vDot = _mm_and_ps(vDot, g_Mask3);
			vDot = _mm_hadd_ps(vDot, vDot);
			vDot = _mm_hadd_ps(vDot, vDot);
			vDot = _mm_rsqrt_ps(vDot);
			vDot = _mm_mul_ps(vDot, v);
			
			return vDot;

#elif defined(_XM_SSE_INTRINSICS_)
			// Perform the dot product
			VECTOR vDot = _mm_mul_ps(v, v);
			// x=Dot.y, y=Dot.z
			VECTOR vTemp = PERMUTE_PS(vDot, _MM_SHUFFLE(2, 1, 2, 1));
			// Result.x = x+y
			vDot = _mm_add_ss(vDot, vTemp);
			// x=Dot.z
			vTemp = PERMUTE_PS(vTemp, _MM_SHUFFLE(1, 1, 1, 1));
			// Result.x = (x+y)+z
			vDot = _mm_add_ss(vDot, vTemp);
			// Splat x
			vDot = PERMUTE_PS(vDot, _MM_SHUFFLE(0, 0, 0, 0));
			// Get the reciprocal
			vDot = _mm_rsqrt_ps(vDot);
			// Perform the normalization
			vDot = _mm_mul_ps(vDot, v);
			
			return vDot;
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
			VECTOR vLengthSq = _mm_dp_ps(v, v, 0x7f);
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
			vResult = _mm_or_ps(vTemp1, vTemp2);
			
			return vResult;

#elif defined(_SSE3_INTRINSICS_)
			// Perform the dot product on x,y and z only
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			vLengthSq = _mm_and_ps(vLengthSq, g_Mask3);
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
			vResult = _mm_div_ps(V, vResult);
			// Any that are infinity, set to zero
			vResult = _mm_and_ps(vResult, vZeroMask);
			// Select qnan or result based on infinite length
			VECTOR vTemp1 = _mm_andnot_ps(vLengthSq, g_QNaN);
			VECTOR vTemp2 = _mm_and_ps(vResult, vLengthSq);
			vResult = _mm_or_ps(vTemp1, vTemp2);
			
			return vResult;

#elif defined(_SSE2_INTRINSICS_)
			// Perform the dot product on x,y and z only
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			VECTOR vTemp = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(2, 1, 2, 1));
			vLengthSq = _mm_add_ss(vLengthSq, vTemp);
			vTemp = PERMUTE_PS(vTemp, _MM_SHUFFLE(1, 1, 1, 1));
			vLengthSq = _mm_add_ss(vLengthSq, vTemp);
			vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(0, 0, 0, 0));
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
			vResult = _mm_or_ps(vTemp1, vTemp2);
			
			return vResult;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV ClampLength(A_VECTOR v, float lengthMin, float lengthMax)
		{
			VECTOR clampMax = VEC::Replicate(lengthMax);
			VECTOR clampMin = VEC::Replicate(lengthMin);

			return ClampLengthV(v, clampMin, clampMax);
		}

		FORCE_INLINE VECTOR VEC_CALLCONV ClampLengthV(A_VECTOR v, A_VECTOR lengthMin, A_VECTOR lengthMax) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert((VEC::GetY(lengthMin) == VEC::GetX(lengthMin)) && (VEC::GetZ(lengthMin) == VEC::GetX(lengthMin)));
			assert((VEC::GetY(lengthMax) == VEC::GetX(lengthMax)) && (VEC::GetZ(lengthMax) == VEC::GetX(lengthMax)));
			assert(GreaterOrEqual(lengthMin, VEC::Zero()));
			assert(GreaterOrEqual(lengthMax, VEC::Zero()));
			assert(GreaterOrEqual(lengthMax, lengthMin));
#endif

			VECTOR lengthSq = LengthSq(v);

			const VECTOR zero = VEC::Zero();

			VECTOR rcpLength = VEC::ReciprocalSqrt(lengthSq);

			VECTOR infiniteLength = VEC::EqualInt(lengthSq, g_Infinity.v);
			VECTOR zeroLength = VEC::Equal(lengthSq, zero);

			VECTOR normal = VEC::Multiply(v, rcpLength);

			VECTOR length = VEC::Multiply(lengthSq, rcpLength);

			VECTOR select = VEC::EqualInt(infiniteLength, zeroLength);
			length = VEC::Select(lengthSq, length, select);
			normal = VEC::Select(lengthSq, normal, select);

			VECTOR controlMax = VEC::Greater(length, lengthMax);
			VECTOR controlMin = VEC::Less(length, lengthMin);

			VECTOR clampLength = VEC::Select(length, lengthMax, controlMax);
			clampLength = VEC::Select(clampLength, lengthMin, controlMin);

			VECTOR Result = VEC::Multiply(normal, clampLength);

			// Preserve the original vector (with no precision loss) if the length falls within the given range
			VECTOR control = VEC::EqualInt(controlMax, controlMin);
			Result = VEC::Select(Result, v, control);

			return Result;
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Reflect(A_VECTOR incident, A_VECTOR normal) noexcept
		{
			// Result = incident - (2 * dot(incident, normal)) * normal
			VECTOR Result = Dot(incident, normal);
			Result = VEC::Add(Result, Result);
			Result = VEC::NegativeMultiplySubtract(Result, normal, incident);

			return Result;
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Refract(A_VECTOR incident, A_VECTOR normal, float refractionIndex) noexcept
		{
			VECTOR index = VEC::Replicate(refractionIndex);

			return RefractV(incident, normal, index);
		}

		FORCE_INLINE VECTOR VEC_CALLCONV RefractV(A_VECTOR incident, A_VECTOR normal, A_VECTOR refractionIndex) noexcept
		{
			// Result = refractionIndex * incident - normal * (refractionIndex * dot(incident, normal) + 
			// sqrt(1 - refractionIndex * refractionIndex * (1 - dot(incident, normal) * dot(incident, normal))))
#if defined(_NO_INTRINSICS_)
			const VECTOR zero = VEC::Zero();

			VECTOR IDotN = Dot(incident, normal);

			// R = 1.0f - refractionIndex * refractionIndex * (1.0f - IDotN * IDotN)
			VECTOR R = VEC::NegativeMultiplySubtract(IDotN, IDotN, g_One.v);
			R = VEC::Multiply(R, refractionIndex);
			R = VEC::NegativeMultiplySubtract(R, refractionIndex, g_One.v);

			if (VEC4::LessOrEqual(R, zero))
			{
				// Total internal reflection
				return zero;
			}
			else
			{
				// R = refractionIndex * IDotN + sqrt(R)
				R = VEC::Sqrt(R);
				R = VEC::MultiplyAdd(refractionIndex, IDotN, R);

				// Result = refractionIndex * incident - normal * R
				VECTOR Result = VEC::Multiply(refractionIndex, incident);
				Result = VEC::NegativeMultiplySubtract(normal, R, Result);

				return Result;
			}

#elif defined(_SSE2_INTRINSICS_)
			// Result = refractionIndex * incident - normal * (refractionIndex * dot(incident, normal) +
			// sqrt(1 - refractionIndex * refractionIndex * (1 - dot(incident, normal) * dot(incident, normal))))
			VECTOR IDotN = Dot(incident, normal);
			// R = 1.0f - refractionIndex * refractionIndex * (1.0f - IDotN * IDotN)
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
				// R = refractionIndex * IDotN + sqrt(R)
				R = _mm_sqrt_ps(R);
				R = FMADD_PS(refractionIndex, IDotN, R);
				// Result = refractionIndex * incident - normal * R
				vResult = _mm_mul_ps(refractionIndex, incident);
				vResult = FNMADD_PS(R, normal, vResult);
			}

			return vResult;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Orthogonal(A_VECTOR v) noexcept
		{
			VECTOR zero = VEC::Zero();
			VECTOR z = VEC::SplatZ(v);
			VECTOR yzyy = VEC::Swizzle<SWIZZLE_Y, SWIZZLE_Z, SWIZZLE_Y, SWIZZLE_Y>(v);

			VECTOR negativeV = VEC::Subtract(zero, v);

			VECTOR zIsNegative = VEC::Less(z, zero);
			VECTOR yzyyIsNegative = VEC::Less(yzyy, zero);

			VECTOR s = VEC::Add(yzyy, z);
			VECTOR d = VEC::Subtract(yzyy, z);

			VECTOR select = VEC::EqualInt(zIsNegative, yzyyIsNegative);

			VECTOR R0 = VEC::Permute<PERMUTE_1X, PERMUTE_0X, PERMUTE_0X, PERMUTE_0X>(negativeV, s);
			VECTOR R1 = VEC::Permute<PERMUTE_1X, PERMUTE_0X, PERMUTE_0X, PERMUTE_0X>(v, d);

			return VEC::Select(R1, R0, select);
		}

		FORCE_INLINE VECTOR VEC_CALLCONV AngleBetweenNormalsEst(A_VECTOR N1, A_VECTOR N2) noexcept
		{
			VECTOR Result = Dot(N1, N2);
			Result = VEC::Clamp(Result, g_NegativeOne.v, g_One.v);
			Result = VEC::ACosEst(Result);

			return Result;
		}

		FORCE_INLINE VECTOR VEC_CALLCONV AngleBetweenNormals(A_VECTOR N1, A_VECTOR N2) noexcept
		{
			VECTOR Result = Dot(N1, N2);
			Result = VEC::Clamp(Result, g_NegativeOne.v, g_One.v);
			Result = VEC::ACos(Result);

			return Result;
		}

		FORCE_INLINE VECTOR VEC_CALLCONV AngleBetweenVectors(A_VECTOR V1, A_VECTOR V2) noexcept
		{
			VECTOR L1 = ReciprocalLength(V1);
			VECTOR L2 = ReciprocalLength(V2);

			VECTOR dot = Dot(V1, V2);

			L1 = VEC::Multiply(L1, L2);

			VECTOR cosAngle = VEC::Multiply(dot, L1);
			cosAngle = VEC::Clamp(cosAngle, g_NegativeOne.v, g_One.v);

			return VEC::ACos(cosAngle);
		}

		FORCE_INLINE VECTOR VEC_CALLCONV LinePointDistance(A_VECTOR linePoint1, A_VECTOR linePoint2, A_VECTOR point) noexcept
		{
			// Given a vector, pointVector, from linePoint1 to point,
			// and a vector, lineVector, from linePoint1 to linePoint2, 
			// the scaled distance pointProjectionScale from linePoint1 
			// to the perpendicular projection of pointVector onto the
			// line is defined as:
			//
			//	pointProjectionScale = dot(pointVector, lineVector) / LengthSq(lineVector)

			VECTOR pointVector = VEC::Subtract(point, linePoint1);
			VECTOR lineVector = VEC::Subtract(linePoint2, linePoint1);

			VECTOR lengthSq = LengthSq(lineVector);

			VECTOR pointProjectionScale = Dot(pointVector, lineVector);
			pointProjectionScale = VEC::Divide(pointProjectionScale, lengthSq);

			VECTOR distanceVector = VEC::Multiply(lineVector, pointProjectionScale);
			distanceVector = VEC::Subtract(pointVector, distanceVector);

			return Length(distanceVector);
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV ComponentsFromNormal(VECTOR* pParallel, VECTOR* pPerpendicular, A_VECTOR v, A_VECTOR normal) noexcept
		{
#if defined(DEBUG) || (_DEBUG)
			assert(pParallel !- nullptr);
			assert(pPerpendicular != nullptr);
#endif

			VECTOR scale = Dot(v, normal);

			VECTOR parallel = VEC::Multiply(normal, scale);

			*pParallel = parallel;
			*pPerpendicular = VEC::Subtract(v, parallel);
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Rotate(A_VECTOR v, A_VECTOR rotationQuaternion) noexcept
		{
			VECTOR A = VEC::Select(g_Select1110.v, v, g_Select1110.v);
			VECTOR Q = QUAT::Conjugate(rotationQuaternion);
			VECTOR Result = QUAT::Multiply(Q, A);

			return QUAT::Multiply(Result, rotationQuaternion);
		}

		FORCE_INLINE VECTOR VEC_CALLCONV InverseRotate(A_VECTOR v, A_VECTOR rotationQuaternion) noexcept
		{
			VECTOR A = VEC::Select(g_Select1110.v, v, g_Select1110.v);
			VECTOR Result = QUAT::Multiply(rotationQuaternion, A);
			VECTOR Q = QUAT::Conjugate(rotationQuaternion);

			return QUAT::Multiply(Result, Q);
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Transform(A_VECTOR v, A_MATRIX m) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR Z = VEC::SplatZ(v);
			VECTOR Y = VEC::SplatY(v);
			VECTOR X = VEC::SplatX(v);

			VECTOR Result = VEC::MultiplyAdd(Z, m.r[2], m.r[3]);
			Result = VEC::MultiplyAdd(Y, m.r[1], Result);
			Result = VEC::MultiplyAdd(X, m.r[0], Result);

			return Result;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vResult = PERMUTE_PS(v, _MM_SHUFFLE(2, 2, 2, 2));
			vResult = FMADD_PS(vResult, m.r[2], m.r[3]);
			VECTOR vTemp = PERMUTE_PS(v, _MM_SHUFFLE(1, 1, 1, 1)); // Y
			vResult = FMADD_PS(vTemp, m.r[1], vResult);
			vTemp = PERMUTE_PS(v, _MM_SHUFFLE(0, 0, 0, 0)); // X
			vResult = FMADD_PS(vTemp, m.r[0], vResult);
			
			return vResult;
#endif
		}

#ifdef _PREFAST_
#pragma prefast(push)
#pragma prefast(disable : 26015 26019, "PREfast noise: Esp:1307" )
#endif

		_Use_decl_annotations_
		FORCE_INLINE Float4* VEC_CALLCONV TransformStream(Float4* pOutputStream, size_t outputStride, const Float3* pInputStream, size_t inputStride, size_t vectorCount, A_MATRIX m) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pOutputStream != nullptr);
			assert(pInputStream != nullptr);

			assert(inputStride >= sizeof(Float3));
			_Analysis_assume_(inputStride >= sizeof(Float3));

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
				VECTOR V = VEC::LoadFloat3(reinterpret_cast<const Float3*>(pInputVector));
				VECTOR Z = VEC::SplatZ(V);
				VECTOR Y = VEC::SplatY(V);
				VECTOR X = VEC::SplatX(V);

				VECTOR Result = VEC::MultiplyAdd(Z, row2, row3);
				Result = VEC::MultiplyAdd(Y, row1, Result);
				Result = VEC::MultiplyAdd(X, row0, Result);

				VEC::StoreFloat4(reinterpret_cast<Float4*>(pOutputVector), Result);

				pInputVector += inputStride;
				pOutputVector += outputStride;
			}

			return pOutputStream;

#elif defined(_SSE2_INTRINSICS_)
			const uint8_t* pInputVector = reinterpret_cast<const uint8_t*>(pInputStream);
			uint8_t* pOutputVector = reinterpret_cast<uint8_t*>(pOutputStream);

			const VECTOR row0 = m.r[0];
			const VECTOR row1 = m.r[1];
			const VECTOR row2 = m.r[2];
			const VECTOR row3 = m.r[3];

			size_t i = 0;
			size_t four = vectorCount >> 2;
			if (four > 0)
			{
				if (inputStride == sizeof(Float3))
				{
					if (!(reinterpret_cast<uintptr_t>(pOutputStream) & 0xF) && !(outputStride & 0xF))
					{
						// Packed input, aligned output
						for (size_t j = 0; j < four; ++j)
						{
							__m128 V1 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
							__m128 L2 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 16));
							__m128 L3 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 32));
							pInputVector += sizeof(Float3) * 4;

							// Unpack the 4 vectors (.w components are junk)
							3UNPACK3INTO4(V1, L2, L3);

							// Result 1
							VECTOR Z = PERMUTE_PS(V1, _MM_SHUFFLE(2, 2, 2, 2));
							VECTOR Y = PERMUTE_PS(V1, _MM_SHUFFLE(1, 1, 1, 1));
							VECTOR X = PERMUTE_PS(V1, _MM_SHUFFLE(0, 0, 0, 0));

							VECTOR vTemp = FMADD_PS(Z, row2, row3);
							VECTOR vTemp2 = _mm_mul_ps(Y, row1);
							VECTOR vTemp3 = _mm_mul_ps(X, row0);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);
							STREAM_PS(reinterpret_cast<float*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							// Result 2
							Z = PERMUTE_PS(V2, _MM_SHUFFLE(2, 2, 2, 2));
							Y = PERMUTE_PS(V2, _MM_SHUFFLE(1, 1, 1, 1));
							X = PERMUTE_PS(V2, _MM_SHUFFLE(0, 0, 0, 0));

							vTemp = FMADD_PS(Z, row2, row3);
							vTemp2 = _mm_mul_ps(Y, row1);
							vTemp3 = _mm_mul_ps(X, row0);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);
							STREAM_PS(reinterpret_cast<float*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							// Result 3
							Z = PERMUTE_PS(V3, _MM_SHUFFLE(2, 2, 2, 2));
							Y = PERMUTE_PS(V3, _MM_SHUFFLE(1, 1, 1, 1));
							X = PERMUTE_PS(V3, _MM_SHUFFLE(0, 0, 0, 0));

							vTemp = FMADD_PS(Z, row2, row3);
							vTemp2 = _mm_mul_ps(Y, row1);
							vTemp3 = _mm_mul_ps(X, row0);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);
							STREAM_PS(reinterpret_cast<float*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							// Result 4
							Z = PERMUTE_PS(V4, _MM_SHUFFLE(2, 2, 2, 2));
							Y = PERMUTE_PS(V4, _MM_SHUFFLE(1, 1, 1, 1));
							X = PERMUTE_PS(V4, _MM_SHUFFLE(0, 0, 0, 0));

							vTemp = FMADD_PS(Z, row2, row3);
							vTemp2 = _mm_mul_ps(Y, row1);
							vTemp3 = _mm_mul_ps(X, row0);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);
							XM_STREAM_PS(reinterpret_cast<float*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							i += 4;
						}
					}
					else
					{
						// Packed input, unaligned output
						for (size_t j = 0; j < four; ++j)
						{
							__m128 V1 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
							__m128 L2 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 16));
							__m128 L3 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 32));
							pInputVector += sizeof(Float3) * 4;

							// Unpack the 4 vectors (.w components are junk)
							3UNPACK3INTO4(V1, L2, L3);

							// Result 1
							VECTOR Z = PERMUTE_PS(V1, _MM_SHUFFLE(2, 2, 2, 2));
							VECTOR Y = PERMUTE_PS(V1, _MM_SHUFFLE(1, 1, 1, 1));
							VECTOR X = PERMUTE_PS(V1, _MM_SHUFFLE(0, 0, 0, 0));

							VECTOR vTemp = FMADD_PS(Z, row2, row3);
							VECTOR vTemp2 = _mm_mul_ps(Y, row1);
							VECTOR vTemp3 = _mm_mul_ps(X, row0);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);
							_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							// Result 2
							Z = PERMUTE_PS(V2, _MM_SHUFFLE(2, 2, 2, 2));
							Y = PERMUTE_PS(V2, _MM_SHUFFLE(1, 1, 1, 1));
							X = PERMUTE_PS(V2, _MM_SHUFFLE(0, 0, 0, 0));

							vTemp = FMADD_PS(Z, row2, row3);
							vTemp2 = _mm_mul_ps(Y, row1);
							vTemp3 = _mm_mul_ps(X, row0);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);
							_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							// Result 3
							Z = PERMUTE_PS(V3, _MM_SHUFFLE(2, 2, 2, 2));
							Y = PERMUTE_PS(V3, _MM_SHUFFLE(1, 1, 1, 1));
							X = PERMUTE_PS(V3, _MM_SHUFFLE(0, 0, 0, 0));

							vTemp = FMADD_PS(Z, row2, row3);
							vTemp2 = _mm_mul_ps(Y, row1);
							vTemp3 = _mm_mul_ps(X, row0);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);
							_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							// Result 4
							Z = PERMUTE_PS(V4, _MM_SHUFFLE(2, 2, 2, 2));
							Y = PERMUTE_PS(V4, _MM_SHUFFLE(1, 1, 1, 1));
							X = PERMUTE_PS(V4, _MM_SHUFFLE(0, 0, 0, 0));

							vTemp = FMADD_PS(Z, row2, row3);
							vTemp2 = _mm_mul_ps(Y, row1);
							vTemp3 = _mm_mul_ps(X, row0);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);
							_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							i += 4;
						}
					}
				}
			}

			if (!(reinterpret_cast<uintptr_t>(pOutputStream) & 0xF) && !(outputStride & 0xF))
			{
				// Aligned output
				for (; i < vectorCount; ++i)
				{
					VECTOR V = VEC::LoadFloat3(reinterpret_cast<const Float3*>(pInputVector));
					pInputVector += inputStride;

					VECTOR Z = PERMUTE_PS(v, _MM_SHUFFLE(2, 2, 2, 2));
					VECTOR Y = PERMUTE_PS(V, _MM_SHUFFLE(1, 1, 1, 1));
					VECTOR X = PERMUTE_PS(V, _MM_SHUFFLE(0, 0, 0, 0));

					VECTOR vTemp = FMADD_PS(Z, row2, row3);
					VECTOR vTemp2 = _mm_mul_ps(Y, row1);
					VECTOR vTemp3 = _mm_mul_ps(X, row0);
					vTemp = _mm_add_ps(vTemp, vTemp2);
					vTemp = _mm_add_ps(vTemp, vTemp3);

					STREAM_PS(reinterpret_cast<float*>(pOutputVector), vTemp);
					pOutputVector += outputStride;
				}
			}
			else
			{
				// Unaligned output
				for (; i < vectorCount; ++i)
				{
					VECTOR V = VEC::LoadFloat3(reinterpret_cast<const Float3*>(pInputVector));
					pInputVector += inputStride;

					VECTOR Z = PERMUTE_PS(V, _MM_SHUFFLE(2, 2, 2, 2));
					VECTOR Y = PERMUTE_PS(V, _MM_SHUFFLE(1, 1, 1, 1));
					VECTOR X = PERMUTE_PS(V, _MM_SHUFFLE(0, 0, 0, 0));

					VECTOR vTemp = FMADD_PS(Z, row2, row3);
					VECTOR vTemp2 = _mm_mul_ps(Y, row1);
					VECTOR vTemp3 = _mm_mul_ps(X, row0);
					vTemp = _mm_add_ps(vTemp, vTemp2);
					vTemp = _mm_add_ps(vTemp, vTemp3);

					_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), vTemp);
					pOutputVector += outputStride;
				}
			}

			SFENCE();

			return pOutputStream;
#endif
		}

#ifdef _PREFAST_
#pragma prefast(pop)
#endif

		FORCE_INLINE VECTOR VEC_CALLCONV TransformCoord(A_VECTOR v, A_MATRIX m) noexcept
		{
			VECTOR Z = VEC::SplatZ(v);
			VECTOR Y = VEC::SplatY(v);
			VECTOR X = VEC::SplatX(v);

			VECTOR Result = VEC::MultiplyAdd((Z, m.r[2], m.r[3]));
			Result = VEC::MultiplyAdd(Y, m.r[1], Result);
			Result = VEC::MultiplyAdd(X, m.r[0], Result);

			VECTOR W = VEC::SplatW(Result);

			return VEC::Divide(Result, W);
		}

#ifdef _PREFAST_
#pragma prefast(push)
#pragma prefast(disable : 26015 26019, "PREfast noise: Esp:1307" )
#endif

		_Use_decl_annotations_
		FORCE_INLINE Float3* VEC_CALLCONV TransformCoordStream(Float3* pOutputStream, size_t outputStride, const Float3* pInputStream, size_t inputStride, size_t vectorCount, A_MATRIX m) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pOutputStream != nullptr);
			assert(pInputStream != nullptr);

			assert(inputStride >= sizeof(Float3));
			_Analysis_assume_(inputStride >= sizeof(Float3));

			assert(outputStride >= sizeof(Float3));
			_Analysis_assume_(outputStride >= sizeof(Float3));
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
				VECTOR V = VEC::LoadFloat3(reinterpret_cast<const Float3*>(pInputVector));
				VECTOR Z = VEC::SplatZ(V);
				VECTOR Y = VEC::SplatY(V);
				VECTOR X = VEC::SplatX(V);

				VECTOR Result = VEC::MultiplyAdd(Z, row2, row3);
				Result = VEC::MultiplyAdd(Y, row1, Result);
				Result = VEC::MultiplyAdd(X, row0, Result);

				VECTOR W = VEC::SplatW(Result);

				Result = VEC::Divide(Result, W);

				VEC::StoreFloat3(reinterpret_cast<Float3*>(pOutputVector), Result);

				pInputVector += inputStride;
				pOutputVector += outputStride;
			}

			return pOutputStream;

#elif defined(_SSE2_INTRINSICS_)
			auto pInputVector = reinterpret_cast<const uint8_t*>(pInputStream);
			auto pOutputVector = reinterpret_cast<uint8_t*>(pOutputStream);

			const VECTOR row0 = m.r[0];
			const VECTOR row1 = m.r[1];
			const VECTOR row2 = m.r[2];
			const VECTOR row3 = m.r[3];

			size_t i = 0;
			size_t four = vectorCount >> 2;
			if (four > 0)
			{
				if (inputStride == sizeof(Float3))
				{
					if (outputStride == sizeof(Float3))
					{
						if (!(reinterpret_cast<uintptr_t>(pOutputStream) & 0xF))
						{
							// Packed input, aligned & packed output
							for (size_t j = 0; j < four; ++j)
							{
								__m128 V1 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
								__m128 L2 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 16));
								__m128 L3 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 32));
								pInputVector += sizeof(Float3) * 4;

								// Unpack the 4 vectors (.w components are junk)
								3UNPACK3INTO4(V1, L2, L3);

								// Result 1
								VECTOR Z = PERMUTE_PS(V1, _MM_SHUFFLE(2, 2, 2, 2));
								VECTOR Y = PERMUTE_PS(V1, _MM_SHUFFLE(1, 1, 1, 1));
								VECTOR X = PERMUTE_PS(V1, _MM_SHUFFLE(0, 0, 0, 0));

								VECTOR vTemp = FMADD_PS(Z, row2, row3);
								VECTOR vTemp2 = _mm_mul_ps(Y, row1);
								VECTOR vTemp3 = _mm_mul_ps(X, row0);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								VECTOR W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));

								V1 = _mm_div_ps(vTemp, W);

								// Result 2
								Z = PERMUTE_PS(V2, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V2, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V2, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = FMADD_PS(Z, row2, row3);
								vTemp2 = _mm_mul_ps(Y, row1);
								vTemp3 = _mm_mul_ps(X, row0);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));

								V2 = _mm_div_ps(vTemp, W);

								// Result 3
								Z = PERMUTE_PS(V3, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V3, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V3, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = FMADD_PS(Z, row2, row3);
								vTemp2 = _mm_mul_ps(Y, row1);
								vTemp3 = _mm_mul_ps(X, row0);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));

								V3 = _mm_div_ps(vTemp, W);

								// Result 4
								Z = PERMUTE_PS(V4, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V4, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V4, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = FMADD_PS(Z, row2, row3);
								vTemp2 = _mm_mul_ps(Y, row1);
								vTemp3 = _mm_mul_ps(X, row0);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));

								V4 = _mm_div_ps(vTemp, W);

								// Pack and store the vectors
								3PACK4INTO3(vTemp);
								STREAM_PS(reinterpret_cast<float*>(pOutputVector), V1);
								STREAM_PS(reinterpret_cast<float*>(pOutputVector + 16), vTemp);
								STREAM_PS(reinterpret_cast<float*>(pOutputVector + 32), V3);
								pOutputVector += sizeof(Float3) * 4;
								i += 4;
							}
						}
						else
						{
							// Packed input, unaligned & packed output
							for (size_t j = 0; j < four; ++j)
							{
								__m128 V1 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
								__m128 L2 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 16));
								__m128 L3 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 32));
								pInputVector += sizeof(Float3) * 4;

								// Unpack the 4 vectors (.w components are junk)
								3UNPACK3INTO4(V1, L2, L3);

								// Result 1
								VECTOR Z = PERMUTE_PS(V1, _MM_SHUFFLE(2, 2, 2, 2));
								VECTOR Y = PERMUTE_PS(V1, _MM_SHUFFLE(1, 1, 1, 1));
								VECTOR X = PERMUTE_PS(V1, _MM_SHUFFLE(0, 0, 0, 0));

								VECTOR vTemp = FMADD_PS(Z, row2, row3);
								VECTOR vTemp2 = _mm_mul_ps(Y, row1);
								VECTOR vTemp3 = _mm_mul_ps(X, row0);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								VECTOR W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));

								V1 = _mm_div_ps(vTemp, W);

								// Result 2
								Z = PERMUTE_PS(V2, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V2, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V2, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = FMADD_PS(Z, row2, row3);
								vTemp2 = _mm_mul_ps(Y, row1);
								vTemp3 = _mm_mul_ps(X, row0);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));

								V2 = _mm_div_ps(vTemp, W);

								// Result 3
								Z = PERMUTE_PS(V3, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V3, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V3, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = FMADD_PS(Z, row2, row3);
								vTemp2 = _mm_mul_ps(Y, row1);
								vTemp3 = _mm_mul_ps(X, row0);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));

								V3 = _mm_div_ps(vTemp, W);

								// Result 4
								Z = PERMUTE_PS(V4, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V4, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V4, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = FMADD_PS(Z, row2, row3);
								vTemp2 = _mm_mul_ps(Y, row1);
								vTemp3 = _mm_mul_ps(X, row0);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));

								V4 = _mm_div_ps(vTemp, W);

								// Pack and store the vectors
								3PACK4INTO3(vTemp);
								_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), V1);
								_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector + 16), vTemp);
								_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector + 32), V3);
								pOutputVector += sizeof(Float3) * 4;
								i += 4;
							}
						}
					}
					else
					{
						// Packed input, unpacked output
						for (size_t j = 0; j < four; ++j)
						{
							__m128 V1 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
							__m128 L2 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 16));
							__m128 L3 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 32));
							pInputVector += sizeof(Float3) * 4;

							// Unpack the 4 vectors (.w components are junk)
							3UNPACK3INTO4(V1, L2, L3);

							// Result 1
							VECTOR Z = PERMUTE_PS(V1, _MM_SHUFFLE(2, 2, 2, 2));
							VECTOR Y = PERMUTE_PS(V1, _MM_SHUFFLE(1, 1, 1, 1));
							VECTOR X = PERMUTE_PS(V1, _MM_SHUFFLE(0, 0, 0, 0));

							VECTOR vTemp = FMADD_PS(Z, row2, row3);
							VECTOR vTemp2 = _mm_mul_ps(Y, row1);
							VECTOR vTemp3 = _mm_mul_ps(X, row0);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);

							VECTOR W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));

							vTemp = _mm_div_ps(vTemp, W);
							VEC::StoreFloat3(reinterpret_cast<Float3*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							// Result 2
							Z = PERMUTE_PS(V2, _MM_SHUFFLE(2, 2, 2, 2));
							Y = PERMUTE_PS(V2, _MM_SHUFFLE(1, 1, 1, 1));
							X = PERMUTE_PS(V2, _MM_SHUFFLE(0, 0, 0, 0));

							vTemp = FMADD_PS(Z, row2, row3);
							vTemp2 = _mm_mul_ps(Y, row1);
							vTemp3 = _mm_mul_ps(X, row0);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);

							W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));

							vTemp = _mm_div_ps(vTemp, W);
							VEC::StoreFloat3(reinterpret_cast<Float3*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							// Result 3
							Z = PERMUTE_PS(V3, _MM_SHUFFLE(2, 2, 2, 2));
							Y = PERMUTE_PS(V3, _MM_SHUFFLE(1, 1, 1, 1));
							X = PERMUTE_PS(V3, _MM_SHUFFLE(0, 0, 0, 0));

							vTemp = FMADD_PS(Z, row2, row3);
							vTemp2 = _mm_mul_ps(Y, row1);
							vTemp3 = _mm_mul_ps(X, row0);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);

							W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));

							vTemp = _mm_div_ps(vTemp, W);
							VEC::StoreFloat3(reinterpret_cast<Float3*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							// Result 4
							Z = PERMUTE_PS(V4, _MM_SHUFFLE(2, 2, 2, 2));
							Y = PERMUTE_PS(V4, _MM_SHUFFLE(1, 1, 1, 1));
							X = PERMUTE_PS(V4, _MM_SHUFFLE(0, 0, 0, 0));

							vTemp = FMADD_PS(Z, row2, row3);
							vTemp2 = _mm_mul_ps(Y, row1);
							vTemp3 = _mm_mul_ps(X, row0);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);

							W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));

							vTemp = _mm_div_ps(vTemp, W);
							VEC::StoreFloat3(reinterpret_cast<FLoat3*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							i += 4;
						}
					}
				}
			}

			for (; i < vectorCount; i++)
			{
				VECTOR V = VEC::LoadFloat3(reinterpret_cast<const Float3*>(pInputVector));
				pInputVector += inputStride;

				VECTOR Z = PERMUTE_PS(V, _MM_SHUFFLE(2, 2, 2, 2));
				VECTOR Y = PERMUTE_PS(V, _MM_SHUFFLE(1, 1, 1, 1));
				VECTOR X = PERMUTE_PS(V, _MM_SHUFFLE(0, 0, 0, 0));

				VECTOR vTemp = FMADD_PS(Z, row2, row3);
				VECTOR vTemp2 = _mm_mul_ps(Y, row1);
				VECTOR vTemp3 = _mm_mul_ps(X, row0);
				vTemp = _mm_add_ps(vTemp, vTemp2);
				vTemp = _mm_add_ps(vTemp, vTemp3);

				VECTOR W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));

				vTemp = _mm_div_ps(vTemp, W);

				VEC::StoreFloat3(reinterpret_cast<Float3*>(pOutputVector), vTemp);
				pOutputVector += outputStride;
			}

			SFENCE();

			return pOutputStream;
#endif
		}

#ifdef _PREFAST_
#pragma prefast(pop)
#endif

		FORCE_INLINE VECTOR VEC_CALLCONV TransformNormal(A_VECTOR v, A_MATRIX m) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR Z = VEC::SplatZ(v);
			VECTOR Y = VEC::SplatY(v);
			VECTOR X = VEC::SplatX(v);

			VECTOR Result = VEC::Multiply(Z, m.r[2]);
			Result = VEC::MultiplyAdd(Y, m.r[1], Result);
			Result = VEC::MultiplyAdd(X, m.r[0], Result);

			return Result;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vResult = PERMUTE_PS(v, _MM_SHUFFLE(2, 2, 2, 2));
			vResult = _mm_mul_ps(vResult, m.r[2]);
			VECTOR vTemp = PERMUTE_PS(v, _MM_SHUFFLE(1, 1, 1, 1)); // Y
			vResult = FMADD_PS(vTemp, m.r[1], vResult);
			vTemp = PERMUTE_PS(v, _MM_SHUFFLE(0, 0, 0, 0)); // X
			vResult = FMADD_PS(vTemp, m.r[0], vResult);
			
			return vResult;
#endif
		}

#ifdef _PREFAST_
#pragma prefast(push)
#pragma prefast(disable : 26015 26019, "PREfast noise: Esp:1307" )
#endif

		_Use_decl_annotations_
		FORCE_INLINE Float3* VEC_CALLCONV TransformNormalStream(Float3* pOutputStream, size_t outputStride, const Float3* pInputStream, size_t inputStride, size_t vectorCount, A_MATRIX m) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pOutputStream != nullptr);
			assert(pInputStream != nullptr);

			assert(inputStride >= sizeof(Float3));
			_Analysis_assume_(inputStride >= sizeof(Float3));

			assert(outputStride >= sizeof(Float3));
			_Analysis_assume_(outputStride >= sizeof(Float3));
#endif

#if defined(_NO_INTRINSICS_)
			const uint8_t* pInputVector = reinterpret_cast<const uint8_t*>(pInputStream);
			uint8_t* pOutputVector = reinterpret_cast<uint8_t*>(pOutputStream);

			const VECTOR row0 = m.r[0];
			const VECTOR row1 = m.r[1];
			const VECTOR row2 = m.r[2];

			for (size_t i = 0; i < vectorCount; i++)
			{
				VECTOR V = VEC::LoadFloat3(reinterpret_cast<const Float3*>(pInputVector));
				VECTOR Z = VEC::SplatZ(V);
				VECTOR Y = VEC::SplatY(V);
				VECTOR X = VEC::SplatX(V);

				VECTOR Result = VEC::Multiply(Z, row2);
				Result = VEC::MultiplyAdd(Y, row1, Result);
				Result = VEC::MultiplyAdd(X, row0, Result);

				VEC::StoreFloat3(reinterpret_cast<Float3*>(pOutputVector), Result);

				pInputVector += inputStride;
				pOutputVector += outputStride;
			}

			return pOutputStream;

#elif defined(_XM_SSE_INTRINSICS_)
			const uint8_t* pInputVector = reinterpret_cast<const uint8_t*>(pInputStream);
			uint8_t* pOutputVector = reinterpret_cast<uint8_t*>(pOutputStream);

			const VECTOR row0 = m.r[0];
			const VECTOR row1 = m.r[1];
			const VECTOR row2 = m.r[2];

			size_t i = 0;
			size_t four = vectorCount >> 2;
			if (four > 0)
			{
				if (inputStride == sizeof(Float3))
				{
					if (outputStride == sizeof(Float3))
					{
						if (!(reinterpret_cast<uintptr_t>(pOutputStream) & 0xF))
						{
							// Packed input, aligned & packed output
							for (size_t j = 0; j < four; ++j)
							{
								__m128 V1 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
								__m128 L2 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 16));
								__m128 L3 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 32));
								pInputVector += sizeof(Float3) * 4;

								// Unpack the 4 vectors (.w components are junk)
								3UNPACK3INTO4(V1, L2, L3);

								// Result 1
								VECTOR Z = PERMUTE_PS(V1, _MM_SHUFFLE(2, 2, 2, 2));
								VECTOR Y = PERMUTE_PS(V1, _MM_SHUFFLE(1, 1, 1, 1));
								VECTOR X = PERMUTE_PS(V1, _MM_SHUFFLE(0, 0, 0, 0));

								VECTOR vTemp = _mm_mul_ps(Z, row2);
								VECTOR vTemp2 = _mm_mul_ps(Y, row1);
								VECTOR vTemp3 = _mm_mul_ps(X, row0);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								V1 = _mm_add_ps(vTemp, vTemp3);

								// Result 2
								Z = PERMUTE_PS(V2, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V2, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V2, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = _mm_mul_ps(Z, row2);
								vTemp2 = _mm_mul_ps(Y, row1);
								vTemp3 = _mm_mul_ps(X, row0);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								V2 = _mm_add_ps(vTemp, vTemp3);

								// Result 3
								Z = PERMUTE_PS(V3, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V3, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V3, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = _mm_mul_ps(Z, row2);
								vTemp2 = _mm_mul_ps(Y, row1);
								vTemp3 = _mm_mul_ps(X, row0);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								V3 = _mm_add_ps(vTemp, vTemp3);

								// Result 4
								Z = PERMUTE_PS(V4, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V4, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V4, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = _mm_mul_ps(Z, row2);
								vTemp2 = _mm_mul_ps(Y, row1);
								vTemp3 = _mm_mul_ps(X, row0);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								V4 = _mm_add_ps(vTemp, vTemp3);

								// Pack and store the vectors
								3PACK4INTO3(vTemp);
								STREAM_PS(reinterpret_cast<float*>(pOutputVector), V1);
								STREAM_PS(reinterpret_cast<float*>(pOutputVector + 16), vTemp);
								STREAM_PS(reinterpret_cast<float*>(pOutputVector + 32), V3);
								pOutputVector += sizeof(Float3) * 4;
								i += 4;
							}
						}
						else
						{
							// Packed input, unaligned & packed output
							for (size_t j = 0; j < four; ++j)
							{
								__m128 V1 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
								__m128 L2 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 16));
								__m128 L3 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 32));
								pInputVector += sizeof(Float3) * 4;

								// Unpack the 4 vectors (.w components are junk)
								3UNPACK3INTO4(V1, L2, L3);

								// Result 1
								VECTOR Z = PERMUTE_PS(V1, _MM_SHUFFLE(2, 2, 2, 2));
								VECTOR Y = PERMUTE_PS(V1, _MM_SHUFFLE(1, 1, 1, 1));
								VECTOR X = PERMUTE_PS(V1, _MM_SHUFFLE(0, 0, 0, 0));

								VECTOR vTemp = _mm_mul_ps(Z, row2);
								VECTOR vTemp2 = _mm_mul_ps(Y, row1);
								VECTOR vTemp3 = _mm_mul_ps(X, row0);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								V1 = _mm_add_ps(vTemp, vTemp3);

								// Result 2
								Z = PERMUTE_PS(V2, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V2, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V2, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = _mm_mul_ps(Z, row2);
								vTemp2 = _mm_mul_ps(Y, row1);
								vTemp3 = _mm_mul_ps(X, row0);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								V2 = _mm_add_ps(vTemp, vTemp3);

								// Result 3
								Z = PERMUTE_PS(V3, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V3, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V3, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = _mm_mul_ps(Z, row2);
								vTemp2 = _mm_mul_ps(Y, row1);
								vTemp3 = _mm_mul_ps(X, row0);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								V3 = _mm_add_ps(vTemp, vTemp3);

								// Result 4
								Z = PERMUTE_PS(V4, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V4, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V4, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = _mm_mul_ps(Z, row2);
								vTemp2 = _mm_mul_ps(Y, row1);
								vTemp3 = _mm_mul_ps(X, row0);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								V4 = _mm_add_ps(vTemp, vTemp3);

								// Pack and store the vectors
								3PACK4INTO3(vTemp);
								_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), V1);
								_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector + 16), vTemp);
								_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector + 32), V3);
								pOutputVector += sizeof(Float3) * 4;
								i += 4;
							}
						}
					}
					else
					{
						// Packed input, unpacked output
						for (size_t j = 0; j < four; ++j)
						{
							__m128 V1 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
							__m128 L2 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 16));
							__m128 L3 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 32));
							pInputVector += sizeof(Float3) * 4;

							// Unpack the 4 vectors (.w components are junk)
							3UNPACK3INTO4(V1, L2, L3);

							// Result 1
							VECTOR Z = PERMUTE_PS(V1, _MM_SHUFFLE(2, 2, 2, 2));
							VECTOR Y = PERMUTE_PS(V1, _MM_SHUFFLE(1, 1, 1, 1));
							VECTOR X = PERMUTE_PS(V1, _MM_SHUFFLE(0, 0, 0, 0));

							VECTOR vTemp = _mm_mul_ps(Z, row2);
							VECTOR vTemp2 = _mm_mul_ps(Y, row1);
							VECTOR vTemp3 = _mm_mul_ps(X, row0);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);

							VEC::StoreFloat3(reinterpret_cast<Float3*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							// Result 2
							Z = PERMUTE_PS(V2, _MM_SHUFFLE(2, 2, 2, 2));
							Y = PERMUTE_PS(V2, _MM_SHUFFLE(1, 1, 1, 1));
							X = PERMUTE_PS(V2, _MM_SHUFFLE(0, 0, 0, 0));

							vTemp = _mm_mul_ps(Z, row2);
							vTemp2 = _mm_mul_ps(Y, row1);
							vTemp3 = _mm_mul_ps(X, row0);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);

							VEC::StoreFloat3(reinterpret_cast<Float3*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							// Result 3
							Z = PERMUTE_PS(V3, _MM_SHUFFLE(2, 2, 2, 2));
							Y = PERMUTE_PS(V3, _MM_SHUFFLE(1, 1, 1, 1));
							X = PERMUTE_PS(V3, _MM_SHUFFLE(0, 0, 0, 0));

							vTemp = _mm_mul_ps(Z, row2);
							vTemp2 = _mm_mul_ps(Y, row1);
							vTemp3 = _mm_mul_ps(X, row0);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);

							VEC::StoreFloat3(reinterpret_cast<Float3*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							// Result 4
							Z = PERMUTE_PS(V4, _MM_SHUFFLE(2, 2, 2, 2));
							Y = PERMUTE_PS(V4, _MM_SHUFFLE(1, 1, 1, 1));
							X = PERMUTE_PS(V4, _MM_SHUFFLE(0, 0, 0, 0));

							vTemp = _mm_mul_ps(Z, row2);
							vTemp2 = _mm_mul_ps(Y, row1);
							vTemp3 = _mm_mul_ps(X, row0);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);

							VEC::StoreFloat3(reinterpret_cast<Float3*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							i += 4;
						}
					}
				}
			}

			for (; i < vectorCount; i++)
			{
				VECTOR V = VEC::LoadFloat3(reinterpret_cast<const Float3*>(pInputVector));
				pInputVector += inputStride;

				VECTOR Z = PERMUTE_PS(V, _MM_SHUFFLE(2, 2, 2, 2));
				VECTOR Y = PERMUTE_PS(V, _MM_SHUFFLE(1, 1, 1, 1));
				VECTOR X = PERMUTE_PS(V, _MM_SHUFFLE(0, 0, 0, 0));

				VECTOR vTemp = _mm_mul_ps(Z, row2);
				VECTOR vTemp2 = _mm_mul_ps(Y, row1);
				VECTOR vTemp3 = _mm_mul_ps(X, row0);
				vTemp = _mm_add_ps(vTemp, vTemp2);
				vTemp = _mm_add_ps(vTemp, vTemp3);

				VEC::StoreFloat3(reinterpret_cast<Float3*>(pOutputVector), vTemp);
				pOutputVector += outputStride;
			}

			SFENCE();

			return pOutputStream;
#endif
		}

#ifdef _PREFAST_
#pragma prefast(pop)
#endif

		FORCE_INLINE VECTOR VEC_CALLCONV Project(A_VECTOR v, 
			float viewportX, float viewportY, 
			float viewportWidth, float viewportHeight, 
			float viewportMinZ, float viewportMaxZ, 
			A_MATRIX projection, B_MATRIX view, B_MATRIX world) noexcept
		{
			const float halfViewportWidth = viewportWidth * 0.5f;
			const float halfViewportHeight = viewportHeight * 0.5f;

			VECTOR scale = VEC::Set(halfViewportWidth, -halfViewportHeight, viewportMaxZ - viewportMinZ, 0.0f);
			VECTOR offset = VEC::Set(viewportX + halfViewportWidth, viewportY + halfViewportHeight, viewportMinZ, 0.0f);

			MATRIX transform = MAT::Multiply(world, view);
			transform = MAT::Multiply(transform, projection);

			VECTOR Result = TransformCoord(v, transform);

			Result = VEC::MultiplyAdd(Result, scale, offset);

			return Result;
		}

#ifdef _PREFAST_
#pragma prefast(push)
#pragma prefast(disable : 26015 26019, "PREfast noise: Esp:1307" )
#endif

		_Use_decl_annotations_
		FORCE_INLINE Float3* VEC_CALLCONV ProjectStream(Float3* pOutput, 
			size_t outputStride, 
			const Float3* pInputStream, 
			size_t inputStride, size_t vectorCount, 
			float viewportX, float viewportY, 
			float viewportWidth, float viewportHeight, 
			float viewportMinZ, float viewportMaxZ, 
			A_MATRIX projection, B_MATRIX view, B_MATRIX world) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pOutputStream != nullptr);
			assert(pInputStream != nullptr);

			assert(inputStride >= sizeof(Float3));
			_Analysis_assume_(inputStride >= sizeof(Float3));

			assert(outputStride >= sizeof(Float3));
			_Analysis_assume_(outputStride >= sizeof(Float3));
#endif

#if defined(_NO_INTRINSICS_)
			const float halfViewportWidth = viewportWidth * 0.5f;
			const float halfViewportHeight = viewportHeight * 0.5f;

			VECTOR scale = VEC::Set(halfViewportWidth, -halfViewportHeight, viewportMaxZ - viewportMinZ, 1.0f);
			VECTOR offset = VEC::Set(viewportX + halfViewportWidth, viewportY + halfViewportHeight, viewportMinZ, 0.0f);

			MATRIX transform = MAT::Multiply(world, view);
			transform = MAT::Multiply(transform, projection);

			const uint8_t* pInputVector = reinterpret_cast<const uint8_t*>(pInputStream);
			uint8_t* pOutputVector = reinterpret_cast<uint8_t*>(pOutputStream);

			for (size_t i = 0; i < vectorCount; i++)
			{
				VECTOR V = VEC::LoadFloat3(reinterpret_cast<const Float3*>(pInputVector));

				VECTOR Result = TransformCoord(V, transform);
				Result = VEC::MultiplyAdd(Result, scale, offset);

				VEC::StoreFloat3(reinterpret_cast<Float3*>(pOutputVector), Result);

				pInputVector += inputStride;
				pOutputVector += outputStride;
			}

			return pOutputStream;

#elif defined(_SSE2_INTRINSICS_)
			const float halfViewportWidth = viewportWidth * 0.5f;
			const float halfViewportHeight = viewportHeight * 0.5f;

			VECTOR scale = VEC::Set(halfViewportWidth, -halfViewportHeight, viewportMaxZ - viewportMinZ, 1.0f);
			VECTOR offset = VEC::Set(viewportX + halfViewportWidth, viewportY + halfViewportHeight, viewportMinZ, 0.0f);

			MATRIX transform = MAT::Multiply(world, view);
			transform = MAT::Multiply(transform, projection);

			const uint8_t* pInputVector = reinterpret_cast<const uint8_t*>(pInputStream);
			uint8_t* pOutputVector = reinterpret_cast<uint8_t*>(pOutputStream);

			size_t i = 0;
			size_t four = vectorCount >> 2;
			if (four > 0)
			{
				if (inputStride == sizeof(Float3))
				{
					if (outputStride == sizeof(Float3))
					{
						if (!(reinterpret_cast<uintptr_t>(pOutputStream) & 0xF))
						{
							// Packed input, aligned & packed output
							for (size_t j = 0; j < four; ++j)
							{
								__m128 V1 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
								__m128 L2 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 16));
								__m128 L3 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 32));
								pInputVector += sizeof(Float3) * 4;

								// Unpack the 4 vectors (.w components are junk)
								3UNPACK3INTO4(V1, L2, L3);

								// Result 1
								VECTOR Z = PERMUTE_PS(V1, _MM_SHUFFLE(2, 2, 2, 2));
								VECTOR Y = PERMUTE_PS(V1, _MM_SHUFFLE(1, 1, 1, 1));
								VECTOR X = PERMUTE_PS(V1, _MM_SHUFFLE(0, 0, 0, 0));

								VECTOR vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
								VECTOR vTemp2 = _mm_mul_ps(Y, transform.r[1]);
								VECTOR vTemp3 = _mm_mul_ps(X, transform.r[0]);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								VECTOR W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
								vTemp = _mm_div_ps(vTemp, W);
								V1 = FMADD_PS(vTemp, scale,offset);

								// Result 2
								Z = PERMUTE_PS(V2, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V2, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V2, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
								vTemp2 = _mm_mul_ps(Y, transform.r[1]);
								vTemp3 = _mm_mul_ps(X, transform.r[0]);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
								vTemp = _mm_div_ps(vTemp, W);
								V2 = FMADD_PS(vTemp, scale, offset);

								// Result 3
								Z = PERMUTE_PS(V3, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V3, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V3, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
								vTemp2 = _mm_mul_ps(Y, transform.r[1]);
								vTemp3 = _mm_mul_ps(X, transform.r[0]);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
								vTemp = _mm_div_ps(vTemp, W);
								V3 = FMADD_PS(vTemp, scale, offset);

								// Result 4
								Z = PERMUTE_PS(V4, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V4, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V4, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
								vTemp2 = _mm_mul_ps(Y, transform.r[1]);
								vTemp3 = _mm_mul_ps(X, transform.r[0]);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
								vTemp = _mm_div_ps(vTemp, W);
								V4 = FMADD_PS(vTemp, scale, offset);

								// Pack and store the vectors
								3PACK4INTO3(vTemp);
								STREAM_PS(reinterpret_cast<float*>(pOutputVector), V1);
								STREAM_PS(reinterpret_cast<float*>(pOutputVector + 16), vTemp);
								STREAM_PS(reinterpret_cast<float*>(pOutputVector + 32), V3);
								pOutputVector += sizeof(Float3) * 4;
								i += 4;
							}
						}
						else
						{
							// Packed input, unaligned & packed output
							for (size_t j = 0; j < four; ++j)
							{
								__m128 V1 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
								__m128 L2 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 16));
								__m128 L3 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 32));
								pInputVector += sizeof(Float3) * 4;

								// Unpack the 4 vectors (.w components are junk)
								3UNPACK3INTO4(V1, L2, L3);

								// Result 1
								VECTOR Z = PERMUTE_PS(V1, _MM_SHUFFLE(2, 2, 2, 2));
								VECTOR Y = PERMUTE_PS(V1, _MM_SHUFFLE(1, 1, 1, 1));
								VECTOR X = PERMUTE_PS(V1, _MM_SHUFFLE(0, 0, 0, 0));

								VECTOR vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
								VECTOR vTemp2 = _mm_mul_ps(Y, transform.r[1]);
								VECTOR vTemp3 = _mm_mul_ps(X, transform.r[0]);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								VECTOR W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
								vTemp = _mm_div_ps(vTemp, W);
								V1 = FMADD_PS(vTemp, scale, offset);

								// Result 2
								Z = PERMUTE_PS(V2, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V2, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V2, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
								vTemp2 = _mm_mul_ps(Y, transform.r[1]);
								vTemp3 = _mm_mul_ps(X, transform.r[0]);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
								vTemp = _mm_div_ps(vTemp, W);
								V2 = FMADD_PS(vTemp, scale, offset);

								// Result 3
								Z = PERMUTE_PS(V3, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V3, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V3, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
								vTemp2 = _mm_mul_ps(Y, transform.r[1]);
								vTemp3 = _mm_mul_ps(X, transform.r[0]);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
								vTemp = _mm_div_ps(vTemp, W);
								V3 = FMADD_PS(vTemp, scale, offset);

								// Result 4
								Z = PERMUTE_PS(V4, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V4, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V4, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
								vTemp2 = _mm_mul_ps(Y, transform.r[1]);
								vTemp3 = _mm_mul_ps(X, transform.r[0]);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
								vTemp = _mm_div_ps(vTemp, W);
								V4 = FMADD_PS(vTemp, scale, offset);

								// Pack and store the vectors
								3PACK4INTO3(vTemp);
								_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), V1);
								_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector + 16), vTemp);
								_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector + 32), V3);
								pOutputVector += sizeof(Float3) * 4;
								i += 4;
							}
						}
					}
					else
					{
						// Packed input, unpacked output
						for (size_t j = 0; j < four; ++j)
						{
							__m128 V1 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
							__m128 L2 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 16));
							__m128 L3 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 32));
							pInputVector += sizeof(Float3) * 4;

							// Unpack the 4 vectors (.w components are junk)
							3UNPACK3INTO4(V1, L2, L3);

							// Result 1
							VECTOR Z = PERMUTE_PS(V1, _MM_SHUFFLE(2, 2, 2, 2));
							VECTOR Y = PERMUTE_PS(V1, _MM_SHUFFLE(1, 1, 1, 1));
							VECTOR X = PERMUTE_PS(V1, _MM_SHUFFLE(0, 0, 0, 0));

							VECTOR vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
							VECTOR vTemp2 = _mm_mul_ps(Y, transform.r[1]);
							VECTOR vTemp3 = _mm_mul_ps(X, transform.r[0]);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);

							VECTOR W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
							vTemp = _mm_div_ps(vTemp, W);
							vTemp = FMADD_PS(vTemp, scale, offset);

							StoreFloat3(reinterpret_cast<Float3*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							// Result 2
							Z = PERMUTE_PS(V2, _MM_SHUFFLE(2, 2, 2, 2));
							Y = PERMUTE_PS(V2, _MM_SHUFFLE(1, 1, 1, 1));
							X = PERMUTE_PS(V2, _MM_SHUFFLE(0, 0, 0, 0));

							vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
							vTemp2 = _mm_mul_ps(Y, transform.r[1]);
							vTemp3 = _mm_mul_ps(X, transform.r[0]);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);

							W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
							vTemp = _mm_div_ps(vTemp, W);
							vTemp = FMADD_PS(vTemp, scale, offset);

							VEC::StoreFloat3(reinterpret_cast<Float3*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							// Result 3
							Z = PERMUTE_PS(V3, _MM_SHUFFLE(2, 2, 2, 2));
							Y = PERMUTE_PS(V3, _MM_SHUFFLE(1, 1, 1, 1));
							X = PERMUTE_PS(V3, _MM_SHUFFLE(0, 0, 0, 0));

							vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
							vTemp2 = _mm_mul_ps(Y, transform.r[1]);
							vTemp3 = _mm_mul_ps(X, transform.r[0]);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);

							W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
							vTemp = _mm_div_ps(vTemp, W);
							vTemp = FMADD_PS(vTemp, scale, offset);

							VEC::StoreFloat3(reinterpret_cast<Float3*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							// Result 4
							Z = PERMUTE_PS(V4, _MM_SHUFFLE(2, 2, 2, 2));
							Y = PERMUTE_PS(V4, _MM_SHUFFLE(1, 1, 1, 1));
							X = PERMUTE_PS(V4, _MM_SHUFFLE(0, 0, 0, 0));

							vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
							vTemp2 = _mm_mul_ps(Y, transform.r[1]);
							vTemp3 = _mm_mul_ps(X, transform.r[0]);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);

							W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
							vTemp = _mm_div_ps(vTemp, W);
							vTemp = FMADD_PS(vTemp, scale, offset);

							VEC::StoreFloat3(reinterpret_cast<Float3*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							i += 4;
						}
					}
				}
			}

			for (; i < vectorCount; i++)
			{
				VECTOR V = VEC::LoadFloat3(reinterpret_cast<const Float3*>(pInputVector));
				pInputVector += inputStride;

				VECTOR Z = PERMUTE_PS(V, _MM_SHUFFLE(2, 2, 2, 2));
				VECTOR Y = PERMUTE_PS(V, _MM_SHUFFLE(1, 1, 1, 1));
				VECTOR X = PERMUTE_PS(V, _MM_SHUFFLE(0, 0, 0, 0));

				VECTOR vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
				VECTOR vTemp2 = _mm_mul_ps(Y, transform.r[1]);
				VECTOR vTemp3 = _mm_mul_ps(X, transform.r[0]);
				vTemp = _mm_add_ps(vTemp, vTemp2);
				vTemp = _mm_add_ps(vTemp, vTemp3);

				VECTOR W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
				vTemp = _mm_div_ps(vTemp, W);
				vTemp = FMADD_PS(vTemp, scale, offset);

				VEC::StoreFloat3(reinterpret_cast<Float3*>(pOutputVector), vTemp);
				pOutputVector += outputStride;
			}

			SFENCE();

			return pOutputStream;
#endif
		}

#ifdef _PREFAST_
#pragma prefast(pop)
#endif

		FORCE_INLINE VECTOR VEC_CALLCONV Unproject(A_VECTOR v, 
			float viewportX, viewportY, 
			float viewportWidth, float viewportHeight, 
			float viewportMinZ, float viewportMaxZ, 
			A_MATRIX projection, B_MATRIX view, B_MATRIX world) noexcept
		{
			static const VECTOR_F32 D = { { { -1.0f, 1.0f, 0.0f, 0.0f } } };

			VECTOR scale = VEC::Set(viewportWidth * 0.5f, -viewportHeight * 0.5f, viewportMaxZ - viewportMinZ, 1.0f);
			scale = VEC::Reciprocal(scale);

			VECTOR offset = VEC::Set(-viewportX, -viewportY, -viewportMinZ, 0.0f);
			offset = VEC::MultiplyAdd(scale, offset, D.v);

			MATRIX transform = MAT::Multiply(world, view);
			transform = MAT::Multiply(transform, projection);
			transform = MAT::Inverse(nullptr, transform);

			VECTOR Result = VEC::MultiplyAdd(v, scale, offset);

			return TransformCoord(Result, transform);
		}

#ifdef _PREFAST_
#pragma prefast(push)
#pragma prefast(disable : 26015 26019, "PREfast noise: Esp:1307" )
#endif

		_Use_decl_annotations_
		FORCE_INLINE Float3* VEC_CALLCONV UnprojectStream(Float3* pOutputStream, 
			size_t outputStride, 
			const Float3* pInputStream, 
			size_t inputStride, 
			size_t vectorCount, 
			float viewportX, float viewportY, 
			float viewportWidth, float viewportHeight, 
			float viewportMinZ, float viewportMaxZ, 
			A_MATRIX projection, B_MATRIX view, B_MATRIX world) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pOutputStream != nullptr);
			assert(pInputStream != nullptr);

			assert(inputStride >= sizeof(Float3));
			_Analysis_assume_(inputStride >= sizeof(Float3));

			assert(outputStride >= sizeof(Float3));
			_Analysis_assume_(outputStride >= sizeof(Float3));
#endif

#if defined(_NO_INTRINSICS_)
			static const VECTOR_F32 D = { { { -1.0f, 1.0f, 0.0f, 0.0f } } };

			VECTOR scale = VEC::Set(viewportWidth * 0.5f, -viewportHeight * 0.5f, viewportMaxZ - viewportMinZ, 1.0f);
			scale = VEC::Reciprocal(scale);

			VECTOR offset = VEC::Set(-viewportX, -viewportY, -viewportMinZ, 0.0f);
			offset = VEC::MultiplyAdd(scale, offset, D.v);

			MATRIX transform = MAT::Multiply(world, view);
			transform = MAT::Multiply(transform, projection);
			transform = MAT::Inverse(nullptr, transform);

			const uint8_t* pInputVector = reinterpret_cast<const uint8_t*>(pInputStream);
			uint8_t* pOutputVector = reinterpret_cast<uint8_t*>(pOutputStream);

			for (size_t i = 0; i < vectorCount; i++)
			{
				VECTOR V = VEC::LoadFloat3(reinterpret_cast<const Float3*>(pInputVector));

				VECTOR Result = VEC::MultiplyAdd(V, scale, offset);

				Result = TransformCoord(Result, transform);

				VEC::StoreFloat3(reinterpret_cast<Float3*>(pOutputVector), Result);

				pInputVector += inputStride;
				pOutputVector += outputStride;
			}

			return pOutputStream;

#elif defined(_SSE2_INTRINSICS_)
			static const VECTOR_F32 D = { { { -1.0f, 1.0f, 0.0f, 0.0f } } };

			VECTOR scale = VEC::Set(viewportWidth * 0.5f, -viewportHeight * 0.5f, viewportMaxZ - viewportMinZ, 1.0f);
			scale = VEC::Reciprocal(scale);

			VECTOR offset = VEC::Set(-viewportX, -viewportY, -viewportMinZ, 0.0f);
			offset = _mm_mul_ps(scale, offset);
			offset = _mm_add_ps(offset, D);

			MATRIX transform = MAT::Multiply(world, view);
			transform = MAT::Multiply(transform, projection);
			transform = MAT::Inverse(nullptr, transform);

			auto pInputVector = reinterpret_cast<const uint8_t*>(pInputStream);
			auto pOutputVector = reinterpret_cast<uint8_t*>(pOutputStream);

			size_t i = 0;
			size_t four = vectorCount >> 2;
			if (four > 0)
			{
				if (inputStride == sizeof(Float3))
				{
					if (outputStride == sizeof(Float3))
					{
						if (!(reinterpret_cast<uintptr_t>(pOutputStream) & 0xF))
						{
							// Packed input, aligned & packed output
							for (size_t j = 0; j < four; ++j)
							{
								__m128 V1 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
								__m128 L2 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 16));
								__m128 L3 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 32));
								pInputVector += sizeof(Float3) * 4;

								// Unpack the 4 vectors (.w components are junk)
								3UNPACK3INTO4(V1, L2, L3);

								// Result 1
								V1 = FMADD_PS(V1, scale, offset);

								VECTOR Z = PERMUTE_PS(V1, _MM_SHUFFLE(2, 2, 2, 2));
								VECTOR Y = PERMUTE_PS(V1, _MM_SHUFFLE(1, 1, 1, 1));
								VECTOR X = PERMUTE_PS(V1, _MM_SHUFFLE(0, 0, 0, 0));

								VECTOR vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
								VECTOR vTemp2 = _mm_mul_ps(Y, transform.r[1]);
								VECTOR vTemp3 = _mm_mul_ps(X, transform.r[0]);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								VECTOR W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
								V1 = _mm_div_ps(vTemp, W);

								// Result 2
								V2 = FMADD_PS(V2, scale, offset);

								Z = PERMUTE_PS(V2, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V2, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V2, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
								vTemp2 = _mm_mul_ps(Y, transform.r[1]);
								vTemp3 = _mm_mul_ps(X, transform.r[0]);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
								V2 = _mm_div_ps(vTemp, W);

								// Result 3
								V3 = FMADD_PS(V3, scale, offset);

								Z = PERMUTE_PS(V3, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V3, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V3, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
								vTemp2 = _mm_mul_ps(Y, transform.r[1]);
								vTemp3 = _mm_mul_ps(X, transform.r[0]);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
								V3 = _mm_div_ps(vTemp, W);

								// Result 4
								V4 = FMADD_PS(V4, scale, offset);

								Z = PERMUTE_PS(V4, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V4, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V4, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
								vTemp2 = _mm_mul_ps(Y, transform.r[1]);
								vTemp3 = _mm_mul_ps(X, transform.r[0]);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
								V4 = _mm_div_ps(vTemp, W);

								// Pack and store the vectors
								3PACK4INTO3(vTemp);
								STREAM_PS(reinterpret_cast<float*>(pOutputVector), V1);
								STREAM_PS(reinterpret_cast<float*>(pOutputVector + 16), vTemp);
								STREAM_PS(reinterpret_cast<float*>(pOutputVector + 32), V3);
								pOutputVector += sizeof(Float3) * 4;
								i += 4;
							}
						}
						else
						{
							// Packed input, unaligned & packed output
							for (size_t j = 0; j < four; ++j)
							{
								__m128 V1 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
								__m128 L2 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 16));
								__m128 L3 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 32));
								pInputVector += sizeof(Float3) * 4;

								// Unpack the 4 vectors (.w components are junk)
								3UNPACK3INTO4(V1, L2, L3);

								// Result 1
								V1 = FMADD_PS(V1, scale, offset);

								VECTOR Z = PERMUTE_PS(V1, _MM_SHUFFLE(2, 2, 2, 2));
								VECTOR Y = PERMUTE_PS(V1, _MM_SHUFFLE(1, 1, 1, 1));
								VECTOR X = PERMUTE_PS(V1, _MM_SHUFFLE(0, 0, 0, 0));

								VECTOR vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
								VECTOR vTemp2 = _mm_mul_ps(Y, transform.r[1]);
								VECTOR vTemp3 = _mm_mul_ps(X, transform.r[0]);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								VECTOR W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
								V1 = _mm_div_ps(vTemp, W);

								// Result 2
								V2 = FMADD_PS(V2, scale, offset);

								Z = PERMUTE_PS(V2, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V2, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V2, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
								vTemp2 = _mm_mul_ps(Y, transform.r[1]);
								vTemp3 = _mm_mul_ps(X, transform.r[0]);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
								V2 = _mm_div_ps(vTemp, W);

								// Result 3
								V3 = FMADD_PS(V3, scale, offset);

								Z = PERMUTE_PS(V3, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V3, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V3, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
								vTemp2 = _mm_mul_ps(Y, transform.r[1]);
								vTemp3 = _mm_mul_ps(X, transform.r[0]);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
								V3 = _mm_div_ps(vTemp, W);

								// Result 4
								V4 = FMADD_PS(V4, scale, offset);

								Z = PERMUTE_PS(V4, _MM_SHUFFLE(2, 2, 2, 2));
								Y = PERMUTE_PS(V4, _MM_SHUFFLE(1, 1, 1, 1));
								X = PERMUTE_PS(V4, _MM_SHUFFLE(0, 0, 0, 0));

								vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
								vTemp2 = _mm_mul_ps(Y, transform.r[1]);
								vTemp3 = _mm_mul_ps(X, transform.r[0]);
								vTemp = _mm_add_ps(vTemp, vTemp2);
								vTemp = _mm_add_ps(vTemp, vTemp3);

								W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
								V4 = _mm_div_ps(vTemp, W);

								// Pack and store the vectors
								3PACK4INTO3(vTemp);
								_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), V1);
								_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector + 16), vTemp);
								_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector + 32), V3);
								pOutputVector += sizeof(Float3) * 4;
								i += 4;
							}
						}
					}
					else
					{
						// Packed input, unpacked output
						for (size_t j = 0; j < four; ++j)
						{
							__m128 V1 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
							__m128 L2 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 16));
							__m128 L3 = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector + 32));
							pInputVector += sizeof(Float3) * 4;

							// Unpack the 4 vectors (.w components are junk)
							3UNPACK3INTO4(V1, L2, L3);

							// Result 1
							V1 = FMADD_PS(V1, scale, offset);

							VECTOR Z = PERMUTE_PS(V1, _MM_SHUFFLE(2, 2, 2, 2));
							VECTOR Y = PERMUTE_PS(V1, _MM_SHUFFLE(1, 1, 1, 1));
							VECTOR X = PERMUTE_PS(V1, _MM_SHUFFLE(0, 0, 0, 0));

							VECTOR vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
							VECTOR vTemp2 = _mm_mul_ps(Y, transform.r[1]);
							VECTOR vTemp3 = _mm_mul_ps(X, transform.r[0]);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);

							VECTOR W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
							vTemp = _mm_div_ps(vTemp, W);

							VEC::StoreFloat3(reinterpret_cast<Float3*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							// Result 2
							V2 = FMADD_PS(V2, scale, offset);

							Z = PERMUTE_PS(V2, _MM_SHUFFLE(2, 2, 2, 2));
							Y = PERMUTE_PS(V2, _MM_SHUFFLE(1, 1, 1, 1));
							X = PERMUTE_PS(V2, _MM_SHUFFLE(0, 0, 0, 0));

							vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
							vTemp2 = _mm_mul_ps(Y, transform.r[1]);
							vTemp3 = _mm_mul_ps(X, transform.r[0]);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);

							W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
							vTemp = _mm_div_ps(vTemp, W);

							VEC::StoreFloat3(reinterpret_cast<Float3*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							// Result 3
							V3 = FMADD_PS(V3, scale, offset);

							Z = PERMUTE_PS(V3, _MM_SHUFFLE(2, 2, 2, 2));
							Y = PERMUTE_PS(V3, _MM_SHUFFLE(1, 1, 1, 1));
							X = PERMUTE_PS(V3, _MM_SHUFFLE(0, 0, 0, 0));

							vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
							vTemp2 = _mm_mul_ps(Y, transform.r[1]);
							vTemp3 = _mm_mul_ps(X, transform.r[0]);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);

							W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
							vTemp = _mm_div_ps(vTemp, W);

							VEC::StoreFloat3(reinterpret_cast<Float3*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							// Result 4
							V4 = FMADD_PS(V4, scale, offset);

							Z = PERMUTE_PS(V4, _MM_SHUFFLE(2, 2, 2, 2));
							Y = PERMUTE_PS(V4, _MM_SHUFFLE(1, 1, 1, 1));
							X = PERMUTE_PS(V4, _MM_SHUFFLE(0, 0, 0, 0));

							vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
							vTemp2 = _mm_mul_ps(Y, transform.r[1]);
							vTemp3 = _mm_mul_ps(X, transform.r[0]);
							vTemp = _mm_add_ps(vTemp, vTemp2);
							vTemp = _mm_add_ps(vTemp, vTemp3);

							W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
							vTemp = _mm_div_ps(vTemp, W);

							VEC::StoreFloat3(reinterpret_cast<Float3*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							i += 4;
						}
					}
				}
			}

			for (; i < vectorCount; i++)
			{
				VECTOR V = VEC::LoadFloat3(reinterpret_cast<const Float3*>(pInputVector));
				pInputVector += inputStride;

				V = _mm_mul_ps(V, scale);
				V = _mm_add_ps(V, offset);

				VECTOR Z = PERMUTE_PS(V, _MM_SHUFFLE(2, 2, 2, 2));
				VECTOR Y = PERMUTE_PS(V, _MM_SHUFFLE(1, 1, 1, 1));
				VECTOR X = PERMUTE_PS(V, _MM_SHUFFLE(0, 0, 0, 0));

				VECTOR vTemp = FMADD_PS(Z, transform.r[2], transform.r[3]);
				VECTOR vTemp2 = _mm_mul_ps(Y, transform.r[1]);
				VECTOR vTemp3 = _mm_mul_ps(X, transform.r[0]);
				vTemp = _mm_add_ps(vTemp, vTemp2);
				vTemp = _mm_add_ps(vTemp, vTemp3);

				VECTOR W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
				vTemp = _mm_div_ps(vTemp, W);

				VEC::StoreFloat3(reinterpret_cast<Float3*>(pOutputVector), vTemp);
				pOutputVector += outputStride;
			}

			SFENCE();

			return pOutputStream;
#endif
		}

#ifdef _PREFAST_
#pragma prefast(pop)
#endif
	}
}

#endif // !ULTREALITY_MATH_FLOAT3_INL
