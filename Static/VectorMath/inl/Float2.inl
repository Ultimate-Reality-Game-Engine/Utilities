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
		Vector::StoreFloat2(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE Float2& VEC_CALLCONV Float2::operator=(A_VECTOR v) noexcept
	{
		Vector::StoreFloat2(this, v);

		return *this;
	}

	FORCE_INLINE VECTOR VEC_CALLCONV Float2::Load() noexcept
	{
		return Vector::LoadFloat2(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV Float2::Store(A_VECTOR v) noexcept
	{
		Vector::StoreFloat2(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE AFloat2::AFloat2(A_VECTOR v) noexcept
	{
		Vector::StoreAFloat2(this, v);
	}

	_Use_decl_annotations_
	FORCE_INLINE AFloat2& VEC_CALLCONV AFloat2::operator=(A_VECTOR v) noexcept
	{
		Vector::StoreAFloat2(this, v);

		return *this;
	}

	FORCE_INLINE VECTOR VEC_CALLCONV AFloat2::Load() noexcept
	{
		return Vector::LoadAFloat2(this);
	}

	_Use_decl_annotations_
	FORCE_INLINE void VEC_CALLCONV AFloat2::Store(A_VECTOR v) noexcept
	{
		Vector::StoreAFloat2(this, v);
	}

	namespace Vector
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

	namespace Vector2
	{
		FORCE_INLINE bool VEC_CALLCONV Equal(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return (((V1.vector4_f32[0] == V2.vector4_f32[0]) && (V1.vector4_f32[1] == V2.vector4_f32[1])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmpeq_ps(V1, V2);
			
			// Ignore z and w
			return (((_mm_movemask_ps(vTemp) & 3) == 3) != 0);
#endif
		}

		FORCE_INLINE uint32_t VEC_CALLCONV EqualR(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			uint32_t CR = 0;
			if ((V1.vector4_f32[0] == V2.vector4_f32[0]) && 
				(V1.vector4_f32[1] == V2.vector4_f32[1]))
			{
				CR = CRMASK_CR6TRUE;
			}
			else if ((V1.vector4_f32[0] != V2.vector4_f32[0]) && 
				(V1.vector4_f32[1] != V2.vector4_f32[1]))
			{
				CR = CRMASK_CR6FALSE;
			}

			return CR;
#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmpeq_ps(V1, V2);

			// Ignore z and w
			int iTest = _mm_movemask_ps(vTemp) & 3;
			uint32_t CR = 0;
			if (iTest == 3)
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
			float dx = fabsf(V1.vector4_f32[0] - V2.vector4_f32[0]);
			float dy = fabsf(V1.vector4_f32[1] - V2.vector4_f32[1]);
			return ((dx <= epsilon.vector4_f32[0]) &&
				(dy <= epsilon.vector4_f32[1]));

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vDelta = _mm_sub_ps(V1, V2);

			VECTOR vTemp = _mm_setzero_ps();
			vTemp = _mm_sub_ps(vTemp, vDelta);
			vTemp = _mm_max_ps(vTemp, vDelta);
			vTemp = _mm_cmple_ps(vTemp, epsilon);

			// Ignore z and w
			return (((_mm_movemask_ps(vTemp) & 3) == 0x3) != 0);
#endif
		}

		FORCE_INLINE bool VEC_CALLCONV NotEqual(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return (((V1.vector4_f32[0] != V2.vector4_f32[0]) || (V1.vector4_f32[1] != V2.vector4_f32[1])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmpeq_ps(V!, V2);

			// Ignore z and w
			return (((_mm_movemask_ps(vTemp) & 3) != 3) != 0);
#endif
		}

		FORCE_INLINE bool VEC_CALLCONV Greater(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return (((V1.vector4_f32[0] > V2.vector4_f32[0]) && (V1.vector4_f32[1] > V2.vector4_f32[1])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmpgt_ps(V1, V2);

			// Ignore z and w
			return (((_mm_movemask_ps(vTemp) & 3) == 3) != 0);
#endif
		}

		FORCE_INLINE uint32_t VEC_CALLCONV GreaterR(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSiCS_)
			uint32_t CR = 0;
			if ((V1.vector4_f32[0] > V2.vector4_f32[0]) &&
				(V1.vector4_f32[1] > V2.vector4_f32[1]))
			{
				CR = CRMASK_CR6TRUE;
			}
			else if ((V1.vector4_f32[0] <= V2.vector4_f32[0]) &&
				(V1.vector4_f32[1] <= V2.vector4_f32[1]))
			{
				CR = CRMASK_CR6FALSE;
			}
			
			return CR;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmpgt_ps(V1, V2);
			int iTest = _mm_movemask_ps(vTemp) & 3;
			if (iTest == 3)
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
			return (((V1.vector4_f32[0] >= V2.vector4_f32[0]) && (V1.vector4_f32[1] >= V2.vector4_f32[1])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmpge_ps(V1, V2);

			return (((_mm_movemask_ps(vTemp) & 3) == 3) != 0);
#endif
		}

		FORCE_INLINE uint32_t VEC_CALLCONV GreaterOrEqualR(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			uint32_t CR = 0;
			if ((V1.vector4_f32[0] >= V2.vector4_f32[0]) &&
				(V1.vector4_f32[1] >= V2.vector4_f32[1]))
			{
				CR = CRMASK_CR6TRUE;
			}
			else if ((V1.vector4_f32[0] < V2.vector4_f32[0]) &&
				(V1.vector4_f32[1] < V2.vector4_f32[1]))
			{
				CR = CRMASK_CR6FALSE;
			}

			return CR;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmpge_ps(V1, V2);
			int iTest = _mm_movemask_ps(vTemp) & 3;
			uint32_t CR = 0;
			if (iTest == 3)
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
			return (((V1.vector4_f32[0] < V2.vector4_f32[0]) && (V1.vector4_f32[1] < V2.vector4_f32[1])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmplt_ps(V1, V2);

			return (((_mm_movemask_ps(vTemp) & 3) == 3) != 0);
#endif
		}

		FORCE_INLINE bool VEC_CALLCONV LessOrEqual(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return (((V1.vector4_f32[0] <= V2.vector4_f32[0]) && (V1.vector4_f32[1] <= V2.vector4_f32[1])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmple_ps(V1, V2);

			return (((_mm_movemask_ps(vTemp) & 3) == 3) != 0);
#endif
		}

		FORCE_INLINE bool VEC_CALLCONV InBounds(A_VECTOR v, A_VECTOR bounds) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return (((V.vector4_f32[0] <= bounds.vector4_f32[0] && v.vector4_f32[0] >= -bounds.vector4_f32[0]) &&
        			 (V.vector4_f32[1] <= bounds.vector4_f32[1] && v.vector4_f32[1] >= -bounds.vector4_f32[1])) != 0);

#elif defined(_SSE2_INTRINSICS_)
			// Test if v is less than or equal to the bounds
			VECTOR vTemp1 = _mm_cmple_ps(v, bounds);
			// Negate the bounds
			VECTOR vTemp2 = _mm_mul_ps(bounds, g_NegativeOne);
			// Test if negated bounds are greater or eqal to v
			vTemp2 = _mm_cmple_ps(vTemp2, v);
			// Blend the answers
			vTemp1 = _mm_and_ps(vTemp1, vTemp2);

			// return if x and y are in bounds, ignore z and w
			return(((_mm_movemask_ps(vTemp1) & 0x3) == 0x3) != 0);
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
        			ISNAN(v.vector4_f32[1]));

#elif defined(_SSE2_INTRINSICS_)
#if defined(__clang__) && defined(__FINITE_MATH_ONLY__)
			ALIGNED(16) float tmp[4];
			_mm_store_ps(tmp, v);
			
			return isnan(tmp[0]) || isnan(tmp[1]);
#else
			// Test against itself. NaN is always not equal
			VECTOR vTempNan = _mm_cmpneq_ps(v, v);
			// If x or y are NaN, the mask is non-zero
			return ((_mm_movemask_ps(vTempNan) & 3) != 0);
#endif
#endif
		}

#if !defined(_NO_INTRINSICS_) && defined(_MSC_VER) && !defined(__INTEL_COMPILER)
#pragma float_control(pop)
#endif

		FORCE_INLINE bool VEC_CALLCONV IsInfinite(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return (ISINF(v.vector4_f32[0]) || ISINF(v.vector4_f32[1]));

#elif defined(_SSE2_INTRINSICS_)
			// Mask off the sign bit
			__m128 vTemp = _mm_and_ps(v, g_AbsMask);
			
			// Compare to infinity
			vTemp = _mm_cmpeq_ps(vTemp, g_Infinity);

			// Check if x or y are infinite, ignore z and w
			return ((_mm_movemask_ps(vTemp) & 3) != 0);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Dot(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 Result;
			Result.f[0] =
				Result.f[1] =
				Result.f[2] =
				Result.f[3] = V1.vector4_f32[0] * V2.vector4_f32[0] + V1.vector4_f32[1] * V2.vector4_f32[1];
			
			return Result.v;

#elif defined(_SSE2_INTRINSICS_)
			// Perform dot product on x and y
			VECTOR vLengthSq = _mm_mul_ps(V1, V2);

			// Initialize vTemp with y splatt
			VECTOR vTemp = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(1, 1, 1, 1));

			// x + y
			vLengthSq = _mm_add_ps(vLengthSq, vTemp);
			vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(0, 0, 0, 0));

			return vLengthSq;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Cross(A_VECTOR V1, A_VECTOR V2) noexcept
		{
			// [V1.x * V2.y - V1.y * V2.x, V1.x * V2.y - V1.y * V2.x]
#if defined(_NO_INTRINSICS_)
			float fCross = (V1.vector4_f32[0] * V2.vector4_f32[1]) - (V1.vector4_f32[1] * V2.vector4_f32[0]);
			VECTOR_F32 vResult;
			vResult.f[0] = 
			vResult.f[1] = 
			vResult.f[2] = 
			vResult.f[3] = fCross;

			return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
			// Swap x and y
			VECTOR vResult = PERMUTE_PS(V2, _MM_SHUFFLE(0, 1, 0, 1));

			// Perform multiplications
			vResult = _mm_mul_ps(vResult, V1);

			// Splat y
			VECTOR vTemp = PERMUTE_PS(vResult, _MM_SHUFFLE(1, 1, 1, 1));

			// Perform the subtractions
			vResult = _mm_sub_ps(vResult, vTemp);

			// Splat the cross product
			vResult = PERMUTE_PS(vResult, _MM_SHUFFLE(0, 0, 0, 0));

			return vResult;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV LengthSq(A_VECTOR v) noexcept
		{
			return Dot(v, v);
		}

		FORCE_INLINE VECTOR VEC_CALLCONV ReciprocalLengthEst(A_VECTOR v) noexcept
		{
# if defined(_NO_INTRINSICS_)
			VECTOR Result;
			Result = LengthSq(v);
			Result = Vector::ReciprocalSqrtEst(Result);

			return Result;

#elif defined(_SSE4_INTRINSICS_)
			VECTOR vTemp = _mm_dp_ps(v, v, 0x3f);

			return _mm_rsqrt_ps(vTemp);

#elif defined(_SSE3_INTRINSICS_)
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			VECTOR vTemp = _mm_hadd_ps(vLengthSq, vLengthSq);
			vLengthSq = _mm_rsqrt_ss(vTemp);
			vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(0, 0, 0, 0));

			return vLengthSq;

#elif defined(_SSE2_INTRINSICS_)
			// Perform dot on x and y
			VECTOR vLengthSq = _mm_mul_ps(v, v);

			// Initialize vTemp to y splat
			VECTOR vTemp = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(1, 1, 1, 1));

			// x + y
			vLengthSq = _mm_add_ss(vLengthSq, vTemp);
			vLengthSq = _mm_rsqrt_ss(vLengthSq);
			vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(0, 0, 0, 0));

			return vLengthSq;
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
			VECTOR vTemp = _mm_dp_ds(v, v, 0x3f);
			VECTOR vLengthSq = _mm_sqrt_ps(vTemp);

			return _mm_div_ps(g_One, vLengthSq);

#elif defined(_SSE3_INTRINSICS_)
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			VECTOR vTemp = _mm_hadd_ps(vLengthSq, vLengthSq);
			vLengthSq = _mm_sqrt_ss(vTemp);
			vLengthSq = _mm_div_ss(g_One, vLengthSq);
			vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(0, 0, 0, 0));

			return vLengthSq;

#elif defined(_SSE2_INTRINSICS_)
			// Perform dot on x and y
			VECTOR vLengthSq = _mm_mul_ps(v, v);

			// Initialize vTemp with y splat
			VECTOR vTemp = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(1, 1, 1, 1));

			// x + y
			vLengthSq = _mm_add_ss(vLengthSq, vTemp);
			vLengthSq = _mm_sqrt_ss(vLengthSq);
			vLengthSq = _mm_div_ss(g_One, vLengthSq);
			vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(0, 0, 0, 0));

			return vLengthSq;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV LengthEst(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR Result;
			Result = LengthSq(v);
			Result = Vector::SqrtEst(Result);

			return Result;

#elif defined(_SSE4_INTRINSICS_)
			VECTOR vTemp = _mm_dp_ps(v, v, 0x3f);

			return _mm_sqrt_ps(vTemp);

#elif defined(_SSE3_INTRINSICS)
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			VECTOR vTemp = _mm_hadd_ps(vLengthSq, vLengthSq);
			vLengthSq = _mm_sqrt_ss(vTemp);
			vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(0, 0, 0, 0));

			return vLengthSq;

#elif defined(_SSE2_INTRINSICS_)
			// Perform the dot product on x and y
			VECTOR vLengthSq = _mm_mul_ps(v, v);

			// Initialize vTemp with y splat
			VECTOR vTemp = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(1, 1, 1, 1));

			// x + y
			vLengthSq = _mm_add_ss(vLengthSq, vTemp);
			vLengthSq = _mm_sqrt_ss(vLengthSq);
			vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(0, 0, 0, 0));

			return vLengthSq;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Length(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR Result;
			Result = LengthSq(v);
			Result = Vector::Sqrt(Result);

			return Result;

#elif defined(_SSE4_INTRINSICS_)
			VECTOR vTemp = _mm_dp_ps(v, v, 0x3f);

			return _mm_sqrt_ps(vTemp);

#elif defined(_SSE3_INTRINSICS_)
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			VECTOR vTemp = _mm_hadd_ps(vLengthSq, vLengthSq);
			vLengthSq = _mm_sqrt_ss(vTemp);
			vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(0, 0, 0, 0));

			return vLengthSq;

#elif defined(_SSE2_INTRINSICS_)
			// Perform dot product on x and y
			VECTOR vLengthSq = _mm_mul_ps(v, v);

			// Initialize vTemp with y splat
			VECTOR vTemp = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(1, 1, 1, 1));

			// x + y
			vLengthSq = _mm_add_ss(vLengthSq, vTemp);
			vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(0, 0, 0, 0));
			vLengthSq = _mm_sqrt_ps(vLengthSq);

			return vLengthSq;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV NormalizeEst(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR Result;
			Result = ReciprocalLength(v);
			Result = Vector::Multiply(v, Result);

			return Result;

#elif defined(_SSE4_INTRINSICS_)
			VECTOR vTemp = _mm_dp_ps(v, v, 0x3f);
			VECTOR vResult = _mm_rsqrt_ps(vTemp);

			return _mm_mul_ps(vResult, v);

#elif defined(_SSE3_INTRINSICS_)
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
			vLengthSq = _mm_rsqrt_ss(vLengthSq);
			vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(0, 0, 0, 0));
			vLengthSq = _mm_mul_ps(vLengthSq, v);

			return vLengthSq;

#elif defined(_SSE2_INTRINSICS_)
			// Perform dot product on x and y
			VECTOR vLengthSq = _mm_mul_ps(v, v);

			// Initialize vTemp with y splat
			VECTOR vTemp = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(1, 1, 1, 1));

			// x + y
			vLengthSq = _mm_add_ss(vLengthSq, vTemp);
			vLengthSq = _mm_rsqrt_ss(vLengthSq);
			vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(0, 0, 0, 0));
			vLengthSq = _mm_mul_ps(vLengthSq, v);

			return vLengthSq;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Normalize(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR Result = Length(v);
			float fLength = Result.vector4_f32[0];

			// Prevent divide by zero
			if (fLength > 0)
			{
				fLength = 1.0f / fLength;
			}

			Result.vector4_f32[0] = v.vector4_f32[0] * fLength;
			Result.vector4_f32[1] = v.vector4_f32[1] * fLength;
			Result.vector4_f32[2] = v.vector4_f32[2] * fLength;
			Result.vector4_f32[3] = v.vector4_f32[3] * fLength;

			return Result;

#elif defined(_SSE4_INTRINSICS_)
			VECTOR vLengthSq = _mm_dp_ps(v, v, 0x3f);
			// Prepare for the division
			VECTOR vResult = _mm_sqrt_ps(vLengthSq);
			// Create zero with a single instruction
			VECTOR vZeroMask = _mm_setzero_ps();
			// Test for a divide by zero (Must be FP to detect -0.0)
			vZeroMask = _mm_cmpneq_ps(vZeroMask, vResult);
			// Failsafe on zero (Or epsilon) length planes
			// If the length is infinity, set the elements to zero
			vLengthSq = _mm_cmpneq_ps(vLengthSq, g_Infinity);
			// Reciprocal mul to perform the normalization
			vResult = _mm_div_ps(V, vResult);
			// Any that are infinity, set to zero
			vResult = _mm_and_ps(vResult, vZeroMask);
			// Select qnan or result based on infinite length
			VECTOR vTemp1 = _mm_andnot_ps(vLengthSq, g_QNaN);
			VECTOR vTemp2 = _mm_and_ps(vResult, vLengthSq);
			vResult = _mm_or_ps(vTemp1, vTemp2);
			
			return vResult;

#elif defined(_SSE3_INTRINSICS_)
			// Perform the dot product on x and y only
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
			vLengthSq = _mm_moveldup_ps(vLengthSq);
			// Prepare for the division
			VECTOR vResult = _mm_sqrt_ps(vLengthSq);
			// Create zero with a single instruction
			VECTOR vZeroMask = _mm_setzero_ps();
			// Test for a divide by zero (Must be FP to detect -0.0)
			vZeroMask = _mm_cmpneq_ps(vZeroMask, vResult);
			// Failsafe on zero (Or epsilon) length planes
			// If the length is infinity, set the elements to zero
			vLengthSq = _mm_cmpneq_ps(vLengthSq, g_Infinity);
			// Reciprocal mul to perform the normalization
			vResult = _mm_div_ps(v, vResult);
			// Any that are infinity, set to zero
			vResult = _mm_and_ps(vResult, vZeroMask);
			// Select qnan or result based on infinite length
			VECTOR vTemp1 = _mm_andnot_ps(vLengthSq, g_QNaN);
			VECTOR vTemp2 = _mm_and_ps(vResult, vLengthSq);
			vResult = _mm_or_ps(vTemp1, vTemp2);
			
			return vResult;

#elif defined(_SSE2_INTRINSICS_)
			// Perform the dot product on x and y only
			VECTOR vLengthSq = _mm_mul_ps(v, v);
			VECTOR vTemp = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(1, 1, 1, 1));
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
			// Reciprocal mul to perform the normalization
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

		FORCE_INLINE VECTOR VEC_CALLCONV ClampLength(A_VECTOR v, float lengthMin, float lengthMax) noexcept
		{
			VECTOR clampMax = Vector::Replicate(lengthMax);
			VECTOR clampMin = Vector::Replicate(lengthMin);

			return ClampLengthV(v, clampMin, clampMax);
		}

		FORCE_INLINE VECTOR VEC_CALLCONV ClampLengthV(A_VECTOR v, A_VECTOR lengthMin, A_VECTOR lengthMax) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert((Vector::GetY(lengthMin) == Vector::GetX(lengthMin)));
			assert((Vector::GETY(lengthMax) == Vector::GetX(lengthMax)));
			assert(GreaterOrEqual(lengthMin, g_Zero));
			assert(GreaterOrEqual(lengthMax, g_Zero));
			assert(GreaterOrEqual(lengthMax, lengthMin));
#endif

			VECTOR lengthSq = LengthSq(v);

			const VECTOR zero = Vector::Zero();

			VECTOR rcpLength = Vector::ReciprocalSqrt(lengthSq);

			VECTOR infiniteLength = Vector::EqualInt(lengthSq, g_Infinity.v);
			VECTOR zeroLength = Vector::Equal(lengthSq, zero);

			VECTOR length = Vector::Multiply(lengthSq, rcpLength);

			VECTOR normal = Vector::Multiply(v, rcpLength);

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
			// Result = incident - (2 * dot(incident, normal)) * normal
			VECTOR Result;
			Result = Dot(incident, normal);
			Result = Vector::Add(Result, Result);
			Result = Vector::NegativeMultiplySubtract(Result, normal, incident);

			return Result;
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Refract(A_VECTOR incident, A_VECTOR normal, float refractionIndex) noexcept
		{
			VECTOR index = Vector::Replicate(refractionIndex);

			return RefractV(incident, normal, index);
		}

		FORCE_INLINE VECTOR VEC_CALLCONV RefractV(A_VECTOR incident, A_VECTOR normal, A_VECTOR refractionIndex) noexcept
		{
			// Result = refractionIndex * incident - normal * (refractionIndex * dot(incident, normal) + 
			// sqrt(1 - refractionIndex * refractionIndex * (1 - dot(incident, normal) * dot(incident, normal))))
#if defined(_NO_INTRINSICS_)
			float IDotN = (incident.vector4_f32[0] * normal.vector4_f32[0]) + (incident.vector4_f32[1] * normal.vector4_f32[1]);
			// R = 1.0f - refractionIndex * refractionIndex * (1.0f - IDotN * IDotN)
			float RY = 1.0f - (IDotN * IDotN);
			float RX = 1.0f - (RY * refractionIndex.vector4_f32[0] * refractionIndex.vector4_f32[0]);
			RY = 1.0f - (RY * refractionIndex.vector4_f32[1] * refractionIndex.vector4_f32[1]);
			if (RX >= 0.0f)
			{
				RX = (refractionIndex.vector4_f32[0] * incident.vector4_f32[0]) - (normal.vector4_f32[0] * ((refractionIndex.vector4_f32[0] * IDotN) + sqrtf(RX)));
			}
			else
			{
				RX = 0.0f;
			}
			if (RY >= 0.0f)
			{
				RY = (refractionIndex.vector4_f32[1] * incident.vector4_f32[1]) - (normal.vector4_f32[1] * ((refractionIndex.vector4_f32[1] * IDotN) + sqrtf(RY)));
			}
			else
			{
				RY = 0.0f;
			}

			VECTOR vResult;
			vResult.vector4_f32[0] = RX;
			vResult.vector4_f32[1] = RY;
			vResult.vector4_f32[2] = 0.0f;
			vResult.vector4_f32[3] = 0.0f;
			
			return vResult;

#elif defined(_SSE2_INTRINSICS_)
			// Result = refractionIndex * incident - normal * (refractionIndex * dot(incident, normal) +
    		// sqrt(1 - refractionIndex * refractionIndex * (1 - dot(incident, normal) * dot(incident, normal))))
    		// Get the 2D Dot product of incident-normal
			VECTOR IDotN = Dot(incident, normal);

			// vTemp 1.0f - refractionIndex * refractionIndex * (1.0f - IDotN * IDotN)
			VECTOR vTemp = FNMADD_PS(IDotN, IDotN, g_One);
			vTemp = _mm_mul_ps(vTemp, refractionIndex);
			vTemp = FNMADD_PS(vTemp, refractionIndex, g_One);

			// If any terms are <= 0, sqrt() will fail, punt to zero
			VECTOR vMask = _mm_cmpgt_ps(vTemp, g_Zero);

			// R = refractionIndex * IDotN + sqrt(R)
			vTemp = _mm_sqrt_ps(vTemp);
			vTemp = FMADD_PS(refractionIndex, IDotN, vTemp);

			// Result = refractionIndex * incident - normal * R
			VECTOR vResult = _mm_mul_ps(refractionIndex, incident);
			vResult = FNMADD_PS(vTemp, normal, vResult);
			
			return _mm_and_ps(vResult, vMask);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Orthogonal(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 Result = { { {
				-v.vector4_f32[1], 
				v.vector4_f32[0], 
				0.0f, 
				0.0f
			} } };

			return Result.v;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vResult = PERMUTE_PS(v, _MM_SHUFFLE(3, 2, 0, 1));
			
			return _mm_mul_ps(vResult, g_NegateX);
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

			VECTOR Dot = Dot(V1, V2);

			L1 = Vector::Multiply(L1, L2);

			VECTOR CosAngle = Vector::Multiply(Dot, L1);
			CosAngle = Vector::Clamp(CosAngle, g_NegativeOne.v, g_One.v);

			return Vector::ACos(CosAngle);
		}

		FORCE_INLINE VECTOR VEC_CALLCONV LinePointDistance(A_VECTOR linePoint1, A_VECTOR linePoint2, A_VECTOR point) noexcept
		{
			// Given a vector PointVector from LinePoint1 to Point and a vector
			// LineVector from LinePoint1 to LinePoint2, the scaled distance
			// PointProjectionScale from LinePoint1 to the perpendicular projection
			// of PointVector onto the line is defined as:
			//
			//     PointProjectionScale = dot(PointVector, LineVector) / LengthSq(LineVector)

			VECTOR pointVector = Vector::Subtract(point, linePoint1);
			VECTOR lineVector = Vector::Subtract(linePoint2, linePoint1);

			VECTOR lengthSq = LengthSq(lineVector);

			VECTOR pointProjectionScale = Dot(pointVector, lineVector);
			pointProjectionScale = Vector::Divide(pointProjectionScale, lengthSq);

			VECTOR distanceVector = Vector::Multiply(lineVector, pointProjectionScale);
			distanceVector = Vector::Subtract(pointVector, distanceVector);

			return Length(distanceVector);
		}

		FORCE_INLINE VECTOR VEC_CALLCONV IntersectLine(A_VECTOR line1Point1, A_VECTOR line1Point2, A_VECTOR line2Point1, B_VECTOR line2Point2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR V1 = Vector::Subtract(line1Point2, line1Point1);
			VECTOR V2 = Vector::Subtract(line2Point2, line2Point1);
			VECTOR V3 = Vector::Subtract(line1Point1, line2Point1);

			VECTOR C1 = Cross(V1, V2);
			VECTOR C2 = Cross(V2, V3);

			VECTOR Result;
			const VECTOR zero = Vector::Zero();
			if (NearEqual(C1, zero, g_Epsilon.v))
			{
				if (NearEqual(C2, zero, g_Epsilon.v))
				{
					// Coincident
					Result = g_Infinity.v;
				}
				else
				{
					// Parallel
					Result = g_GNaN.v;
				}
			}
			else
			{
				// Intersection point = line1Point1 + V1 * (C2 / C1)
				VECTOR scale = Vector::Reciprocal(C1);
				scale = Vector::Multiply(C2, scale);
				Result = Vector::MultiplyAdd(V1, scale, line1Point1);
			}

			return Result;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR V1 = _mm_sub_ps(line1Point2, line1Point1);
			VECTOR V2 = _mm_sub_ps(line2Point2, line2Point1);
			VECTOR V3 = _mm_sub_ps(line1Point1, line2Point1);

			// Generate the cross product
			VECTOR C1 = Cross(V1, V2);
			VECTOR C2 = Cross(V2, V3);

			// If C1 is not close to epsilon, use the calculated value
			VECTOR vResultMask = _mm_setzero_ps();
			vResultMask = _mm_sub_ps(vResultMask, C1);
			vResultMask = _mm_max_ps(vResultMask, C1);

			// 0xFFFFFFFF is the calculated value is to be used
			vResultMask = _mm_cmpgt_ps(vResultMask, g_Epsilon);

			// If C1 is close to epsilon, which fail type is? Infinity or NaN
			VECTOR vFailMask = _mm_setzero_ps();
			vFailMask = _mm_sub_ps(vFailMask, C2);
			vFailMask = _mm_max_ps(vFailMask, C2);
			vFailMask = _mm_cmple_ps(vFailMask, g_Epsilon);

			VECTOR vFail = _mm_and_ps(vFailMask, g_Infinity);
			vFailMask = _mm_andnot_ps(vFailMask, g_QNaN);

			// vFail is NaN or Infinite
			vFail = _mm_or_ps(vFail, vFailMask);

			// Intersection point = line1Point1 + V1 * (C2 / C1)
			VECTOR vResult = _mm_div_ps(C2, C1);
			vResult = FMADD_PS(vResult, V1, line1Point1);

			// Use result, or failure value
			vResult = _mm_and_ps(vResult, vResultMask);
			vResultMask = _mm_andnot_ps(vResultMask, vFail);

			return _mm_or_ps(vResult, vResultMask);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Transform(A_VECTOR v, A_VECTOR m) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR y = Vector::SplatY(v);
			VECTOR x = Vector::SplatX(v);

			VECTOR Result = Vector::MultiplyAdd(y, m.r[1], m.r[3]);
			Result = Vector::MultiplyAdd(x, m.r[0], Result);

			return Result;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vResult = PERMUTE_PS(v, _MM_SHUFFLE(1, 1, 1, 1)); // y
			vResult = FMADD_PS(vResult, m.r[1], m.r[3]);
			VECTOR vTemp = PERMUTE_PS(v, _MM_SHUFFLE(0, 0, 0, 0)); // x
			vResult = FMADD_PS(vTemp, m.r[0], vResult);

			return vResult;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE Float4* VEC_CALLCONV TransformStream(Float4* pOutputStream, size_t outputStride, const Float2* pInputStream, size_t inputStride, size_t vectorCount, A_MATRIX m) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pOutputStream != nullptr);
			assert(pInputStream != nullptr);

			assert(inputStride >= sizeof(Float2));
			_Analysis_assume_(inputStride <= sizeof(Float2));

			assert(outputStride >= sizeof(Float4));
			_Analysis_assume_(outputStride >= sizeof(Float4));
#endif

//#if defined(_NO_INTRINSICS_)
			const uint8_t* pInputVector = reinterpret_cast<const uint8_t*>(pInputStream);
			uint8_t* pOutputVector = reinterpret_cast<uint8_t*>(pOutputStream);

			const VECTOR row0 = m.r[0];
			const VECTOR row1 = m.r[1];
			const VECTOR row3 = m.r[3];

			for (size_t i = 0; i < vectorCountl i++)
			{
				VECTOR v = Vector::LoadFloat2(reinterpret_cast<const Float2*>(pInputVector));
				VECTOR y = Vector::SplatY(v);
				VECTOR x = Vector::SplatX(v);

				VECTOR Result = Vector::MultiplyAdd(y, row1, row3);
				Result = Vector::MultiplyAdd(x, row0, Result);

#ifdef _PREFAST_
#pragma prefast(push)
#pragma prefast(disable : 26015, "PREfast noise: Esp:1307" )
#endif

				Vector::StoreFloat4(reinterpret_cast<XMFLOAT4*>(pOutputVector), Result);

#ifdef _PREFAST_
#pragma prefast(pop)
#endif

				pInputVector += inputStride;
				pOutputVector += outputStride;
			}

			return pOutputStream;

//#elif defined(_AVX2_INTRINSICS_)
			const uint8_t* pInputVector = reinterpret_cast<const uint8_t*>(pInputStream);
			uint8_t* pOutputVector = reinterpret_cast<uint8_t*>(pOutputStream);

			size_t i = 0;
			size_t four = vectorCount >> 2;
			if (four > 0)
			{
				__m256 row0 = _mm256_broadcast_ps(&m.r[0]);
				__m256 row1 = _mm256_broadcast_ps(&m.r[1]);
				__m256 row3 = _mm256_broadcast_ps(&m.r[3]);

				if (inputStride == sizeof(Float2))
				{
					if (outputStride == sizeof(Float4))
					{
						if (!(reinterpret_cast<uintptr_t>(pOutputStream) & 0x1F))
						{
							// Packed input, aligned & packed output
							for (size_t j = 0; j < four; ++j)
							{
								__m256 VV = _mm256_loadu_ps(reinterpret_cast<const float*>(pInputVector));
								pInputVector += sizeof(Float2) * 4;

								__m256 Y2 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(3, 3, 3, 3));
								__m256 X2 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(2, 2, 2, 2));
								__m256 Y1 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(1, 1, 1, 1));
								__m256 X1 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(0, 0, 0, 0));

								__m256 vTempB = _mm256_fmadd_ps(Y1, row1, row3);
								__m256 vTempB2 = _mm256_fmadd_ps(Y2, row1, row3);
								__m256 vTempA = _mm256_mul_ps(X1, row0);
								__m256 vTempA2 = _mm256_mul_ps(X2, row0);
								vTempA = _mm256_add_ps(vTempA, vTempB);
								vTempA2 = _mm256_add_ps(vTempA2, vTempB2);

								X1 = _mm256_insertf128_ps(vTempA, _mm256_castps256_ps128(vTempA2), 1);
								_256_STREAM_PS(reinterpret_cast<float*>(pOutputVector), X1);
								pOutputVector += sizeof(Float4) * 2;

								X2 = _mm256_insertf128_ps(vTempA2, _mm256_extractf128_ps(vTempA, 1), 0);
								_256_STREAM_PS(reinterpret_cast<float*>(pOutputVector), X2);
								pOutputVector += sizeof(Float4) * 2;

								i += 4;
							}
						}
						else
						{
							// Packed input, packed output
							for (size_t j = 0; j < four; ++j)
							{
								__m256 VV = _mm256_loadu_ps(reinterpret_cast<const float*>(pInputVector));
								pInputVector += sizeof(Float2) * 4;

								__m256 Y2 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(3, 3, 3, 3));
								__m256 X2 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(2, 2, 2, 2));
								__m256 Y1 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(1, 1, 1, 1));
								__m256 X1 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(0, 0, 0, 0));

								__m256 vTempB = _mm256_fmadd_ps(Y1, row1, row3);
								__m256 vTempB2 = _mm256_fmadd_ps(Y2, row1, row3);
								__m256 vTempA = _mm256_mul_ps(X1, row0);
								__m256 vTempA2 = _mm256_mul_ps(X2, row0);
								vTempA = _mm256_add_ps(vTempA, vTempB);
								vTempA2 = _mm256_add_ps(vTempA2, vTempB2);

								X1 = _mm256_insertf128_ps(vTempA, _mm256_castps256_ps128(vTempA2), 1);
								_mm256_storeu_ps(reinterpret_cast<float*>(pOutputVector), X1);
								pOutputVector += sizeof(Float4) * 2;

								X2 = _mm256_insertf128_ps(vTempA2, _mm256_extractf128_ps(vTempA, 1), 0);
								_mm256_storeu_ps(reinterpret_cast<float*>(pOutputVector), X2);
								pOutputVector += sizeof(Float4) * 2;

								i += 4;
							}
						}
					}
					else
					{
						// Packed input, unpacked output
						for (size_t j = 0; j < four; ++j)
						{
							__m256 VV = _mm256_loadu_ps(reinterpret_cast<const float*>(pInputVector));
							pInputVector += sizeof(Float2) * 4;

							__m256 Y2 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(3, 3, 3, 3));
							__m256 X2 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(2, 2, 2, 2));
							__m256 Y1 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(1, 1, 1, 1));
							__m256 X1 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(0, 0, 0, 0));

							__m256 vTempB = _mm256_fmadd_ps(Y1, row1, row3);
							__m256 vTempB2 = _mm256_fmadd_ps(Y2, row1, row3);
							__m256 vTempA = _mm256_mul_ps(X1, row0);
							__m256 vTempA2 = _mm256_mul_ps(X2, row0);
							vTempA = _mm256_add_ps(vTempA, vTempB);
							vTempA2 = _mm256_add_ps(vTempA2, vTempB2);

							_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), _mm256_castps256_ps128(vTempA));
							pOutputVector += outputStride;

							_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), _mm256_castps256_ps128(vTempA2));
							pOutputVector += outputStride;

							_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), _mm256_extractf128_ps(vTempA, 1));
							pOutputVector += outputStride;

							_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), _mm256_extractf128_ps(vTempA2, 1));
							pOutputVector += outputStride;

							i += 4;
						}
					}
				}
			}

			if (i < vectorCount)
			{
				const VECTOR row0 = m.r[0];
				const VECTOR row1 = m.r[1];
				const VECTOR row3 = m.r[3];

				for (; i < vectorCount; i++)
				{
					__m128 xy = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(pInputVector)));
					pInputVector += inputStride;

					VECTOR Y = PERMUTE_PS(xy, _MM_SHUFFLE(1, 1, 1, 1));
					VECTOR X = PERMUTE_PS(xy, _MM_SHUFFLE(0, 0, 0, 0));

					VECTOR vTemp = XM_FMADD_PS(Y, row1, row3);
					VECTOR vTemp2 = _mm_mul_ps(X, row0);
					vTemp = _mm_add_ps(vTemp, vTemp2);

					_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), vTemp);
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
			const VECTOR row3 = m.r[3];

			size_t i = 0;
			size_t two = vectorCount >> 1;
			if (two > 0)
			{
				if (inputStride == sizeof(Float2))
				{
					if (!(reinterpret_cast<uintptr_t>(pOutputStream) & 0xF) && !(outputStride & 0xF))
					{
						// Packed input, aligned output
						for (size_t j = 0; j < two; ++j)
						{
							VECTOR V = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
							pInputVector += sizeof(Float2) * 2;

							VECTOR Y = PERMUTE_PS(V, _MM_SHUFFLE(1, 1, 1, 1));
							VECTOR X = PERMUTE_PS(V, _MM_SHUFFLE(0, 0, 0, 0));

							VECTOR vTemp = FMADD_PS(Y, row1, row3);
							VECTOR vTemp2 = _mm_mul_ps(X, row0);
							vTemp = _mm_add_ps(vTemp, vTemp2);

							STREAM_PS(reinterpret_cast<float*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							Y = PERMUTE_PS(V, _MM_SHUFFLE(3, 3, 3, 3));
							X = PERMUTE_PS(V, _MM_SHUFFLE(2, 2, 2, 2));

							vTemp = FMADD_PS(Y, row1, row3);
							vTemp2 = _mm_mul_ps(X, row0);
							vTemp = _mm_add_ps(vTemp, vTemp2);

							STREAM_PS(reinterpret_cast<float*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							i += 2;
						}
					}
					else
					{
						// Packed input, unaligned output
						for (size_t j = 0; j < two; ++j)
						{
							VECTOR V = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
							pInputVector += sizeof(Float2) * 2;

							VECTOR Y = PERMUTE_PS(V, _MM_SHUFFLE(1, 1, 1, 1));
							VECTOR X = PERMUTE_PS(V, _MM_SHUFFLE(0, 0, 0, 0));

							VECTOR vTemp = FMADD_PS(Y, row1, row3);
							VECTOR vTemp2 = _mm_mul_ps(X, row0);
							vTemp = _mm_add_ps(vTemp, vTemp2);

							_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							Y = PERMUTE_PS(V, _MM_SHUFFLE(3, 3, 3, 3));
							X = PERMUTE_PS(V, _MM_SHUFFLE(2, 2, 2, 2));

							vTemp = FMADD_PS(Y, row1, row3);
							vTemp2 = _mm_mul_ps(X, row0);
							vTemp = _mm_add_ps(vTemp, vTemp2);

							_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), vTemp);
							pOutputVector += outputStride;

							i += 2;
						}
					}
				}
			}

			if (!(reinterpret_cast<uintptr_t>(pInputVector) & 0xF) && !(InputStride & 0xF))
			{
				if (!(reinterpret_cast<uintptr_t>(pOutputStream) & 0xF) && !(outputStride & 0xF))
				{
					// Aligned input, aligned output
					for (; i < VectorCount; i++)
					{
						VECTOR V = _mm_castsi128_ps(_mm_loadl_epi64(reinterpret_cast<const __m128i*>(pInputVector)));
						pInputVector += inputStride;

						VECTOR Y = PERMUTE_PS(V, _MM_SHUFFLE(1, 1, 1, 1));
						VECTOR X = PERMUTE_PS(V, _MM_SHUFFLE(0, 0, 0, 0));

						VECTOR vTemp = FMADD_PS(Y, row1, row3);
						VECTOR vTemp2 = _mm_mul_ps(X, row0);
						vTemp = _mm_add_ps(vTemp, vTemp2);

						STREAM_PS(reinterpret_cast<float*>(pOutputVector), vTemp);
						pOutputVector += outputStride;
					}
				}
				else
				{
					// Aligned input, unaligned output
					for (; i < vectorCount; i++)
					{
						VECTOR V = _mm_castsi128_ps(_mm_loadl_epi64(reinterpret_cast<const __m128i*>(pInputVector)));
						pInputVector += inputStride;

						VECTOR Y = PERMUTE_PS(V, _MM_SHUFFLE(1, 1, 1, 1));
						VECTOR X = PERMUTE_PS(V, _MM_SHUFFLE(0, 0, 0, 0));

						VECTOR vTemp = FMADD_PS(Y, row1, row3);
						VECTOR vTemp2 = _mm_mul_ps(X, row0);
						vTemp = _mm_add_ps(vTemp, vTemp2);

						_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), vTemp);
						pOutputVector += outputStride;
					}
				}
			}
			else
			{
				// Unaligned input
				for (; i < vectorCount; i++)
				{
					__m128 xy = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(pInputVector)));
					pInputVector += inputStride;

					VECTOR Y = PERMUTE_PS(xy, _MM_SHUFFLE(1, 1, 1, 1));
					VECTOR X = PERMUTE_PS(xy, _MM_SHUFFLE(0, 0, 0, 0));

					VECTOR vTemp = FMADD_PS(Y, row1, row3);
					VECTOR vTemp2 = _mm_mul_ps(X, row0);
					vTemp = _mm_add_ps(vTemp, vTemp2);

					_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), vTemp);
					pOutputVector += outputStride;
				}
			}

			SFENCE();

			return pOutputStream;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV TransformCoord(A_VECTOR v, A_VECTOR m) noexcept
		{
			VECTOR y = Vector::SplatY(v);
			VECTOR x = Vector::SplatX(v);

			VECTOR Result = Vector::MultiplyAdd(y, m.r[1], m.r[3]);
			Result = Vector::MultiplyAdd(x, m.r[0], Result);

			VECTOR w = Vector::SplatW(Result);

			return Vector::Divide(Result, w);
		}

		_Use_decl_annotations_
		FORCE_INLINE Float2* VEC_CALLCONV TransformCoordStream(Float2* pOutputStream, size_t outputStride, const Float2* pInputStream, size_t inputStride, size_t vectorCount, A_MATRIX m) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pOutputStream != nullptr);
			assert(pInputStream != nullptr);

			assert(inputStride >= sizeof(Float2));
			_Analysis_assume_(inputStride >= sizeof(Float2));

			assert(outputStride >= sizeof(Float2));
			_Analysis_assume_(outputStride >= sizeof(Float2));
#endif

//#if defined(_NO_INTRINSICS_)
			const uint8_t* pInputVector = reinterpret_cast<const uint8_t*>(pInputStream);
			uint8_t* pOutputVector = reinterpret_cast<uint8_t*>(pOutputStream);

			const VECTOR row0 = m.r[0];
			const VECTOR row1 = m.r[1];
			const VECTOR row3 = m.r[3];

			for (size_t i = 0; i < vectorCount; i++)
			{
				VECTOR V = Vector::LoadFloat2(reinterpret_cast<const Float2*>(pInputVector));
				VECTOR Y = Vector::SplatY(V);
				VECTOR X = Vector::SplatX(V);

				VECTOR Result = Vector::MultiplyAdd(Y, row1, row3);
				Result = Vector::MultiplyAdd(X, row0, Result);

				VECTOR W = Vector::SplatW(Result);

				Result = Vector::Divide(Result, W);

#ifdef _PREFAST_
#pragma prefast(push)
#pragma prefast(disable : 26015, "PREfast noise: Esp:1307" )
#endif

				Vector::StoreFloat2(reinterpret_cast<Float2*>(pOutputVector), Result);

#ifdef _PREFAST_
#pragma prefast(pop)
#endif

				pInputVector += inputStride;
				pOutputVector += outputStride;
			}

			return pOutputStream;

//#elif defined(_AVX2_INTRINSICS_)
			const uint8_t* pInputVector = reinterpret_cast<const uint8_t*>(pInputStream);
			uint8_t* pOutputVector = reinterpret_cast<uint8_t*>(pOutputStream);

			size_t i = 0;
			size_t four = vectorCount >> 2;
			if (four > 0)
			{
				__m256 row0 = _mm256_broadcast_ps(&m.r[0]);
				__m256 row1 = _mm256_broadcast_ps(&m.r[1]);
				__m256 row3 = _mm256_broadcast_ps(&m.r[3]);

				if (inputStride == sizeof(Float2))
				{
					if (outputStride == sizeof(Float2))
					{
						if (!(reinterpret_cast<uintptr_t>(pOutputStream) & 0x1F))
						{
							// Packed input, aligned & packed output
							for (size_t j = 0; j < four; ++j)
							{
								__m256 VV = _mm256_loadu_ps(reinterpret_cast<const float*>(pInputVector));
								pInputVector += sizeof(Float2) * 4;

								__m256 Y2 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(3, 3, 3, 3));
								__m256 X2 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(2, 2, 2, 2));
								__m256 Y1 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(1, 1, 1, 1));
								__m256 X1 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(0, 0, 0, 0));

								__m256 vTempB = _mm256_fmadd_ps(Y1, row1, row3);
								__m256 vTempB2 = _mm256_fmadd_ps(Y2, row1, row3);
								__m256 vTempA = _mm256_mul_ps(X1, row0);
								__m256 vTempA2 = _mm256_mul_ps(X2, row0);
								vTempA = _mm256_add_ps(vTempA, vTempB);
								vTempA2 = _mm256_add_ps(vTempA2, vTempB2);

								__m256 W = _mm256_shuffle_ps(vTempA, vTempA, _MM_SHUFFLE(3, 3, 3, 3));
								vTempA = _mm256_div_ps(vTempA, W);

								W = _mm256_shuffle_ps(vTempA2, vTempA2, _MM_SHUFFLE(3, 3, 3, 3));
								vTempA2 = _mm256_div_ps(vTempA2, W);

								X1 = _mm256_shuffle_ps(vTempA, vTempA2, 0x44);
								_256_STREAM_PS(reinterpret_cast<float*>(pOutputVector), X1);
								pOutputVector += sizeof(Float2) * 4;

								i += 4;
							}
						}
						else
						{
							// Packed input, packed output
							for (size_t j = 0; j < four; ++j)
							{
								__m256 VV = _mm256_loadu_ps(reinterpret_cast<const float*>(pInputVector));
								pInputVector += sizeof(Float2) * 4;

								__m256 Y2 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(3, 3, 3, 3));
								__m256 X2 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(2, 2, 2, 2));
								__m256 Y1 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(1, 1, 1, 1));
								__m256 X1 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(0, 0, 0, 0));

								__m256 vTempB = _mm256_fmadd_ps(Y1, row1, row3);
								__m256 vTempB2 = _mm256_fmadd_ps(Y2, row1, row3);
								__m256 vTempA = _mm256_mul_ps(X1, row0);
								__m256 vTempA2 = _mm256_mul_ps(X2, row0);
								vTempA = _mm256_add_ps(vTempA, vTempB);
								vTempA2 = _mm256_add_ps(vTempA2, vTempB2);

								__m256 W = _mm256_shuffle_ps(vTempA, vTempA, _MM_SHUFFLE(3, 3, 3, 3));
								vTempA = _mm256_div_ps(vTempA, W);

								W = _mm256_shuffle_ps(vTempA2, vTempA2, _MM_SHUFFLE(3, 3, 3, 3));
								vTempA2 = _mm256_div_ps(vTempA2, W);

								X1 = _mm256_shuffle_ps(vTempA, vTempA2, 0x44);
								_mm256_storeu_ps(reinterpret_cast<float*>(pOutputVector), X1);
								pOutputVector += sizeof(Float2) * 4;

								i += 4;
							}
						}
					}
					else
					{
						// Packed input, unpacked output
						for (size_t j = 0; j < four; ++j)
						{
							__m256 VV = _mm256_loadu_ps(reinterpret_cast<const float*>(pInputVector));
							pInputVector += sizeof(Float2) * 4;

							__m256 Y2 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(3, 3, 3, 3));
							__m256 X2 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(2, 2, 2, 2));
							__m256 Y1 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(1, 1, 1, 1));
							__m256 X1 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(0, 0, 0, 0));

							__m256 vTempB = _mm256_fmadd_ps(Y1, row1, row3);
							__m256 vTempB2 = _mm256_fmadd_ps(Y2, row1, row3);
							__m256 vTempA = _mm256_mul_ps(X1, row0);
							__m256 vTempA2 = _mm256_mul_ps(X2, row0);
							vTempA = _mm256_add_ps(vTempA, vTempB);
							vTempA2 = _mm256_add_ps(vTempA2, vTempB2);

							__m256 W = _mm256_shuffle_ps(vTempA, vTempA, _MM_SHUFFLE(3, 3, 3, 3));
							vTempA = _mm256_div_ps(vTempA, W);

							W = _mm256_shuffle_ps(vTempA2, vTempA2, _MM_SHUFFLE(3, 3, 3, 3));
							vTempA2 = _mm256_div_ps(vTempA2, W);

							_mm_store_sd(reinterpret_cast<double*>(pOutputVector),
								_mm_castps_pd(_mm256_castps256_ps128(vTempA)));
							pOutputVector += outputStride;

							_mm_store_sd(reinterpret_cast<double*>(pOutputVector),
								_mm_castps_pd(_mm256_castps256_ps128(vTempA2)));
							pOutputVector += outputStride;

							_mm_store_sd(reinterpret_cast<double*>(pOutputVector),
								_mm_castps_pd(_mm256_extractf128_ps(vTempA, 1)));
							pOutputVector += outputStride;

							_mm_store_sd(reinterpret_cast<double*>(pOutputVector),
								_mm_castps_pd(_mm256_extractf128_ps(vTempA2, 1)));
							pOutputVector += outputStride;

							i += 4;
						}
					}
				}
			}

			if (i < vectorCount)
			{
				const VECTOR row0 = m.r[0];
				const VECTOR row1 = m.r[1];
				const VECTOR row3 = m.r[3];

				for (; i < vectorCount; i++)
				{
					__m128 xy = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(pInputVector)));
					pInputVector += inputStride;

					VECTOR Y = PERMUTE_PS(xy, _MM_SHUFFLE(1, 1, 1, 1));
					VECTOR X = PERMUTE_PS(xy, _MM_SHUFFLE(0, 0, 0, 0));

					VECTOR vTemp = FMADD_PS(Y, row1, row3);
					VECTOR vTemp2 = _mm_mul_ps(X, row0);
					vTemp = _mm_add_ps(vTemp, vTemp2);

					VECTOR W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));
					vTemp = _mm_div_ps(vTemp, W);

					_mm_store_sd(reinterpret_cast<double*>(pOutputVector), _mm_castps_pd(vTemp));
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
			const VECTOR row3 = m.r[3];

			size_t i = 0;
			size_t two = vectorCount >> 1;
			if (two > 0)
			{
				if (inputStride == sizeof(Float2))
				{
					if (outputStride == sizeof(Float2))
					{
						if (!(reinterpret_cast<uintptr_t>(pOutputStream) & 0xF))
						{
							// Packed input, aligned & packed output
							for (size_t j = 0; j < two; ++j)
							{
								VECTOR V = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
								pInputVector += sizeof(Float2) * 2;

								// Result 1
								VECTOR Y = PERMUTE_PS(V, _MM_SHUFFLE(1, 1, 1, 1));
								VECTOR X = PERMUTE_PS(V, _MM_SHUFFLE(0, 0, 0, 0));

								VECTOR vTemp = XM_FMADD_PS(Y, row1, row3);
								VECTOR vTemp2 = _mm_mul_ps(X, row0);
								vTemp = _mm_add_ps(vTemp, vTemp2);

								VECTOR W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));

								VECTOR V1 = _mm_div_ps(vTemp, W);

								// Result 2
								Y = PERMUTE_PS(V, _MM_SHUFFLE(3, 3, 3, 3));
								X = PERMUTE_PS(V, _MM_SHUFFLE(2, 2, 2, 2));

								vTemp = FMADD_PS(Y, row1, row3);
								vTemp2 = _mm_mul_ps(X, row0);
								vTemp = _mm_add_ps(vTemp, vTemp2);

								W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));

								VECTOR V2 = _mm_div_ps(vTemp, W);

								vTemp = _mm_movelh_ps(V1, V2);

								STREAM_PS(reinterpret_cast<float*>(pOutputVector), vTemp);
								pOutputVector += sizeof(Float2) * 2;

								i += 2;
							}
						}
						else
						{
							// Packed input, unaligned & packed output
							for (size_t j = 0; j < two; ++j)
							{
								VECTOR V = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
								pInputVector += sizeof(Float2) * 2;

								// Result 1
								VECTOR Y = PERMUTE_PS(V, _MM_SHUFFLE(1, 1, 1, 1));
								VECTOR X = PERMUTE_PS(V, _MM_SHUFFLE(0, 0, 0, 0));

								VECTOR vTemp = FMADD_PS(Y, row1, row3);
								VECTOR vTemp2 = _mm_mul_ps(X, row0);
								vTemp = _mm_add_ps(vTemp, vTemp2);

								VECTOR W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));

								VECTOR V1 = _mm_div_ps(vTemp, W);

								// Result 2
								Y = PERMUTE_PS(V, _MM_SHUFFLE(3, 3, 3, 3));
								X = PERMUTE_PS(V, _MM_SHUFFLE(2, 2, 2, 2));

								vTemp = FMADD_PS(Y, row1, row3);
								vTemp2 = _mm_mul_ps(X, row0);
								vTemp = _mm_add_ps(vTemp, vTemp2);

								W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));

								VECTOR V2 = _mm_div_ps(vTemp, W);

								vTemp = _mm_movelh_ps(V1, V2);

								_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), vTemp);
								pOutputVector += sizeof(Float2) * 2;

								i += 2;
							}
						}
					}
					else
					{
						// Packed input, unpacked output
						for (size_t j = 0; j < two; ++j)
						{
							VECTOR V = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
							pInputVector += sizeof(Float2) * 2;

							// Result 1
							VECTOR Y = PERMUTE_PS(V, _MM_SHUFFLE(1, 1, 1, 1));
							VECTOR X = PERMUTE_PS(V, _MM_SHUFFLE(0, 0, 0, 0));

							VECTOR vTemp = FMADD_PS(Y, row1, row3);
							VECTOR vTemp2 = _mm_mul_ps(X, row0);
							vTemp = _mm_add_ps(vTemp, vTemp2);

							VECTOR W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));

							vTemp = _mm_div_ps(vTemp, W);

							_mm_store_sd(reinterpret_cast<double*>(pOutputVector), _mm_castps_pd(vTemp));
							pOutputVector += outputStride;

							// Result 2
							Y = PERMUTE_PS(V, _MM_SHUFFLE(3, 3, 3, 3));
							X = PERMUTE_PS(V, _MM_SHUFFLE(2, 2, 2, 2));

							vTemp = FMADD_PS(Y, row1, row3);
							vTemp2 = _mm_mul_ps(X, row0);
							vTemp = _mm_add_ps(vTemp, vTemp2);

							W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));

							vTemp = _mm_div_ps(vTemp, W);

							_mm_store_sd(reinterpret_cast<double*>(pOutputVector), _mm_castps_pd(vTemp));
							pOutputVector += outputStride;

							i += 2;
						}
					}
				}
			}

			if (!(reinterpret_cast<uintptr_t>(pInputVector) & 0xF) && !(inputStride & 0xF))
			{
				// Aligned input
				for (; i < vectorCount; i++)
				{
					VECTOR V = _mm_castsi128_ps(_mm_loadl_epi64(reinterpret_cast<const __m128i*>(pInputVector)));
					pInputVector += inputStride;

					VECTOR Y = PERMUTE_PS(V, _MM_SHUFFLE(1, 1, 1, 1));
					VECTOR X = PERMUTE_PS(V, _MM_SHUFFLE(0, 0, 0, 0));

					VECTOR vTemp = FMADD_PS(Y, row1, row3);
					VECTOR vTemp2 = _mm_mul_ps(X, row0);
					vTemp = _mm_add_ps(vTemp, vTemp2);

					VECTOR W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));

					vTemp = _mm_div_ps(vTemp, W);

					_mm_store_sd(reinterpret_cast<double*>(pOutputVector), _mm_castps_pd(vTemp));
					pOutputVector += outputStride;
				}
			}
			else
			{
				// Unaligned input
				for (; i < vectorCount; i++)
				{
					__m128 xy = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(pInputVector)));
					pInputVector += inputStride;

					VECTOR Y = PERMUTE_PS(xy, _MM_SHUFFLE(1, 1, 1, 1));
					VECTOR X = PERMUTE_PS(xy, _MM_SHUFFLE(0, 0, 0, 0));

					VECTOR vTemp = FMADD_PS(Y, row1, row3);
					VECTOR vTemp2 = _mm_mul_ps(X, row0);
					vTemp = _mm_add_ps(vTemp, vTemp2);

					VECTOR W = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 3, 3, 3));

					vTemp = _mm_div_ps(vTemp, W);

					_mm_store_sd(reinterpret_cast<double*>(pOutputVector), _mm_castps_pd(vTemp));
					pOutputVector += outputStride;
				}
			}

			SFENCE();

			return pOutputStream;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV TransformNormal(A_VECTOR v, A_MATRIX m) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR Y = Vector::SplatY(v);
			VECTOR X = Vector::SplatX(v);

			VECTOR Result = Vector::Multiply(Y, m.r[1]);
			Result = Vector::MultiplyAdd(X, m.r[0], Result);

			return Result;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vResult = PERMUTE_PS(v, _MM_SHUFFLE(1, 1, 1, 1)); // y
			vResult = _mm_mul_ps(vResult, m.r[1]);
			
			VECTOR vTemp = PERMUTE_PS(v, _MM_SHUFFLE(0, 0, 0, 0)); // x
			vResult = FMADD_PS(vTemp. m.r[0], vResult);
			
			return vResult;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE Float2* VEC_CALLCONV TransformNormalStream(Float2* pOutputStream, size_t outputStride, const Float2* pInputStream, size_t inputStride, size_t vectorCount, A_MATRIX m) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pOutputStream != nullptr);
			assert(pInputStream != nullptr);

			assert(inputStride >= sizeof(Float2));
			_Analysis_assume_(inputStride >= sizeof(Float2));

			assert(outputStride >= sizeof(Float2));
			_Analysis_assume_(outputStride >= sizeof(Float2));
#endif

//#if defined(_NO_INTRINSICS_)
			const uint8_t* pInputVector = reinterpret_cast<const uint8_t*>(pInputStream);
			uint8_t* pOutputVector = reinterpret_cast<uint8_t*>(pOutputStream);

			const VECTOR row0 = m.r[0];
			const VECTOR row1 = m.r[1];

			for (size_t i = 0; i < vectorCount; i++)
			{
				VECTOR V = Vector::Float2(reinterpret_cast<const Float2*>(pInputVector));
				VECTOR Y = Vector::SplatY(v);
				VECTOR X = Vector::SplatX(v);

				VECTOR Result = Vector::Multiply(Y, row1);
				Result = Vector::MultiplyAdd(X, row0, Result);

#ifdef _PREFAST_
#pragma prefast(push)
#pragma prefast(disable : 26015, "PREfast noise: Esp:1307" )
#endif

				Vector::StoreFloat2(reinterpret_cast<Float2*>(pOutputVector), Result);

#ifdef _PREFAST_
#pragma prefast(pop)
#endif

				pInputVector += InputStride;
				pOutputVector += OutputStride;
			}

			return pOutputStream;

#elif defined(_AVX2_INTRINSICS_)
			const uint8_t* pInputVector = reinterpret_cast<const uint8_t*>(pInputStream);
			uint8_t* pOutputVector = reinterpret_cast<uint8_t*>(pOutputStream);

			size_t i = 0;
			size_t four = vectorCount >> 2;
			if (four > 0)
			{
				__m256 row0 = _mm256_broadcast_ps(&m.r[0]);
				__m256 row1 = _mm256_broadcast_ps(&m.r[1]);

				if (inputStride == sizeof(Float2))
				{
					if (outputStride == sizeof(Float2))
					{
						if (!(reinterpret_cast<uintptr_t>(pOutputStream) & 0x1F))
						{
							// Packed input, aligned & packed output
							for (size_t j = 0; j < four; ++j)
							{
								__m256 VV = _mm256_loadu_ps(reinterpret_cast<const float*>(pInputVector));
								pInputVector += sizeof(Float2) * 4;

								__m256 Y2 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(3, 3, 3, 3));
								__m256 X2 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(2, 2, 2, 2));
								__m256 Y1 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(1, 1, 1, 1));
								__m256 X1 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(0, 0, 0, 0));

								__m256 vTempA = _mm256_mul_ps(Y1, row1);
								__m256 vTempB = _mm256_mul_ps(Y2, row1);
								vTempA = _mm256_fmadd_ps(X1, row0, vTempA);
								vTempB = _mm256_fmadd_ps(X2, row0, vTempB);

								X1 = _mm256_shuffle_ps(vTempA, vTempB, 0x44);
								_256_STREAM_PS(reinterpret_cast<float*>(pOutputVector), X1);
								pOutputVector += sizeof(Float2) * 4;

								i += 4;
							}
						}
						else
						{
							// Packed input, packed output
							for (size_t j = 0; j < four; ++j)
							{
								__m256 VV = _mm256_loadu_ps(reinterpret_cast<const float*>(pInputVector));
								pInputVector += sizeof(Float2) * 4;

								__m256 Y2 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(3, 3, 3, 3));
								__m256 X2 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(2, 2, 2, 2));
								__m256 Y1 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(1, 1, 1, 1));
								__m256 X1 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(0, 0, 0, 0));

								__m256 vTempA = _mm256_mul_ps(Y1, row1);
								__m256 vTempB = _mm256_mul_ps(Y2, row1);
								vTempA = _mm256_fmadd_ps(X1, row0, vTempA);
								vTempB = _mm256_fmadd_ps(X2, row0, vTempB);

								X1 = _mm256_shuffle_ps(vTempA, vTempB, 0x44);
								_mm256_storeu_ps(reinterpret_cast<float*>(pOutputVector), X1);
								pOutputVector += sizeof(Float2) * 4;

								i += 4;
							}
						}
					}
					else
					{
						// Packed input, unpacked output
						for (size_t j = 0; j < four; ++j)
						{
							__m256 VV = _mm256_loadu_ps(reinterpret_cast<const float*>(pInputVector));
							pInputVector += sizeof(Float2) * 4;

							__m256 Y2 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(3, 3, 3, 3));
							__m256 X2 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(2, 2, 2, 2));
							__m256 Y1 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(1, 1, 1, 1));
							__m256 X1 = _mm256_shuffle_ps(VV, VV, _MM_SHUFFLE(0, 0, 0, 0));

							__m256 vTempA = _mm256_mul_ps(Y1, row1);
							__m256 vTempB = _mm256_mul_ps(Y2, row1);
							vTempA = _mm256_fmadd_ps(X1, row0, vTempA);
							vTempB = _mm256_fmadd_ps(X2, row0, vTempB);

							_mm_store_sd(reinterpret_cast<double*>(pOutputVector),
								_mm_castps_pd(_mm256_castps256_ps128(vTempA)));
							pOutputVector += outputStride;

							_mm_store_sd(reinterpret_cast<double*>(pOutputVector),
								_mm_castps_pd(_mm256_castps256_ps128(vTempB)));
							pOutputVector += outputStride;

							_mm_store_sd(reinterpret_cast<double*>(pOutputVector),
								_mm_castps_pd(_mm256_extractf128_ps(vTempA, 1)));
							pOutputVector += outputStride;

							_mm_store_sd(reinterpret_cast<double*>(pOutputVector),
								_mm_castps_pd(_mm256_extractf128_ps(vTempB, 1)));
							pOutputVector += outputStride;

							i += 4;
						}
					}
				}
			}

			if (i < vectorCount)
			{
				const VECTOR row0 = m.r[0];
				const VECTOR row1 = m.r[1];

				for (; i < vectorCount; i++)
				{
					__m128 xy = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(pInputVector)));
					pInputVector += inputStride;

					VECTOR Y = PERMUTE_PS(xy, _MM_SHUFFLE(1, 1, 1, 1));
					VECTOR X = PERMUTE_PS(xy, _MM_SHUFFLE(0, 0, 0, 0));

					VECTOR vTemp = _mm_mul_ps(Y, row1);
					vTemp = FMADD_PS(X, row0, vTemp);

					_mm_store_sd(reinterpret_cast<double*>(pOutputVector), _mm_castps_pd(vTemp));
					pOutputVector += outputStride;
				}
			}

			SFENCE();

			return pOutputStream;
#elif defined(_SSE2_INTRINSICS_)
			auto pInputVector = reinterpret_cast<const uint8_t*>(pInputStream);
			auto pOutputVector = reinterpret_cast<uint8_t*>(pOutputStream);

			const VECTOR row0 = m.r[0];
			const VECTOR row1 = m.r[1];

			size_t i = 0;
			size_t two = vectorCount >> 1;
			if (two > 0)
			{
				if (inputStride == sizeof(Float2))
				{
					if (outputStride == sizeof(Float2))
					{
						if (!(reinterpret_cast<uintptr_t>(pOutputStream) & 0xF))
						{
							// Packed input, aligned & packed output
							for (size_t j = 0; j < two; ++j)
							{
								VECTOR V = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
								pInputVector += sizeof(Float2) * 2;

								// Result 1
								VECTOR Y = PERMUTE_PS(V, _MM_SHUFFLE(1, 1, 1, 1));
								VECTOR X = PERMUTE_PS(V, _MM_SHUFFLE(0, 0, 0, 0));

								VECTOR vTemp = _mm_mul_ps(Y, row1);
								VECTOR V1 = XM_FMADD_PS(X, row0, vTemp);

								// Result 2
								Y = PERMUTE_PS(V, _MM_SHUFFLE(3, 3, 3, 3));
								X = PERMUTE_PS(V, _MM_SHUFFLE(2, 2, 2, 2));

								vTemp = _mm_mul_ps(Y, row1);
								VECTOR V2 = FMADD_PS(X, row0, vTemp);

								vTemp = _mm_movelh_ps(V1, V2);

								STREAM_PS(reinterpret_cast<float*>(pOutputVector), vTemp);
								pOutputVector += sizeof(Float2) * 2;

								i += 2;
							}
						}
						else
						{
							// Packed input, unaligned & packed output
							for (size_t j = 0; j < two; ++j)
							{
								VECTOR V = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
								pInputVector += sizeof(Float2) * 2;

								// Result 1
								VECTOR Y = PERMUTE_PS(V, _MM_SHUFFLE(1, 1, 1, 1));
								VECTOR X = PERMUTE_PS(V, _MM_SHUFFLE(0, 0, 0, 0));

								VECTOR vTemp = _mm_mul_ps(Y, row1);
								VECTOR V1 = FMADD_PS(X, row0, vTemp);

								// Result 2
								Y = PERMUTE_PS(V, _MM_SHUFFLE(3, 3, 3, 3));
								X = PERMUTE_PS(V, _MM_SHUFFLE(2, 2, 2, 2));

								vTemp = _mm_mul_ps(Y, row1);
								VECTOR V2 = FMADD_PS(X, row0, vTemp);

								vTemp = _mm_movelh_ps(V1, V2);

								_mm_storeu_ps(reinterpret_cast<float*>(pOutputVector), vTemp);
								pOutputVector += sizeof(Float2) * 2;

								i += 2;
							}
						}
					}
					else
					{
						// Packed input, unpacked output
						for (size_t j = 0; j < two; ++j)
						{
							VECTOR V = _mm_loadu_ps(reinterpret_cast<const float*>(pInputVector));
							pInputVector += sizeof(Float2) * 2;

							// Result 1
							VECTOR Y = PERMUTE_PS(V, _MM_SHUFFLE(1, 1, 1, 1));
							VECTOR X = PERMUTE_PS(V, _MM_SHUFFLE(0, 0, 0, 0));

							VECTOR vTemp = _mm_mul_ps(Y, row1);
							vTemp = XM_FMADD_PS(X, row0, vTemp);

							_mm_store_sd(reinterpret_cast<double*>(pOutputVector), _mm_castps_pd(vTemp));
							pOutputVector += outputStride;

							// Result 2
							Y = PERMUTE_PS(V, _MM_SHUFFLE(3, 3, 3, 3));
							X = PERMUTE_PS(V, _MM_SHUFFLE(2, 2, 2, 2));

							vTemp = _mm_mul_ps(Y, row1);
							vTemp = FMADD_PS(X, row0, vTemp);

							_mm_store_sd(reinterpret_cast<double*>(pOutputVector), _mm_castps_pd(vTemp));
							pOutputVector += outputStride;

							i += 2;
						}
					}
				}
			}

			if (!(reinterpret_cast<uintptr_t>(pInputVector) & 0xF) && !(inputStride & 0xF))
			{
				// Aligned input
				for (; i < vectorCount; i++)
				{
					VECTOR V = _mm_castsi128_ps(_mm_loadl_epi64(reinterpret_cast<const __m128i*>(pInputVector)));
					pInputVector += inputStride;

					VECTOR Y = PERMUTE_PS(V, _MM_SHUFFLE(1, 1, 1, 1));
					VECTOR X = PERMUTE_PS(V, _MM_SHUFFLE(0, 0, 0, 0));

					VECTOR vTemp = _mm_mul_ps(Y, row1);
					vTemp = FMADD_PS(X, row0, vTemp);

					_mm_store_sd(reinterpret_cast<double*>(pOutputVector), _mm_castps_pd(vTemp));
					pOutputVector += outputStride;
				}
			}
			else
			{
				// Unaligned input
				for (; i < vectorCount; i++)
				{
					__m128 xy = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(pInputVector)));
					pInputVector += inputStride;

					VECTOR Y = PERMUTE_PS(xy, _MM_SHUFFLE(1, 1, 1, 1));
					VECTOR X = PERMUTE_PS(xy, _MM_SHUFFLE(0, 0, 0, 0));

					VECTOR vTemp = _mm_mul_ps(Y, row1);
					vTemp = FMADD_PS(X, row0, vTemp);

					_mm_store_sd(reinterpret_cast<double*>(pOutputVector), _mm_castps_pd(vTemp));
					pOutputVector += outputStride;
				}
			}

			SFENCE();

			return pOutputStream;
#endif
		}
	}
}

#endif // !ULTREALITY_MATH_FLOAT2_INL
