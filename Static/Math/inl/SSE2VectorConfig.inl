#ifndef ULTREALITY_MATH_SSE2_CONFIG_INL
#define ULTREALITY_MATH_SSE2_CONFIG_INL

#if defined(__GNUC__) or defined(__clang__)
#define FORCE_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline
#endif

namespace UltReality::Math
{
	// VECTOR_F32 ****************************************************************
	// ***************************************************************************
	FORCE_INLINE VECTOR_F32::operator VECTOR() const noexcept
	{
		return v;
	}

	FORCE_INLINE VECTOR_F32::operator const float* () const noexcept
	{
		return f;
	}

#if defined(_SSE2_INTRINSICS_)
	FORCE_INLINE VECTOR_F32::operator __m128i() const noexcept
	{
		return _mm_castps_si128(v);
	}

	FORCE_INLINE VECTOR_F32::operator __m128d() const noexcept
	{
		return _mm_castps_pd(v);
	}
#endif

	// End VECTOR_F32 ************************************************************
	// ***************************************************************************

	// VECTOR_I32 ****************************************************************
	// ***************************************************************************

	FORCE_INLINE VECTOR_I32::operator VECTOR() const noexcept
	{
		return v;
	}

#if defined(_SSE2_INTRINSICS_)
	FORCE_INLINE VECTOR_I32::operator __m128i() const noexcept
	{
		return _mm_castps_si128(v);
	}

	FORCE_INLINE VECTOR_I32::operator __m128d() const noexcept
	{
		return _mm_castps_pd(v);
	}
#endif

	// End VECTOR_I32 ************************************************************
	// ***************************************************************************

	// VECTOR_U32 ****************************************************************
	// ***************************************************************************

	FORCE_INLINE VECTOR_U32::operator VECTOR() const noexcept
	{
		return v;
	}

#if defined(_SSE2_INTRINSICS_)
	FORCE_INLINE VECTOR_U32::operator __m128i() const noexcept
	{
		return _mm_castps_si128(v);
	}

	FORCE_INLINE VECTOR_U32::operator __m128d() const noexcept
	{
		return _mm_castps_pd(v);
	}
#endif

	// End VECTOR_U32 ************************************************************
	// ***************************************************************************

	// VECTOR_U8 *****************************************************************
	// ***************************************************************************

	FORCE_INLINE VECTOR_U8::operator VECTOR() const noexcept
	{
		return v;
	}

#if defined(_SSE2_INTRINSICS_)
	FORCE_INLINE VECTOR_U8::operator __m128i() const noexcept
	{
		return _mm_castps_si128(v);
	}

	FORCE_INLINE VECTOR_U8::operator __m128d() const noexcept
	{
		return _mm_castps_pd(v);
	}
#endif

	// End VECTOR_U8 *************************************************************
	// ***************************************************************************

	// VECTOR Overloads/Operators ************************************************
	// ***************************************************************************

#if !defined(_NO_VECTOR_OVERLOADS_)
	FORCE_INLINE VECTOR VEC_CALLCONV operator+ (A_VECTOR v) noexcept
	{
		return v;
	}

	FORCE_INLINE VECTOR VEC_CALLCONV operator- (A_VECTOR v) noexcept
	{
		return VEC::Negate(v);
	}

	FORCE_INLINE VECTOR& VEC_CALLCONV operator+= (VECTOR& v1, A_VECTOR v2) noexcept
	{
		v1 = VEC::Add(v1, v2);
		return v1;
	}

	FORCE_INLINE VECTOR& VEC_CALLCONV operator-= (VECTOR& v1, A_VECTOR v2) noexcept
	{
		v1 = VEC::Subtract(v1, v2);
		return v1;
	}

	FORCE_INLINE VECTOR& VEC_CALLCONV operator*= (VECTOR& v1, A_VECTOR v2) noexcept
	{
		v1 = VEC::Multiply(v1, v2);
		return v1;
	}

	FORCE_INLINE VECTOR& VEC_CALLCONV operator/= (VECTOR& v1, A_VECTOR v2) noexcept
	{
		v1 = VEC::Divide(v1, v2);
		return v1;
	}

	FORCE_INLINE VECTOR& operator*= (VECTOR& v, const float s) noexcept
	{
		v = VEC::Scale(v, s);
		return v;
	}

	FORCE_INLINE VECTOR& operator/= (VECTOR& v, const float s) noexcept
	{
		VECTOR v_s = VEC::Replicate(s);
		v = VEC::Divide(v, v_s);
		return v;
	}

	FORCE_INLINE VECTOR VEC_CALLCONV operator+ (A_VECTOR v1, A_VECTOR v2) noexcept
	{
		return VEC::Add(v1, v2);
	}

	FORCE_INLINE VECTOR VEC_CALLCONV operator- (A_VECTOR v1, A_VECTOR v2) noexcept
	{
		return VEC::Subtract(v1, v2);
	}

	FORCE_INLINE VECTOR VEC_CALLCONV operator* (A_VECTOR v1, A_VECTOR v2) noexcept
	{
		return VEC::Multiply(v1, v2);
	}

	FORCE_INLINE VECTOR VEC_CALLCONV operator/ (A_VECTOR v1, A_VECTOR v2) noexcept
	{
		return VEC::Divide(v1, v2);
	}

	FORCE_INLINE VECTOR VEC_CALLCONV operator* (A_VECTOR v, const float s) noexcept
	{
		return VEC::Scale(v, s);
	}

	FORCE_INLINE VECTOR VEC_CALLCONV operator/ (A_VECTOR v, const float s) noexcept
	{
		VECTOR v_s = VEC::Replicate(s);
		return VEC::Divide(v, v_s);
	}

	FORCE_INLINE VECTOR VEC_CALLCONV operator* (const float s, A_VECTOR v) noexcept
	{
		return VEC::Scale(v, s);
	}
#endif

	// End VECTOR Overloads/Operators ********************************************
	// ***************************************************************************

	namespace VEC
	{
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4701)
		// C4701: false positives
#endif
		FORCE_INLINE VECTOR VEC_CALLCONV ConvertIntToFloat(A_VECTOR vInt, uint32_t divExponent) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(divExponent < 32);
#endif

#if defined(_NO_INTRINSICS_)
			float fScale = 1.0f / static_cast<float>(1U << divExponent);
			uint32_t elementIndex = 0;
			VECTOR Result;
			do
			{
				auto iTemp = static_cast<int32_t>(vInt.vector4_u32[elementIndex]);
				Result.vector4_f32[elementIndex] = static_cast<float>(iTemp) * fScale;
			} while (++elementIndex < 4);

			return Result;
#elif defined(_SSE2_INTRINSICS_)
			// Convert to floats
			VECTOR vResult = _mm_cvtepi32_ps(_mm_castps_si128(vInt));

			// Convert divExponent into 1.0f/(1<<divExponent)
			uint32_t uScale = 0x3F800000U - (divExponent << 23);

			// Splat the scalar value
			__m128i vScale = _mm_set1_epi32(static_cast<int>(uScale));
			vResult = _mm_mul_ps(vResult, _mm_castsi128_ps(vScale));
			
			return vResult;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV ConvertFloatToInt(A_VECTOR vFloat, uint32_t mulExponent) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(mulExponent < 32);
#endif

#if defined(_NO_INTRINSICS_)
			// Get the scalar factor.
			auto fScale = static_cast<float>(1U << mulExponent);
			uint32_t elementIndex = 0;
			VECTOR Result;
			do 
			{
				int32_t iResult;
				float fTemp = vFloat.vector4_f32[elementIndex] * fScale;
				if (fTemp <= -(65536.0f * 32768.0f))
				{
					iResult = (-0x7FFFFFFF) - 1;
				}
				else if (fTemp > (65536.0f * 32768.0f) - 128.0f)
				{
					iResult = 0x7FFFFFFF;
				}
				else {
					iResult = static_cast<int32_t>(fTemp);
				}
				Result.vector4_u32[elementIndex] = static_cast<uint32_t>(iResult);
			} while (++elementIndex < 4);
			
			return Result;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vResult = _mm_set_ps1(static_cast<float>(1U << mulExponent));
			vResult = _mm_mul_ps(vResult, vFloat);
			
			// In case of positive overflow, detect it
			VECTOR vOverflow = _mm_cmpgt_ps(vResult, g_MaxInt);
			
			// Float to int conversion
			__m128i vResulti = _mm_cvttps_epi32(vResult);
			
			// If there was positive overflow, set to 0x7FFFFFFF
			vResult = _mm_and_ps(vOverflow, g_AbsMask);
			vOverflow = _mm_andnot_ps(vOverflow, _mm_castsi128_ps(vResulti));
			vOverflow = _mm_or_ps(vOverflow, vResult);
			return vOverflow;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV ConvertUIntToFloat(A_VECTOR vUInt, uint32_t divExponent) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(divExponent < 32);
#endif

#if defined(_NO_INTRINSICS_)
			float fScale = 1.0f / static_cast<float>(1U << divExponent);
			uint32_t elementIndex = 0;
			VECTOR Result;
			do 
			{
				Result.vector4_f32[elementIndex] = static_cast<float>(vUInt.vector4_u32[elementIndex]) * fScale;
			} while (++elementIndex < 4);
			
			return Result;

#elif defined(_SSE2_INTRINSICS_)
			// For the values that are higher than 0x7FFFFFFF, a fixup is needed
			// Determine which ones need the fix.
			VECTOR vMask = _mm_and_ps(vUInt, g_NegativeZero);
			// Force all values positive
			VECTOR vResult = _mm_xor_ps(vUInt, vMask);
			// Convert to floats
			vResult = _mm_cvtepi32_ps(_mm_castps_si128(vResult));
			// Convert 0x80000000 -> 0xFFFFFFFF
			__m128i iMask = _mm_srai_epi32(_mm_castps_si128(vMask), 31);
			// For only the ones that are too big, add the fixup
			vMask = _mm_and_ps(_mm_castsi128_ps(iMask), g_FixUnsigned);
			vResult = _mm_add_ps(vResult, vMask);
			// Convert DivExponent into 1.0f/(1<<DivExponent)
			uint32_t uScale = 0x3F800000U - (divExponent << 23);
			// Splat
			iMask = _mm_set1_epi32(static_cast<int>(uScale));
			vResult = _mm_mul_ps(vResult, _mm_castsi128_ps(iMask));
			return vResult;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV ConvertFloatToUInt(A_VECTOR vFloat, uint32_t mulExponent) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(mulExponent < 32);
#endif

#if defined(_NO_INTRINSICS_)
			// Get the scalar factor.
			auto fScale = static_cast<float>(1U << mulExponent);
			uint32_t elementIndex = 0;
			VECTOR Result;
			do 
			{
				uint32_t uResult;
				float fTemp = vFloat.vector4_f32[elementIndex] * fScale;
				if (fTemp <= 0.0f)
				{
					uResult = 0;
				}
				else if (fTemp >= (65536.0f * 65536.0f))
				{
					uResult = 0xFFFFFFFFU;
				}
				else {
					uResult = static_cast<uint32_t>(fTemp);
				}
				Result.vector4_u32[elementIndex] = uResult;
			} while (++elementIndex < 4);
			
			return Result;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vResult = _mm_set_ps1(static_cast<float>(1U << mulExponent));
			vResult = _mm_mul_ps(vResult, vFloat);
			// Clamp to >=0
			vResult = _mm_max_ps(vResult, g_Zero);
			// Any numbers that are too big, set to 0xFFFFFFFFU
			VECTOR vOverflow = _mm_cmpgt_ps(vResult, g_MaxUInt);
			VECTOR vValue = g_UnsignedFix;
			// Too large for a signed integer?
			VECTOR vMask = _mm_cmpge_ps(vResult, vValue);
			// Zero for number's lower than 0x80000000, 32768.0f*65536.0f otherwise
			vValue = _mm_and_ps(vValue, vMask);
			// Perform fixup only on numbers too large (Keeps low bit precision)
			vResult = _mm_sub_ps(vResult, vValue);
			__m128i vResulti = _mm_cvttps_epi32(vResult);
			// Convert from signed to unsigned pnly if greater than 0x80000000
			vMask = _mm_and_ps(vMask, g_NegativeZero);
			vResult = _mm_xor_ps(_mm_castsi128_ps(vResulti), vMask);
			// On those that are too large, set to 0xFFFFFFFF
			vResult = _mm_or_ps(vResult, vOverflow);
			return vResult;
#endif
		}
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4068 4214 4204 4365 4616 4640 6001 6101)
		// C4068/4616: ignore unknown pragmas
		// C4214/4204: nonstandard extension used
		// C4365/4640: Off by default noise
		// C6001/6101: False positives
#endif

#ifdef _PREFAST_
#pragma prefast(push)
#pragma prefast(disable : 25000, "FXMVECTOR is 16 bytes")
#pragma prefast(disable : 26495, "Union initialization confuses /analyze")
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-equal"
#pragma clang diagnostic ignored "-Wundefined-reinterpret-cast"
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif

		FORCE_INLINE VECTOR VEC_CALLCONV SetBinaryConstant(uint32_t C0, uint32_t C1, uint32_t C2, uint32_t C3) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 vResult;
			vResult.u[0] = (0 - (C0 & 1)) & 0x3F800000;
			vResult.u[0] = (0 - (C1 & 1)) & 0x3F800000;
			vResult.u[0] = (0 - (C2 & 1)) & 0x3F800000;
			vResult.u[0] = (0 - (C3 & 1)) & 0x3F800000;

			return vResult;
#elif defined(_SSE2_INTRINSICS_)
			static const VECTOR_U32 g_vMask1 = { { {1,1,1,1} } };
			// Move params to a vector
			__m128i vTemp = _mm_set_epi32(static_cast<int>(C3), static_cast<int>(C2), static_cast<int>(C1), static_cast<int>(C0));
			// Mask off the low bits
			vTemp = _mm_and_si128(vTemp, g_vMask1);
			// 0xFFFFFFFF on true bits
			vTemp = _mm_cmpeq_epi32(vTemp, g_vMask1);
			// 0xFFFFFFFF -> 1.0f, 0x00000000 -> 0.0f
			vTemp = _mm_and_si128(vTemp, g_One);

			return _mm_castsi128_ps(vTemp);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SplatConstant(int32_t intConstant, uint32_t divExponent) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(intConstant >= -16 && intConstant <= 15);
			assert(divExponent < 32);
#endif

#if defined(_NO_INTRINSICS_)
			using VECTOR::ConvertIntToFloat;

			VECTOR_I32 V = { { { intConstant, intConstant, intConstant, intConstant } } };

			return ConvertIntToFloat(V.v, divExponent);
#elif defined(_SSE2_INTRINSICS_)
			// Splat the int
			__m128i vScale = _mm_set1_epi32(intConstant);
			// Convert to a float
			VECTOR vResult = _mm_cvtepi32_ps(vScale);
			// Convert divExponent into 1.0f/(1<<divExponent)
			uint32_t uScale = 0x3F800000U - (divExponent << 23);
			// Splat the scalar value (it's really a float)
			vScale = _mm_set1_epi32(static_cast<int>(uScale));
			// Multiply by the reciprocal (perform a right shift by divExponent)
			vResult = _mm_mul_ps(vResult, _mm_castsi128_ps(vScale));

			return vResult;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SplatConstantInt(int32_t intConstant) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(intConstant >= -16 && intConstant <= 15);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR_I32 V = { { { intConstant, intConstant, intConstant, intConstant } } };
			return V.v
#elif defined(_SSE2_INTRINSICS_)
			__m128i V = _mm_set1_epi32(intConstant);
			return _mm_castsi128_ps(V);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadInt(const uint32_t* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR V;
			V.vector4_u32[0] = *pSource;
			V.vector4_u32[1] = 0;
			V.vector4_u32[2] = 0;
			V.vector4_u32[3] = 0;

			return V;
#elif defined(_SSE2_INTRINSICS_)
			return _mm_load_ss(reinterpret_cast<const float*>(pSource));
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadFloat(const float* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR V;
			V.vector4_f32[0] = *pSource;
			V.vector4_f32[1] = 0;
			V.vector4_f32[2] = 0;
			V.vector4_f32[3] = 0;

			return V;
#elif defined(_SSE2_INTRINSICS_)
			return _mm_load_ss(pSource);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadInt2(const uint32_t* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR V;
			V.vector4_u32[0] = pSource[0];
			V.vector4_u32[1] = pSource[1];
			V.vector4_u32[2] = 0;
			V.vector4_u32[3] = 0;

			return V;
#elif defined(_SSE2_INTRINSICS_)
			return _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(pSource)));
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadAInt2(const uint32_t* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
			assert((reinterpret_cast<uintptr_t>(pSource) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR V;
			V.vector4_u32[0] = pSource[0];
			V.vector4_u32[1] = pSource[1];
			V.vector4_u32[2] = 0;
			V.vector4_u32[3] = 0;

			return V;
#elif defined(_SSE2_INTRINSICS_)
			return _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(pSource)));
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadInt3(const uint32_t* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR V;
			V.vector4_u32[0] = pSource[0];
			V.vector4_u32[1] = pSource[1];
			V.vector4_u32[2] = pSource[2];
			V.vector4_u32[3] = 0;

			return V;
#elif defined(_SSE2_INTRINSICS_)
			__m128 xy = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(pSource)));
			__m128 z = _mm_load_ss(reinterpret_cast<const float*>(pSource + 2));
			return _mm_movelh_ps(xy, z);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadAInt3(const uint32_t* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
			assert((reinterpret_cast<uintptr_t>(pSource) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR V;
			V.vector4_u32[0] = pSource[0];
			V.vector4_u32[1] = pSource[1];
			V.vector4_u32[2] = pSource[2];
			V.vector4_u32[3] = 0;

			return V;
#elif defined(_SSE2_INTRINSICS_)
			// Load the first two elements (aligned load)
			__m128 xy = _mm_load_ps(reinterpret_cast<const float*>(pSource));

			// Mask out z and set w to 0
			__m128 z = _mm_load_ss(reinterpret_cast<const float*>(pSource + 2));
			return _mm_movelh_ps(xy, z);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadInt4(const uint32_t* pSource) noexcept
		{
#if defined(_DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR V;
			V.vector4_u32[0] = pSource[0];
			V.vector4_u32[1] = pSource[1];
			V.vector4_u32[2] = pSource[2];
			V.vector4_u32[3] = pSource[3];

			return V;
#elif defined(_SSE2_INTRINSICS_)
			__m128i V = _mm_loadu_si128(reinterpret_cast<const __m128i*>(pSource));
			return _mm_castsi128_ps(V);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV LoadAInt4(const uint32_t* pSource) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pSource != nullptr);
			assert((reinterpret_cast<uintptr_t>(pSource) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR V;
			V.vector4_u32[0] = pSource[0];
			V.vector4_u32[1] = pSource[1];
			V.vector4_u32[2] = pSource[2];
			V.vector4_u32[3] = pSource[3];

			return V;
#elif defined(_SSE2_INTRINSICS_)
			__m128i V = _mm_load_si128(reinterpret_cast<const __m128i*>(pSource));
			return _mm_castsi128_ps(V);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreInt(uint32_t* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			*pDestination = GetIntX(v);

#elif defined(_SSE2_INTRINSICS_)
			_mm_store_ss(reinterpret_cast<float*>(pDestination), v);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreFloat(float* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			*pDestination = GetX(v);

#elif defined(_SSE2_INTRINSICS_)
			_mm_store_ss(pDestination, v);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreInt2(uint32_t* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination[0] = v.vector4_u32[0];
			pDestination[1] = v.vector4_u32[1];

#elif defined(_SSE2_INTRINSICS_)
			_mm_store_sd(reinterpret_cast<double*>(pDestination), _mm_castps_pd(v));
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreAInt2(uint32_t* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
			assert((reinterpret_cast<uintptr_t>(pDestination) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination[0] = v.vector4_u32[0];
			pDestination[1] = v.vector4_u32[1];

#elif defined(_SSE2_INTRINSICS_)
			_mm_store_sd(reinterpret_cast<double*>(pDestination), _mm_castps_pd(v));
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreInt3(uint32_t* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination[0] = v.vector4_u32[0];
			pDestination[1] = v.vector4_u32[1];
			pDestination[2] = v.vector4_u32[2];

#elif defined(_SSE2_INTRINSICS_)
			_mm_store_sd(reinterpret_cast<double*>(pDestination), _mm_castps_pd(v));
			__m128 z = PERMUTE_PS(v, _MM_SHUFFLE(2, 2, 2, 2));
			_mm_store_ss(reinterpret_cast<float*>(&pDestination[2]), z);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreAInt3(uint32_t* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
			assert((reinterpret_cast<uintptr_t>(pDestination) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination[0] = v.vector4_u32[0];
			pDestination[1] = v.vector4_u32[1];
			pDestination[2] = v.vector4_u32[2];

#elif defined(_SSE2_INTRINSICS_)
			_mm_store_sd(reinterpret_cast<double*>(pDestination), _mm_castps_pd(v));
			__m128 z = _mm_movehl_ps(v, v);
			_mm_store_ss(reinterpret_cast<float*>(&pDestination[2]), z);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreInt4(uint32_t* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination[0] = v.vector4_u32[0];
			pDestination[1] = v.vector4_u32[1];
			pDestination[2] = v.vector4_u32[2];
			pDestination[3] = v.vector4_u32[3];

#elif defined(_SSE2_INTRINSICS_)
			_mm_storeu_si128(reinterpret_cast<__m128i*>(pDestination), _mm_castps_si128(v));
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV StoreAInt4(uint32_t* pDestination, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pDestination != nullptr);
			assert((reinterpret_cast<uintptr_t>(pDestination) & 0xF) == 0);
#endif

#if defined(_NO_INTRINSICS_)
			pDestination[0] = v.vector4_u32[0];
			pDestination[1] = v.vector4_u32[1];
			pDestination[2] = v.vector4_u32[2];
			pDestination[3] = v.vector4_u32[3];

#elif defined(_SSE2_INTRINSICS_)
			_mm_store_si128(reinterpret_cast<__m128i*>(pDestination), _mm_castps_si128(v));
#endif
		}

#if defined(_NO_INTRINSICS_)
#define ISNAN(x) isnan(x)
#define ISINF(x) isinf(x)
#endif

#if defined(_SSE2_INTRINSICS_)
#define UNPACK3INTO4(l1, l2, l3) \
	VECTOR V3 = _mm_shuffle_ps(l2, l3, _MM_SHUFFLE(0, 0, 3, 2));\
	VECTOR V2 = _mm_shuffle_ps(l2, l1, _MM_SHUFFLE(3, 3, 1, 0));\
	V2 = PERMUTE_PS(V2, _MM_SHUFFLE(1, 1, 0, 2));\
	VECTOR V4 = _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(l3), 32 / 8))

#define PACK4INTO3(v2x) \
	v2x = _mm_shuffle_ps(V2, V3, _MM_SHUFFLE(1, 0, 2, 1));\
	V2 = _mm_shuffle_ps(V2, V1, _MM_SHUFFLE(2, 2, 0, 0));\
	V1 = _mm_shuffle_ps(V1, V2, _MM_SHUFFLE(0, 2, 1, 0));\
	V3 = _mm_shuffle_ps(V3, V4, _MM_SHUFFLE(0, 0, 2, 2));\
	V3 = _mm_shuffle_ps(V3, V4, _MM_SHUFFLE(2, 1, 2, 0))
#endif

		FORCE_INLINE VECTOR VEC_CALLCONV Zero() noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 vResult = { { { 0.0f, 0.0f, 0.0f, 0.0f } } };
			return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_setzero_ps();
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Set(float x, float y, float z, float w) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 vResult = { { { x, y, z, w } } };
			return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_set_ps(w, z, y, x);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SetInt(uint32_t x, uint32_t y, uint32_t z, uint32_t w) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 vResult = { { { x, y, z, w } } };
			return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
			__m128i v = _mm_set_epi32(static_cast<int>(w), static_cast<int>(z), static_cast<int>(y), static_cast<int>(x));
			return _mm_castsi128_ps(v);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Replicate(float value) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 vResult;
			vResult.f[0] = vResult[1] = vResult[2] = vResult[3] = value;
			return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_set_ps1(value);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV ReplicatePtr(const float* pValue) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 vResult;
			vResult.f[0] = vResult[1] = vResult[2] = vResult[3] = *pValue;
			return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_load_ps1(pValue);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV ReplicateInt(uint32_t value) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 vResult;
			vResult.u[0] = vResult.u[1] = vResult.u[2] = vResult.u[3] = value;
			return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
			__m128i vTemp = _mm_set1_epi32(static_cast<int>(value));
			return _mm_castsi128_ps(vTemp);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV ReplicateIntPtr(const uint32_t* pValue) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 vResult;
			vResult.u[0] = vResult.u[1] = vResult.u[2] = vResult.u[3] = *pValue;
			return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_load_ps1(reinterpret_cast<const float*>(pValue));
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV TrueInt() noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 vResult = { { { 0xFFFFFFFFU, 0xFFFFFFFFU, 0xFFFFFFFFU, 0xFFFFFFFFU } } };
			return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
			__m128i v = _mm_set1_epi32(-1);
			return _mm_castsi128_ps(v);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV FalseInt() noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 vResult = { { { 0.0f, 0.0f, 0.0f, 0.0f } } };
			return vResult;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_setzero_ps();
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SplatX(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 vResult;
			vResult.f[0] = vResult.f[1] = vResult.f[2] = vResult.f[3] = v.vector4_f32[0];
			return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
			return PERMUTE_PS(v, _MM_SHUFFLE(0, 0, 0, 0));
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SplatY(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 vResult;
			vResult.f[0] = vResult.f[1] = vResult.f[2] = vResult.f[3] = v.vector4_f32[1];
			return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
			return PERMUTE_PS(v, _MM_SHUFFLE(1, 1, 1, 1));
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SplatZ(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 vResult;
			vResult.f[0] = vResult.f[1] = vResult.f[2] = vResult.f[3] = v.vector4_f32[2];
			return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
			return PERMUTE_PS(v, _MM_SHUFFLE(2, 2, 2, 2));
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SplatW(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 vResult;
			vResult.f[0] = vResult.f[1] = vResult.f[2] = vResult.f[3] = vResult.vector4_f32[3];
			return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
			return PERMUTE_PS(v, _MM_SHUFFLE(3, 3, 3, 3));
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SplatOne() noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 vResult;
			vResult.f[0] = vResult.f[1] = vResult.f[2] = vResult.f[3] = 1.0f;
			return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
			return g_One;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SplatInfinity() noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 vResult;
			vResult.u[0] = vResult.u[1] = vResult.u[2] = vResult.u[3] = 0x7F800000;
			return vResult;

#elif defined(_SSE2_INTRINSICS_)
			return g_Infinity;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SplatQNaN() noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 vResult;
			vResult.u[0] = vResult.u[1] = vResult.u[2] = vResult.u[3] = 0x7FC00000;
			return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
			return g_QNaN;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SplatEpsilon() noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 vResult;
			vResult.u[0] = vResult.u[1] = vResult.u[2] = vResult.u[3] = 0x34000000;
			return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
			return g_Epsilon;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SplatSignMask() noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 vResult;
			vResult.u[0] = vResult.u[1] = vResult.u[2] = vResult.u[3] = 0x80000000U;
			return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
			__m128i v = _mm_set1_epi32(static_cast<int>(0x80000000));
			return _mm_castsi128_ps(v);
#endif
		}

		FORCE_INLINE float VEC_CALLCONV GetByIndex(A_VECTOR v, size_t index) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(index < 4);
			_Analysis_assume_(index < 4);
#endif

#if defined(_NO_INTRINSICS_)
			return v.vector4_f32[index];

#else
			VECTOR_F32 U;
			U.v = v;
			return U.f[index];
#endif
		}

		FORCE_INLINE float VEC_CALLCONV GetX(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return v.vector4_f32[0];

#elif defined(_SSE2_INTRINSICS_)
			return _mm_cvtss_f32(v);
#endif
		}

		FORCE_INLINE float VEC_CALLCONV GetY(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return v.vector4_f32[1];

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = PERMUTE_PS(v, _MM_SHUFFLE(1, 1, 1, 1));
			return _mm_cvtss_f32(vTemp);
#endif
		}

		FORCE_INLINE float VEC_CALLCONV GetZ(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return v.vector4_f32[2];

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = PERMUTE_PS(v, _MM_SHUFFLE(2, 2, 2, 2));
			return _mm_cvtss_f32(vTemp);
#endif
		}

		FORCE_INLINE float VEC_CALLCONV GetW(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return v.vector4_f32[3];

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = PERMUTE_PS(v, _MM_SHUFFLE(3, 3, 3, 3));
			return _mm_cvtss_f32(vTemp);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV GetByIndexPtr(float* f, A_VECTOR v, size_t index) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(f != nullptr);
			assert(index < 4);
			_Analysis_assume_(index < 4);
#endif

#if defined(_NO_INTRINSICS_)
			*f = v.vector4_f32[index];

#else
			VECTOR_F32 U;
			U.v = v;
			*f = U.f[index];
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV GetXPtr(float* x, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(x != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			*x = v.vector4_f32[0];

#elif defined(_SSE2_INTRINSICS_)
			_mm_store_ss(x, v);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV GetYPtr(float* y, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(y != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			*y = v.vector4_f32[1];

#elif defined(_SSE4_INTRINSICS_)
			*(reinterpret_cast<int*>(y)) = _mm_extract_ps(v, 1);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = PERMUTE_PS(v, _MM_SHUFFLE(1, 1, 1, 1));
			_mm_store_ss(y, vTemp);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV GetZPtr(float* z, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(z != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			*z = v.vector4_f32[2];

#elif defined(_SSE4_INTRINSICS_)
			*(reinterpret_cast<int*>(z)) = _mm_extract_ps(v, 2);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = PERMUTE_PS(v, _MM_SHUFFLE(2, 2, 2, 2));
			_mm_store_ss(z, vTemp);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV GetWPtr(float* w, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(w != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			*w = v.vector4_f32[3];

#elif defined(_SSE4_INTRINSICS_)
			*(reinterpret_cast<int*>(w)) = _mm_extract_ps(v, 3);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = PERMUTE_PS(v, _MM_SHUFFLE(3, 3, 3, 3));
			_mm_store_ss(w, vTemp);
#endif
		}

		FORCE_INLINE uint32_t VEC_CALLCONV GetIntByIndex(A_VECTOR v, size_t index) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(index < 4);
			_Analysis_assume_(index < 4);
#endif

#if defined(_NO_INTRINSICS_)
			return v.vector4_u32[index];

#else
			VECTOR_U32 U;
			U.v = v;
			return U.u[index];
#endif
		}

		FORCE_INLINE uint32_t VEC_CALLCONV GetIntX(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return v.vector4_u32[0];

#elif defined(_SSE2_INTRINSICS_)
			return static_cast<uint32_t>(_mm_cvtsi128_si32(_mm_castps_si128(v)));
#endif
		}

		FORCE_INLINE uint32_t VEC_CALLCONV GetIntY(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return v.vector4_u32[1];

#elif defined(_SSE4_INTRINSICS_)
			__m128i V1 = _mm_castps_si128(v);
			return static_cast<uint32_t>(_mm_extract_epi32(V1, 1));

#elif defined(_SSE2_INTRINSICS_)
			__m128i vResulti = _mm_shuffle_epi32(_mm_castps_si128(v), _MM_SHUFFLE(1, 1, 1, 1));
			return static_cast<uint32_t>(_mm_cvtsi128_si32(vResulti));
#endif
		}

		FORCE_INLINE uint32_t VEC_CALLCONV GetIntZ(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return v.vector4_u32[2];

#elif defined(_SSE4_INTRINSICS_)
			__m128i V1 = _mm_castps_si128(v);
			return static_cast<uint32_t>(_mm_extract_epi32(V1, 2));

#elif defined(_SSE2_INTRINSICS_)
			__m128i vResulti = _mm_shuffle_epi32(_mm_castps_si128(v), _MM_SHUFFLE(2, 2, 2, 2));
			return static_cast<uint32_t>(_mm_cvtsi128_si32(vResulti));
#endif
		}

		FORCE_INLINE uint32_t VEC_CALLCONV GetIntW(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			return v.vector4_u32[3];

#elif defined(_SSE4_INTRINSICS_)
			__m128i V1 = _mm_castps_si128(v);
			return static_cast<uint32_t>(_mm_extract_epi32(V1, 3));

#elif defined(_SSE2_INTRINSICS_)
			__m128i vResulti = _mm_shuffle_epi32(_mm_castps_si128(v), _MM_SHUFFLE(3, 3, 3, 3));
			return static_cast<uint32_t>(_mm_cvtsi128_si32(vResulti));
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV GetIntByIndexPtr(uint32_t* i, A_VECTOR v, size_t index) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(i != nullptr);
			assert(index < 4);
			_Analysis_assume_(index < 4);
#endif

#if defined(_NO_INTRINSICS_)
			*i = v.vector4_u32[index];

#else
			VECTOR_U32 U;
			U.v = v;
			*i = U.u[index];
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV GetIntXPtr(uint32_t* x, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(x != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			*x = v.vector4_u32[0];

#elif defined(_SSE2_INTRINSICS_)
			_mm_store_ss(reinterpret_cast<float*>(x), v);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV GetIntYPtr(uint32_t* y, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(y != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			*y = v.vector4_u32[1];

#elif defined(_SSE4_INTRINSICS_)
			__m128i V1 = _mm_castps_si128(v);
			*y = static_cast<uint32_t>(_mm_extract_epi32(V1, 1));

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vResult = PERMUTE_PS(v, _MM_SHUFFLE(1, 1, 1, 1));
			_mm_store_ss(reinterpret_cast<float*>(y), vResult);
#endif
		}
		
		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV GetIntZPtr(uint32_t* z, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(z != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			*z = v.vector4_u32[2];

#elif defined(_SSE4_INTRINSICS_)
			__m128i V1 = _mm_castps_si128(v);
			*z = static_cast<uint32_t>(_mm_extract_epi32(V1, 2));

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vResult = PERMUTE_PS(v, _MM_SHUFFLE(2, 2, 2, 2));
			_mm_store_ss(reinterpret_cast<float*>(z), vResult);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE void VEC_CALLCONV GetIntWPtr(uint32_t* w, A_VECTOR v) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(w != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			*w = v.vector4_u32[3];

#elif defined(_SSE4_INTRINSICS_)
			__m128i V1 = _mm_castps_si128(v);
			*w = static_cast<uint32_t>(_mm_extract_epi32(V1, 3));

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vResult = PERMUTE_PS(v, _MM_SHUFFLE(3, 3, 3, 3));
			_mm_store_ss(reinterpret_cast<float*>(w), vResult);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SetByIndex(A_VECTOR v, float f, size_t index) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(index < 4);
			_Analysis_assume_(index < 4);
#endif

			VECTOR_F32 U;
			U.v = v;
			U.f[index] = f;
			return U.v;
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SetX(A_VECTOR v, float x) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 U = { { { x, v.vector4_f32[1], v.vector4_f32[2], v.vector_f32[3] } } };
			return U.v;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vResult = _mm_set_ss(x);
			vResult = _mm_move_ss(v, vResult);
			return vResult;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SetY(A_VECTOR v, float y) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 U = { { { v.vector4_f32[0], y, v.vector4_f32[2], v.vector_f32[3] } } };
			return U.v;

#elif defined(_SSE4_INTRINSICS_)
			VECTOR vResult = _mm_set_ss(y);
			vResult = _mm_insert_ps(v, vResult, 0x10);
			return vResult;

#elif defined(_SSE2_INTRINSICS_)
			// Swap y and x
			VECTOR vResult = PERMUTE_PS(v, _MM_SHUFFLE(3, 2, 0, 1));

			// Convert input to vector
			VECTOR vTemp = _mm_set_ss(y);

			// Replace the x component (y)
			vResult = _mm_move_ss(vResult, vTemp);

			// Swap y and x again
			vResult = PERMUTE_PS(vResult, _MM_SHUFFLE(3, 2, 0, 1));

			return vResult;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SetZ(A_VECTOR v, float z) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 U = { { { v.vector4_f32[0], v.vector4_f32[1], z, v.vector_f32[3] } } };
			return U.v;

#elif defined(_SSE4_INTRINSICS_)
			VECTOR vResult = _mm_set_ss(z);
			vResult = _mm_insert_ps(v, vResult, 0x20);
			return vResult;

#elif defined(_SSE2_INTRINSICS_)
			// Swap z and x
			VECTOR vResult = PERMUTE_PS(v, _MM_SHUFFLE(3, 0, 1, 2));

			// Convert input to vector
			VECTOR vTemp = _mm_set_ss(z);

			// Replace the x component (z)
			vResult = _mm_move_ss(vResult, vTemp);

			// Swap z and x again
			vResult = PERMUTE_PS(vResult, _MM_SHUFFLE(3, 0, 1, 2));

			return vResult;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SetW(A_VECTOR v, float w) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 U = { { { v.vector4_f32[0], v.vector4_f32[1], v.vector4_f32[2], w}}};
			return U.v;

#elif defined(_SSE4_INTRINSICS_)
			VECTOR vResult = _mm_set_ss(w);
			vResult = _mm_insert_ps(v, vResult, 0x30);
			return vResult;

#elif defined(_SSE2_INTRINSICS_)
			// Swap w and x
			VECTOR vResult = PERMUTE_PS(v, _MM_SHUFFLE(0, 2, 1, 3));

			// Convert input to vector
			VECTOR vTemp = _mm_set_ss(w);

			// Replace the x component (w)
			vResult = _mm_move_ss(vResult, vTemp);

			// Swap w and x again
			vResult = PERMUTE_PS(vResult, _MM_SHUFFLE(0, 2, 1, 3));

			return vResult;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV SetByIndexPtr(A_VECTOR v, const float* f, size_t index) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(f != nullptr);
			assert(index < 4);
			_Analysis_assume_(index < 4);
#endif

			VECTOR_F32 U;
			U.v = v;
			U.f[index] = *f;
			return U.v;
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV SetXPtr(A_VECTOR v, const float* x) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(x != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR_F32 U = { { { *x, v.vector4_f32[1], v.vector4_f32[2], v.vector_f32[3] } } };
			return U.v;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vResult = _mm_load_ss(x);
			vResult = _mm_move_ss(v, vResult);
			return vResult;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV SetYPtr(A_VECTOR v, const float* y) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(y != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR_F32 U = { { { v.vector4_f32[0], *y, v.vector4_f32[2], v.vector_f32[3] } } };
			return U.v;

#elif defined(_SSE2_INTRINSICS_)
			// Swap y and x
			VECTOR vResult = PERMUTE_PS(v, _MM_SHUFFLE(3, 2, 0, 1));

			// Convert input to vector
			VECTOR vTemp = _mm_load_ss(y);

			// Replace the x component (y)
			vResult = _mm_move_ss(vResult, vTemp);

			// Swap y and x again
			vResult = PERMUTE_PS(vResult, _MM_SHUFFLE(3, 2, 0, 1));

			return vResult;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV SetZPtr(A_VECTOR v, const float* z) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(z != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR_F32 U = { { { v.vector4_f32[0], v.vector4_f32[1], *z, v.vector_f32[3] } } };
			return U.v;

#elif defined(_SSE2_INTRINSICS_)
			// Swap z and x
			VECTOR vResult = PERMUTE_PS(v, _MM_SHUFFLE(3, 0, 1, 2));

			// Convert input to vector
			VECTOR vTemp = _mm_load_ss(z);

			// Replace the x component (z)
			vResult = _mm_move_ss(vResult, vTemp);

			// Swap z and x again
			vResult = PERMUTE_PS(vResult, _MM_SHUFFLE(3, 0, 1, 2));

			return vResult;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV SetWPtr(A_VECTOR v, const float* w) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(w != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR_F32 U = { { { v.vector4_f32[0], v.vector4_f32[1], v.vector4_f32[2], *w}} };
			return U.v;

#elif defined(_SSE2_INTRINSICS_)
			// Swap w and x
			VECTOR vResult = PERMUTE_PS(v, _MM_SHUFFLE(0, 2, 1, 3));

			// Convert input to vector
			VECTOR vTemp = _mm_load_ss(w);

			// Replace the x component (w)
			vResult = _mm_move_ss(vResult, vTemp);

			// Swap w and x again
			vResult = PERMUTE_PS(vResult, _MM_SHUFFLE(0, 2, 1, 3));

			return vResult;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SetIntX(A_VECTOR v, uint32_t x) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 U = { { { x, v.vector4_u32[1], v.vector4_u32[2], v.vector_u32[3] } } };
			return U.v;

#elif defined(_SSE2_INTRINSICS_)
			__m128i vTemp = _mm_cvtsi32_si128(static_cast<int>(x));
			VECTOR vResult = _mm_move_ss(v, _mm_castsi128_ps(vTemp));
			return vResult;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SetIntY(A_VECTOR v, uint32_t y) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 U = { { { v.vector4_u32[0], y, v.vector4_u32[2], v.vector_u32[3] } } };
			return U.v;

#elif defined(_SSE4_INTRINSICS_)
			__m128i vResult = _mm_castps_si128(v);
			vResult = _mm_insert_epi32(vResult, static_cast<int>(y), 1);
			return _mm_castsi128_ps(vResult);

#elif defined(_SSE2_INTRINSICS_)
			// Swap y and x
			VECTOR vResult = PERMUTE_PS(v, _MM_SHUFFLE(3, 2, 0, 1));

			// Convert input to vector
			__m128i vTemp = _mm_cvtsi32_si128(static_cast<int>(y));

			// Replace the x component (y)
			vResult = _mm_move_ss(vResult, _mm_castsi128_ps(vTemp));

			// Swap y and x again
			vResult = PERMUTE_PS(vResult, _MM_SHUFFLE(3, 2, 0, 1));
			return vResult;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SetIntZ(A_VECTOR v, uint32_t z) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 U = { { { v.vector4_u32[0], v.vector4_u32[1], z, v.vector_u32[3] } } };
			return U.v;

#elif defined(_SSE4_INTRINSICS_)
			__m128i vResult = _mm_castps_si128(v);
			vResult = _mm_insert_epi32(vResult, static_cast<int>(z), 2);
			return _mm_castsi128_ps(vResult);

#elif defined(_SSE2_INTRINSICS_)
			// Swap z and x
			VECTOR vResult = PERMUTE_PS(v, _MM_SHUFFLE(3, 0, 1, 2));

			// Convert input to vector
			__m128i vTemp = _mm_cvtsi32_si128(static_cast<int>(z));

			// Replace the x component (z)
			vResult = _mm_move_ss(vResult, _mm_castsi128_ps(vTemp));

			// Swap z and x again
			vResult = PERMUTE_PS(vResult, _MM_SHUFFLE(3, 0, 1, 2));
			return vResult;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SetIntW(A_VECTOR v, uint32_t w) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 U = { { { v.vector4_u32[0], v.vector4_u32[1], v.vector_u32[2], w } } };
			return U.v;

#elif defined(_SSE4_INTRINSICS_)
			__m128i vResult = _mm_castps_si128(v);
			vResult = _mm_insert_epi32(vResult, static_cast<int>(w), 3);
			return _mm_castsi128_ps(vResult);

#elif defined(_SSE2_INTRINSICS_)
			// Swap w and x
			VECTOR vResult = PERMUTE_PS(v, _MM_SHUFFLE(0, 2, 1, 3));

			// Convert input to vector
			__m128i vTemp = _mm_cvtsi32_si128(static_cast<int>(w));

			// Replace the x component (w)
			vResult = _mm_move_ss(vResult, _mm_castsi128_ps(vTemp));

			// Swap w and x again
			vResult = PERMUTE_PS(vResult, _MM_SHUFFLE(0, 2, 1, 3));
			return vResult;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV SetIntByIndexPtr(A_VECTOR v, const uint32_t* i, size_t index) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(i != nullptr);
			assert(index < 4);
			_Analysis_assume_(index < 4);
#endif

			VECTOR_U32 vTemp;
			vTemp.v = v;
			vTemp.u[index] = *i;
			return vTemp;
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV SetIntXPtr(A_VECTOR v, const uint32_t* x) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(x != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR_U32 U = { { { *x, v.vector4_u32[1], v.vector4_u32[2], v.vector_u32[3] } } };
			return U.v;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_load_ss(reinterpret_cast<const float*>(x));
			VECTOR vResult = _mm_move_ss(v, vTemp);
			return vResult;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV SetIntYPtr(A_VECTOR v, const uint32_t* y) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(y != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR_U32 U = { { { v.vector4_u32[0], *y, v.vector4_u32[2], v.vector_u32[3] } } };
			return U.v;

#elif defined(_SSE2_INTRINSICS_)
			// Swap y and x
			VECTOR vResult = PERMUTE_PS(v, _MM_SHUFFLE(3, 2, 0, 1));

			// Convert input to vector
			VECTOR vTemp = _mm_load_ss(reinterpret_cast<const float*>(y));

			// Replace the x component (y)
			vResult = _mm_move_ss(vResult, vTemp);

			// Swap y and x again
			vResult = PERMUTE_PS(vResult, _MM_SHUFFLE(3, 2, 0, 1));
			return vResult;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV SetIntZPtr(A_VECTOR v, const uint32_t* z) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(z != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR_U32 U = { { { v.vector4_u32[0], v.vector4_u32[1], *z, v.vector_u32[3] } } };
			return U.v;

#elif defined(_SSE2_INTRINSICS_)
			// Swap z and x
			VECTOR vResult = PERMUTE_PS(v, _MM_SHUFFLE(3, 0, 1, 2));

			// Convert input to vector
			VECTOR vTemp = _mm_load_ss(reinterpret_cast<const float*>(z));

			// Replace the x component (z)
			vResult = _mm_move_ss(vResult, vTemp);

			// Swap z and x again
			vResult = PERMUTE_PS(vResult, _MM_SHUFFLE(3, 0, 1, 2));
			return vResult;
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV SetIntWPtr(A_VECTOR v, const uint32_t* w) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(w != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR_U32 U = { { { v.vector4_u32[0], v.vector4_u32[1], v.vector_u32[2], *w } } };
			return U.v;

#elif defined(_SSE2_INTRINSICS_)
			// Swap w and x
			VECTOR vResult = PERMUTE_PS(v, _MM_SHUFFLE(0, 2, 1, 3));

			// Convert input to vector
			VECTOR vTemp = _mm_load_ss(reinterpret_cast<const float*>(w));

			// Replace the x component (w)
			vResult = _mm_move_ss(vResult, vTemp);

			// Swap w and x again
			vResult = PERMUTE_PS(vResult, _MM_SHUFFLE(0, 2, 1, 3));
			return vResult;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Swizzle(A_VECTOR v, uint32_t E0, uint32_t E1, uint32_t E2, uint32_t E3) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert((E0 < 4) && (E1 < 4) && (E2 < 4) && (E3 < 4));
			_Analysis_assume_((E0 < 4) && (E1 < 4) && (E2 < 4) && (E3 < 4));
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR_F32 Result = { { { v.vector4_f32[E0], v.vector4_f32[E1], v.vector4_f32[E2], v.vector4_f32[E3] } } };
			return Result;

#elif defined(_SSE2_INTRINSICS_)
			unsigned int elem[4] = { E0, E1, E2, E3 };
			__m128i vControl = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&elem[0]));
			return _mm_permutevar_ps(v, vControl);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Permute(A_VECTOR V1, A_VECTOR V2, uint32_t permuteX, uint32_t permuteY, uint32_t permuteZ, uint32_t permuteW) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(permuteX <= 7 && permuteY <= 7 && permuteZ <= 7 && permuteW <= 7);
			_Analysis_assume_(permuteX <= 7 && permuteY <= 7 && permuteZ <= 7 && permuteW <= 7);
#endif

#if defined(_NO_INTRINSICS_)
			const uint32_t* aPtr[2];
			aPtr[0] = reinterpret_cast<const uint32_t*>(&V1);
			aPtr[1] = reinterpret_cast<const uint32_t*>(&V2);

			VECTOR Result;
			auto pWork = reinterpret_cast<uint32_t*>(&Result);

			const uint32_t i0 = permuteX & 3;
			const uint32_t vi0 = permuteX >> 2;
			pWork[0] = aPtr[vi0][i0];

			const uint32_t i1 = permuteY & 3;
			const uint32_t vi1 = permuteY >> 2;
			pWork[1] = aPtr[vi1][i1];

			const uint32_t i2 = permuteZ & 3;
			const uint32_t vi2 = permuteZ >> 2;
			pWork[2] = aPtr[vi2][i2];

			const uint32_t i3 = permuteW & 3;
			const uint32_t vi3 = permuteW >> 2;
			pWork[3] = aPtr[vi3][i3];

			return Result;

#elif defined(_SSE2_INTRINSICS_)
			static const VECTOR_U32 three = { { { 3, 3, 3, 3 } } };

			ALIGNED(16) unsigned int elem[4] = { permuteX, permuteY, permuteZ, permuteW };
			__m128i vControl = _mm_load_si128(reinterpret_cast<const __m128i*>(&elem[0]));

			__m128i vSelect = _mm_cmpgt_epi32(vControl, three);
			vControl = _mm_castps_si128(_mm_and_ps(_mm_castsi128_ps(vControl), three));

			__m128 shuffled1 = _mm_permutevar_ps(V1, vControl);
			__m128 shuffled2 = _mm_permutevar_ps(V2, vControl);

			__m128 masked1 = _mm_andnot_ps(_mm_castsi128_ps(vSelect), shuffled1);
			__m128 masked2 = _mm_and_ps(_mm_castsi128_ps(vSelect), shuffled2);

			return _mm_or_ps(masked1, masked2);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SelectControl(uint32_t vectorIndex0, uint32_t vectorIndex1, uint32_t vectorIndex2, uint32_t vectorIndex3) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(vectorIndex0 < 2);
			assert(vectorIndex1 < 2);
			assert(vectorIndex2 < 2);
			assert(vectorIndex3 < 2);
			_Analysis_assume_(vectorIndex0 < 2);
			_Analysis_assume_(vectorIndex1 < 2);
			_Analysis_assume_(vectorIndex2 < 2);
			_Analysis_assume_(vectorIndex3 < 2);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR control;
			const uint32_t controlElement[] = { SELECT_NONE, SELECT_ALL };

			control.vector4_u32[0] = controlElement[vectorIndex0];
			control.vector4_u32[1] = controlElement[vectorIndex1];
			control.vector4_u32[2] = controlElement[vectorIndex2];
			control.vector4_u32[3] = controlElement[vectorIndex3];

			return control;

#elif defined(_SSE2_INTRINSICS_)
			// x = index0, y = index1, z = index2, w = index3
			__m128i vTemp = _mm_set_epi32(static_cast<int>(vectorIndex3), static_cast<int>(vectorIndex2), static_cast<int>(vectorIndex1), static_cast<int>(vectorIndex0));

			// Non zero entries set to 0xFFFFFFFF
			vTemp = _mm_cmpgt_epi32(vTemp, g_Zero);

			return _mm_castsi128_ps(vTemp);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Select(A_VECTOR V1, A_VECTOR V2, A_VECTOR control) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 result = 
			{ { {
					(V1.vector4_u32[0] & ~control.vector4_u32[0]) | (V2.vector4_u32[0] & control.vector4_u32[0]),
					(V1.vector4_u32[1] & ~control.vector4_u32[1]) | (V2.vector4_u32[1] & control.vector4_u32[1]),
					(V1.vector4_u32[2] & ~control.vector4_u32[2]) | (V2.vector4_u32[2] & control.vector4_u32[2]),
					(V1.vector4_u32[3] & ~control.vector4_u32[3]) | (V2.vector4_u32[3] & control.vector4_u32[3])
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp1 = _mm_andnot_ps(control, V1);
			VECTOR vTemp2 = _mm_and_ps(V2, control);

			return _mm_or_ps(vTemp1, vTemp2);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV MergeXY(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 result =
			{ { {
					V1.vector4_u32[0],
					V2.vector4_u32[0],
					V1.vector4_u32[1],
					V1.vector4_u32[1]
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_unpacklo_ps(V1, V2);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV MergeZW(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 result =
			{ { {
					V1.vector4_u32[2],
					V2.vector4_u32[2],
					V1.vector4_u32[3],
					V2.vector4_u32[3]
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_unpackhi_ps(V1, V2);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV ShiftLeft(A_VECTOR V1, A_VECTOR V2, uint32_t elements) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(elements < 4);
			_Analysis_assume_(elements < 4);
#endif

			return Permute(V1, V2, elements, (elements + 1), (elements + 2), (elements + 3));
		}

		FORCE_INLINE VECTOR VEC_CALLCONV RotateLeft(A_VECTOR v, uint32_t elements) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(elements < 4);
			_Analysis_assume_(elements < 4);
#endif

			return Swizzle(v, elements & 3, (elements + 1) & 3, (elements + 2) & 3, (elements + 3) & 3);
		}

		FORCE_INLINE VECTOR VEC_CALLCONV RotateRight(A_VECTOR v, uint32_t elements) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(elements < 4);
			_Analysis_assume_(elements < 4);
#endif

			return Swizzle(v, (4 - elements) & 3, (5 - elements) & 3, (6 - elements) & 3, (7 - elements) & 3);
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Insert(
			A_VECTOR vDestination,
			A_VECTOR vSource,
			uint32_t VSLeftRotateElements,
			uint32_t select0, uint32_t select1, uint32_t select2, uint32_t select3) noexcept
		{
			VECTOR control = SelectControl(select0 & 1, select1 & 1, select2 & 1, select3 & 1);
			return Select(vDestination, RotateLeft(vSource, VSLeftRotateElements), control);
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Equal(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 control = { { {
				(V1.vector4_f32[0] == V2.vector4_f32[0]) ? 0xFFFFFFFF : 0,
				(V1.vector4_f32[1] == V2.vector4_f32[1]) ? 0xFFFFFFFF : 0,
				(V1.vector4_f32[2] == V2.vector4_f32[2]) ? 0xFFFFFFFF : 0,
				(V1.vector4_f32[3] == V2.vector4_f32[3]) ? 0xFFFFFFFF : 0
			} } };

			return control.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_cmpeq_ps(V1, V2);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV EqualR(uint32_t* pCR, A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pCR != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			uint32_t ux = (V1.vector4_f32[0] == V2.vector4_f32[0]) ? 0xFFFFFFFFU : 0;
			uint32_t uy = (V1.vector4_f32[1] == V2.vector4_f32[1]) ? 0xFFFFFFFFU : 0;
			uint32_t uz = (V1.vector4_f32[2] == V2.vector4_f32[2]) ? 0xFFFFFFFFU : 0;
			uint32_t uw = (V1.vector4_f32[3] == V2.vector4_f32[3]) ? 0xFFFFFFFFU : 0;
			uint32_t CR = 0;

			if (ux & uy & uz & uw) // All elements are equal
			{
				CR = CRMASK_CR6TRUE;
			}
			else if (!(ux | uy | uz | uw)) // None of the elements are equal
			{
				CR = CRMASK_CR6FALSE;
			}
			*pCR = CR;

			VECTOR_U32 control = { { { ux, uy, uz, uw } } };
			return control;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmpeq_ps(V1, V2);
			uint32_t CR = 0;
			int iTest = _mm_movemask_ps(vTemp);
			if (iTest == 0xF) // All elements are equal
			{
				CR = CRMASK_CR6TRUE;
			}
			else if (!iTest) // None of the elements are equal
			{
				CR = CRMASK_CR6FALSE;
			}
			*pCR = CR;

			return vTemp;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV EqualInt(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 control = { { {
				(V1.vector4_u32[0] == V2.vector4_u32[0]) ? 0xFFFFFFFF : 0,
				(V1.vector4_u32[1] == V2.vector4_u32[1]) ? 0xFFFFFFFF : 0,
				(V1.vector4_u32[2] == V2.vector4_u32[2]) ? 0xFFFFFFFF : 0,
				(V1.vector4_u32[3] == V2.vector4_u32[3]) ? 0xFFFFFFFF : 0
			} } };

			return control.v;

#elif defined(_SSE2_INTRINSICS_)
			__m128i v = _mm_cmpeq_epi32(_mm_castps_si128(V1), _mm_castps_si128(V2));
			return _mm_castsi128_ps(v);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV EqualIntR(uint32_t* pCR, A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pCR != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR control = EqualInt(V1, V2);

			*pCR = 0;
			if (4EqualInt(control, TrueInt())) // All elements are equal
			{
				*pCR |= CRMASK_CR6TRUE;
			}
			else if (4EqualInt(control, FalseInt())) // None of the elements are equal
			{
				*pCR |= CRMASK_CR6FALSE;
			}

			return control;

#elif defined(_SSE2_INTRINSICS_)
			__m128i v = _mm_cmpeq_epi32(_mm_castps_si128(V1), _mm_castps_si128(V2));
			int iTemp = _mm_movemask_ps(_mm_castsi128_ps(v));
			uint32_t CR = 0;
			if (iTemp == 0xF) // All elements are equal
			{
				CR = CRMASK_CR6TRUE;
			}
			else if (!iTemp) // None of the elements are equal
			{
				CR = CRMASK_CR6FALSE;
			}
			*pCR = CR;

			return _mm_castsi128_ps(v);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV NearEqual(A_VECTOR V1, A_VECTOR V2, A_VECTOR epsilon) noexcept
		{
#if defined(_NO_INTRINSICS_)
			float deltaX = V1.vector4_f32[0] - V2.vector4_f32[0];
			float deltaY = V1.vector4_f32[1] - V2.vector4_f32[1];
			float deltaZ = V1.vector4_f32[2] - V2.vector4_f32[2];
			float deltaW = V1.vector4_f32[3] - V2.vector4_f32[3];

			deltaX = fabsf(deltaX);
			deltaY = fabsf(deltaY);
			deltaZ = fabsf(deltaZ);
			deltaW = fabsf(deltaW);

			VECTOR_U32 control = { { {
				(deltaX <= epsilon.vector4_f32[0]) ? 0xFFFFFFFFU : 0,
				(deltaY <= epsilon.vector4_f32[1]) ? 0xFFFFFFFFU : 0,
				(deltaZ <= epsilon.vector4_f32[2]) ? 0xFFFFFFFFU : 0,
				(deltaW <= epsilon.vector4_f32[3]) ? 0xFFFFFFFFU : 0
			} } };

			return control.v;

#elif defined(_SSE2_INTRINSICS_)
			// Get difference vector
			VECTOR vDelta = _mm_sub_ps(V1, V2);

			// Get absolute value of difference
			VECTOR vTemp = _mm_setzero_ps();
			vTemp = _mm_sub_ps(vTemp, vDelta);
			vTemp = _mm_max_ps(vTemp, vDelta);
			vTemp = _mm_cmple_ps(vTemp, epsilon);
			
			return vTemp;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV NotEqual(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 control = { { {
				(V1.vector4_f32[0] != V2.vector4_f32[0]) ? 0xFFFFFFFF : 0,
				(V1.vector4_f32[1] != V2.vector4_f32[1]) ? 0xFFFFFFFF : 0,
				(V1.vector4_f32[2] != V2.vector4_f32[2]) ? 0xFFFFFFFF : 0,
				(V1.vector4_f32[3] != V2.vector4_f32[3]) ? 0xFFFFFFFF : 0
			} } };

			return control.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_cmpneq_ps(V1, V2);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV NotEqualInt(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 control = { { {
				(V1.vector4_u32[0] != V2.vector4_u32[0]) ? 0xFFFFFFFFU : 0,
				(V1.vector4_u32[1] != V2.vector4_u32[1]) ? 0xFFFFFFFFU : 0,
				(V1.vector4_u32[2] != V2.vector4_u32[2]) ? 0xFFFFFFFFU : 0,
				(V1.vector4_u32[3] != V2.vector4_u32[3]) ? 0xFFFFFFFFU : 0
			} } };

			return control.v;

#elif defined(_SSE2_INTRINSICS_)
			__m128i v = _mm_cmpeq_epi32(_mm_castps_si128(V1), _mm_castps_si128(V2));
			return _mm_xor_ps(_mm_castsi128_ps(v), g_NegOneMask);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Greater(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 control = { { {
				(V1.vector4_f32[0] > V2.vector4_f32[0]) ? 0xFFFFFFFF : 0,
				(V1.vector4_f32[1] > V2.vector4_f32[1]) ? 0xFFFFFFFF : 0,
				(V1.vector4_f32[2] > V2.vector4_f32[2]) ? 0xFFFFFFFF : 0,
				(V1.vector4_f32[3] > V2.vector4_f32[3]) ? 0xFFFFFFFF : 0
			} } };

			return control.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_cmpgt_ps(V1, V2);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV GreaterR(uint32_t* pCR, A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pCR != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			uint32_t ux = (V1.vector4_f32[0] > V2.vector4_f32[0]) ? 0xFFFFFFFFU : 0;
			uint32_t uy = (V1.vector4_f32[1] > V2.vector4_f32[1]) ? 0xFFFFFFFFU : 0;
			uint32_t uz = (V1.vector4_f32[2] > V2.vector4_f32[2]) ? 0xFFFFFFFFU : 0;
			uint32_t uw = (V1.vector4_f32[3] > V2.vector4_f32[3]) ? 0xFFFFFFFFU : 0;
			uint32_t CR = 0;
			if (ux & uy & uz & uw) // All elements of V1 are greater than those of V2
			{
				CR = CRMASK_CR6TRUE;
			}
			else if (!(ux | uy | uz | uw)) // None of the elements of V1 are greater than those of V2
			{
				CR = CRMASK_CR6FALSE;
			}
			*pCR = CR;

			VECTOR_U32 control = { { { ux, uy, uz, uw } } };
			return control.v;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmpgt_ps(V1, V2);
			uint32_t CR = 0;
			int iTest = _mm_movemask_ps(vTemp);
			if (iTest == 0xF) // All elements of V1 are greater than those of V2
			{
				CR = CRMASK_CR6TRUE;
			}
			else if (!iTest) // None of the elements of V1 are greater than those of V2
			{
				CR = CRMASK_CR6FALSE;
			}
			*pCR = CR;

			return vTemp;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV GreaterOrEqual(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 control = { { {
				(V1.vector4_f32[0] >= V2.vector4_f32[0]) ? 0xFFFFFFFF : 0,
				(V1.vector4_f32[1] >= V2.vector4_f32[1]) ? 0xFFFFFFFF : 0,
				(V1.vector4_f32[2] >= V2.vector4_f32[2]) ? 0xFFFFFFFF : 0,
				(V1.vector4_f32[3] >= V2.vector4_f32[3]) ? 0xFFFFFFFF : 0
			} } };

			return control.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_cmpge_ps(V1, V2);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV GreaterOrEqualR(uint32_t* pCR, A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pCR != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			uint32_t ux = (V1.vector4_f32[0] >= V2.vector4_f32[0]) ? 0xFFFFFFFFU : 0;
			uint32_t uy = (V1.vector4_f32[1] >= V2.vector4_f32[1]) ? 0xFFFFFFFFU : 0;
			uint32_t uz = (V1.vector4_f32[2] >= V2.vector4_f32[2]) ? 0xFFFFFFFFU : 0;
			uint32_t uw = (V1.vector4_f32[3] >= V2.vector4_f32[3]) ? 0xFFFFFFFFU : 0;
			uint32_t CR = 0;
			if (ux & uy & uz & uw) // All elements of V1 are greater than or equal to those of V2
			{
				CR = CRMASK_CR6TRUE;
			}
			else if (!(ux | uy | uz | uw)) // None of the elements of V1 are greater than or equal to those of V2
			{
				CR = CRMASK_CR6FALSE;
			}
			*pCR = CR;

			VECTOR_U32 control = { { { ux, uy, uz, uw } } };
			return control.v;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = _mm_cmpge_ps(V1, V2);
			uint32_t CR = 0;
			int iTest = _mm_movemask_ps(vTemp);
			if (iTest == 0xF) // All the elements of V1 are greater than or equal to those of V2
			{
				CR = CRMASK_CR6TRUE;
			}
			else if (!iTest) // None of the elements of V1 are greater than or equal to those of V2
			{
				CR = CRMASK_CR6FALSE;
			}
			*pCR = CR;

			return vTemp;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Less(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 control = { { {
				(V1.vector4_f32[0] < V2.vector4_f32[0]) ? 0xFFFFFFFF : 0,
				(V1.vector4_f32[1] < V2.vector4_f32[1]) ? 0xFFFFFFFF : 0,
				(V1.vector4_f32[2] < V2.vector4_f32[2]) ? 0xFFFFFFFF : 0,
				(V1.vector4_f32[3] < V2.vector4_f32[3]) ? 0xFFFFFFFF : 0
			} } };

			return control.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_cmplt_ps(V1, V2);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV LessOrEqual(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 control = { { {
				(V1.vector4_f32[0] <= V2.vector4_f32[0]) ? 0xFFFFFFFF : 0,
				(V1.vector4_f32[1] <= V2.vector4_f32[1]) ? 0xFFFFFFFF : 0,
				(V1.vector4_f32[2] <= V2.vector4_f32[2]) ? 0xFFFFFFFF : 0,
				(V1.vector4_f32[3] <= V2.vector4_f32[3]) ? 0xFFFFFFFF : 0
			} } };

			return control.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_cmple_ps(V1, V2);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV InBounds(A_VECTOR v, A_VECTOR bounds) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 control = { { {
				(v.vector4_f32[0] <= bounds.vector4_f32[0] && v.vector4_f32[0] >= -bounds.vector4_f32[0]) ? 0xFFFFFFFF : 0,
				(v.vector4_f32[1] <= bounds.vector4_f32[1] && v.vector4_f32[1] >= -bounds.vector4_f32[1]) ? 0xFFFFFFFF : 0,
				(v.vector4_f32[2] <= bounds.vector4_f32[2] && v.vector4_f32[2] >= -bounds.vector4_f32[2]) ? 0xFFFFFFFF : 0,
				(v.vector4_f32[3] <= bounds.vector4_f32[3] && v.vector4_f32[3] >= -bounds.vector4_f32[3]) ? 0xFFFFFFFF : 0
			} } };

			return control.v;

#elif defined(_SSE2_INTRINSICS_)
			// Test if provided vector is less than or equal to bounds vector
			VECTOR vInUpperBound = _mm_cmple_ps(v, bounds);

			// Negate the bounds and cache in bInLowerBounds for next operation
			VECTOR vInLowerBound = _mm_mul_ps(bounds, g_NegativeOne);

			// Test if provided vector is greater than or equal to negative bounds vector (cached in vInLowerBound)
			vInLowerBound = _mm_cmpge_ps(v, vInLowerBound);

			// Blend the inequalities. A component is within bounds if it satisfies both conditions
			return _mm_and_ps(vInUpperBound, vInLowerBound);
#endif
		}

		_Use_decl_annotations_
		FORCE_INLINE VECTOR VEC_CALLCONV InBoundsR(uint32_t* pCR, A_VECTOR v, A_VECTOR bounds) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			assert(pCR != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
			uint32_t ux = (v.vector4_f32[0] <= bounds.vector4_f32[0] && v.vector4_f32[0] >= -bounds.vector4_f32[0]) ? 0xFFFFFFFFU : 0;
			uint32_t uy = (v.vector4_f32[1] <= bounds.vector4_f32[1] && v.vector4_f32[1] >= -bounds.vector4_f32[1]) ? 0xFFFFFFFFU : 0;
			uint32_t uz = (v.vector4_f32[2] <= bounds.vector4_f32[2] && v.vector4_f32[2] >= -bounds.vector4_f32[2]) ? 0xFFFFFFFFU : 0;
			uint32_t uw = (v.vector4_f32[3] <= bounds.vector4_f32[3] && v.vector4_f32[3] >= -bounds.vector4_f32[3]) ? 0xFFFFFFFFU : 0;
			uint32_t CR = 0;
			if (ux & uy & uz & uw) // All elements of v are within bounds
			{
				CR = CRMASK_CR6BOUNDS;
			}
			*pCR = CR;

			VECTOR_U32 control = { { { ux, uy, uz, uw } } };
			return control.v;

#elif defined(_SSE2_INTRINSICS_)
			// Test if provided vector is less than or equal to bounds vector
			VECTOR vInUpperBound = _mm_cmple_ps(v, bounds);

			// Negate the bounds and cache in bInLowerBounds for next operation
			VECTOR vInLowerBound = _mm_mul_ps(bounds, g_NegativeOne);

			// Test if provided vector is greater than or equal to negative bounds vector (cached in vInLowerBound)
			vInLowerBound = _mm_cmpge_ps(v, vInLowerBound);

			// Blend the inequalities. A component is within bounds if it satisfies both conditions
			vInUpperBound = _mm_and_ps(vInUpperBound, vInLowerBound);

			uint32_t CR = 0;
			if (_mm_movemask_ps(vInUpperBound) == 0xF) // All elements of v are within bounds
			{
				CR = CRMASK_CR6BOUNDS;
			}
			*pCR = CR;

			return vInUpperBound;
#endif
		}

#if !defined(_NO_INTRINSICS_) && defined(_MSC_VER) && !defined(__INTEL_COMPILER)
#pragma float_control(push)
#pragma float_control(precise, on)
#endif

		FORCE_INLINE VECTOR VEC_CALLCONV IsNaN(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 control = { { {
				ISNAN(v.vector4_f32[0]) ? 0xFFFFFFFFU : 0,
				ISNAN(v.vector4_f32[1]) ? 0xFFFFFFFFU : 0,
				ISNAN(v.vector4_f32[2]) ? 0xFFFFFFFFU : 0,
				ISNAN(v.vector4_f32[3]) ? 0xFFFFFFFFU : 0
			} } };

			return control.v;

#elif defined(_SSE2_INTRINSICS_)
#if defined(__clang__) && defined(__FINITE_MATH_ONLY__)
			ALIGNED(16) float tmp[4];
			_mm_store_ps(tmp, v);
			VECTOR_U32 vResult = { { {
				isnan(tmp[0]) ? 0xFFFFFFFF : 0,
				isnan(tmp[1]) ? 0xFFFFFFFF : 0,
				isnan(tmp[2]) ? 0xFFFFFFFF : 0,
				isnan(tmp[3]) ? 0xFFFFFFFF : 0
			} } };

			return vResult.v;

#else
			// Test against itself; NaN is always not equal
			return _mm_cmpneq_ps(v, v);
#endif
#endif
		}

#if !defined(_NO_INTRINSICS_) && defined(_MSC_VER) && !defined(__INTEL_COMPILER)
#pragma float_control(pop)
#endif

		FORCE_INLINE VECTOR VEC_CALLCONV IsInfinite(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 control = { { {
				ISINF(v.vector4_f32[0]) ? 0xFFFFFFFFU : 0,
				ISINF(v.vector4_f32[1]) ? 0xFFFFFFFFU : 0,
				ISINF(v.vector4_f32[2]) ? 0xFFFFFFFFU : 0,
				ISINF(v.vector4_f32[3]) ? 0xFFFFFFFFU : 0
			} } };

			return control.v;

#elif defined(_SSE2_INTRINSICS_)
			// Mask off the sign bit
			__m128 vTemp = _mm_and_ps(v, g_AbsMask);

			// Compare to infinity
			vTemp = _mm_cmpeq_ps(vTemp, g_Infinity);

			return vTemp;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Min(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 result = { { {
				(V1.vector4_f32[0] < V2.vector4_f32[0]) ? V1.vector4_f32[0] : V2.vector4_f32[0],
				(V1.vector4_f32[1] < V2.vector4_f32[1]) ? V1.vector4_f32[1] : V2.vector4_f32[1],
				(V1.vector4_f32[2] < V2.vector4_f32[2]) ? V1.vector4_f32[2] : V2.vector4_f32[2],
				(V1.vector4_f32[3] < V2.vector4_f32[3]) ? V1.vector4_f32[3] : V2.vector4_f32[3]
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_min_ps(V1, V2);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Max(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 result = { { {
				(V1.vector4_f32[0] > V2.vector4_f32[0]) ? V1.vector4_f32[0] : V2.vector4_f32[0],
				(V1.vector4_f32[1] > V2.vector4_f32[1]) ? V1.vector4_f32[1] : V2.vector4_f32[1],
				(V1.vector4_f32[2] > V2.vector4_f32[2]) ? V1.vector4_f32[2] : V2.vector4_f32[2],
				(V1.vector4_f32[3] > V2.vector4_f32[3]) ? V1.vector4_f32[3] : V2.vector4_f32[3]
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_max_ps(V1, V2);
#endif
		}

		namespace
		{
			// Round to nearest (even) a.k.a. banker's rounding
			FORCE_INLINE float round_to_nearest(float x) noexcept
			{
				float i = floorf(x);
				x -= i;
				if (x < 0.5f)
					return i;
				if (x > 0.5f)
					return i + 1.0f;

				float int_part;
				(void)modff(i / 2.0f, &int_part);
				if ((2.0f * int_part) == i)
					return i;

				return i + 1.0f;
			}
		}

#if !defined(_NO_INTRINSICS_) && defined(_MSC_VER) && !defined(__INTEL_COMPILER)
#pragma float_control(push)
#pragma float_control(precise, on)
#endif

		FORCE_INLINE VECTOR VEC_CALLCONV Round(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 result = { { {
				round_to_nearest(v.vector4_f32[0]),
				round_to_nearest(v.vector4_f32[1]),
				round_to_nearest(v.vector4_f32[2]),
				round_to_nearest(v.vector4_f32[3])
			} } };

			return result.v;

#elif defined(_SSE4_INTRINSICS_)
			return _mm_round_ps(v, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);

#elif defined(_SSE2_INTRINSICS_)
			__m128 sign = _mm_and_ps(v, g_NegativeZero);
			__m128 sMagic = _mm_or_ps(g_NoFraction, sign);
			__m128 R1 = _mm_and_ps(v, sMagic);
			R1 = _mm_sub_ps(R1, sMagic);
			__m128 R2 = _mm_and_ps(v, g_AbsMask);
			__m128 mask = _mm_cmple_ps(R2, g_NoFraction);
			R2 = _mm_andnot_ps(mask, v);
			R1 = _mm_and_ps(R1, mask);

			VECTOR vResult = _mm_xor_ps(R1, R2);
			return vResult;
#endif
		}

#if !defined(_NO_INTRINSICS_) && defined(_MSC_VER) && !defined(__INTEL_COMPILER)
#pragma float_control(pop)
#endif

		FORCE_INLINE VECTOR VEC_CALLCONV Truncate(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR result;

			for (uint32_t i = 0; i < 4; i++)
			{
				if (ISNAN(v.vector4_f32[i]))
				{
					result.vector4_u32[i] = 0x7FC00000;
				}
				else if (fabsf(v.vector4_f32[i]) < 8388608.0f)
				{
					result.vector4_f32[i] = static_cast<float>(static_cast<int32_t>(v.vector4_f32[i]));
				}
				else
				{
					result.vector4_f32[i] = v.vector4_f32[i];
				}
			}

			return result

#elif defined(_SSE4_INTRINSICS_)
			return _mm_round_ps(v, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);

#elif defined(_SSE2_INTRINSICS_)
			__m128i vTest = _mm_and_si128(_mm_castps_si128(v), g_AbsMask);

			// Test for elements greater than 8388608 (All floats with no fractional, NAN and INF)
			vTest = _mm_cmplt_epi32(vTest, g_NoFraction);

			// Convert to int and back to float to truncate
			__m128i vInt = _mm_cvttps_epi32(v);

			VECTOR vResult = _mm_cvtepi32_ps(vInt);

			// All elements less than 8388608 will use the result from round to int
			vResult = _mm_and_ps(vResult, _mm_castsi128_ps(vTest));

			// All elements greater than 8388608 will carry on the original value
			vTest = _mm_andnot_si128(vTest, _mm_castps_si128(v));
			vResult = _mm_or_ps(vResult, _mm_castsi128_ps(vTest));

			return vResult;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Floor(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 result = { { {
				floorf(v.vector4_f32[0]),
				floorf(v.vector4_f32[1]),
				floorf(v.vector4_f32[2]),
				floorf(v.vector4_f32[3])
			} } };

			return result.v;

#elif defined(_SSE4_INTRINSICS_)
			return _mm_floor_ps(v);

#elif defined(_SSE2_INTRINSICS_)
			__m128i vTest = _mm_and_si128(_mm_castps_si128(v), g_AbsMask);
			vTest = _mm_cmplt_epi32(vTest, g_NoFraction);

			// Cast elements to int to truncate
			__m128i vInt = _mm_cvttps_epi32(v);

			VECTOR vResult = _mm_cvtepi32_ps(vInt);
			__m128 vLarger = _mm_cmpgt_ps(vResult, v);

			// 0 -> 0, 0xFFFFFFFF -> -1.0f
			vLarger = _mm_cvtepi32_ps(_mm_castps_si128(vLarger));
			vResult = _mm_and_ps(vResult, vLarger);

			// All elements less than 8388608 will use the result from round to int
			vResult = _mm_and_ps(vResult, _mm_castsi128_ps(vTest));

			// All elements greater than 8388608 will carry on the original value
			vTest = _mm_andnot_si128(vTest, _mm_castps_si128(v));
			vResult = _mm_or_ps(vResult, _mm_castsi128_ps(vTest));

			return vResult;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Ceiling(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 result = { { {
				ceilf(v.vector4_f32[0]),
				ceilf(v.vector4_f32[1]),
				ceilf(v.vector4_f32[2]),
				ceilf(v.vector4_f32[3])
			} } };
			return result.v;

#elif defined(_SSE4_INTRINSICS_)
			return _mm_ceil_ps(v);

#elif defined(_SSE2_INTRINSICS_)
			__m128i vTest = _mm_and_si128(_mm_castps_si128(v), g_AbsMask);
			vTest = _mm_cmplt_epi32(vTest, g_NoFraction);

			// Cast elements to int to truncate
			__m128i vInt = _mm_cvttps_epi32(v);

			VECTOR vResult = _mm_cvtepi32_ps(vInt);

			__m128 vSmaller = _mm_cmplt_ps(vResult, v);

			// 0 -> 0, 0xFFFFFFFF -> -1.0f
			vSmaller = _mm_cvtepi32_ps(_mm_castps_si128(vSmaller));
			vResult = _mm_sub_ps(vResult, vSmaller);

			// All elements less than 8388608 will use the result from round to int
			vResult = _mm_and_ps(vResult, _mm_castsi128_ps(vTest));

			// All elements greater than 8388608 will carry on the original value
			vTest = _mm_andnot_si128(vTest, _mm_castps_si128(v));
			vResult = _mm_or_ps(vResult, _mm_castsi128_ps(vTest));

			return vResult;
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Clamp(A_VECTOR v, A_VECTOR min, A_VECTOR max) noexcept
		{
#if defined(DEBUG) || defined(_DEBUG)
			//assert(Vector4LessOrEqual(min, max));
#endif

#if defined(_NO_INTRINSICS_)
			VECTOR result = Max(min, v);
			
			return Min(max, result);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vResult = _mm_max_ps(min, v);
			return _mm_min_ps(max, vResult);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Saturate(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			const VECTOR zero = Zero();
			return Clamp(v, zero, g_One.v);

#elif defined(_SSE2_INTRINSICS_)
			// Set elements < 0 to 0
			VECTOR vResult = _mm_max_ps(v, g_Zero);

			// Set elements > 1 to 1
			return _mm_min_ps(vResult, g_One);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV AndInt(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 result = { { {
				V1.vector4_u32[0] & V2.vector4_u32[0],
				V1.vector4_u32[1] & V2.vector4_u32[1],
				V1.vector4_u32[2] & V2.vector4_u32[2],
				V1.vector4_u32[3] & V2.vector4_u32[3]
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_and_ps(V1, V2);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV AndCInt(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 result = { { {
				V1.vector4_u32[0] & ~V2.vector4_u32[0],
				V1.vector4_u32[1] & ~V2.vector4_u32[1],
				V1.vector4_u32[2] & ~V2.vector4_u32[2],
				V1.vector4_u32[3] & ~V2.vector4_u32[3]
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_andnot_ps(V1, V2);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV OrInt(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 result = { { {
				V1.vector4_u32[0] | V2.vector4_u32[0],
				V1.vector4_u32[1] | V2.vector4_u32[1],
				V1.vector4_u32[2] | V2.vector4_u32[2],
				V1.vector4_u32[3] | V2.vector4_u32[3]
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_or_ps(V1, V2);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV NorInt(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 result = { { {
				~(V1.vector4_u32[0] | V2.vector4_u32[0]),
				~(V1.vector4_u32[1] | V2.vector4_u32[1]),
				~(V1.vector4_u32[2] | V2.vector4_u32[2]),
				~(V1.vector4_u32[3] | V2.vector4_u32[3])
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vResult = _mm_or_ps(V1, V2);
			return _mm_andnot_ps(vResult, g_NegOneMask);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV XorInt(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_U32 result = { { {
				V1.vector4_u32[0] ^ V2.vector4_u32[0],
				V1.vector4_u32[1] ^ V2.vector4_u32[1],
				V1.vector4_u32[2] ^ V2.vector4_u32[2],
				V1.vector4_u32[3] ^ V2.vector4_u32[3]
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_xor_ps(V1, V2);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Negate(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 result = { { {
				-v.vector4_f32[0],
				-v.vector4_f32[1],
				-v.vector4_f32[2],
				-v.vector4_f32[3]
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR z = _mm_setzero_ps();
			return _mm_sub_ps(z, v);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Add(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 result = { { {
				V1.vector4_f32[0] + V2.vector4_f32[0],
				V1.vector4_f32[1] + V2.vector4_f32[1],
				V1.vector4_f32[2] + V2.vector4_f32[2],
				V1.vector4_f32[3] + V2.vector4_f32[3]
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_add_ps(V1, V2);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Sum(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 result;
			result.f[0] =
			result.f[1] =
			result.f[2] =
			result.f[3] = v.vector4_f32[0] + v.vector4_f32[1] + v.vector4_f32[2] + v.vector4_f32[3];

			return result.v;

#elif defined(_SSE3_INTRINSICS_)
			VECTOR vTemp = _mm_hadd_ps(v, v);
			return _mm_hadd_ps(vTemp, vTemp);

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vTemp = PERMUTE_PS(v, _MM_SHUFFLE(2, 3, 0, 1));
			VECTOR vTemp2 = _mm_add_ps(v, vTemp);
			vTemp = PERMUTE_PS(vTemp2, _MM_SHUFFLE(1, 0, 3, 2));
			return _mm_add_ps(vTemp, vTemp2);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV AddAngles(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			const VECTOR zero = Zero();

			VECTOR result = Add(V1, V2);

			VECTOR mask = Less(result, g_NegativePi.v);
			VECTOR offset = Select(offset, g_NegativeTwoPi.v, mask);

			return Add(result, offset);

#elif defined(_SSE2_INTRINSICS_)
			// Adjust the angles
			VECTOR vResult = _mm_and_ps(V1, V2);

			// Determine which elements are less than Pi
			VECTOR vOffset = _mm_cmplt_ps(vResult, g_NegativePi);
			vOffset = _mm_and_ps(vOffset, g_TwoPi);

			// Add 2Pi to all entries that are less than -Pi
			vResult = _mm_add_ps(vResult, vOffset);

			// Determine which elements are greater than or equal to Pi
			vOffset = _mm_cmpge_ps(vResult, g_Pi);
			vOffset = _mm_and_ps(vOffset, g_TwoPi);

			// Subtract 2Pi from all entries that are greater than Pi and return the result
			return _mm_sub_ps(vResult, vOffset);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Subtract(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 result = { { {
				V1.vector4_f32[0] - V2.vector4_f32[0],
				V1.vector4_f32[1] - V2.vector4_f32[1],
				V1.vector4_f32[2] - V2.vector4_f32[2],
				V1.vector4_f32[3] - V2.vector4_f32[3]
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_sub_ps(V1, V2);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SubtractAngles(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			const VECTOR zero = Zero();

			VECTOR result = Subtract(V1, V2);

			VECTOR mask = Less(result, g_NegativePi.v);
			VECTOR offset = Select(zero, g_TwoPi.v, mask);

			mask = GreaterOrEqual(result, g_Pi.v);
			offset = Select(offset, g_NegativeTwoPi.v, mask);

			return Add(result, offset);

#elif defined(_SSE2_INTRINSICS_)
			// Adjust the angles
			VECTOR vResult = _mm_sub_ps(V1, V2);

			// Determine which elements are less than Pi
			VECTOR vOffset = _mm_cmplt_ps(vResult, g_NegativePi);
			vOffset = _mm_and_ps(vOffset, g_TwoPi);

			// Add 2Pi to all entries that are less than -Pi
			vResult = _mm_and_ps(vResult, vOffset);

			// Determine which elements are greater than or equal to Pi
			vOffset = _mm_cmpge_ps(vResult, g_Pi);
			vOffset = _mm_and_ps(vOffset, g_TwoPi);

			// Subtract 2Pi from all entries that are grater than Pi and return result
			return _mm_sub_ps(vResult, vOffset);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Multiply(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 result = { { {
				V1.vector4_f32[0] * V2.vector4_f32[0],
				V1.vector4_f32[1] * V2.vector4_f32[1],
				V1.vector4_f32[2] * V2.vector4_f32[2],
				V1.vector4_f32[3] * V2.vector4_f32[3]
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_mul_ps(V1, V2);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV MultiplyAdd(A_VECTOR V1, A_VECTOR V2, A_VECTOR V3) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 result = { { {
				V1.vector4_f32[0] * V2.vector4_f32[0] + V3.vector4_f32[0],
				V1.vector4_f32[1] * V2.vector4_f32[1] + V3.vector4_f32[1],
				V1.vector4_f32[2] * V2.vector4_f32[2] + V3.vector4_f32[2],
				V1.vector4_f32[3] * V2.vector4_f32[3] + V3.vector4_f32[3]
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			FMADD_PS(V1, V2, V3);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Divide(A_VECTOR V1, A_VECTOR V2) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 result = { { {
				V1.vector4_f32[0] / V2.vector4_f32[0],
				V1.vector4_f32[1] / V2.vector4_f32[1],
				V1.vector4_f32[2] / V2.vector4_f32[2],
				V1.vector4_f32[3] / V2.vector4_f32[3]
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_div_ps(V1, V2);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV NegativeMultiplySubtract(A_VECTOR V1, A_VECTOR V2, A_VECTOR V3) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 result = { { {
				V3.vector4_f32[0] - (V1.vector4_f32[0] * V2.vector4_f32[0]),
				V3.vector4_f32[1] - (V1.vector4_f32[1] * V2.vector4_f32[1]),
				V3.vector4_f32[2] - (V1.vector4_f32[2] * V2.vector4_f32[2]),
				V3.vector4_f32[3] - (V1.vector4_f32[3] * V2.vector4_f32[3])
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			FNMADD_PS(V1, V2, V3);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Scale(A_VECTOR v, float scaleFactor) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 result = { { {
				v.vector4_f32[0] * scaleFactor,
				v.vector4_f32[1] * scaleFactor,
				v.vector4_f32[2] * scaleFactor,
				v.vector4_f32[3] * scaleFactor
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vResult = _mm_set_ps1(scaleFactor);
			return _mm_mul_ps(vResult, v);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV ReciprocalEst(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 result = { { {
				1.0f / v.vector4_f32[0],
				1.0f / v.vector4_f32[1],
				1.0f / v.vector4_f32[2],
				1.0f / v.vector4_f32[3]
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_rcp_ps(v);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Reciprocal(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 result = { { {
				1.0f / v.vector4_f32[0],
				1.0f / v.vector4_f32[1],
				1.0f / v.vector4_f32[2],
				1.0f / v.vector4_f32[3]
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_div_ps(g_One, v);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV SqrtEst(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 result = { { {
				sqrtf(v.vector4_f32[0]),
				sqrtf(v.vector4_f32[1]),
				sqrtf(v.vector4_f32[2]),
				sqrtf(v.vector4_f32[3])
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_sqrt_ps(v);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Sqrt(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 result = { { {
				sqrtf(v.vector4_f32[0]),
				sqrtf(v.vector4_f32[1]),
				sqrtf(v.vector4_f32[2]),
				sqrtf(v.vector4_f32[3])
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_sqrt_ps(v);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV ReciprocalSqrtEst(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 result = { { {
				1.0f / sqrtf(v.vector4_f32[0]),
				1.0f / sqrtf(v.vector4_f32[1]),
				1.0f / sqrtf(v.vector4_f32[2]),
				1.0f / sqrtf(v.vector4_f32[3])
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			return _mm_rsqrt_ps(v);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV ReciprocalSqrt(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 result = { { {
				1.0f / sqrtf(v.vector4_f32[0]),
				1.0f / sqrtf(v.vector4_f32[1]),
				1.0f / sqrtf(v.vector4_f32[2]),
				1.0f / sqrtf(v.vector4_f32[3])
			} } };

			return result.v;

#elif defined(_SSE2_INTRINSICS_)
			VECTOR vResult = _mm_sqrt_ps(v);
			return _mm_div_ps(g_One, vResult);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Exp2(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 result = { { {
				exp2f(V.vector4_f32[0]),
				exp2f(V.vector4_f32[1]),
				exp2f(V.vector4_f32[2]),
				exp2f(V.vector4_f32[3])
			} } };

			return result.v;

#elif defined(_SVML_INTRINSICS_)
			return _mm_exp2_ps(v);

#elif defined(_SSE2_INTRINSICS_)
			__m128i vTrunci = _mm_cvttps_epi32(v);
			__m128 vTrunc = _mm_cvtepi32_ps(vTrunci);
			__m128 vY = _mm_sub_ps(v, vTrunc);

			__m128 vPoly = FMADD_PS(g_ExpEst7, vY, g_ExpEst6);
			vPoly = FMADD_PS(vPoly, vY, g_ExpEst5);
			vPoly = FMADD_PS(vPoly, vY, g_ExpEst4);
			vPoly = FMADD_PS(vPoly, vY, g_ExpEst3);
			vPoly = FMADD_PS(vPoly, vY, g_ExpEst2);
			vPoly = FMADD_PS(vPoly, vY, g_ExpEst1);
			vPoly = FMADD_PS(vPoly, vY, g_One);

			__m128i vBiased = _mm_add_epi32(vTrunci, g_ExponentBias);
			vBiased = _mm_slli_epi32(vBiased, 23);
			__m128 vResult0 = _mm_div_ps(_mm_castsi128_ps(vBiased), vPoly);

			vBiased = _mm_add_epi32(vTrunci, g_253);
			vBiased = _mm_slli_epi32(vBiased, 23);
			__m128 vResult1 = _mm_div_ps(_mm_castsi128_ps(vBiased), vPoly);
			vResult1 = _mm_mul_ps(g_MinNormal.v, vResult1);

			// Use selection to handle the cases
			// if (v is NaN) -> QNaN;
			// else if (v sign bit set)
			//		if(v > -150)
			//			if(v.exponent < -126) -> vResult1
			//			else -> vResult0
			//		else -> +0
			// else
			//		if(v < 128) -> vResult0
			//		else -> +inf

			__m128i vInputi = _mm_castps_si128(v);

			__m128i vComp = _mm_cmplt_epi32(vInputi, g_Bin128);
			__m128i vSelect0 = _mm_and_si128(vComp, _mm_castps_si128(vResult0));
			__m128i vSelect1 = _mm_andnot_si128(vComp, g_Infinity);
			__m128i vResult2 = _mm_or_si128(vSelect0, vSelect1);

			vComp = _mm_cmplt_epi32(vTrunci, g_SubnormalExponent);
			vSelect1 = _mm_and_si128(vComp, _mm_castps_si128(vResult1));
			vSelect0 = _mm_andnot_si128(vComp, _mm_castps_si128(vResult0));
			__m128i vResult3 = _mm_or_si128(vSelect0, vSelect1);

			vComp = _mm_cmplt_epi32(vInputi, g_BinNeg150);
			vSelect0 = _mm_and_si128(vComp, vResult3);
			vSelect1 = _mm_andnot_si128(vComp, g_Zero);
			__m128i vResult4 = _mm_or_si128(vSelect0, vSelect1);

			__m128i vSign = _mm_and_si128(vInputi, g_NegativeZero);
			vComp = _mm_cmpeq_epi32(vSign, g_NegativeZero);
			vSelect0 = _mm_and_si128(vComp, vResult4);
			vSelect1 = _mm_andnot_si128(vComp, vResult2);
			__m128i vResult5 = _mm_or_si128(vSelect0, vSelect1);

			__m128i vT0 = _mm_and_si128(vInputi, g_QNaNTest);
			__m128i vT1 = _mm_and_si128(vInputi, g_Infinity);
			vT0 = _mm_cmpeq_epi32(vT0, g_Zero);
			vT1 = _mm_cmpeq_epi32(vT1, g_Infinity);
			__m128i vIsNaN = _mm_andnot_si128(vT0, vT1);

			vSelect0 = _mm_and_si128(vIsNaN, g_QNaN);
			vSelect1 = _mm_andnot_si128(vIsNaN, vResult5);
			__m128i vResult = _mm_or_si128(vSelect0, vSelect1);

			return _mm_castsi128_ps(vResult);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Exp10(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 result = { { {
				powf(10.0f, v.vector4_f32[0]),
				powf(10.0f, v.vector4_f32[1]),
				powf(10.0f, v.vector4_f32[2]),
				powf(10.0f, v.vector4_f32[3])
			} } };

			return result.v;

#elif defined(_SVML_INTRINSICS_)
			return _mm_exp10_ps(v);

#else
			VECTOR vTen = Multiply(g_Lg10, v);
			return Exp2(vTen);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV ExpE(A_VECTOR v) noexcept
		{
#if defined(_NO_INTRINSICS_)
			VECTOR_F32 result = { { {
				expf(v.vector4_f32[0]),
				expf(v.vector4_f32[1]),
				expf(v.vector4_f32[2]),
				expf(v.vector4_f32[3])
			} } };

			return result.v;

#elif defined(_SVML_INTRINSICS_)
			return _mm_exp_ps(v);

#else
			// exp(v) = exp2(vin*log2(e))
			VECTOR vE = Multiply(g_LgE, v);
			return Exp2(vE);
#endif
		}

		FORCE_INLINE VECTOR VEC_CALLCONV Exp(A_VECTOR v) noexcept
		{
			return Exp2(v);
		}

#if defined(_SSE2_INTRINSICS_)
		namespace
		{
			FORCE_INLINE __m128i multi_sll_epi32(__m128i value, __m128i count) noexcept
			{
				__m128i v = _mm_shuffle_epi32(value, _MM_SHUFFLE(0, 0, 0, 0));
				__m128i c = _mm_shuffle_epi32(count, _MM_SHUFFLE(0, 0, 0, 0));
				c = _mm_and_si128(c, g_MaskX);
				__m128i r0 = _mm_sll_epi32(v, c);

				v = _mm_shuffle_epi32(value, _MM_SHUFFLE(1, 1, 1, 1));
				c = _mm_shuffle_epi32(count, _MM_SHUFFLE(1, 1, 1, 1));
				c = _mm_and_si128(c, g_MaskX);
				__m128i r1 = _mm_sll_epi32(v, c);

				v = _mm_shuffle_epi32(value, _MM_SHUFFLE(2, 2, 2, 2));
				c = _mm_shuffle_epi32(count, _MM_SHUFFLE(2, 2, 2, 2));
				c = _mm_and_si128(c, g_MaskX);
				__m128i r2 = _mm_sll_epi32(v, c);

				v = _mm_shuffle_epi32(value, _MM_SHUFFLE(3, 3, 3, 3));
				c = _mm_shuffle_epi32(count, _MM_SHUFFLE(3, 3, 3, 3));
				c = _mm_and_si128(c, g_MaskX);
				__m128i r3 = _mm_sll_epi32(v, c);

				// (r0,r0,r1,r1)
				__m128 r01 = _mm_shuffle_ps(_mm_castsi128_ps(r0), _mm_castsi128_ps(r1), _MM_SHUFFLE(0, 0, 0, 0));
				// (r2,r2,r3,r3)
				__m128 r23 = _mm_shuffle_ps(_mm_castsi128_ps(r2), _mm_castsi128_ps(r3), _MM_SHUFFLE(0, 0, 0, 0));
				// (r0,r1,r2,r3)
				__m128 result = _mm_shuffle_ps(r01, r23, _MM_SHUFFLE(2, 0, 2, 0));
				return _mm_castps_si128(result);
			}

			FORCE_INLINE __m128i multi_srl_epi32(__m128i value, __m128i count) noexcept
			{
				__m128i v = _mm_shuffle_epi32(value, _MM_SHUFFLE(0, 0, 0, 0));
				__m128i c = _mm_shuffle_epi32(count, _MM_SHUFFLE(0, 0, 0, 0));
				c = _mm_and_si128(c, g_MaskX);
				__m128i r0 = _mm_srl_epi32(v, c);

				v = _mm_shuffle_epi32(value, _MM_SHUFFLE(1, 1, 1, 1));
				c = _mm_shuffle_epi32(count, _MM_SHUFFLE(1, 1, 1, 1));
				c = _mm_and_si128(c, g_MaskX);
				__m128i r1 = _mm_srl_epi32(v, c);

				v = _mm_shuffle_epi32(value, _MM_SHUFFLE(2, 2, 2, 2));
				c = _mm_shuffle_epi32(count, _MM_SHUFFLE(2, 2, 2, 2));
				c = _mm_and_si128(c, g_MaskX);
				__m128i r2 = _mm_srl_epi32(v, c);

				v = _mm_shuffle_epi32(value, _MM_SHUFFLE(3, 3, 3, 3));
				c = _mm_shuffle_epi32(count, _MM_SHUFFLE(3, 3, 3, 3));
				c = _mm_and_si128(c, g_MaskX);
				__m128i r3 = _mm_srl_epi32(v, c);

				// (r0,r0,r1,r1)
				__m128 r01 = _mm_shuffle_ps(_mm_castsi128_ps(r0), _mm_castsi128_ps(r1), _MM_SHUFFLE(0, 0, 0, 0));
				// (r2,r2,r3,r3)
				__m128 r23 = _mm_shuffle_ps(_mm_castsi128_ps(r2), _mm_castsi128_ps(r3), _MM_SHUFFLE(0, 0, 0, 0));
				// (r0,r1,r2,r3)
				__m128 result = _mm_shuffle_ps(r01, r23, _MM_SHUFFLE(2, 0, 2, 0));
				return _mm_castps_si128(result);
			}

			FORCE_INLINE __m128i GetLeadingBit(const __m128i value) noexcept
			{
				static const VECTOR_I32 g_XM0000FFFF = { { { 0x0000FFFF, 0x0000FFFF, 0x0000FFFF, 0x0000FFFF } } };
				static const VECTOR_I32 g_XM000000FF = { { { 0x000000FF, 0x000000FF, 0x000000FF, 0x000000FF } } };
				static const VECTOR_I32 g_XM0000000F = { { { 0x0000000F, 0x0000000F, 0x0000000F, 0x0000000F } } };
				static const VECTOR_I32 g_XM00000003 = { { { 0x00000003, 0x00000003, 0x00000003, 0x00000003 } } };

				__m128i v = value, r, c, b, s;

				c = _mm_cmpgt_epi32(v, g_XM0000FFFF);   // c = (v > 0xFFFF)
				b = _mm_srli_epi32(c, 31);              // b = (c ? 1 : 0)
				r = _mm_slli_epi32(b, 4);               // r = (b << 4)
				v = multi_srl_epi32(v, r);              // v = (v >> r)

				c = _mm_cmpgt_epi32(v, g_XM000000FF);   // c = (v > 0xFF)
				b = _mm_srli_epi32(c, 31);              // b = (c ? 1 : 0)
				s = _mm_slli_epi32(b, 3);               // s = (b << 3)
				v = multi_srl_epi32(v, s);              // v = (v >> s)
				r = _mm_or_si128(r, s);                 // r = (r | s)

				c = _mm_cmpgt_epi32(v, g_XM0000000F);   // c = (v > 0xF)
				b = _mm_srli_epi32(c, 31);              // b = (c ? 1 : 0)
				s = _mm_slli_epi32(b, 2);               // s = (b << 2)
				v = multi_srl_epi32(v, s);              // v = (v >> s)
				r = _mm_or_si128(r, s);                 // r = (r | s)

				c = _mm_cmpgt_epi32(v, g_XM00000003);   // c = (v > 0x3)
				b = _mm_srli_epi32(c, 31);              // b = (c ? 1 : 0)
				s = _mm_slli_epi32(b, 1);               // s = (b << 1)
				v = multi_srl_epi32(v, s);              // v = (v >> s)
				r = _mm_or_si128(r, s);                 // r = (r | s)

				s = _mm_srli_epi32(v, 1);
				r = _mm_or_si128(r, s);
				return r;
			}
		}
#endif // _SSE2_INTRINSICS_

		FORCE_INLINE VECTOR VEC_CALLCONV Log2(A_VECTOR v) noexcept
		{

		}

#ifdef __clang__
#pragma clang diagnostic pop
#endif
#ifdef _PREFAST_
#pragma prefast(pop)
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif
	}
}

#endif // !ULTREALITY_MATH_SSE2_CONFIG_INL
