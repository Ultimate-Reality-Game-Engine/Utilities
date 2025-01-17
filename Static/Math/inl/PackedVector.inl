#ifndef ULTREALITY_MATH_PACKED_VECTOR_INL
#define ULTREALITY_MATH_PACKED_VECTOR_INL

#if defined(__GNUC__) or defined(__clang__)
#define FORCE_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline
#endif

namespace UltReality::Math
{
    namespace PackedVector
    {
        FORCE_INLINE float ConvertHalfToFloat(HALF value) noexcept
        {
#if defined(_F16C_INTRINSICS_) && !defined(_NO_INTRINSICS_)
            __m128i V1 = _mm_cvtsi32_si128(static_cast<int>(value));
            __m128 V2 = _mm_cvtph_ps(V1);
            return _mm_cvtss_f32(V2);

#else
            uint32_t Mantissa = static_cast<uint32_t>(value & 0x03FF);

            uint32_t Exponent = (value & 0x7C00);
            if (Exponent == 0x7C00) // INF/NAN
            {
                Exponent = 0x8f;
            }
            else if (Exponent != 0) // The value is normalized
            {
                Exponent = static_cast<uint32_t>((static_cast<int>(value) >> 10) & 0x1F);
            }
            else if (Mantissa != 0) // The value is denormalized
            {
                // Normalize the value in the resulting float
                Exponent = 1;

                do
                {
                    Exponent--;
                    Mantissa <<= 1;
                } while ((Mantissa & 0x0400) == 0);

                Mantissa &= 0x03FF;
            }
            else                        // The value is zero
            {
                Exponent = static_cast<uint32_t>(-112);
            }

            uint32_t Result =
                ((static_cast<uint32_t>(value) & 0x8000) << 16) // Sign
                | ((Exponent + 112) << 23)                      // Exponent
                | (Mantissa << 13);                             // Mantissa

            return reinterpret_cast<float*>(&Result)[0];
#endif // !_XM_F16C_INTRINSICS_
        }

#ifdef _PREFAST_
#pragma prefast(push)
#pragma prefast(disable : 26015 26019, "PREfast noise: Esp:1307" )
#endif

        _Use_decl_annotations_
        FORCE_INLINE float* ConvertHalfToFloatStream(
            float* pOutputStream, 
            size_t outputStride, 
            const HALF* pInputStream, 
            size_t inputStride, 
            size_t halfCount
        ) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pOutputStream != nullptr);
            assert(pInputStream != nullptr);

            assert(inputStride >= sizeof(HALF));
            _Analysis_assume_(inputStride >= sizeof(HALF));

            assert(outputStride >= sizeof(float));
            _Analysis_assume_(outputStride >= sizeof(float));
#endif

#if defined(_F16C_INTRINSICS_) && !defined(_NO_INTRINSICS_)
            const uint8_t* pHalf = reinterpret_cast<const uint8_t*>(pInputStream);
            uint8_t* pFloat = reinterpret_cast<uint8_t*>(pOutputStream);

            size_t i = 0;
            size_t four = halfCount >> 2;
            if (four > 0)
            {
                if (inputStride == sizeof(HALF))
                {
                    if (outputStride == sizeof(float))
                    {
                        if ((reinterpret_cast<uintptr_t>(pFloat) & 0xF) == 0)
                        {
                            // Packed input, aligned & packed output
                            for (size_t j = 0; j < four; ++j)
                            {
                                __m128i HV = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(pHalf));
                                pHalf += inputStride * 4;

                                __m128 FV = _mm_cvtph_ps(HV);

                                STREAM_PS(reinterpret_cast<float*>(pFloat), FV);
                                pFloat += outputStride * 4;
                                i += 4;
                            }
                        }
                        else
                        {
                            // Packed input, packed output
                            for (size_t j = 0; j < four; ++j)
                            {
                                __m128i HV = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(pHalf));
                                pHalf += inputStride * 4;

                                __m128 FV = _mm_cvtph_ps(HV);

                                _mm_storeu_ps(reinterpret_cast<float*>(pFloat), FV);
                                pFloat += outputStride * 4;
                                i += 4;
                            }
                        }
                    }
                    else
                    {
                        // Packed input, scattered output
                        for (size_t j = 0; j < four; ++j)
                        {
                            __m128i HV = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(pHalf));
                            pHalf += inputStride * 4;

                            __m128 FV = _mm_cvtph_ps(HV);

                            _mm_store_ss(reinterpret_cast<float*>(pFloat), FV);
                            pFloat += outputStride;
                            *reinterpret_cast<int*>(pFloat) = _mm_extract_ps(FV, 1);
                            pFloat += outputStride;
                            *reinterpret_cast<int*>(pFloat) = _mm_extract_ps(FV, 2);
                            pFloat += outputStride;
                            *reinterpret_cast<int*>(pFloat) = _mm_extract_ps(FV, 3);
                            pFloat += outputStride;
                            i += 4;
                        }
                    }
                }
                else if (outputStride == sizeof(float))
                {
                    if ((reinterpret_cast<uintptr_t>(pFloat) & 0xF) == 0)
                    {
                        // Scattered input, aligned & packed output
                        for (size_t j = 0; j < four; ++j)
                        {
                            uint16_t H1 = *reinterpret_cast<const HALF*>(pHalf);
                            pHalf += inputStride;
                            uint16_t H2 = *reinterpret_cast<const HALF*>(pHalf);
                            pHalf += inputStride;
                            uint16_t H3 = *reinterpret_cast<const HALF*>(pHalf);
                            pHalf += inputStride;
                            uint16_t H4 = *reinterpret_cast<const HALF*>(pHalf);
                            pHalf += inputStride;

                            __m128i HV = _mm_setzero_si128();
                            HV = _mm_insert_epi16(HV, H1, 0);
                            HV = _mm_insert_epi16(HV, H2, 1);
                            HV = _mm_insert_epi16(HV, H3, 2);
                            HV = _mm_insert_epi16(HV, H4, 3);
                            __m128 FV = _mm_cvtph_ps(HV);

                            STREAM_PS(reinterpret_cast<float*>(pFloat), FV);
                            pFloat += outputStride * 4;
                            i += 4;
                        }
                    }
                    else
                    {
                        // Scattered input, packed output
                        for (size_t j = 0; j < four; ++j)
                        {
                            uint16_t H1 = *reinterpret_cast<const HALF*>(pHalf);
                            pHalf += inputStride;
                            uint16_t H2 = *reinterpret_cast<const HALF*>(pHalf);
                            pHalf += inputStride;
                            uint16_t H3 = *reinterpret_cast<const HALF*>(pHalf);
                            pHalf += inputStride;
                            uint16_t H4 = *reinterpret_cast<const HALF*>(pHalf);
                            pHalf += inputStride;

                            __m128i HV = _mm_setzero_si128();
                            HV = _mm_insert_epi16(HV, H1, 0);
                            HV = _mm_insert_epi16(HV, H2, 1);
                            HV = _mm_insert_epi16(HV, H3, 2);
                            HV = _mm_insert_epi16(HV, H4, 3);
                            __m128 FV = _mm_cvtph_ps(HV);

                            _mm_storeu_ps(reinterpret_cast<float*>(pFloat), FV);
                            pFloat += outputStride * 4;
                            i += 4;
                        }
                    }
                }
                else
                {
                    // Scattered input, scattered output
                    for (size_t j = 0; j < four; ++j)
                    {
                        uint16_t H1 = *reinterpret_cast<const HALF*>(pHalf);
                        pHalf += inputStride;
                        uint16_t H2 = *reinterpret_cast<const HALF*>(pHalf);
                        pHalf += inputStride;
                        uint16_t H3 = *reinterpret_cast<const HALF*>(pHalf);
                        pHalf += inputStride;
                        uint16_t H4 = *reinterpret_cast<const HALF*>(pHalf);
                        pHalf += inputStride;

                        __m128i HV = _mm_setzero_si128();
                        HV = _mm_insert_epi16(HV, H1, 0);
                        HV = _mm_insert_epi16(HV, H2, 1);
                        HV = _mm_insert_epi16(HV, H3, 2);
                        HV = _mm_insert_epi16(HV, H4, 3);
                        __m128 FV = _mm_cvtph_ps(HV);

                        _mm_store_ss(reinterpret_cast<float*>(pFloat), FV);
                        pFloat += outputStride;
                        *reinterpret_cast<int*>(pFloat) = _mm_extract_ps(FV, 1);
                        pFloat += outputStride;
                        *reinterpret_cast<int*>(pFloat) = _mm_extract_ps(FV, 2);
                        pFloat += outputStride;
                        *reinterpret_cast<int*>(pFloat) = _mm_extract_ps(FV, 3);
                        pFloat += outputStride;
                        i += 4;
                    }
                }
            }

            for (; i < halfCount; ++i)
            {
                *reinterpret_cast<float*>(pFloat) = ConvertHalfToFloat(reinterpret_cast<const HALF*>(pHalf)[0]);
                pHalf += inputStride;
                pFloat += outputStride;
            }

            SFENCE();

            return pOutputStream;

#else
            const uint8_t* pHalf = reinterpret_cast<const uint8_t*>(pInputStream);
            uint8_t* pFloat = reinterpret_cast<uint8_t*>(pOutputStream);

            for (size_t i = 0; i < halfCount; i++)
            {
                *reinterpret_cast<float*>(pFloat) = ConvertHalfToFloat(reinterpret_cast<const HALF*>(pHalf)[0]);
                pHalf += inputStride;
                pFloat += outputStride;
            }

            return pOutputStream;
#endif // !_F16C_INTRINSICS_
        }

        FORCE_INLINE HALF ConvertFloatToHalf(float value) noexcept
        {
#if defined(_F16C_INTRINSICS_) && !defined(_NO_INTRINSICS_)
            __m128 V1 = _mm_set_ss(value);
            __m128i V2 = _mm_cvtps_ph(V1, _MM_FROUND_TO_NEAREST_INT);
            return static_cast<HALF>(_mm_extract_epi16(V2, 0)); 

#else
            uint32_t Result;

            uint32_t* IValue = reinterpret_cast<uint32_t*>(&value)[0];
            uint32_t sign = (IValue & 0x80000000U) >> 16U;
            IValue = IValue & 0x7FFFFFFFU;      // Hack off the sign
            if (IValue >= 0x47800000 /*e+16*/)
            {
                // The number is too large to be represented as a half. Return infinity or NaN
                Result = 0x7C00U | ((IValue > 0x7F800000) ? (0x200 | ((IValue >> 13U) & 0x3FFU)) : 0U);
            }
            else if (IValue <= 0x33000000U /*e-25*/)
            {
                Result = 0;
            }
            else if (IValue < 0x38800000U /*e-14*/)
            {
                // The number is too small to be represented as a normalized half.
                // Convert it to a denormalized value.
                uint32_t shift = 125U - (IValue >> 23U);
                IValue = 0x800000U | (IValue & 0x7FFFFFU);
                Result = IValue >> (shift + 1);
                uint32_t s = (IValue & ((1U << Shift) - 1)) != 0;
                Result += (Result | s) & ((IValue >> shift) & 1U);
            }
            else
            {
                // Rebias the exponent to represent the value as a normalized half.
                IValue += 0xC8000000U;
                Result = ((IValue + 0x0FFFU + ((IValue >> 13U) & 1U)) >> 13U) & 0x7FFFU;
            }
            return static_cast<HALF>(Result | sign);
#endif // !_F16C_INTRINSICS_
        }

        _Use_decl_annotations_
        FORCE_INLINE HALF* ConvertFloatToHalfStream(
            HALF* pOutputStream, 
            size_t outputStride, 
            const float* pInputStream, 
            size_t inputStride, 
            size_t floatCount
        ) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pOutputStream != nullptr);
            assert(pInputStream != nullptr);

            assert(inputStride >= sizeof(float));
            _Analysis_assume_(inputStride >= sizeof(float));

            assert(outputStride >= sizeof(HALF));
            _Analysis_assume_(outputStride >= sizeof(HALF));
#endif // DEBUG

#if defined(_F16C_INTRINSICS_) && !defined(_NO_INTRINSICS_)
            const uint8_t* pFloat = reinterpret_cast<const uint8_t*>(pInputStream);
            uint8_t* pHalf = reinterpret_cast<uint8_t*>(pOutputStream);

            size_t i = 0;
            size_t four = floatCount >> 2;
            if (four > 0)
            {
                if (inputStride == sizeof(float))
                {
                    if (outputStride == sizeof(HALF))
                    {
                        if ((reinterpret_cast<uintptr_t>(pFloat) & 0xF) == 0)
                        {
                            // Aligned and packed input, packed output
                            for (size_t j = 0; j < four; ++j)
                            {
                                __m128 FV = _mm_load_ps(reinterpret_cast<const float*>(pFloat));
                                pFloat += inputStride * 4;

                                __m128i HV = _mm_cvtps_ph(FV, _MM_FROUND_TO_NEAREST_INT);

                                _mm_storel_epi64(reinterpret_cast<__m128i*>(pHalf), HV);
                                pHalf += outputStride * 4;
                                i += 4;
                            }
                        }
                        else
                        {
                            // Packed input, packed output
                            for (size_t j = 0; j < four; ++j)
                            {
                                __m128 FV = _mm_loadu_ps(reinterpret_cast<const float*>(pFloat));
                                pFloat += inputStride * 4;

                                __m128i HV = _mm_cvtps_ph(FV, _MM_FROUND_TO_NEAREST_INT);

                                _mm_storel_epi64(reinterpret_cast<__m128i*>(pHalf), HV);
                                pHalf += outputStride * 4;
                                i += 4;
                            }
                        }
                    }
                    else
                    {
                        if ((reinterpret_cast<uintptr_t>(pFloat) & 0xF) == 0)
                        {
                            // Aligned & packed input, scattered output
                            for (size_t j = 0; j < four; ++j)
                            {
                                __m128 FV = _mm_load_ps(reinterpret_cast<const float*>(pFloat));
                                pFloat += inputStride * 4;

                                __m128i HV = _mm_cvtps_ph(FV, _MM_FROUND_TO_NEAREST_INT);

                                *reinterpret_cast<HALF*>(pHalf) = static_cast<HALF>(_mm_extract_epi16(HV, 0));
                                pHalf += outputStride;
                                *reinterpret_cast<HALF*>(pHalf) = static_cast<HALF>(_mm_extract_epi16(HV, 1));
                                pHalf += outputStride;
                                *reinterpret_cast<HALF*>(pHalf) = static_cast<HALF>(_mm_extract_epi16(HV, 2));
                                pHalf += outputStride;
                                *reinterpret_cast<HALF*>(pHalf) = static_cast<HALF>(_mm_extract_epi16(HV, 3));
                                pHalf += outputStride;
                                i += 4;
                            }
                        }
                        else
                        {
                            // Packed input, scattered output
                            for (size_t j = 0; j < four; ++j)
                            {
                                __m128 FV = _mm_loadu_ps(reinterpret_cast<const float*>(pFloat));
                                pFloat += inputStride * 4;

                                __m128i HV = _mm_cvtps_ph(FV, _MM_FROUND_TO_NEAREST_INT);

                                *reinterpret_cast<HALF*>(pHalf) = static_cast<HALF>(_mm_extract_epi16(HV, 0));
                                pHalf += outputStride;
                                *reinterpret_cast<HALF*>(pHalf) = static_cast<HALF>(_mm_extract_epi16(HV, 1));
                                pHalf += outputStride;
                                *reinterpret_cast<HALF*>(pHalf) = static_cast<HALF>(_mm_extract_epi16(HV, 2));
                                pHalf += outputStride;
                                *reinterpret_cast<HALF*>(pHalf) = static_cast<HALF>(_mm_extract_epi16(HV, 3));
                                pHalf += outputStride;
                                i += 4;
                            }
                        }
                    }
                }
                else if (outputStride == sizeof(HALF))
                {
                    // Scattered input, packed output
                    for (size_t j = 0; j < four; ++j)
                    {
                        __m128 FV1 = _mm_load_ss(reinterpret_cast<const float*>(pFloat));
                        pFloat += inputStride;

                        __m128 FV2 = _mm_broadcast_ss(reinterpret_cast<const float*>(pFloat));
                        pFloat += inputStride;

                        __m128 FV3 = _mm_broadcast_ss(reinterpret_cast<const float*>(pFloat));
                        pFloat += inputStride;

                        __m128 FV4 = _mm_broadcast_ss(reinterpret_cast<const float*>(pFloat));
                        pFloat += inputStride;

                        __m128 FV = _mm_blend_ps(FV1, FV2, 0x2);
                        __m128 FT = _mm_blend_ps(FV3, FV4, 0x8);
                        FV = _mm_blend_ps(FV, FT, 0xC);

                        __m128i HV = _mm_cvtps_ph(FV, _MM_FROUND_TO_NEAREST_INT);

                        _mm_storel_epi64(reinterpret_cast<__m128i*>(pHalf), HV);
                        pHalf += outputStride * 4;
                        i += 4;
                    }
                }
                else
                {
                    // Scattered input, scattered output
                    for (size_t j = 0; j < four; ++j)
                    {
                        __m128 FV1 = _mm_load_ss(reinterpret_cast<const float*>(pFloat));
                        pFloat += inputStride;

                        __m128 FV2 = _mm_broadcast_ss(reinterpret_cast<const float*>(pFloat));
                        pFloat += inputStride;

                        __m128 FV3 = _mm_broadcast_ss(reinterpret_cast<const float*>(pFloat));
                        pFloat += inputStride;

                        __m128 FV4 = _mm_broadcast_ss(reinterpret_cast<const float*>(pFloat));
                        pFloat += inputStride;

                        __m128 FV = _mm_blend_ps(FV1, FV2, 0x2);
                        __m128 FT = _mm_blend_ps(FV3, FV4, 0x8);
                        FV = _mm_blend_ps(FV, FT, 0xC);

                        __m128i HV = _mm_cvtps_ph(FV, _MM_FROUND_TO_NEAREST_INT);

                        *reinterpret_cast<HALF*>(pHalf) = static_cast<HALF>(_mm_extract_epi16(HV, 0));
                        pHalf += outputStride;
                        *reinterpret_cast<HALF*>(pHalf) = static_cast<HALF>(_mm_extract_epi16(HV, 1));
                        pHalf += outputStride;
                        *reinterpret_cast<HALF*>(pHalf) = static_cast<HALF>(_mm_extract_epi16(HV, 2));
                        pHalf += outputStride;
                        *reinterpret_cast<HALF*>(pHalf) = static_cast<HALF>(_mm_extract_epi16(HV, 3));
                        pHalf += outputStride;
                        i += 4;
                    }
                }
            }

            for (; i < floatCount; ++i)
            {
                *reinterpret_cast<HALF*>(pHalf) = ConvertFloatToHalf(reinterpret_cast<const float*>(pFloat)[0]);
                pFloat += inputStride;
                pHalf += outputStride;
            }

            return pOutputStream;

#else
            const uint8_t* pFloat = reinterpret_cast<const uint8_t*>(pInputStream);
            uint8_t* pHalf = reinterpret_cast<uint8_t*>(pOutputStream);

            for (size_t i = 0; i < floatCount; i++)
            {
                *reinterpret_cast<HALF*>(pHalf) = ConvertFloatToHalf(reinterpret_cast<const float*>(pFloat)[0]);
                pFloat += inputStride;
                pHalf += outputStride;
            }
            return pOutputStream;
#endif // !_F16C_INTRINSICS_
        }

#ifdef _PREFAST_
#pragma prefast(pop)
#endif

#ifdef _PREFAST_
#pragma prefast(push)
#pragma prefast(disable:28931, "PREfast noise: Esp:1266")
#endif

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadColor(const COLOR* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            // int32_t -> Float conversions are done in one instruction.
            // uint32_t -> Float calls a runtime function. Keep in int32_t
            int32_t iColor = static_cast<int32_t>(pSource->c);
            VECTOR_F32 vColor = { { {
                    static_cast<float>((iColor >> 16) & 0xFF)* (1.0f / 255.0f),
                    static_cast<float>((iColor >> 8) & 0xFF)* (1.0f / 255.0f),
                    static_cast<float>(iColor & 0xFF)* (1.0f / 255.0f),
                    static_cast<float>((iColor >> 24) & 0xFF)* (1.0f / 255.0f)
                } } };
            return vColor.v;

#elif defined(_SSE2_INTRINSICS_)
            // Splat the color in all four entries
            __m128i vInt = _mm_set1_epi32(static_cast<int>(pSource->c));
            // Shift R&0xFF0000, G&0xFF00, B&0xFF, A&0xFF000000
            vInt = _mm_and_si128(vInt, g_MaskA8R8G8B8);
            // a is unsigned! Flip the bit to convert the order to signed
            vInt = _mm_xor_si128(vInt, g_FlipA8R8G8B8);
            // Convert to floating point numbers
            VECTOR vTemp = _mm_cvtepi32_ps(vInt);
            // RGB + 0, A + 0x80000000.f to undo the signed order.
            vTemp = _mm_add_ps(vTemp, g_FixAA8R8G8B8);
            // Convert 0-255 to 0.0f-1.0f
            return _mm_mul_ps(vTemp, g_NormalizeA8R8G8B8);
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadHalf2(const HALF2* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_F16C_INTRINSICS_) && !defined(_NO_INTRINSICS_)
            __m128 V = _mm_load_ss(reinterpret_cast<const float*>(pSource));
            return _mm_cvtph_ps(_mm_castps_si128(V));
#else
            VECTOR_F32 vResult = { { {
                    ConvertHalfToFloat(pSource->x),
                    ConvertHalfToFloat(pSource->y),
                    0.0f,
                    0.0f
                } } };
            return vResult.v;
#endif // !_F16C_INTRINSICS_
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadShortN2(const SHORTN2* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR_F32 vResult = { { {
                    (pSource->x == -32768) ? -1.f : (static_cast<float>(pSource->x)* (1.0f / 32767.0f)),
                    (pSource->y == -32768) ? -1.f : (static_cast<float>(pSource->y)* (1.0f / 32767.0f)),
                    0.0f,
                    0.0f
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            // Splat the two shorts in all four entries (WORD alignment okay,
            // DWORD alignment preferred)
            __m128 vTemp = _mm_load_ps1(reinterpret_cast<const float*>(&pSource->x));
            // Mask x&0xFFFF, y&0xFFFF0000,z&0,w&0
            vTemp = _mm_and_ps(vTemp, g_MaskX16Y16);
            // x needs to be sign extended
            vTemp = _mm_xor_ps(vTemp, g_FlipX16Y16);
            // Convert to floating point numbers
            vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
            // x - 0x8000 to undo the signed order.
            vTemp = _mm_add_ps(vTemp, g_FixX16Y16);
            // Convert -1.0f - 1.0f
            vTemp = _mm_mul_ps(vTemp, g_NormalizeX16Y16);
            // Clamp result (for case of -32768)
            return _mm_max_ps(vTemp, g_NegativeOne);
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadShort2(const SHORT2* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
            VECTOR_F32 vResult = { { {
                    static_cast<float>(pSource->x),
                    static_cast<float>(pSource->y),
                    0.f,
                    0.f
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            // Splat the two shorts in all four entries (WORD alignment okay,
            // DWORD alignment preferred)
            __m128 vTemp = _mm_load_ps1(reinterpret_cast<const float*>(&pSource->x));
            // Mask x&0xFFFF, y&0xFFFF0000,z&0,w&0
            vTemp = _mm_and_ps(vTemp, g_MaskX16Y16);
            // x needs to be sign extended
            vTemp = _mm_xor_ps(vTemp, g_FlipX16Y16);
            // Convert to floating point numbers
            vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
            // x - 0x8000 to undo the signed order.
            vTemp = _mm_add_ps(vTemp, g_FixX16Y16);
            // Y is 65536 too large
            return _mm_mul_ps(vTemp, g_FixupY16);
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadUShortN2(const USHORTN2* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR_F32 vResult = { { {
                    static_cast<float>(pSource->x) / 65535.0f,
                    static_cast<float>(pSource->y) / 65535.0f,
                    0.f,
                    0.f
                } } };
            return vResult.v;

#elif defined(_XM_SSE_INTRINSICS_)
            static const VECTOR_F32 FixupY16 = { { { 1.0f / 65535.0f, 1.0f / (65535.0f * 65536.0f), 0.0f, 0.0f } } };
            static const VECTOR_F32 FixaddY16 = { { { 0, 32768.0f * 65536.0f, 0, 0 } } };
            // Splat the two shorts in all four entries (WORD alignment okay,
            // DWORD alignment preferred)
            __m128 vTemp = _mm_load_ps1(reinterpret_cast<const float*>(&pSource->x));
            // Mask x&0xFFFF, y&0xFFFF0000,z&0,w&0
            vTemp = _mm_and_ps(vTemp, g_MaskX16Y16);
            // y needs to be sign flipped
            vTemp = _mm_xor_ps(vTemp, g_FlipY);
            // Convert to floating point numbers
            vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
            // y + 0x8000 to undo the signed order.
            vTemp = _mm_add_ps(vTemp, FixaddY16);
            // Y is 65536 times too large
            vTemp = _mm_mul_ps(vTemp, FixupY16);
            return vTemp;
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadUShort2(const USHORT2* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR_F32 vResult = { { {
                    static_cast<float>(pSource->x),
                    static_cast<float>(pSource->y),
                    0.f,
                    0.f
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 FixaddY16 = { { { 0, 32768.0f, 0, 0 } } };
            // Splat the two shorts in all four entries (WORD alignment okay,
            // DWORD alignment preferred)
            __m128 vTemp = _mm_load_ps1(reinterpret_cast<const float*>(&pSource->x));
            // Mask x&0xFFFF, y&0xFFFF0000,z&0,w&0
            vTemp = _mm_and_ps(vTemp, g_MaskX16Y16);
            // y needs to be sign flipped
            vTemp = _mm_xor_ps(vTemp, g_FlipY);
            // Convert to floating point numbers
            vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
            // Y is 65536 times too large
            vTemp = _mm_mul_ps(vTemp, g_FixupY16);
            // y + 0x8000 to undo the signed order.
            vTemp = _mm_add_ps(vTemp, FixaddY16);
            return vTemp;
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadByteN2(const BYTEN2* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR_F32 vResult = { { {
                    (pSource->x == -128) ? -1.f : (static_cast<float>(pSource->x)* (1.0f / 127.0f)),
                    (pSource->y == -128) ? -1.f : (static_cast<float>(pSource->y)* (1.0f / 127.0f)),
                    0.0f,
                    0.0f
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 Scale = { { { 1.0f / 127.0f, 1.0f / (127.0f * 256.0f), 0, 0 } } };
            static const VECTOR_U32 Mask = { { { 0xFF, 0xFF00, 0, 0 } } };
            // Splat the color in all four entries (x,z,y,w)
            __m128i vInt = LOADU_SI16(&pSource->v);
            XMVECTOR vTemp = PERMUTE_PS(_mm_castsi128_ps(vInt), _MM_SHUFFLE(0, 0, 0, 0));
            // Mask
            vTemp = _mm_and_ps(vTemp, Mask);
            // x,y and z are unsigned! Flip the bits to convert the order to signed
            vTemp = _mm_xor_ps(vTemp, g_XorByte4);
            // Convert to floating point numbers
            vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
            // x, y and z - 0x80 to complete the conversion
            vTemp = _mm_add_ps(vTemp, g_AddByte4);
            // Fix y, z and w because they are too large
            vTemp = _mm_mul_ps(vTemp, Scale);
            // Clamp result (for case of -128)
            return _mm_max_ps(vTemp, g_NegativeOne);
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadByte2(const BYTE2* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR_F32 vResult = { { {
                    static_cast<float>(pSource->x),
                    static_cast<float>(pSource->y),
                    0.0f,
                    0.0f
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 Scale = { { { 1.0f, 1.0f / 256.0f, 1.0f / 65536.0f, 1.0f / (65536.0f * 256.0f) } } };
            static const VECTOR_U32 Mask = { { { 0xFF, 0xFF00, 0, 0 } } };
            // Splat the color in all four entries (x,z,y,w)
            __m128i vInt = LOADU_SI16(&pSource->v);
            VECTOR vTemp = PERMUTE_PS(_mm_castsi128_ps(vInt), _MM_SHUFFLE(0, 0, 0, 0));
            // Mask
            vTemp = _mm_and_ps(vTemp, Mask);
            // x,y and z are unsigned! Flip the bits to convert the order to signed
            vTemp = _mm_xor_ps(vTemp, g_XorByte4);
            // Convert to floating point numbers
            vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
            // x, y and z - 0x80 to complete the conversion
            vTemp = _mm_add_ps(vTemp, g_AddByte4);
            // Fix y, z and w because they are too large
            return _mm_mul_ps(vTemp, Scale);
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadUByteN2(const UBYTEN2* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR_F32 vResult = { { {
                    static_cast<float>(pSource->x)* (1.0f / 255.0f),
                    static_cast<float>(pSource->y)* (1.0f / 255.0f),
                    0.0f,
                    0.0f
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 Scale = { { { 1.0f / 255.0f, 1.0f / (255.0f * 256.0f), 0, 0 } } };
            static const VECTOR_U32 Mask = { { { 0xFF, 0xFF00, 0, 0 } } };
            // Splat the color in all four entries (x,z,y,w)
            __m128i vInt = LOADU_SI16(&pSource->v);
            VECTOR vTemp = PERMUTE_PS(_mm_castsi128_ps(vInt), _MM_SHUFFLE(0, 0, 0, 0));
            // Mask
            vTemp = _mm_and_ps(vTemp, Mask);
            // w is signed! Flip the bits to convert the order to unsigned
            vTemp = _mm_xor_ps(vTemp, g_FlipW);
            // Convert to floating point numbers
            vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
            // w + 0x80 to complete the conversion
            vTemp = _mm_add_ps(vTemp, g_AddUDec4);
            // Fix y, z and w because they are too large
            return _mm_mul_ps(vTemp, Scale);
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadUByte2(const UBYTE2* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR_F32 vResult = { { {
                    static_cast<float>(pSource->x),
                    static_cast<float>(pSource->y),
                    0.0f,
                    0.0f
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 Scale = { { { 1.0f, 1.0f / 256.0f, 0, 0 } } };
            static const VECTOR_U32 Mask = { { { 0xFF, 0xFF00, 0, 0 } } };
            // Splat the color in all four entries (x,z,y,w)
            __m128i vInt = LOADU_SI16(&pSource->v);
            VECTOR vTemp = PERMUTE_PS(_mm_castsi128_ps(vInt), _MM_SHUFFLE(0, 0, 0, 0));
            // Mask
            vTemp = _mm_and_ps(vTemp, Mask);
            // w is signed! Flip the bits to convert the order to unsigned
            vTemp = _mm_xor_ps(vTemp, g_FlipW);
            // Convert to floating point numbers
            vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
            // w + 0x80 to complete the conversion
            vTemp = _mm_add_ps(vTemp, g_AddUDec4);
            // Fix y, z and w because they are too large
            return _mm_mul_ps(vTemp, Scale);
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadU565(const U565* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR_F32 vResult = { { {
                    float(pSource->v & 0x1F),
                    float((pSource->v >> 5) & 0x3F),
                    float((pSource->v >> 11) & 0x1F),
                    0.f,
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_I32 U565And = { { { 0x1F, 0x3F << 5, 0x1F << 11, 0 } } };
            static const VECTOR_F32 U565Mul = { { { 1.0f, 1.0f / 32.0f, 1.0f / 2048.f, 0 } } };
            // Get the 16 bit value and splat it
            __m128i vInt = LOADU_SI16(&pSource->v);
            VECTOR vResult = PERMUTE_PS(_mm_castsi128_ps(vInt), _MM_SHUFFLE(0, 0, 0, 0));
            // Mask off x, y and z
            vResult = _mm_and_ps(vResult, U565And);
            // Convert to float
            vResult = _mm_cvtepi32_ps(_mm_castps_si128(vResult));
            // Normalize x, y, and z
            vResult = _mm_mul_ps(vResult, U565Mul);
            return vResult;
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadFloat3PK(const FLOAT3PK* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

            ALIGNED(16) uint32_t Result[4];
            uint32_t Mantissa;
            uint32_t Exponent;

            // X Channel (6-bit mantissa)
            Mantissa = pSource->xm;

            if (pSource->xe == 0x1f) // INF or NAN
            {
                Result[0] = static_cast<uint32_t>(0x7f800000 | (static_cast<int>(pSource->xm) << 17));
            }
            else
            {
                if (pSource->xe != 0) // The value is normalized
                {
                    Exponent = pSource->xe;
                }
                else if (Mantissa != 0) // The value is denormalized
                {
                    // Normalize the value in the resulting float
                    Exponent = 1;

                    do
                    {
                        Exponent--;
                        Mantissa <<= 1;
                    } while ((Mantissa & 0x40) == 0);

                    Mantissa &= 0x3F;
                }
                else // The value is zero
                {
                    Exponent = static_cast<uint32_t>(-112);
                }

                Result[0] = ((Exponent + 112) << 23) | (Mantissa << 17);
            }

            // Y Channel (6-bit mantissa)
            Mantissa = pSource->ym;

            if (pSource->ye == 0x1f) // INF or NAN
            {
                Result[1] = static_cast<uint32_t>(0x7f800000 | (static_cast<int>(pSource->ym) << 17));
            }
            else
            {
                if (pSource->ye != 0) // The value is normalized
                {
                    Exponent = pSource->ye;
                }
                else if (Mantissa != 0) // The value is denormalized
                {
                    // Normalize the value in the resulting float
                    Exponent = 1;

                    do
                    {
                        Exponent--;
                        Mantissa <<= 1;
                    } while ((Mantissa & 0x40) == 0);

                    Mantissa &= 0x3F;
                }
                else // The value is zero
                {
                    Exponent = static_cast<uint32_t>(-112);
                }

                Result[1] = ((Exponent + 112) << 23) | (Mantissa << 17);
            }

            // Z Channel (5-bit mantissa)
            Mantissa = pSource->zm;

            if (pSource->ze == 0x1f) // INF or NAN
            {
                Result[2] = static_cast<uint32_t>(0x7f800000 | (static_cast<int>(pSource->zm) << 17));
            }
            else
            {
                if (pSource->ze != 0) // The value is normalized
                {
                    Exponent = pSource->ze;
                }
                else if (Mantissa != 0) // The value is denormalized
                {
                    // Normalize the value in the resulting float
                    Exponent = 1;

                    do
                    {
                        Exponent--;
                        Mantissa <<= 1;
                    } while ((Mantissa & 0x20) == 0);

                    Mantissa &= 0x1F;
                }
                else // The value is zero
                {
                    Exponent = static_cast<uint32_t>(-112);
                }

                Result[2] = ((Exponent + 112) << 23) | (Mantissa << 18);
            }

            return LoadFloat3A(reinterpret_cast<const FLOAT3A*>(&Result));
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadFloat3SE(const FLOAT3SE* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

            union { float f; int32_t i; } fi;
            fi.i = 0x33800000 + (pSource->e << 23);
            float Scale = fi.f;

            VECTOR_F32 v = { { {
                    Scale * float(pSource->xm),
                    Scale * float(pSource->ym),
                    Scale * float(pSource->zm),
                    1.0f } } };
            return v;
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadHalf4(const HALF4* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_F16C_INTRINSICS_) && !defined(_NO_INTRINSICS_)
            __m128i V = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(pSource));
            return _mm_cvtph_ps(V);
#else
            VECTOR_F32 vResult = { { {
                    ConvertHalfToFloat(pSource->x),
                    ConvertHalfToFloat(pSource->y),
                    ConvertHalfToFloat(pSource->z),
                    ConvertHalfToFloat(pSource->w)
                } } };
            return vResult.v;
#endif // !_F16C_INTRINSICS_
        }

        
    } // namespace PackedVector
    
} // namespace UltReality::Math

#endif // !ULTREALITY_MATH_PACKED_VECTOR_INL
