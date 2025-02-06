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

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadShortN4(const SHORTN4* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR_F32 vResult = { { {
                    (pSource->x == -32768) ? -1.f : (static_cast<float>(pSource->x)* (1.0f / 32767.0f)),
                    (pSource->y == -32768) ? -1.f : (static_cast<float>(pSource->y)* (1.0f / 32767.0f)),
                    (pSource->z == -32768) ? -1.f : (static_cast<float>(pSource->z)* (1.0f / 32767.0f)),
                    (pSource->w == -32768) ? -1.f : (static_cast<float>(pSource->w)* (1.0f / 32767.0f))
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            // Splat the color in all four entries (x,z,y,w)
            __m128d vIntd = _mm_load1_pd(reinterpret_cast<const double*>(&pSource->x));
            // Shift x&0ffff,z&0xffff,y&0xffff0000,w&0xffff0000
            __m128 vTemp = _mm_and_ps(_mm_castpd_ps(vIntd), g_MaskX16Y16Z16W16);
            // x and z are unsigned! Flip the bits to convert the order to signed
            vTemp = _mm_xor_ps(vTemp, g_FlipX16Y16Z16W16);
            // Convert to floating point numbers
            vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
            // x and z - 0x8000 to complete the conversion
            vTemp = _mm_add_ps(vTemp, g_FixX16Y16Z16W16);
            // Convert to -1.0f - 1.0f
            vTemp = _mm_mul_ps(vTemp, g_NormalizeX16Y16Z16W16);
            // Very important! The entries are x,z,y,w, flip it to x,y,z,w
            vTemp = PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 1, 2, 0));
            // Clamp result (for case of -32768)
            return _mm_max_ps(vTemp, g_NegativeOne);
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadShort4(const SHORT4* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR_F32 vResult = { { {
                    static_cast<float>(pSource->x),
                    static_cast<float>(pSource->y),
                    static_cast<float>(pSource->z),
                    static_cast<float>(pSource->w)
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            // Splat the color in all four entries (x,z,y,w)
            __m128d vIntd = _mm_load1_pd(reinterpret_cast<const double*>(&pSource->x));
            // Shift x&0ffff,z&0xffff,y&0xffff0000,w&0xffff0000
            __m128 vTemp = _mm_and_ps(_mm_castpd_ps(vIntd), g_MaskX16Y16Z16W16);
            // x and z are unsigned! Flip the bits to convert the order to signed
            vTemp = _mm_xor_ps(vTemp, g_FlipX16Y16Z16W16);
            // Convert to floating point numbers
            vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
            // x and z - 0x8000 to complete the conversion
            vTemp = _mm_add_ps(vTemp, g_FixX16Y16Z16W16);
            // Fix y and w because they are 65536 too large
            vTemp = _mm_mul_ps(vTemp, g_FixupY16W16);
            // Very important! The entries are x,z,y,w, flip it to x,y,z,w
            return PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 1, 2, 0));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadUShortN4(const USHORTN4* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR_F32 vResult = { { {
                    static_cast<float>(pSource->x) / 65535.0f,
                    static_cast<float>(pSource->y) / 65535.0f,
                    static_cast<float>(pSource->z) / 65535.0f,
                    static_cast<float>(pSource->w) / 65535.0f
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 FixupY16W16 = { { { 1.0f / 65535.0f, 1.0f / 65535.0f, 1.0f / (65535.0f * 65536.0f), 1.0f / (65535.0f * 65536.0f) } } };
            static const VECTOR_F32 FixaddY16W16 = { { { 0, 0, 32768.0f * 65536.0f, 32768.0f * 65536.0f } } };
            // Splat the color in all four entries (x,z,y,w)
            __m128d vIntd = _mm_load1_pd(reinterpret_cast<const double*>(&pSource->x));
            // Shift x&0ffff,z&0xffff,y&0xffff0000,w&0xffff0000
            __m128 vTemp = _mm_and_ps(_mm_castpd_ps(vIntd), g_MaskX16Y16Z16W16);
            // y and w are signed! Flip the bits to convert the order to unsigned
            vTemp = _mm_xor_ps(vTemp, g_FlipZW);
            // Convert to floating point numbers
            vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
            // y and w + 0x8000 to complete the conversion
            vTemp = _mm_add_ps(vTemp, FixaddY16W16);
            // Fix y and w because they are 65536 too large
            vTemp = _mm_mul_ps(vTemp, FixupY16W16);
            // Very important! The entries are x,z,y,w, flip it to x,y,z,w
            return PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 1, 2, 0));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadUShort4(const USHORT4* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR_F32 vResult = { { {
                    static_cast<float>(pSource->x),
                    static_cast<float>(pSource->y),
                    static_cast<float>(pSource->z),
                    static_cast<float>(pSource->w)
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 FixaddY16W16 = { { { 0, 0, 32768.0f, 32768.0f } } };
            // Splat the color in all four entries (x,z,y,w)
            __m128d vIntd = _mm_load1_pd(reinterpret_cast<const double*>(&pSource->x));
            // Shift x&0ffff,z&0xffff,y&0xffff0000,w&0xffff0000
            __m128 vTemp = _mm_and_ps(_mm_castpd_ps(vIntd), g_MaskX16Y16Z16W16);
            // y and w are signed! Flip the bits to convert the order to unsigned
            vTemp = _mm_xor_ps(vTemp, g_FlipZW);
            // Convert to floating point numbers
            vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
            // Fix y and w because they are 65536 too large
            vTemp = _mm_mul_ps(vTemp, g_FixupY16W16);
            // y and w + 0x8000 to complete the conversion
            vTemp = _mm_add_ps(vTemp, FixaddY16W16);
            // Very important! The entries are x,z,y,w, flip it to x,y,z,w
            return PERMUTE_PS(vTemp, _MM_SHUFFLE(3, 1, 2, 0));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadXDecN4(const XDECN4* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            static const uint32_t SignExtend[] = { 0x00000000, 0xFFFFFC00 };

            uint32_t ElementX = pSource->v & 0x3FF;
            uint32_t ElementY = (pSource->v >> 10) & 0x3FF;
            uint32_t ElementZ = (pSource->v >> 20) & 0x3FF;

            VECTOR_F32 vResult = { { {
                    (ElementX == 0x200) ? -1.f : (static_cast<float>(static_cast<int16_t>(ElementX | SignExtend[ElementX >> 9])) / 511.0f),
                    (ElementY == 0x200) ? -1.f : (static_cast<float>(static_cast<int16_t>(ElementY | SignExtend[ElementY >> 9])) / 511.0f),
                    (ElementZ == 0x200) ? -1.f : (static_cast<float>(static_cast<int16_t>(ElementZ | SignExtend[ElementZ >> 9])) / 511.0f),
                    static_cast<float>(pSource->v >> 30) / 3.0f
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            // Splat the color in all four entries
            __m128 vTemp = _mm_load_ps1(reinterpret_cast<const float*>(&pSource->v));
            // Shift R&0xFF0000, G&0xFF00, B&0xFF, A&0xFF000000
            vTemp = _mm_and_ps(vTemp, g_MaskA2B10G10R10);
            // a is unsigned! Flip the bit to convert the order to signed
            vTemp = _mm_xor_ps(vTemp, g_FlipA2B10G10R10);
            // Convert to floating point numbers
            vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
            // RGB + 0, A + 0x80000000.f to undo the signed order.
            vTemp = _mm_add_ps(vTemp, g_FixAA2B10G10R10);
            // Convert 0-255 to 0.0f-1.0f
            vTemp = _mm_mul_ps(vTemp, g_NormalizeA2B10G10R10);
            // Clamp result (for case of -512)
            return _mm_max_ps(vTemp, g_NegativeOne);
#endif
        }

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4996)
// C4996: ignore deprecation warning
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadXDec4(const XDEC4* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            static const uint32_t SignExtend[] = { 0x00000000, 0xFFFFFC00 };

            uint32_t ElementX = pSource->v & 0x3FF;
            uint32_t ElementY = (pSource->v >> 10) & 0x3FF;
            uint32_t ElementZ = (pSource->v >> 20) & 0x3FF;

            VECTOR_F32 vResult = { { {
                    static_cast<float>(static_cast<int16_t>(ElementX | SignExtend[ElementX >> 9])),
                    static_cast<float>(static_cast<int16_t>(ElementY | SignExtend[ElementY >> 9])),
                    static_cast<float>(static_cast<int16_t>(ElementZ | SignExtend[ElementZ >> 9])),
                    static_cast<float>(pSource->v >> 30)
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_U32 XDec4Xor = { { { 0x200, 0x200 << 10, 0x200 << 20, 0x80000000 } } };
            static const VECTOR_F32 XDec4Add = { { { -512.0f, -512.0f * 1024.0f, -512.0f * 1024.0f * 1024.0f, 32768 * 65536.0f } } };
            // Splat the color in all four entries
            VECTOR vTemp = _mm_load_ps1(reinterpret_cast<const float*>(&pSource->v));
            // Shift R&0xFF0000, G&0xFF00, B&0xFF, A&0xFF000000
            vTemp = _mm_and_ps(vTemp, g_MaskDec4);
            // a is unsigned! Flip the bit to convert the order to signed
            vTemp = _mm_xor_ps(vTemp, XDec4Xor);
            // Convert to floating point numbers
            vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
            // RGB + 0, A + 0x80000000.f to undo the signed order.
            vTemp = _mm_add_ps(vTemp, XDec4Add);
            // Convert 0-255 to 0.0f-1.0f
            vTemp = _mm_mul_ps(vTemp, g_MulDec4);
            return vTemp;
#endif
        }

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadUDecN4(const UDECN4* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)

            uint32_t ElementX = pSource->v & 0x3FF;
            uint32_t ElementY = (pSource->v >> 10) & 0x3FF;
            uint32_t ElementZ = (pSource->v >> 20) & 0x3FF;

            VECTOR_F32 vResult = { { {
                    static_cast<float>(ElementX) / 1023.0f,
                    static_cast<float>(ElementY) / 1023.0f,
                    static_cast<float>(ElementZ) / 1023.0f,
                    static_cast<float>(pSource->v >> 30) / 3.0f
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 UDecN4Mul = { { { 1.0f / 1023.0f, 1.0f / (1023.0f * 1024.0f), 1.0f / (1023.0f * 1024.0f * 1024.0f), 1.0f / (3.0f * 1024.0f * 1024.0f * 1024.0f) } } };
            // Splat the color in all four entries
            VECTOR vTemp = _mm_load_ps1(reinterpret_cast<const float*>(&pSource->v));
            // Shift R&0xFF0000, G&0xFF00, B&0xFF, A&0xFF000000
            vTemp = _mm_and_ps(vTemp, g_MaskDec4);
            // a is unsigned! Flip the bit to convert the order to signed
            vTemp = _mm_xor_ps(vTemp, g_FlipW);
            // Convert to floating point numbers
            vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
            // RGB + 0, A + 0x80000000.f to undo the signed order.
            vTemp = _mm_add_ps(vTemp, g_AddUDec4);
            // Convert 0-255 to 0.0f-1.0f
            vTemp = _mm_mul_ps(vTemp, UDecN4Mul);
            return vTemp;
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadUDecN4_XR(const UDECN4* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif

#if defined(_NO_INTRINSICS_)

            int32_t ElementX = pSource->v & 0x3FF;
            int32_t ElementY = (pSource->v >> 10) & 0x3FF;
            int32_t ElementZ = (pSource->v >> 20) & 0x3FF;

            VECTOR_F32 vResult = { { {
                    static_cast<float>(ElementX - 0x180) / 510.0f,
                    static_cast<float>(ElementY - 0x180) / 510.0f,
                    static_cast<float>(ElementZ - 0x180) / 510.0f,
                    static_cast<float>(pSource->v >> 30) / 3.0f
                } } };

            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 XRMul = { { { 1.0f / 510.0f, 1.0f / (510.0f * 1024.0f), 1.0f / (510.0f * 1024.0f * 1024.0f), 1.0f / (3.0f * 1024.0f * 1024.0f * 1024.0f) } } };
            static const VECTOR_I32 XRBias = { { { 0x180, 0x180 * 1024, 0x180 * 1024 * 1024, 0 } } };
            // Splat the color in all four entries
            VECTOR vTemp = _mm_load_ps1(reinterpret_cast<const float*>(&pSource->v));
            // Mask channels
            vTemp = _mm_and_ps(vTemp, g_MaskDec4);
            // Subtract bias
            vTemp = _mm_castsi128_ps(_mm_sub_epi32(_mm_castps_si128(vTemp), XRBias));
            // a is unsigned! Flip the bit to convert the order to signed
            vTemp = _mm_xor_ps(vTemp, g_FlipW);
            // Convert to floating point numbers
            vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
            // RGB + 0, A + 0x80000000.f to undo the signed order.
            vTemp = _mm_add_ps(vTemp, g_AddUDec4);
            // Convert to 0.0f-1.0f
            return _mm_mul_ps(vTemp, XRMul);
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadUDec4(const UDEC4* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            uint32_t ElementX = pSource->v & 0x3FF;
            uint32_t ElementY = (pSource->v >> 10) & 0x3FF;
            uint32_t ElementZ = (pSource->v >> 20) & 0x3FF;

            VECTOR_F32 vResult = { { {
                    static_cast<float>(ElementX),
                    static_cast<float>(ElementY),
                    static_cast<float>(ElementZ),
                    static_cast<float>(pSource->v >> 30)
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            // Splat the color in all four entries
            VECTOR vTemp = _mm_load_ps1(reinterpret_cast<const float*>(&pSource->v));
            // Shift R&0xFF0000, G&0xFF00, B&0xFF, A&0xFF000000
            vTemp = _mm_and_ps(vTemp, g_MaskDec4);
            // a is unsigned! Flip the bit to convert the order to signed
            vTemp = _mm_xor_ps(vTemp, g_FlipW);
            // Convert to floating point numbers
            vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
            // RGB + 0, A + 0x80000000.f to undo the signed order.
            vTemp = _mm_add_ps(vTemp, g_AddUDec4);
            // Convert 0-255 to 0.0f-1.0f
            vTemp = _mm_mul_ps(vTemp, g_MulDec4);
            return vTemp;
#endif
        }

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4996)
// C4996: ignore deprecation warning
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadDecN4(const DECN4* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            static const uint32_t SignExtend[] = { 0x00000000, 0xFFFFFC00 };
            static const uint32_t SignExtendW[] = { 0x00000000, 0xFFFFFFFC };

            uint32_t ElementX = pSource->v & 0x3FF;
            uint32_t ElementY = (pSource->v >> 10) & 0x3FF;
            uint32_t ElementZ = (pSource->v >> 20) & 0x3FF;
            uint32_t ElementW = pSource->v >> 30;

            VECTOR_F32 vResult = { { {
                    (ElementX == 0x200) ? -1.f : (static_cast<float>(static_cast<int16_t>(ElementX | SignExtend[ElementX >> 9])) / 511.0f),
                    (ElementY == 0x200) ? -1.f : (static_cast<float>(static_cast<int16_t>(ElementY | SignExtend[ElementY >> 9])) / 511.0f),
                    (ElementZ == 0x200) ? -1.f : (static_cast<float>(static_cast<int16_t>(ElementZ | SignExtend[ElementZ >> 9])) / 511.0f),
                    (ElementW == 0x2) ? -1.f : static_cast<float>(static_cast<int16_t>(ElementW | SignExtendW[(ElementW >> 1) & 1]))
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 DecN4Mul = { { { 1.0f / 511.0f, 1.0f / (511.0f * 1024.0f), 1.0f / (511.0f * 1024.0f * 1024.0f), 1.0f / (1024.0f * 1024.0f * 1024.0f) } } };
            // Splat the color in all four entries
            VECTOR vTemp = _mm_load_ps1(reinterpret_cast<const float*>(&pSource->v));
            // Shift R&0xFF0000, G&0xFF00, B&0xFF, A&0xFF000000
            vTemp = _mm_and_ps(vTemp, g_MaskDec4);
            // a is unsigned! Flip the bit to convert the order to signed
            vTemp = _mm_xor_ps(vTemp, g_XorDec4);
            // Convert to floating point numbers
            vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
            // RGB + 0, A + 0x80000000.f to undo the signed order.
            vTemp = _mm_add_ps(vTemp, g_AddDec4);
            // Convert 0-255 to 0.0f-1.0f
            vTemp = _mm_mul_ps(vTemp, DecN4Mul);
            // Clamp result (for case of -512/-1)
            return _mm_max_ps(vTemp, g_NegativeOne);
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadDec4(const DEC4* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            static const uint32_t SignExtend[] = { 0x00000000, 0xFFFFFC00 };
            static const uint32_t SignExtendW[] = { 0x00000000, 0xFFFFFFFC };

            uint32_t ElementX = pSource->v & 0x3FF;
            uint32_t ElementY = (pSource->v >> 10) & 0x3FF;
            uint32_t ElementZ = (pSource->v >> 20) & 0x3FF;
            uint32_t ElementW = pSource->v >> 30;

            VECTOR_F32 vResult = { { {
                    static_cast<float>(static_cast<int16_t>(ElementX | SignExtend[ElementX >> 9])),
                    static_cast<float>(static_cast<int16_t>(ElementY | SignExtend[ElementY >> 9])),
                    static_cast<float>(static_cast<int16_t>(ElementZ | SignExtend[ElementZ >> 9])),
                    static_cast<float>(static_cast<int16_t>(ElementW | SignExtendW[ElementW >> 1]))
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            // Splat the color in all four entries
            VECTOR vTemp = _mm_load_ps1(reinterpret_cast<const float*>(&pSource->v));
            // Shift R&0xFF0000, G&0xFF00, B&0xFF, A&0xFF000000
            vTemp = _mm_and_ps(vTemp, g_MaskDec4);
            // a is unsigned! Flip the bit to convert the order to signed
            vTemp = _mm_xor_ps(vTemp, g_XorDec4);
            // Convert to floating point numbers
            vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
            // RGB + 0, A + 0x80000000.f to undo the signed order.
            vTemp = _mm_add_ps(vTemp, g_AddDec4);
            // Convert 0-255 to 0.0f-1.0f
            vTemp = _mm_mul_ps(vTemp, g_MulDec4);
            return vTemp;
#endif
        }

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadUByteN4(const UBYTEN4* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR_F32 vResult = { { {
                    static_cast<float>(pSource->x) / 255.0f,
                    static_cast<float>(pSource->y) / 255.0f,
                    static_cast<float>(pSource->z) / 255.0f,
                    static_cast<float>(pSource->w) / 255.0f
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 LoadUByteN4Mul = { { { 1.0f / 255.0f, 1.0f / (255.0f * 256.0f), 1.0f / (255.0f * 65536.0f), 1.0f / (255.0f * 65536.0f * 256.0f) } } };
            // Splat the color in all four entries (x,z,y,w)
            VECTOR vTemp = _mm_load1_ps(reinterpret_cast<const float*>(&pSource->x));
            // Mask x&0ff,y&0xff00,z&0xff0000,w&0xff000000
            vTemp = _mm_and_ps(vTemp, g_MaskByte4);
            // w is signed! Flip the bits to convert the order to unsigned
            vTemp = _mm_xor_ps(vTemp, g_FlipW);
            // Convert to floating point numbers
            vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
            // w + 0x80 to complete the conversion
            vTemp = _mm_add_ps(vTemp, g_AddUDec4);
            // Fix y, z and w because they are too large
            vTemp = _mm_mul_ps(vTemp, LoadUByteN4Mul);
            return vTemp;
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadUByte4(const UBYTE4* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR_F32 vResult = { { {
                    static_cast<float>(pSource->x),
                    static_cast<float>(pSource->y),
                    static_cast<float>(pSource->z),
                    static_cast<float>(pSource->w)
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 LoadUByte4Mul = { { { 1.0f, 1.0f / 256.0f, 1.0f / 65536.0f, 1.0f / (65536.0f * 256.0f) } } };
            // Splat the color in all four entries (x,z,y,w)
            VECTOR vTemp = _mm_load1_ps(reinterpret_cast<const float*>(&pSource->x));
            // Mask x&0ff,y&0xff00,z&0xff0000,w&0xff000000
            vTemp = _mm_and_ps(vTemp, g_MaskByte4);
            // w is signed! Flip the bits to convert the order to unsigned
            vTemp = _mm_xor_ps(vTemp, g_FlipW);
            // Convert to floating point numbers
            vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
            // w + 0x80 to complete the conversion
            vTemp = _mm_add_ps(vTemp, g_AddUDec4);
            // Fix y, z and w because they are too large
            vTemp = _mm_mul_ps(vTemp, LoadUByte4Mul);
            return vTemp;
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadByteN4(const BYTEN4* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR_F32 vResult = { { {
                    (pSource->x == -128) ? -1.f : (static_cast<float>(pSource->x) / 127.0f),
                    (pSource->y == -128) ? -1.f : (static_cast<float>(pSource->y) / 127.0f),
                    (pSource->z == -128) ? -1.f : (static_cast<float>(pSource->z) / 127.0f),
                    (pSource->w == -128) ? -1.f : (static_cast<float>(pSource->w) / 127.0f)
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 LoadByteN4Mul = { { { 1.0f / 127.0f, 1.0f / (127.0f * 256.0f), 1.0f / (127.0f * 65536.0f), 1.0f / (127.0f * 65536.0f * 256.0f) } } };
            // Splat the color in all four entries (x,z,y,w)
            VECTOR vTemp = _mm_load1_ps(reinterpret_cast<const float*>(&pSource->x));
            // Mask x&0ff,y&0xff00,z&0xff0000,w&0xff000000
            vTemp = _mm_and_ps(vTemp, g_MaskByte4);
            // x,y and z are unsigned! Flip the bits to convert the order to signed
            vTemp = _mm_xor_ps(vTemp, g_XorByte4);
            // Convert to floating point numbers
            vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
            // x, y and z - 0x80 to complete the conversion
            vTemp = _mm_add_ps(vTemp, g_AddByte4);
            // Fix y, z and w because they are too large
            vTemp = _mm_mul_ps(vTemp, LoadByteN4Mul);
            // Clamp result (for case of -128)
            return _mm_max_ps(vTemp, g_NegativeOne);
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadByte4(const BYTE4* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR_F32 vResult = { { {
                    static_cast<float>(pSource->x),
                    static_cast<float>(pSource->y),
                    static_cast<float>(pSource->z),
                    static_cast<float>(pSource->w)
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 LoadByte4Mul = { { { 1.0f, 1.0f / 256.0f, 1.0f / 65536.0f, 1.0f / (65536.0f * 256.0f) } } };
            // Splat the color in all four entries (x,z,y,w)
            VECTOR vTemp = _mm_load1_ps(reinterpret_cast<const float*>(&pSource->x));
            // Mask x&0ff,y&0xff00,z&0xff0000,w&0xff000000
            vTemp = _mm_and_ps(vTemp, g_MaskByte4);
            // x,y and z are unsigned! Flip the bits to convert the order to signed
            vTemp = _mm_xor_ps(vTemp, g_XorByte4);
            // Convert to floating point numbers
            vTemp = _mm_cvtepi32_ps(_mm_castps_si128(vTemp));
            // x, y and z - 0x80 to complete the conversion
            vTemp = _mm_add_ps(vTemp, g_AddByte4);
            // Fix y, z and w because they are too large
            vTemp = _mm_mul_ps(vTemp, LoadByte4Mul);
            return vTemp;
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadUNibble4(const UNIBBLE4* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR_F32 vResult = { { {
                    float(pSource->v & 0xF),
                    float((pSource->v >> 4) & 0xF),
                    float((pSource->v >> 8) & 0xF),
                    float((pSource->v >> 12) & 0xF)
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_I32 UNibble4And = { { { 0xF, 0xF0, 0xF00, 0xF000 } } };
            static const VECTOR_F32 UNibble4Mul = { { { 1.0f, 1.0f / 16.f, 1.0f / 256.f, 1.0f / 4096.f } } };
            // Get the 16 bit value and splat it
            __m128i vInt = LOADU_SI16(&pSource->v);
            VECTOR vResult = PERMUTE_PS(_mm_castsi128_ps(vInt), _MM_SHUFFLE(0,0,0,0));
            // Mask off x, y and z
            vResult = _mm_and_ps(vResult, UNibble4And);
            // Convert to float
            vResult = _mm_cvtepi32_ps(_mm_castps_si128(vResult));
            // Normalize x, y, and z
            vResult = _mm_mul_ps(vResult, UNibble4Mul);
            return vResult;
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE VECTOR VEC_CALLCONV LoadU555(const U555* pSource) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSource != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR_F32 vResult = { { {
                    float(pSource->v & 0x1F),
                    float((pSource->v >> 5) & 0x1F),
                    float((pSource->v >> 10) & 0x1F),
                    float((pSource->v >> 15) & 0x1)
                } } };
            return vResult.v;

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_I32 U555And = { { { 0x1F, 0x1F << 5, 0x1F << 10, 0x8000 } } };
            static const VECTOR_F32 U555Mul = { { { 1.0f, 1.0f / 32.f, 1.0f / 1024.f, 1.0f / 32768.f } } };
            // Get the 16bit value and splat it
            __m128i vInt = LOADU_SI16(&pSource->v);
            VECTOR vResult = PERMUTE_PS(_mm_castsi128_ps(vInt), _MM_SHUFFLE(0, 0, 0, 0));
            // Mask off x, y and z
            vResult = _mm_and_ps(vResult, U555And);
            // Convert to float
            vResult = _mm_cvtepi32_ps(_mm_castps_si128(vResult));
            // Normalize x, y, and z
            vResult = _mm_mul_ps(vResult, U555Mul);
            return vResult;
#endif
        }

#ifdef _PREFAST_
#pragma prefast(pop)
#endif

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreColor(COLOR* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Saturate(v);
            N = VEC::Multiply(N, g_UByteMax);
            N = VEC::Round(N);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->c = (static_cast<uint32_t>(tmp.w) << 24) |
                (static_cast<uint32_t>(tmp.x) << 16) |
                (static_cast<uint32_t>(tmp.y) << 8) |
                static_cast<uint32_t>(tmp.z);

#elif defined(_SSE2_INTRINSICS_)
            // Set <0 to 0
            VECTOR vResult = _mm_max_ps(v, g_Zero);
            // Set>1 to 1
            vResult = _mm_min_ps(vResult, g_One);
            // Convert to 0-255
            vResult = _mm_mul_ps(vResult, g_UByteMax);
            // Shuffle RGBA to ARGB
            vResult = PERMUTE_PS(vResult, _MM_SHUFFLE(3, 0, 1, 2));
            // Convert to int
            __m128i vInt = _mm_cvtps_epi32(vResult);
            // Mash to shorts
            vInt = _mm_packs_epi32(vInt, vInt);
            // Mash to bytes
            vInt = _mm_packus_epi16(vInt, vInt);
            // Store the color
            _mm_store_ss(reinterpret_cast<float*>(&pDestination->c), _mm_castsi128_ps(vInt));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreHalf2(HALF2* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

#if defined(_F16C_INTRINSICS_) && !defined(_NO_INTRINSICS_)
            __m128i V1 = _mm_cvtps_ph(v, _MM_FROUND_TO_NEAREST_INT);
            _mm_store_ss(reinterpret_cast<float*>(pDestination), _mm_castsi128_ps(V1));
#else
            pDestination->x = ConvertFloatToHalf(VEC::GetX(v));
            pDestination->y = ConvertFloatToHalf(VEC::GetY(v));
#endif // !_F16C_INTRINSICS_
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreShortN2(SHORTN2* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Clamp(v, g_NegativeOne.v, g_One.v);
            N = VEC::Multiply(N, g_ShortMax);
            N = VEC::Round(N);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->x = static_cast<int16_t>(tmp.x);
            pDestination->y = static_cast<int16_t>(tmp.y);

#elif defined(_SSE2_INTRINSICS_)
            VECTOR vResult = _mm_max_ps(v, g_NegativeOne);
            vResult = _mm_min_ps(vResult, g_One);
            vResult = _mm_mul_ps(vResult, g_ShortMax);
            __m128i vResulti = _mm_cvtps_epi32(vResult);
            vResulti = _mm_packs_epi32(vResulti, vResulti);
            _mm_store_ss(reinterpret_cast<float*>(&pDestination->x), _mm_castsi128_ps(vResulti));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreShort2(SHORT2* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Clamp(v, g_ShortMin, g_ShortMax);
            N = VEC::Round(N);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->x = static_cast<int16_t>(tmp.x);
            pDestination->y = static_cast<int16_t>(tmp.y);

#elif defined(_SSE2_INTRINSICS_)
            // Bounds check
            VECTOR vResult = _mm_max_ps(v, g_ShortMin);
            vResult = _mm_min_ps(vResult, g_ShortMax);
            // Convert to int with rounding
            __m128i vInt = _mm_cvtps_epi32(vResult);
            // Pack the ints into shorts
            vInt = _mm_packs_epi32(vInt, vInt);
            _mm_store_ss(reinterpret_cast<float*>(&pDestination->x), _mm_castsi128_ps(vInt));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreUShortN2(USHORTN2* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Saturate(v);
            N = VEC::MultiplyAdd(N, g_UShortMax, g_OneHalf.v);
            N = VEC::Truncate(N);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->x = static_cast<uint16_t>(tmp.x);
            pDestination->y = static_cast<uint16_t>(tmp.y);

#elif defined(_SSE2_INTRINSICS_)
            // Bounds check
            VECTOR vResult = _mm_max_ps(v, g_Zero);
            vResult = _mm_min_ps(vResult, g_One);
            vResult = _mm_mul_ps(vResult, g_UShortMax);
            vResult = _mm_add_ps(vResult, g_OneHalf);
            // Convert to int
            __m128i vInt = _mm_cvttps_epi32(vResult);
            // Since the SSE pack instruction clamps using signed rules,
            // manually extract the values to store them to memory
            pDestination->x = static_cast<uint16_t>(_mm_extract_epi16(vInt, 0));
            pDestination->y = static_cast<uint16_t>(_mm_extract_epi16(vInt, 2));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreUShort2(USHORT2* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Clamp(v, VEC::Zero(), g_UShortMax);
            N = VEC::Round(N);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->x = static_cast<uint16_t>(tmp.x);
            pDestination->y = static_cast<uint16_t>(tmp.y);

#elif defined(_SSE2_INTRINSICS_)
            // Bounds check
            VECTOR vResult = _mm_max_ps(v, g_Zero);
            vResult = _mm_min_ps(vResult, g_UShortMax);
            // Convert to int with rounding
            __m128i vInt = _mm_cvtps_epi32(vResult);
            // Since the SSE pack instruction clamps using signed rules,
            // manually extract the values to store them to memory
            pDestination->x = static_cast<uint16_t>(_mm_extract_epi16(vInt, 0));
            pDestination->y = static_cast<uint16_t>(_mm_extract_epi16(vInt, 2));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreByteN2(BYTEN2* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Clamp(v, g_NegativeOne.v, g_One.v);
            N = VEC::Multiply(N, g_ByteMax);
            N = VEC::Round(N);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->x = static_cast<int8_t>(tmp.x);
            pDestination->y = static_cast<int8_t>(tmp.y);

#elif defined(_SSE2_INTRINSICS_)
            // Clamp to bounds
            VECTOR vResult = _mm_max_ps(v, g_NegativeOne);
            vResult = _mm_min_ps(vResult, g_One);
            // Scale by multiplication
            vResult = _mm_mul_ps(vResult, g_ByteMax);
            // Convert to int by rounding
            __m128i vInt = _mm_cvtps_epi32(vResult);
            // No SSE operations will write to 16-bit values, so we have to extract them manually
            auto x = static_cast<uint16_t>(_mm_extract_epi16(vInt, 0));
            auto y = static_cast<uint16_t>(_mm_extract_epi16(vInt, 2));
            pDestination->v = static_cast<uint16_t>(((static_cast<int>(y) & 0xFF) << 8) | (static_cast<int>(x) & 0xFF));
#endif
        }
        
        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreByte2(BYTE2* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Clamp(v, g_ByteMin, g_ByteMax);
            N = VEC::Round(N);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->x = static_cast<int8_t>(tmp.x);
            pDestination->y = static_cast<int8_t>(tmp.y);

#elif defined(_SSE2_INTRINSICS_)
            // Clamp to bounds
            VECTOR vResult = _mm_max_ps(v, g_ByteMin);
            vResult = _mm_min_ps(vResult, g_ByteMax);
            // Convert to int by rounding
            __m128i vInt = _mm_cvtps_epi32(vResult);
            // No SSE operations will write to 16-bit values, so we have to extract them manually
            auto x = static_cast<uint16_t>(_mm_extract_epi16(vInt, 0));
            auto y = static_cast<uint16_t>(_mm_extract_epi16(vInt, 2));
            pDestination->v = static_cast<uint16_t>(((static_cast<int>(y) & 0xFF) << 8) | (static_cast<int>(x) & 0xFF));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreUByteN2(UBYTEN2* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Saturate(v);
            N = VEC::MultiplyAdd(N, g_UByteMax, g_OneHalf.v);
            N = VEC::Truncate(N);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->x = static_cast<uint8_t>(tmp.x);
            pDestination->y = static_cast<uint8_t>(tmp.y);

#elif defined(_SSE2_INTRINSICS_)
            // Clamp to bounds
            VECTOR vResult = _mm_max_ps(v, g_Zero);
            vResult = _mm_min_ps(vResult, g_One);
            // Scale by multiplication
            vResult = _mm_mul_ps(vResult, g_UByteMax);
            vResult = _mm_add_ps(vResult, g_OneHalf);
            // Convert to int
            __m128i vInt = _mm_cvttps_epi32(vResult);
            // No SSE operations will write to 16-bit values, so we have to extract them manually
            auto x = static_cast<uint16_t>(_mm_extract_epi16(vInt, 0));
            auto y = static_cast<uint16_t>(_mm_extract_epi16(vInt, 2));
            pDestination->v = static_cast<uint16_t>(((static_cast<int>(y) & 0xFF) << 8) | (static_cast<int>(x) & 0xFF));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreUByte2(UBYTE2* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Clamp(v, VEC::Zero(), g_UByteMax);
            N = VEC::Round(N);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->x = static_cast<uint8_t>(tmp.x);
            pDestination->y = static_cast<uint8_t>(tmp.y);

#elif defined(_SSE2_INTRINSICS_)
            // Clamp to bounds
            VECTOR vResult = _mm_max_ps(v, g_Zero);
            vResult = _mm_min_ps(vResult, g_UByteMax);
            // Convert to int by rounding
            __m128i vInt = _mm_cvtps_epi32(vResult);
            // No SSE operations will write to 16-bit values, so we have to extract them manually
            auto x = static_cast<uint16_t>(_mm_extract_epi16(vInt, 0));
            auto y = static_cast<uint16_t>(_mm_extract_epi16(vInt, 2));
            pDestination->v = static_cast<uint16_t>(((static_cast<int>(y) & 0xFF) << 8) | (static_cast<int>(x) & 0xFF));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreU565(U565* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

            static const VECTOR_F32 Max = { { { 31.0f, 63.0f, 31.0f, 0.0f } } };

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Clamp(v, VEC::Zero(), Max.v);
            N = VEC::Round(N);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->v = static_cast<uint16_t>(
                ((static_cast<int>(tmp.z) & 0x1F) << 11)
                | ((static_cast<int>(tmp.y) & 0x3F) << 5)
                | ((static_cast<int>(tmp.x) & 0x1F)));

#elif defined(_SSE2_INTRINSICS_)
            // Bounds check
            VECTOR vResult = _mm_max_ps(v, g_Zero);
            vResult = _mm_min_ps(vResult, Max);
            // Convert to int with rounding
            __m128i vInt = _mm_cvtps_epi32(vResult);
            // No SSE operations will write to 16-bit values, so we have to extract them manually
            auto x = static_cast<uint16_t>(_mm_extract_epi16(vInt, 0));
            auto y = static_cast<uint16_t>(_mm_extract_epi16(vInt, 2));
            auto z = static_cast<uint16_t>(_mm_extract_epi16(vInt, 4));
            pDestination->v = static_cast<uint16_t>(
                ((static_cast<int>(z) & 0x1F) << 11)
                | ((static_cast<int>(y) & 0x3F) << 5)
                | ((static_cast<int>(x) & 0x1F)));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreFloat3PK(FLOAT3PK* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEDUG

            ALIGNED(16) uint32_t IValue[4];
            StoreFloat3A(reinterpret_cast<FLOAT3A*>(&IValue), v);

            uint32_t Result[3];

            // X & Y Channels (5-bit exponent, 6-bit mantissa)
            for (uint32_t j = 0; j < 2; ++j)
            {
                uint32_t Sign = IValue[j] & 0x80000000;
                uint32_t I = IValue[j] & 0x7FFFFFFF;

                if ((I & 0x7F800000) == 0x7F800000)
                {
                    // INF or NAN
                    Result[j] = 0x7C0U;
                    if ((I & 0x7FFFFF) != 0)
                    {
                        Result[j] = 0x7FFU;
                    }
                    else if (Sign)
                    {
                        // -INF is clamped to 0 since 3PK is positive only
                        Result[j] = 0;
                    }
                }
                else if (Sign || I < 0x35800000)
                {
                    // 3PK is positive only, so clamp to zero
                    Result[j] = 0;
                }
                else if (I > 0x477E0000U)
                {
                    // The number is too large to be represented as a float11, set to max
                    Result[j] = 0x7BFU;
                }
                else
                {
                    if (I < 0x38800000U)
                    {
                        // The number is too small to be represented as a normalized float11
                        // Convert it to a denormalized value.
                        uint32_t Shift = 113U - (I >> 23U);
                        I = (0x800000U | (I & 0x7FFFFFU)) >> Shift;
                    }
                    else
                    {
                        // Rebias the exponent to represent the value as a normalized float11
                        I += 0xC8000000U;
                    }

                    Result[j] = ((I + 0xFFFFU + ((I >> 17U) & 1U)) >> 17U) & 0x7ffU;
                }
            }

            // Z Channel (5-bit exponent, 5-bit mantissa)
            uint32_t Sign = IValue[2] & 0x80000000;
            uint32_t I = IValue[2] & 0x7FFFFFFF;

            if ((I & 0x7F800000) == 0x7F800000)
            {
                // INF or NAN
                Result[2] = 0x3E0U;
                if (I & 0x7FFFFF)
                {
                    Result[2] = 0x3FFU;
                }
                else if (Sign || I < 0x36000000)
                {
                    // -INF is clamped to 0 since 3PK is positive only
                    Result[2] = 0;
                }
            }
            else if (Sign)
            {
                // 3PK is positive only, so clamp to zero
                Result[2] = 0;
            }
            else if (I > 0x477C0000U)
            {
                // The number is too large to be represented as a float10, set to max
                Result[2] = 0x3DFU;
            }
            else
            {
                if (I < 0x38800000U)
                {
                    // The number is too small to be represented as a normalized float10
                    // Convert it to a denormalized value.
                    uint32_t Shift = 113U - (I >> 23U);
                    I = (0x800000U | (I & 0x7FFFFFU)) >> Shift;
                }
                else
                {
                    // Rebias the exponent to represent the value as a normalized float10
                    I += 0xC8000000U;
                }

                Result[2] = ((I + 0x1FFFFU + ((I >> 18U) & 1U)) >> 18U) & 0x3ffU;
            }

            // Pack Result into memory
            pDestination->v = (Result[0] & 0x7ff)
                | ((Result[1] & 0x7ff) << 11)
                | ((Result[2] & 0x3ff) << 22);
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreFloat3SE(FLOAT3SE* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

            FLOAT3A tmp;
            StoreFloat3A(&tmp, v);

            static constexpr float maxf9 = float(0x1FF << 7);
            static constexpr float minf9 = float(1.f / (1 << 16));

            float x = (tmp.x >= 0.f) ? ((tmp.x > maxf9) ? maxf9 : tmp.x) : 0.f;
            float y = (tmp.y >= 0.f) ? ((tmp.y > maxf9) ? maxf9 : tmp.y) : 0.f;
            float z = (tmp.z >= 0.f) ? ((tmp.z > maxf9) ? maxf9 : tmp.z) : 0.f;

            const float max_xy = (x > y) ? x : y;
            const float max_xyz = (max_xy > z) ? max_xy : z;

            const float maxColor = (max_xyz > minf9) ? max_xyz : minf9;

            union { float f; int32_t i; } fi;
            fi.f = maxColor;
            fi.i += 0x00004000; // round up leaving 9 bits in fraction (including assumed 1)

            auto exp = static_cast<uint32_t>(fi.i) >> 23;
            pDestination->e = exp - 0x6f;

            fi.i = static_cast<int32_t>(0x83000000 - (exp << 23));
            float ScaleR = fi.f;

            pDestination->xm = static_cast<uint32_t>(round_to_nearest(x * ScaleR));
            pDestination->ym = static_cast<uint32_t>(round_to_nearest(y * ScaleR));
            pDestination->zm = static_cast<uint32_t>(round_to_nearest(z * ScaleR));
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreHalf4(HALF4* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

#if defined(_F16C_INTRINSICS_) && !defined(_NO_INTRINSICS_)
            __m128i V1 = _mm_cvtps_ph(v, _MM_FROUND_TO_NEAREST_INT);
            _mm_storel_epi64(reinterpret_cast<__m128i*>(pDestination), V1);
#else
            FLOAT4A t;
            StoreFloat4A(&t, v);

            pDestination->x = ConvertFloatToHalf(t.x);
            pDestination->y = ConvertFloatToHalf(t.y);
            pDestination->z = ConvertFloatToHalf(t.z);
            pDestination->w = ConvertFloatToHalf(t.w);
#endif // !_F16C_INTRINSICS_
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreShortN4(SHORTN4* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Clamp(v, g_NegativeOne.v, g_One.v);
            N = VEC::Multiply(N, g_ShortMax);
            N = VEC::Round(N);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->x = static_cast<int16_t>(tmp.x);
            pDestination->y = static_cast<int16_t>(tmp.y);
            pDestination->z = static_cast<int16_t>(tmp.z);
            pDestination->w = static_cast<int16_t>(tmp.w);

#elif defined(_SSE2_INTRINSICS_)
            VECTOR vResult = _mm_max_ps(v, g_NegativeOne);
            vResult = _mm_min_ps(vResult, g_One);
            vResult = _mm_mul_ps(vResult, g_ShortMax);
            __m128i vResulti = _mm_cvtps_epi32(vResult);
            vResulti = _mm_packs_epi32(vResulti, vResulti);
            _mm_store_sd(reinterpret_cast<double*>(&pDestination->x), _mm_castsi128_pd(vResulti));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreShort4(SHORT4* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Clamp(v, g_ShortMin, g_ShortMax);
            N = VEC::Round(N);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->x = static_cast<int16_t>(tmp.x);
            pDestination->y = static_cast<int16_t>(tmp.y);
            pDestination->z = static_cast<int16_t>(tmp.z);
            pDestination->w = static_cast<int16_t>(tmp.w);

#elif defined(_SSE2_INTRINSICS_)
            // Bounds check
            VECTOR vResult = _mm_max_ps(v, g_ShortMin);
            vResult = _mm_min_ps(vResult, g_ShortMax);
            // Convert to int with rounding
            __m128i vInt = _mm_cvtps_epi32(vResult);
            // Pack the ints into shorts
            vInt = _mm_packs_epi32(vInt, vInt);
            _mm_store_sd(reinterpret_cast<double*>(&pDestination->x), _mm_castsi128_pd(vInt));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreUShortN4(USHORTN4* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Saturate(v);
            N = VEC::MultiplyAdd(N, g_UShortMax, g_OneHalf.v);
            N = VEC::Truncate(N);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->x = static_cast<uint16_t>(tmp.x);
            pDestination->y = static_cast<uint16_t>(tmp.y);
            pDestination->z = static_cast<uint16_t>(tmp.z);
            pDestination->w = static_cast<uint16_t>(tmp.w);

#elif defined(_SSE2_INTRINSICS_)
            // Bounds check
            VECTOR vResult = _mm_max_ps(v, g_Zero);
            vResult = _mm_min_ps(vResult, g_One);
            vResult = _mm_mul_ps(vResult, g_UShortMax);
            vResult = _mm_add_ps(vResult, g_OneHalf);
            // Convert to int
            __m128i vInt = _mm_cvttps_epi32(vResult);
            // Since the SSE pack instruction clamps using signed rules,
            // manually extract the values to store them to memory
            pDestination->x = static_cast<uint16_t>(_mm_extract_epi16(vInt, 0));
            pDestination->y = static_cast<uint16_t>(_mm_extract_epi16(vInt, 2));
            pDestination->z = static_cast<uint16_t>(_mm_extract_epi16(vInt, 4));
            pDestination->w = static_cast<uint16_t>(_mm_extract_epi16(vInt, 6));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreUShort4(USHORT4* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Clamp(v, VEC::Zero(), g_UShortMax);
            N = VEC::Round(N);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->x = static_cast<uint16_t>(tmp.x);
            pDestination->y = static_cast<uint16_t>(tmp.y);
            pDestination->z = static_cast<uint16_t>(tmp.z);
            pDestination->w = static_cast<uint16_t>(tmp.w);

#elif defined(_SSE2_INTRINSICS_)
            // Bounds check
            VECTOR vResult = _mm_max_ps(v, g_Zero);
            vResult = _mm_min_ps(vResult, g_UShortMax);
            // Convert to int with rounding
            __m128i vInt = _mm_cvtps_epi32(vResult);
            // Since the SSE pack instruction clamps using signed rules,
            // manually extract the values to store them to memory
            pDestination->x = static_cast<uint16_t>(_mm_extract_epi16(vInt, 0));
            pDestination->y = static_cast<uint16_t>(_mm_extract_epi16(vInt, 2));
            pDestination->z = static_cast<uint16_t>(_mm_extract_epi16(vInt, 4));
            pDestination->w = static_cast<uint16_t>(_mm_extract_epi16(vInt, 6));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreXDecN4(XDECN4* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

            static const VECTOR_F32 Min = { { { -1.0f, -1.0f, -1.0f, 0.0f } } };

#if defined(_NO_INTRINSICS_)

            static const VECTOR_F32  Scale = { { { 511.0f, 511.0f, 511.0f, 3.0f } } };

            VECTOR N = VEC::Clamp(v, Min.v, g_One.v);
            N = VEC::Multiply(N, Scale.v);
            N = VEC::Round(N);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->v = static_cast<uint32_t>(
                (static_cast<int>(tmp.w) << 30)
                | ((static_cast<int>(tmp.z) & 0x3FF) << 20)
                | ((static_cast<int>(tmp.y) & 0x3FF) << 10)
                | (static_cast<int>(tmp.x) & 0x3FF));

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 Scale = { { { 511.0f, 511.0f * 1024.0f, 511.0f * 1048576.0f, 3.0f * 536870912.0f } } };
            static const VECTOR_I32 ScaleMask = { { { 0x3FF, 0x3FF << 10, 0x3FF << 20, 0x3 << 29 } } };
            VECTOR vResult = _mm_max_ps(v, Min);
            vResult = _mm_min_ps(vResult, g_One);
            // Scale by multiplication
            vResult = _mm_mul_ps(vResult, Scale);
            // Convert to int (W is unsigned)
            __m128i vResulti = _mm_cvtps_epi32(vResult);
            // Mask off any fraction
            vResulti = _mm_and_si128(vResulti, ScaleMask);
            // To fix W, add itself to shift it up to <<30 instead of <<29
            __m128i vResultw = _mm_and_si128(vResulti, g_MaskW);
            vResulti = _mm_add_epi32(vResulti, vResultw);
            // Do a horizontal or of all 4 entries
            vResult = PERMUTE_PS(_mm_castsi128_ps(vResulti), _MM_SHUFFLE(0, 3, 2, 1));
            vResulti = _mm_or_si128(vResulti, _mm_castps_si128(vResult));
            vResult = PERMUTE_PS(vResult, _MM_SHUFFLE(0, 3, 2, 1));
            vResulti = _mm_or_si128(vResulti, _mm_castps_si128(vResult));
            vResult = PERMUTE_PS(vResult, _MM_SHUFFLE(0, 3, 2, 1));
            vResulti = _mm_or_si128(vResulti, _mm_castps_si128(vResult));
            _mm_store_ss(reinterpret_cast<float*>(&pDestination->v), _mm_castsi128_ps(vResulti));
#endif
        }

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4996)
// C4996: ignore deprecation warning
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreXDec4(XDEC4* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

            static const VECTOR_F32 MinXDec4 = { { { -511.0f, -511.0f, -511.0f, 0.0f } } };
            static const VECTOR_F32 MaxXDec4 = { { { 511.0f, 511.0f, 511.0f, 3.0f } } };

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Clamp(v, MinXDec4, MaxXDec4);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->v = static_cast<uint32_t>(
                (static_cast<int>(tmp.w) << 30)
                | ((static_cast<int>(tmp.z) & 0x3FF) << 20)
                | ((static_cast<int>(tmp.y) & 0x3FF) << 10)
                | ((static_cast<int>(tmp.x) & 0x3FF)));

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 ScaleXDec4 = { { { 1.0f, 1024.0f / 2.0f, 1024.0f * 1024.0f, 1024.0f * 1024.0f * 1024.0f / 2.0f } } };
            static const VECTOR_I32 MaskXDec4 = { { { 0x3FF, 0x3FF << (10 - 1), 0x3FF << 20, 0x3 << (30 - 1) } } };
            // Clamp to bounds
            VECTOR vResult = _mm_max_ps(v, MinXDec4);
            vResult = _mm_min_ps(vResult, MaxXDec4);
            // Scale by multiplication
            vResult = _mm_mul_ps(vResult, ScaleXDec4);
            // Convert to int
            __m128i vResulti = _mm_cvttps_epi32(vResult);
            // Mask off any fraction
            vResulti = _mm_and_si128(vResulti, MaskXDec4);
            // Do a horizontal or of 4 entries
            __m128i vResulti2 = _mm_shuffle_epi32(vResulti, _MM_SHUFFLE(3, 2, 3, 2));
            // x = x|z, y = y|w
            vResulti = _mm_or_si128(vResulti, vResulti2);
            // Move Z to the x position
            vResulti2 = _mm_shuffle_epi32(vResulti, _MM_SHUFFLE(1, 1, 1, 1));
            // Perform a single bit left shift on y|w
            vResulti2 = _mm_add_epi32(vResulti2, vResulti2);
            // i = x|y|z|w
            vResulti = _mm_or_si128(vResulti, vResulti2);
            _mm_store_ss(reinterpret_cast<float*>(&pDestination->v), _mm_castsi128_ps(vResulti));
#endif
        }

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreUDecN4(UDECN4* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            static const VECTOR_F32  Scale = { { { 1023.0f, 1023.0f, 1023.0f, 3.0f } } };

            VECTOR N = VEC::Saturate(v);
            N = VEC::Multiply(N, Scale.v);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->v = static_cast<uint32_t>(
                (static_cast<int>(tmp.w) << 30)
                | ((static_cast<int>(tmp.z) & 0x3FF) << 20)
                | ((static_cast<int>(tmp.y) & 0x3FF) << 10)
                | ((static_cast<int>(tmp.x) & 0x3FF)));

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 ScaleUDecN4 = { { { 1023.0f, 1023.0f * 1024.0f * 0.5f, 1023.0f * 1024.0f * 1024.0f, 3.0f * 1024.0f * 1024.0f * 1024.0f * 0.5f } } };
            static const VECTOR_I32 MaskUDecN4 = { { { 0x3FF, 0x3FF << (10 - 1), 0x3FF << 20, 0x3 << (30 - 1) } } };
            // Clamp to bounds
            VECTOR vResult = _mm_max_ps(v, g_Zero);
            vResult = _mm_min_ps(vResult, g_One);
            // Scale by multiplication
            vResult = _mm_mul_ps(vResult, ScaleUDecN4);
            // Convert to int
            __m128i vResulti = _mm_cvttps_epi32(vResult);
            // Mask off any fraction
            vResulti = _mm_and_si128(vResulti, MaskUDecN4);
            // Do a horizontal or of 4 entries
            __m128i vResulti2 = _mm_shuffle_epi32(vResulti, _MM_SHUFFLE(3, 2, 3, 2));
            // x = x|z, y = y|w
            vResulti = _mm_or_si128(vResulti, vResulti2);
            // Move Z to the x position
            vResulti2 = _mm_shuffle_epi32(vResulti, _MM_SHUFFLE(1, 1, 1, 1));
            // Perform a left shift by one bit on y|w
            vResulti2 = _mm_add_epi32(vResulti2, vResulti2);
            // i = x|y|z|w
            vResulti = _mm_or_si128(vResulti, vResulti2);
            _mm_store_ss(reinterpret_cast<float*>(&pDestination->v), _mm_castsi128_ps(vResulti));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreUDecN4_XR(UDECN4* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination);
#endif // DEBUG

            static const VECTOR_F32 Scale = { { { 510.0f, 510.0f, 510.0f, 3.0f } } };
            static const VECTOR_F32 Bias = { { { 384.0f, 384.0f, 384.0f, 0.0f } } };
            static const VECTOR_F32 C = { { { 1023.f, 1023.f, 1023.f, 3.f } } };

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::MultiplyAdd(v, Scale, Bias);
            N = VEC::Clamp(N, g_Zero, C);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->v = static_cast<uint32_t>(
                (static_cast<int>(tmp.w) << 30)
                | ((static_cast<int>(tmp.z) & 0x3FF) << 20)
                | ((static_cast<int>(tmp.y) & 0x3FF) << 10)
                | ((static_cast<int>(tmp.x) & 0x3FF)));

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 Shift = { { { 1.0f, 1024.0f * 0.5f, 1024.0f * 1024.0f, 1024.0f * 1024.0f * 1024.0f * 0.5f } } };
            static const VECTOR_U32 MaskUDecN4 = { { { 0x3FF, 0x3FF << (10 - 1), 0x3FF << 20, 0x3 << (30 - 1) } } };
            // Scale & bias
            VECTOR vResult = FMADD_PS(v, Scale, Bias);
            // Clamp to bounds
            vResult = _mm_max_ps(vResult, g_Zero);
            vResult = _mm_min_ps(vResult, C);
            // Scale by shift values
            vResult = _mm_mul_ps(vResult, Shift);
            // Convert to int
            __m128i vResulti = _mm_cvttps_epi32(vResult);
            // Mask off any fraction
            vResulti = _mm_and_si128(vResulti, MaskUDecN4);
            // Do a horizontal or of 4 entries
            __m128i vResulti2 = _mm_shuffle_epi32(vResulti, _MM_SHUFFLE(3, 2, 3, 2));
            // x = x|z, y = y|w
            vResulti = _mm_or_si128(vResulti, vResulti2);
            // Move Z to the x position
            vResulti2 = _mm_shuffle_epi32(vResulti, _MM_SHUFFLE(1, 1, 1, 1));
            // Perform a left shift by one bit on y|w
            vResulti2 = _mm_add_epi32(vResulti2, vResulti2);
            // i = x|y|z|w
            vResulti = _mm_or_si128(vResulti, vResulti2);
            _mm_store_ss(reinterpret_cast<float*>(&pDestination->v), _mm_castsi128_ps(vResulti));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreUDec4(UDEC4* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

            static const VECTOR_F32 MaxUDec4 = { { { 1023.0f, 1023.0f, 1023.0f, 3.0f } } };

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Clamp(v, VEC::Zero(), MaxUDec4);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->v = static_cast<uint32_t>(
                (static_cast<int>(tmp.w) << 30)
                | ((static_cast<int>(tmp.z) & 0x3FF) << 20)
                | ((static_cast<int>(tmp.y) & 0x3FF) << 10)
                | ((static_cast<int>(tmp.x) & 0x3FF)));

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 ScaleUDec4 = { { { 1.0f, 1024.0f / 2.0f, 1024.0f * 1024.0f, 1024.0f * 1024.0f * 1024.0f / 2.0f } } };
            static const VECTOR_I32 MaskUDec4 = { { { 0x3FF, 0x3FF << (10 - 1), 0x3FF << 20, 0x3 << (30 - 1) } } };
            // Clamp to bounds
            VECTOR vResult = _mm_max_ps(v, g_Zero);
            vResult = _mm_min_ps(vResult, MaxUDec4);
            // Scale by multiplication
            vResult = _mm_mul_ps(vResult, ScaleUDec4);
            // Convert to int
            __m128i vResulti = _mm_cvttps_epi32(vResult);
            // Mask off any fraction
            vResulti = _mm_and_si128(vResulti, MaskUDec4);
            // Do a horizontal or of 4 entries
            __m128i vResulti2 = _mm_shuffle_epi32(vResulti, _MM_SHUFFLE(3, 2, 3, 2));
            // x = x|z, y = y|w
            vResulti = _mm_or_si128(vResulti, vResulti2);
            // Move Z to the x position
            vResulti2 = _mm_shuffle_epi32(vResulti, _MM_SHUFFLE(1, 1, 1, 1));
            // Perform a left shift by one bit on y|w
            vResulti2 = _mm_add_epi32(vResulti2, vResulti2);
            // i = x|y|z|w
            vResulti = _mm_or_si128(vResulti, vResulti2);
            _mm_store_ss(reinterpret_cast<float*>(&pDestination->v), _mm_castsi128_ps(vResulti));
#endif
        }

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4996)
// C4996: ignore deprecation warning
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreDecN4(DECN4* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
             assert(pDestination);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            static const VECTOR_F32 Scale = { { { 511.0f, 511.0f, 511.0f, 1.0f } } };

            VECTOR N = VEC::Clamp(v, g_NegativeOne.v, g_One.v);
            N = VEC::Multiply(N, Scale.v);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->v = static_cast<uint32_t>(
                (static_cast<int>(tmp.w) << 30)
                | ((static_cast<int>(tmp.z) & 0x3FF) << 20)
                | ((static_cast<int>(tmp.y) & 0x3FF) << 10)
                | ((static_cast<int>(tmp.x) & 0x3FF)));

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 ScaleDecN4 = { { { 511.0f, 511.0f * 1024.0f, 511.0f * 1024.0f * 1024.0f, 1.0f * 1024.0f * 1024.0f * 1024.0f } } };
            // Clamp to bounds
            VECTOR vResult = _mm_max_ps(v, g_NegativeOne);
            vResult = _mm_min_ps(vResult, g_One);
            // Scale by multiplication
            vResult = _mm_mul_ps(vResult, ScaleDecN4);
            // Convert to int
            __m128i vResulti = _mm_cvttps_epi32(vResult);
            // Mask off any fraction
            vResulti = _mm_and_si128(vResulti, g_MaskDec4);
            // Do a horizontal or of 4 entries
            __m128i vResulti2 = _mm_shuffle_epi32(vResulti, _MM_SHUFFLE(3, 2, 3, 2));
            // x = x|z, y = y|w
            vResulti = _mm_or_si128(vResulti, vResulti2);
            // Move Z to the x position
            vResulti2 = _mm_shuffle_epi32(vResulti, _MM_SHUFFLE(1, 1, 1, 1));
            // i = x|y|z|w
            vResulti = _mm_or_si128(vResulti, vResulti2);
            _mm_store_ss(reinterpret_cast<float*>(&pDestination->v), _mm_castsi128_ps(vResulti));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreDec4(DEC4* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

            static const VECTOR_F32 MinDec4 = { { { -511.0f, -511.0f, -511.0f, -1.0f } } };
            static const VECTOR_F32 MaxDec4 = { { { 511.0f, 511.0f, 511.0f, 1.0f } } };

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Clamp(v, MinDec4, MaxDec4);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->v = static_cast<uint32_t>(
                (static_cast<int>(tmp.w) << 30)
                | ((static_cast<int>(tmp.z) & 0x3FF) << 20)
                | ((static_cast<int>(tmp.y) & 0x3FF) << 10)
                | ((static_cast<int>(tmp.x) & 0x3FF)));

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 ScaleDec4 = { { { 1.0f, 1024.0f, 1024.0f * 1024.0f, 1024.0f * 1024.0f * 1024.0f } } };
            // Clamp to bounds
            VECTOR vResult = _mm_max_ps(v, MinDec4);
            vResult = _mm_min_ps(vResult, MaxDec4);
            // Scale by multiplication
            vResult = _mm_mul_ps(vResult, ScaleDec4);
            // Convert to int
            __m128i vResulti = _mm_cvttps_epi32(vResult);
            // Mask off any fraction
            vResulti = _mm_and_si128(vResulti, g_MaskDec4);
            // Do a horizontal or of 4 entries
            __m128i vResulti2 = _mm_shuffle_epi32(vResulti, _MM_SHUFFLE(3, 2, 3, 2));
            // x = x|z, y = y|w
            vResulti = _mm_or_si128(vResulti, vResulti2);
            // Move Z to the x position
            vResulti2 = _mm_shuffle_epi32(vResulti, _MM_SHUFFLE(1, 1, 1, 1));
            // i = x|y|z|w
            vResulti = _mm_or_si128(vResulti, vResulti2);
            _mm_store_ss(reinterpret_cast<float*>(&pDestination->v), _mm_castsi128_ps(vResulti));
#endif
        }

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreUByteN4(UBYTEN4* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Saturate(v);
            N = VEC::Multiply(N, g_UByteMax);
            N = VEC::Truncate(N);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->x = static_cast<uint8_t>(tmp.x);
            pDestination->y = static_cast<uint8_t>(tmp.y);
            pDestination->z = static_cast<uint8_t>(tmp.z);
            pDestination->w = static_cast<uint8_t>(tmp.w);

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 ScaleUByteN4 = { { { 255.0f, 255.0f * 256.0f * 0.5f, 255.0f * 256.0f * 256.0f, 255.0f * 256.0f * 256.0f * 256.0f * 0.5f } } };
            static const VECTOR_I32 MaskUByteN4 = { { { 0xFF, 0xFF << (8 - 1), 0xFF << 16, 0xFF << (24 - 1) } } };
            // Clamp to bounds
            VECTOR vResult = _mm_max_ps(v, g_Zero);
            vResult = _mm_min_ps(vResult, g_One);
            // Scale by multiplication
            vResult = _mm_mul_ps(vResult, ScaleUByteN4);
            // Convert to int
            __m128i vResulti = _mm_cvttps_epi32(vResult);
            // Mask off any fraction
            vResulti = _mm_and_si128(vResulti, MaskUByteN4);
            // Do a horizontal or of 4 entries
            __m128i vResulti2 = _mm_shuffle_epi32(vResulti, _MM_SHUFFLE(3, 2, 3, 2));
            // x = x|z, y = y|w
            vResulti = _mm_or_si128(vResulti, vResulti2);
            // Move Z to the x position
            vResulti2 = _mm_shuffle_epi32(vResulti, _MM_SHUFFLE(1, 1, 1, 1));
            // Perform a single bit left shift to fix y|w
            vResulti2 = _mm_add_epi32(vResulti2, vResulti2);
            // i = x|y|z|w
            vResulti = _mm_or_si128(vResulti, vResulti2);
            _mm_store_ss(reinterpret_cast<float*>(&pDestination->v), _mm_castsi128_ps(vResulti));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreUByte4(UBYTE4* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Clamp(v, VEC::Zero(), g_UByteMax);
            N = VEC::Round(N);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->x = static_cast<uint8_t>(tmp.x);
            pDestination->y = static_cast<uint8_t>(tmp.y);
            pDestination->z = static_cast<uint8_t>(tmp.z);
            pDestination->w = static_cast<uint8_t>(tmp.w);

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 ScaleUByte4 = { { { 1.0f, 256.0f * 0.5f, 256.0f * 256.0f, 256.0f * 256.0f * 256.0f * 0.5f } } };
            static const VECTOR_I32 MaskUByte4 = { { { 0xFF, 0xFF << (8 - 1), 0xFF << 16, 0xFF << (24 - 1) } } };
            // Clamp to bounds
            VECTOR vResult = _mm_max_ps(v, g_Zero);
            vResult = _mm_min_ps(vResult, g_UByteMax);
            // Scale by multiplication
            vResult = _mm_mul_ps(vResult, ScaleUByte4);
            // Convert to int by rounding
            __m128i vResulti = _mm_cvtps_epi32(vResult);
            // Mask off any fraction
            vResulti = _mm_and_si128(vResulti, MaskUByte4);
            // Do a horizontal or of 4 entries
            __m128i vResulti2 = _mm_shuffle_epi32(vResulti, _MM_SHUFFLE(3, 2, 3, 2));
            // x = x|z, y = y|w
            vResulti = _mm_or_si128(vResulti, vResulti2);
            // Move Z to the x position
            vResulti2 = _mm_shuffle_epi32(vResulti, _MM_SHUFFLE(1, 1, 1, 1));
            // Perform a single bit left shift to fix y|w
            vResulti2 = _mm_add_epi32(vResulti2, vResulti2);
            // i = x|y|z|w
            vResulti = _mm_or_si128(vResulti, vResulti2);
            _mm_store_ss(reinterpret_cast<float*>(&pDestination->v), _mm_castsi128_ps(vResulti));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreByteN4(BYTEN4* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Clamp(v, g_NegativeOne.v, g_One.v);
            N = VEC::Multiply(N, g_ByteMax);
            N = VEC::Truncate(N);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->x = static_cast<int8_t>(tmp.x);
            pDestination->y = static_cast<int8_t>(tmp.y);
            pDestination->z = static_cast<int8_t>(tmp.z);
            pDestination->w = static_cast<int8_t>(tmp.w);

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 ScaleByteN4 = { { { 127.0f, 127.0f * 256.0f, 127.0f * 256.0f * 256.0f, 127.0f * 256.0f * 256.0f * 256.0f } } };
            static const VECTOR_I32 MaskByteN4 = { { { 0xFF, 0xFF << 8, 0xFF << 16, static_cast<int>(0xFF000000) } } };
            // Clamp to bounds
            VECTOR vResult = _mm_max_ps(v, g_NegativeOne);
            vResult = _mm_min_ps(vResult, g_One);
            // Scale by multiplication
            vResult = _mm_mul_ps(vResult, ScaleByteN4);
            // Convert to int
            __m128i vResulti = _mm_cvttps_epi32(vResult);
            // Mask off any fraction
            vResulti = _mm_and_si128(vResulti, MaskByteN4);
            // Do a horizontal or of 4 entries
            __m128i vResulti2 = _mm_shuffle_epi32(vResulti, _MM_SHUFFLE(3, 2, 3, 2));
            // x = x|z, y = y|w
            vResulti = _mm_or_si128(vResulti, vResulti2);
            // Move Z to the x position
            vResulti2 = _mm_shuffle_epi32(vResulti, _MM_SHUFFLE(1, 1, 1, 1));
            // i = x|y|z|w
            vResulti = _mm_or_si128(vResulti, vResulti2);
            _mm_store_ss(reinterpret_cast<float*>(&pDestination->v), _mm_castsi128_ps(vResulti));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreByte4(BYTE4* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Clamp(v, g_ByteMin, g_ByteMax);
            N = VEC::Round(N);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->x = static_cast<int8_t>(tmp.x);
            pDestination->y = static_cast<int8_t>(tmp.y);
            pDestination->z = static_cast<int8_t>(tmp.z);
            pDestination->w = static_cast<int8_t>(tmp.w);

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 ScaleByte4 = { { { 1.0f, 256.0f, 256.0f * 256.0f, 256.0f * 256.0f * 256.0f } } };
            static const VECTOR_I32 MaskByte4 = { { { 0xFF, 0xFF << 8, 0xFF << 16, static_cast<int>(0xFF000000) } } };
            // Clamp to bounds
            VECTOR vResult = _mm_max_ps(v, g_ByteMin);
            vResult = _mm_min_ps(vResult, g_ByteMax);
            // Scale by multiplication
            vResult = _mm_mul_ps(vResult, ScaleByte4);
            // Convert to int by rounding
            __m128i vResulti = _mm_cvtps_epi32(vResult);
            // Mask off any fraction
            vResulti = _mm_and_si128(vResulti, MaskByte4);
            // Do a horizontal or of 4 entries
            __m128i vResulti2 = _mm_shuffle_epi32(vResulti, _MM_SHUFFLE(3, 2, 3, 2));
            // x = x|z, y = y|w
            vResulti = _mm_or_si128(vResulti, vResulti2);
            // Move Z to the x position
            vResulti2 = _mm_shuffle_epi32(vResulti, _MM_SHUFFLE(1, 1, 1, 1));
            // i = x|y|z|w
            vResulti = _mm_or_si128(vResulti, vResulti2);
            _mm_store_ss(reinterpret_cast<float*>(&pDestination->v), _mm_castsi128_ps(vResulti));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreUNibble4(UNIBBLE4* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination != nullptr);
#endif // DEBUG

            static const VECTOR_F32 Max = { { { 15.0f, 15.0f, 15.0f, 15.0f } } };

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Clamp(v, VEC::Zero(), Max.v);
            N = VEC::Round(N);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->v = static_cast<uint16_t>(
                ((static_cast<int>(tmp.w) & 0xF) << 12)
                | ((static_cast<int>(tmp.z) & 0xF) << 8)
                | ((static_cast<int>(tmp.y) & 0xF) << 4)
                | (static_cast<int>(tmp.x) & 0xF));

#elif defined(_SSE2_INTRINSICS_)
            // Bounds check
            VECTOR vResult = _mm_max_ps(v, g_Zero);
            vResult = _mm_min_ps(vResult, Max);
            // Convert to int with rounding
            __m128i vInt = _mm_cvtps_epi32(vResult);
            // No SSE operations will write to 16-bit values, so we have to extract them manually
            auto x = static_cast<uint16_t>(_mm_extract_epi16(vInt, 0));
            auto y = static_cast<uint16_t>(_mm_extract_epi16(vInt, 2));
            auto z = static_cast<uint16_t>(_mm_extract_epi16(vInt, 4));
            auto w = static_cast<uint16_t>(_mm_extract_epi16(vInt, 6));
            pDestination->v = static_cast<uint16_t>(
                ((static_cast<int>(w) & 0xF) << 12)
                | ((static_cast<int>(z) & 0xF) << 8)
                | ((static_cast<int>(y) & 0xF) << 4)
                | ((static_cast<int>(x) & 0xF)));
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV StoreU555(U555* pDestination, A_VECTOR v) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pDestination);
#endif // DEBUG

            static const XMVECTORF32 Max = { { { 31.0f, 31.0f, 31.0f, 1.0f } } };

#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Clamp(v, VEC::Zero(), Max.v);
            N = VEC::Round(N);

            FLOAT4A tmp;
            StoreFloat4A(&tmp, N);

            pDestination->v = static_cast<uint16_t>(
                ((tmp.w > 0.f) ? 0x8000 : 0)
                | ((static_cast<int>(tmp.z) & 0x1F) << 10)
                | ((static_cast<int>(tmp.y) & 0x1F) << 5)
                | (static_cast<int>(tmp.x) & 0x1F));

#elif defined(_SSE2_INTRINSICS_)
            // Bounds check
            VECTOR vResult = _mm_max_ps(v, g_Zero);
            vResult = _mm_min_ps(vResult, Max);
            // Convert to int with rounding
            __m128i vInt = _mm_cvtps_epi32(vResult);
            // No SSE operations will write to 16-bit values, so we have to extract them manually
            auto x = static_cast<uint16_t>(_mm_extract_epi16(vInt, 0));
            auto y = static_cast<uint16_t>(_mm_extract_epi16(vInt, 2));
            auto z = static_cast<uint16_t>(_mm_extract_epi16(vInt, 4));
            auto w = static_cast<uint16_t>(_mm_extract_epi16(vInt, 6));
            pDestination->v = static_cast<uint16_t>(
                (static_cast<int>(w) ? 0x8000 : 0)
                | ((static_cast<int>(z) & 0x1F) << 10)
                | ((static_cast<int>(y) & 0x1F) << 5)
                | ((static_cast<int>(x) & 0x1F)));
#endif
        }

        FORCE_INLINE COLOR::COLOR(float _r, float _g, float _b. float _a) noexcept
        {
            StoreColor(this, VEC::Set(_r, _g, _b, _a));
        }

        _Use_decl_annotations_
        FORCE_INLINE COLOR::COLOR(const float* pArray) noexcept
        {
            StoreColor(this, VEC::LoadFLoat4(reinterpret_cast<const FLOAT4*>(pArray)));
        }

        FORCE_INLINE HALF2::HALF2(float _x, float _y) noexcept
        {
            x = ConvertFloatToHalf(_x);
            y = ConvertFloatToHalf(_y);
        }

        _Use_decl_annotations_
        FORCE_INLINE HALF2::HALF2(cosnt float* pArray) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pArray != nullptr);
#endif // DEBUG

            x = ConvertFloatToHalf(pArray[0]);
            y = ConvertFloatToHalf(pArray[1]);
        }

        FORCE_INLINE SHORTN2::SHORTN2(float _x, float _y) noexcept
        {
            StoreShortN2(this, VEC::Set(_x, _y, 0.0f, 0.0f));
        }

        _Use_decl_annotations_
        FORCE_INLINE SHORTN2::SHORTN2(const float* pArray) noexcept
        {
            StoreShortN2(this, VEC::LoadFloat2(reinterpret_cast<const FLOAT2*>(pArray)));
        }

        FORCE_INLINE SHORT2::SHORT2(float _x, float _y) noexcept
        {
            StoreShort2(this, VEC::Set(_x, _y, 0.0f, 0.0f));
        }

        _Use_decl_annotations_
        FORCE_INLINE SHORT2::SHORT2(const float* pArray) noexcept
        {
            StoreShort2(this, VEC::LoadFloat2(reinterpret_cast<const FLOAT2*>(pArray)));
        }

        FORCE_INLINE USHORTN2::USHORTN2(float _x, float _y) noexcept
        {
            StoreUShortN2(this, VEC::Set(_x, _y, 0.0f, 0.0f));
        }

        _Use_decl_annotations_
        FORCE_INLINE USHORTN2::USHORTN2(const float* pArray) noexcept
        {
            StoreUShortN2(this, VEC::LoadFloat2(reinterpret_cast<const FLOAT2*>(pArray)));
        }

        FORCE_INLINE USHORT2::USHORT2(float _x, float _y) noexcept
        {
            StoreUShort2(this, VEC::Set(_x, _y, 0.0f, 0.0f));
        }

        _Use_decl_annotations_
        FORCE_INLINE USHORT2::USHORT2(const float* pArray) noexcept
        {
            StoreUShort2(this, VEC::LoadFloat2(reinterpret_cast<const FLOAT2*>(pArray)));
        }

        FORCE_INLINE BYTEN2::BYTEN2(float _x, float _y) noexcept
        {
            StoreByteN2(this, VEC::Set(_x, _y, 0.0f, 0.0f));
        }

        _Use_decl_annotations_
        FORCE_INLINE BYTEN2::BYTEN2(const float* pArray) noexcept
        {
            StoreByteN2(this, VEC::LoadFloat2(reinterpret_cast<const FLOAT2*>(pArray)));
        }

        FORCE_INLINE BYTE2::BYTE2(float _x, float _y) noexcept
        {
            StoreByte2(this, VEC::Set(_x, _y, 0.0f, 0.0f));
        }

        _Use_decl_annotations_
        FORCE_INLINE BYTE2::BYTE2(const float* pArray) noexcept
        {
            StoreByte2(this, VEC::LoadFloat2(reinterpret_cast<const FLOAT2*>(pArray)));
        }

        FORCE_INLINE UBYTEN2::UBYTEN2(float _x, float _y) noexcept
        {
            StoreUByteN2(this, VEC::Set(_x, _y, 0.0f, 0.0f));
        }

        _Use_decl_annotations_
        FORCE_INLINE UBYTEN2::UBYTEN2(const float* pArray) noexcept
        {
            StoreUByteN2(this, VEC::LoadFloat2(reinterpret_cast<const FLOAT2*>(pArray)));
        }

        FORCE_INLINE UBYTE2::UBYTE2(float _x, float _y) noexcept
        {
            StoreUByte2(this, VEC::Set(_x, _y, 0.0f, 0.0f));
        }

        _Use_decl_annotations_
        FORCE_INLINE UBYTE2::UBYTE2(const float* pArray) noexcept
        {
            StoreUByte2(this, VEC::LoadFloat2(reinterpret_cast<const FLOAT2*>(pArray)));
        }

        FORCE_INLINE U565::U565(float _x, float _y, float _z) noexcept
        {
            StoreU565(this, VEC::Set(_x, _y, _z, 0.0f));
        }

        _Use_decl_annotations_
        FORCE_INLINE U565::U565(const float* pArray) noexcept
        {
            StoreU565(this, VEC::LoadFloat3(reinterpret_cast<const FLOAT3*>(pArray)));
        }

        FORCE_INLINE FLOAT3PK::FLOAT3PK(float _x, float _y, float _z) noexcept
        {
            StoreFloat3PK(this, VEC::Set(_x, _y, _z, 0.0f));
        }

        _Use_decl_annotations_
        FORCE_INLINE FLOAT3PK::FLOAT3PK(const float* pArray) noexcept
        {
            StoreFloat3PK(this, VEC::LoadFloat3(reinterpret_cast<const FLOAT3*>(pArray)));
        }

        FORCE_INLINE FLOAT3SE::FLOAT3SE(float _x, float _y, float _z) noexcept
        {
            StoreFloat3SE(this, VEC::Set(_x, _y, _z, 0.0f));
        }

        _Use_decl_annotations_
        FORCE_INLINE FLOAT3SE::FLOAT3SE(const float* pArray) noexcept
        {
            StoreFloat3SE(this, VEC::LoadFloat3(reinterpret_cast<const FLOAT3*>(pArray)));
        }

        FORCE_INLINE HALF4::HALF4(float _X, float _y, float _z, float _w) noexcept
        {
            x = ConvertFloatToHalf(_x);
            y = ConvertFloatToHalf(_y);
            z = ConvertFloatToHalf(_z);
            w = ConvertFloatToHalf(_w);
        }

        _Use_decl_annotations_
        FORCE_INLINE HALF4::HALF4(const float* pArray) noexcept
        {
            ConvertFloatToHalfStream(&x, sizeof(HALF), pArray, sizeof(float), 4);
        }

        FORCE_INLINE SHORTN4::SHORTN4(float _x, float _y, float _z, float _w) noexcept
        {
            StoreShortN4(this, VEC::Set(_x, _y, _z, _w));
        }

        _Use_decl_annotations_
        FORCE_INLINE SHORTN4::SHORTN4(const float* pArray) noexcept
        {
            StoreShortN4(this, VEC::LoadFloat4(reinterpret_cast<const FLAOT4*>(pArray)));
        }

        FORCE_INLINE SHORT4::SHORT4(float _x, float _y, float _z, float _w) noexcept
        {
            StoreShort4(this, VEC::Set(_x, _y, _z, _w));
        }

        _Use_decl_annotations_
        FORCE_INLINE SHORT4::SHORT4(const float* pArray) noexcept
        {
            StoreShort4(this, VEC::LoadFloat4(reinterpret_cast<const FLOAT4*>(pArray)));
        }

        FORCE_INLINE USHORTN4::USHORTN4(float _x, float _y, float _z, float _w) noexcept
        {
            StoreUShortN4(this, VEC::Set(_x, _y, _z, _w));
        }

        _Use_decl_annotations_
        FORCE_INLINE USHORTN4::USHORTN4(const float* pArray) noexcept
        {
            StoreUShortN4(this, VEC::LoadFloat4(reinterpret_cast<const FLOAT4*>(pArray)));
        }

        FORCE_INLINE USHORT4::USHORT4(float _x, float _y, float _z, float _w) noexcept
        {
            StoreUShort4(this, VEC::Set(_x, _y, _z, _w));
        }

        _Use_decl_annotations_
        FORCE_INLINE USHORT4::USHORT4(const float* pArray) noexcept
        {
            StoreUShort4(this, VEC::LoadFloat4(reinterpret_cast<const FLOAT4*>(pArray)));
        }

        FORCE_INLINE XDECN4::XDECN4(float _x, float _y, float _z, float _w) noexcept
        {
            StoreXDecN4(this, VEC::Set(_x, _y, _z, _w));
        }

        _Use_decl_annotations_
        FORCE_INLINE XDECN4::XDECN4(const float* pArray) noexcept
        {
            StoreXDecN4(this, VEC::LoadFloat4(reinterpret_cast<const FLOAT4*>(pArray)));
        }

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4996)
 // C4996: ignore deprecation warning
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

        FORCE_INLINE XDEC4::XDEC4(float _x, float _y, float _z, float _w) noexcept
        {
            StoreXDec4(this, VEC::Set(_x, _y, _z, _w));
        }

        _Use_decl_annotations_
        FORCE_INLINE XDEC4::XDEC4(const float* pArray) noexcept
        {
            StoreXDec4(this, VEC::LoadFloat4(reinterpret_cast<const FLOAT4*>(pArray)));
        }

        FORCE_INLINE DECN4::DECN4(float _x, float _y, float _z, float _w) noexcept
        {
            StoreDecN4(this, VEC::Set(_x, _y, _z, _w));
        }

        _Use_decl_annotations_
        FORCE_INLINE DECN4::DECN4(const float* pArray) noexcept
        {
            StoreDecN4(this, VEC::LoadFloat4(reinterpret_cast<const FLOAT4*>(pArray)));
        }

        FORCE_INLINE DEC4::DEC4(float _x, float _y, float _z, float _w) noexcept
        {
            StoreDec4(this, VEC::Set(_x, _y, _z, _w));
        }

        _Use_decl_annotations_
        FORCE_INLINE DEC4::DEC4(const float* pArray) noexcept
        {
            StoreDec4(this, VEC:LoadFloat4(reinterpret_cast<const FLOAT4*>(pArray)));
        }

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif

        FORCE_INLINE UDECN4::UDECN4(float _x, float _y, float _z, float _w) noexcept
        {
            StoreUDecN4(this, VEC::Set(_x, _y, _z, _w));
        }

        _Use_decl_annotations_
        FORCE_INLINE UDECN4::UDECN4(const float* pArray) noexcept
        {
            StoreUDecN4(this, VEC::LoadFloat4(reinterpret_cast<const FLOAT4*>(pArray)));
        }

        FORCE_INLINE UDEC4::UDEC4(float _x, float _y, float _z, float _w) noexcept
        {
            StoreUDec4(this, VEC::Set(_x, _y, _z, _w));
        }

        _Use_decl_annotations_
        FORCE_INLINE UDEC4::UDEC4(const float* pArray) noexcept
        {
            StoreUDec4(this, VEC::LoadFloat4(reinterpret_cast<const FLOAT4*>(pArray)));
        }

        FORCE_INLINE BYTEN4::BYTEN4(float _x, float _y, float _z, float _w) noexcept
        {
            StoreByteN4(this, VEC::Set(_x, _y, _z, _w));
        }

        _Use_decl_annotations_
        FORCE_INLINE BYTEN4::BYTEN4(const float* pArray) noexcept
        {
            StoreByteN4(this, VEC::LoadFloat4(reinterpret_cast<const FLOAT4*>(pArray)));
        }

        FORCE_INLINE BYTE4::BYTE4(float _x, float _y, float _z, float _w) noexcept
        {
            StoreByte4(this, VEC::Set(_x, _y, _z, _w));
        }

        _Use_decl_annotations_
        FORCE_INLINE BYTE4::BYTE4(const float* pArray) noexcept
        {
            StoreByte4(this, VEC::LoadFloat4(reinterpret_cast<const FLOAT4*>(pArray)));
        }

        FORCE_INLINE UBYTEN4::UBYTEN4(float _x, float _y, float _z, float _w) noexcept
        {
            StoreUByteN4(this, VEC::Set(_x, _y, _z, _w));
        }

        _Use_decl_annotations_
        FORCE_INLINE UBYTEN4::UBYTEN4(const float* pArray) noexcept
        {
            StoreUByteN4(this, VEC::LoadFloat4(reinterpret_cast<const FLOAT4*>(pArray)));
        }

        FORCE_INLINE UBYTE4::UBYTE4(float _x, float _y, float _z, float _w) noexcept
        {
            StoreUByte4(this, VEC::Set(_x, _y, _z, _w));
        }

        _Use_decl_annotations_
        FORCE_INLINE UBYTE4::UBYTE4(const float* pArray) noexcept
        {
            StoreUByte4(this, VEC::LoadFloat4(reinterpret_cast<const FLOAT4*>(pArray)));
        }

        FORCE_INLINE UNIBBLE4::UNIBBLE4(float _x, float _y, float _z, float _w) noexcept
        {
            StoreUNibble4(this, VEC::Set(_x, _y, _z, _w));
        }

        _Use_decl_annotations_
        FORCE_INLINE UNIBBLE4:UNIBBLE4(const float* pArray) noexcept
        {
            StoreUNibble4(this, VEC::LoadFloat4(reinterpret_cast<const FLOAT4*>(pArray)));
        }

        FORCE_INLINE U555:U555(float _x, float _y, float _z, bool _w) noexcept
        {
            StoreU555(this, VEC::Set(_x, _y, _z, ((_w) ? 1.0f : 0.0f)));
        }

        _Use_decl_annotations_
        FORCE_INLINE U555:U555(const float* pArray, bool _w) noexcept
        {
            VECTOR v = VEC::LoadFloat3(reinterpret_cast<const FLOAT3*>(pArray));
            StoreU555(this, VEC::SetW(v, ((_w) ? 1.0f : 0.0f)));
        }
    } // namespace PackedVector
} // namespace UltReality::Math

#endif // !ULTREALITY_MATH_PACKED_VECTOR_INL
