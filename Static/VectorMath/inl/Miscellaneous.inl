#ifndef ULTREALITY_MATH_MISCELLANEOUS_INL
#define ULTREALITY_MATH_MISCELLANEOUS_INL

#if defined(__GNUC__) or defined(__clang__)
#define FORCE_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline
#endif

namespace UltReality::Math
{
    namespace Miscellaneous
    {
        FORCE_INLINE bool VerifyCPUSupport() noexcept
        {
#if defined(_SSE2_INTRINSICS_) && !defined(_NO_INTRINSICS_)
            int CPUInfo[4] = { -1 };
#if (defined(__clang__) || defined(__GNUC__)) && defined(__cpuid)
            __cpuid(0, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
#else
            __cpuid(CPUInfo, 0);
#endif

#ifdef __AVX2__
            if (CPUInfo[0] < 7)
                return false;
#else
            if (CPUInfo[0] < 1)
                return false;
#endif

#if (defined(__clang__) || defined(__GNUC__)) && defined(__cpuid)
            __cpuid(1, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
#else
            __cpuid(CPUInfo, 1);
#endif

#if defined(__AVX2__) || defined(_AVX2_INTRINSICS_)
            // The compiler can emit FMA3 instructions even without explicit intrinsics use
            if ((CPUInfo[2] & 0x38081001) != 0x38081001)
                return false; // No F16C/AVX/OSXSAVE/SSE4.1/FMA3/SSE3 support
#elif defined(_FMA3_INTRINSICS_) && defined(_F16C_INTRINSICS_)
            if ((CPUInfo[2] & 0x38081001) != 0x38081001)
                return false; // No F16C/AVX/OSXSAVE/SSE4.1/FMA3/SSE3 support
#elif defined(_FMA3_INTRINSICS_)
            if ((CPUInfo[2] & 0x18081001) != 0x18081001)
                return false; // No AVX/OSXSAVE/SSE4.1/FMA3/SSE3 support
#elif defined(_F16C_INTRINSICS_)
            if ((CPUInfo[2] & 0x38080001) != 0x38080001)
                return false; // No F16C/AVX/OSXSAVE/SSE4.1/SSE3 support
#elif defined(__AVX__) || defined(_AVX_INTRINSICS_)
            if ((CPUInfo[2] & 0x18080001) != 0x18080001)
                return false; // No AVX/OSXSAVE/SSE4.1/SSE3 support
#elif defined(_SSE4_INTRINSICS_)
            if ((CPUInfo[2] & 0x80001) != 0x80001)
                return false; // No SSE3/SSE4.1 support
#elif defined(_SSE3_INTRINSICS_)
            if (!(CPUInfo[2] & 0x1))
                return false; // No SSE3 support
#endif

            // The x64 processor model requires SSE2 support, but no harm in checking
            if ((CPUInfo[3] & 0x6000000) != 0x6000000)
                return false; // No SSE2/SSE support

#if defined(__AVX2__) || defined(_AVX2_INTRINSICS_)
#if defined(__clang__) || defined(__GNUC__)
            __cpuid_count(7, 0, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
#else
            __cpuidex(CPUInfo, 7, 0);
#endif
            if (!(CPUInfo[1] & 0x20))
                return false; // No AVX2 support
#endif

            return true;
#endif
        }

        FORCE_INLINE VECTOR VEC_CALLCONV FresnelTerm(A_VECTOR cosIncidentAngle, A_VECTOR refractionIndex) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(!Vector4::IsInfinite(cosIncidentAngle));
#endif

            // Result = 0.5f * (g - c)^2 / (g + c)^2 * ((c * (g + c) - 1)^2 / (c * (g - c) + 1)^2 + 1) where
            // c = cosIncidentAngle
            // g = sqrt(c^2 + refractionIndex^2 - 1)

#if defined(_NO_INTRINSICS_)
            VECTOR G = Vector::MultiplyAdd(refractionIndex, refractionIndex, g_NegativeOne.v);
            G = Vector::MultiplyAdd(cosIncidentAngle, cosIncidentAngle, G);
            G = Vector::Abs(G);
            G = Vector::Sqrt(G);

            VECTOR S = Vector::Add(G, cosIncidentAngle);
            VECTOR D = Vector::Subtract(G, cosIncidentAngle);

            VECTOR V0 = Vector::Multiply(D, D);
            VECTOR V1 = Vector::Multiply(S, S);
            V1 = Vector::Reciprocal(V1);
            V0 = Vector::Multiply(g_OneHalf.v, V0);
            V0 = Vector::Multiply(V0, V1);

            VECTOR V2 = Vector::MultiplyAdd(cosIncidentAngle, S, g_NegativeOne.v);
            VECTOR V3 = Vector::MultiplyAdd(cosIncidentAngle, D, g_One.v);
            V2 = Vector::Multiply(V2, V2);
            V3 = Vector::Multiply(V3, V3);
            V3 = Vector::Reciprocal(V3);
            V2 = Vector::MultiplyAdd(V2, V3, g_One.v);

            VECTOR Result = Vector::Multiply(V0, V2);

            return Vector::Saturate(Result);

#elif defined(_SSE2_INTRINSICS_)
            // G = sqrt(abs((refractionIndex^2-1) + cosIncidentAngle^2))
            VECTOR G = _mm_mul_ps(refractionIndex, refractionIndex);
            VECTOR vTemp = _mm_mul_ps(cosIncidentAngle, cosIncidentAngle);
            G = _mm_sub_ps(G, g_One);
            vTemp = _mm_add_ps(vTemp, G);
            // max((0-vTemp),vTemp) == abs(vTemp)
            // The abs is needed to deal with refraction and cosine being zero
            G = _mm_setzero_ps();
            G = _mm_sub_ps(G, vTemp);
            G = _mm_max_ps(G, vTemp);
            // Last operation, the sqrt()
            G = _mm_sqrt_ps(G);

            // Calc G-C and G+C
            VECTOR GAddC = _mm_add_ps(G, cosIncidentAngle);
            VECTOR GSubC = _mm_sub_ps(G, cosIncidentAngle);
            // Perform the term (0.5f *(g - c)^2) / (g + c)^2
            VECTOR vResult = _mm_mul_ps(GSubC, GSubC);
            vTemp = _mm_mul_ps(GAddC, GAddC);
            vResult = _mm_mul_ps(vResult, g_OneHalf);
            vResult = _mm_div_ps(vResult, vTemp);
            // Perform the term ((c * (g + c) - 1)^2 / (c * (g - c) + 1)^2 + 1)
            GAddC = _mm_mul_ps(GAddC, cosIncidentAngle);
            GSubC = _mm_mul_ps(GSubC, cosIncidentAngle);
            GAddC = _mm_sub_ps(GAddC, g_One);
            GSubC = _mm_add_ps(GSubC, g_One);
            GAddC = _mm_mul_ps(GAddC, GAddC);
            GSubC = _mm_mul_ps(GSubC, GSubC);
            GAddC = _mm_div_ps(GAddC, GSubC);
            GAddC = _mm_add_ps(GAddC, g_One);
            // Multiply the two term parts
            vResult = _mm_mul_ps(vResult, GAddC);
            // Clamp to 0.0 - 1.0f
            vResult = _mm_max_ps(vResult, g_Zero);
            
            return _mm_min_ps(vResult, g_One);
#endif
        }
    } // namespace Miscellaneous
} // namespace UltReality::Math

#endif // !ULTREALITY_MATH_MISCELLANEOUS_INL