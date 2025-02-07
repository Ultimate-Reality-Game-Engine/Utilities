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
            assert(!VEC4::IsInfinite(cosIncidentAngle));
#endif

            // Result = 0.5f * (g - c)^2 / (g + c)^2 * ((c * (g + c) - 1)^2 / (c * (g - c) + 1)^2 + 1) where
            // c = cosIncidentAngle
            // g = sqrt(c^2 + refractionIndex^2 - 1)

#if defined(_NO_INTRINSICS_)
            VECTOR G = VEC::MultiplyAdd(refractionIndex, refractionIndex, g_NegativeOne.v);
            G = VEC::MultiplyAdd(cosIncidentAngle, cosIncidentAngle, G);
            G = VEC::Abs(G);
            G = VEC::Sqrt(G);

            VECTOR S = VEC::Add(G, cosIncidentAngle);
            VECTOR D = VEC::Subtract(G, cosIncidentAngle);

            VECTOR V0 = VEC::Multiply(D, D);
            VECTOR V1 = VEC::Multiply(S, S);
            V1 = VEC::Reciprocal(V1);
            V0 = VEC::Multiply(g_OneHalf.v, V0);
            V0 = VEC::Multiply(V0, V1);

            VECTOR V2 = VEC::MultiplyAdd(cosIncidentAngle, S, g_NegativeOne.v);
            VECTOR V3 = VEC::MultiplyAdd(cosIncidentAngle, D, g_One.v);
            V2 = VEC::Multiply(V2, V2);
            V3 = VEC::Multiply(V3, V3);
            V3 = VEC::Reciprocal(V3);
            V2 = VEC::MultiplyAdd(V2, V3, g_One.v);

            VECTOR Result = VEC::Multiply(V0, V2);

            return VEC::Saturate(Result);

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

        FORCE_INLINE bool ScalarNearEqual(float S1, float S2, float epsilon) noexcept
        {
            float delta = S1 - S2;

            return (fabsf(delta) <= epsilon);
        }

        FORCE_INLINE float ScalarModAngle(float angle) noexcept
        {
            // Note: The modulo is performed with unsigned math only to work
            // around a precision error on numbers that are close to PI

            // Normalize the range from 0.0f to _2PI
            angle = angle + _PI;
            // Perform the modulo, unsigned
            float fTemp = fabsf(angle);
            fTemp = fTemp - (_2PI * static_cast<float>(static_cast<int32_t>(fTemp / _2PI)));
            // Restore the number to the range of -_PI to _PI-epsilon
            fTemp = fTemp - _PI;
            // If the modulo'd value was negative, restore negation
            if (angle < 0.0f)
            {
                fTemp = -fTemp;
            }
            
            return fTemp;
        }

        FORCE_INLINE float ScalarSine(float value) noexcept
        {
            // Map Value to y in [-pi,pi], x = 2*pi*quotient + remainder.
            float quotient = _1OVER2PI * value;
            if (value >= 0.0f)
            {
                quotient = static_cast<float>(static_cast<int>(quotient + 0.5f));
            }
            else
            {
                quotient = static_cast<float>(static_cast<int>(quotient - 0.5f));
            }
            float y = value - _2PI * quotient;

            // Map y to [-pi/2,pi/2] with sin(y) = sin(Value).
            if (y > _PIOVER2)
            {
                y = _PI - y;
            }
            else if (y < -_PIOVER2)
            {
                y = -_PI - y;
            }

            // 11-degree minimax approximation
            float y2 = y * y;
            return (((((-2.3889859e-08f * y2 + 2.7525562e-06f) * y2 - 0.00019840874f) * y2 + 0.0083333310f) * y2 - 0.16666667f) * y2 + 1.0f) * y;
        }

        FORCE_INLINE float ScalarSineEst(float value) noexcept
        {
            // Map Value to y in [-pi,pi], x = 2*pi*quotient + remainder.
            float quotient = _1OVER2PI * value;
            if (value >= 0.0f)
            {
                quotient = static_cast<float>(static_cast<int>(quotient + 0.5f));
            }
            else
            {
                quotient = static_cast<float>(static_cast<int>(quotient - 0.5f));
            }
            float y = value - _2PI * quotient;

            // Map y to [-pi/2,pi/2] with sin(y) = sin(Value).
            if (y > _PIOVER2)
            {
                y = _PI - y;
            }
            else if (y < -_PIOVER2)
            {
                y = -_PI - y;
            }

            // 7-degree minimax approximation
            float y2 = y * y;
            return (((-0.00018524670f * y2 + 0.0083139502f) * y2 - 0.16665852f) * y2 + 1.0f) * y;
        }

        FORCE_INLINE float ScalarCos(float value) noexcept
        {
            // Map Value to y in [-pi,pi], x = 2*pi*quotient + remainder.
            float quotient = _1OVER2PI * value;
            if (value >= 0.0f)
            {
                quotient = static_cast<float>(static_cast<int>(quotient + 0.5f));
            }
            else
            {
                quotient = static_cast<float>(static_cast<int>(quotient - 0.5f));
            }
            float y = value - _2PI * quotient;

            // Map y to [-pi/2,pi/2] with cos(y) = sign*cos(x).
            float sign;
            if (y > _PIOVER2)
            {
                y = _PI - y;
                sign = -1.0f;
            }
            else if (y < -_PIOVER2)
            {
                y = -_PI - y;
                sign = -1.0f;
            }
            else
            {
                sign = +1.0f;
            }

            // 10-degree minimax approximation
            float y2 = y * y;
            float p = ((((-2.6051615e-07f * y2 + 2.4760495e-05f) * y2 - 0.0013888378f) * y2 + 0.041666638f) * y2 - 0.5f) * y2 + 1.0f;
            return sign * p;
        }

        FORCE_INLINE float ScalarCosEst(float value) noexcept
        {
            // Map Value to y in [-pi,pi], x = 2*pi*quotient + remainder.
            float quotient = _1OVER2PI * value;
            if (value >= 0.0f)
            {
                quotient = static_cast<float>(static_cast<int>(quotient + 0.5f));
            }
            else
            {
                quotient = static_cast<float>(static_cast<int>(quotient - 0.5f));
            }
            float y = value - _2PI * quotient;

            // Map y to [-pi/2,pi/2] with cos(y) = sign*cos(x).
            float sign;
            if (y > _PIOVER2)
            {
                y = _PI - y;
                sign = -1.0f;
            }
            else if (y < -_PIOVER2)
            {
                y = -_PI - y;
                sign = -1.0f;
            }
            else
            {
                sign = +1.0f;
            }

            // 6-degree minimax approximation
            float y2 = y * y;
            float p = ((-0.0012712436f * y2 + 0.041493919f) * y2 - 0.49992746f) * y2 + 1.0f;
            return sign * p;
        }

        _Use_decl_annotations_
        FORCE_INLINE void ScalarSineCos(float* pSine, float* pCos, float value) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSine != nullptr);
            assert(pCos != nullptr);
#endif

            // Map Value to y in [-pi,pi], x = 2*pi*quotient + remainder.
            float quotient = _1OVER2PI * value;
            if (value >= 0.0f)
            {
                quotient = static_cast<float>(static_cast<int>(quotient + 0.5f));
            }
            else
            {
                quotient = static_cast<float>(static_cast<int>(quotient - 0.5f));
            }
            float y = value - _2PI * quotient;

            // Map y to [-pi/2,pi/2] with sin(y) = sin(Value).
            float sign;
            if (y > _PIOVER2)
            {
                y = _PI - y;
                sign = -1.0f;
            }
            else if (y < -_PIOVER2)
            {
                y = -_PI - y;
                sign = -1.0f;
            }
            else
            {
                sign = +1.0f;
            }

            float y2 = y * y;

            // 11-degree minimax approximation
            *pSin = (((((-2.3889859e-08f * y2 + 2.7525562e-06f) * y2 - 0.00019840874f) * y2 + 0.0083333310f) * y2 - 0.16666667f) * y2 + 1.0f) * y;

            // 10-degree minimax approximation
            float p = ((((-2.6051615e-07f * y2 + 2.4760495e-05f) * y2 - 0.0013888378f) * y2 + 0.041666638f) * y2 - 0.5f) * y2 + 1.0f;
            *pCos = sign * p;
        }

        _Use_decl_annotations_
        FORCE_INLINE void ScalarSinCosEst(float* pSine, float* pCos, float value) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pSine != nullptr);
            assert(pCos != nullptr);
#endif

            // Map Value to y in [-pi,pi], x = 2*pi*quotient + remainder.
            float quotient = _1OVER2PI * value;
            if (value >= 0.0f)
            {
                quotient = static_cast<float>(static_cast<int>(quotient + 0.5f));
            }
            else
            {
                quotient = static_cast<float>(static_cast<int>(quotient - 0.5f));
            }
            float y = value - _2PI * quotient;

            // Map y to [-pi/2,pi/2] with sin(y) = sin(Value).
            float sign;
            if (y > _PIOVER2)
            {
                y = _PI - y;
                sign = -1.0f;
            }
            else if (y < -_PIOVER2)
            {
                y = -_PI - y;
                sign = -1.0f;
            }
            else
            {
                sign = +1.0f;
            }

            float y2 = y * y;

            // 7-degree minimax approximation
            *pSin = (((-0.00018524670f * y2 + 0.0083139502f) * y2 - 0.16665852f) * y2 + 1.0f) * y;

            // 6-degree minimax approximation
            float p = ((-0.0012712436f * y2 + 0.041493919f) * y2 - 0.49992746f) * y2 + 1.0f;
            *pCos = sign * p;
        }

        FORCE_INLINE float ScalarASine(float value) noexcept
        {
            // Clamp input to [-1,1].
            bool nonnegative = (value >= 0.0f);
            float x = fabsf(value);
            float omx = 1.0f - x;
            if (omx < 0.0f)
            {
                omx = 0.0f;
            }
            float root = sqrtf(omx);

            // 7-degree minimax approximation
            float result = ((((((-0.0012624911f * x + 0.0066700901f) * x - 0.0170881256f) * x + 0.0308918810f) * x - 0.0501743046f) * x + 0.0889789874f) * x - 0.2145988016f) * x + 1.5707963050f;
            result *= root;  // acos(|x|)

            // acos(x) = pi - acos(-x) when x < 0, asin(x) = pi/2 - acos(x)
            return (nonnegative ? _PIOVER2 - result : result - _PIOVER2);
        }

        FORCE_INLINE float ScalarASineEst(float value) noexcept
        {
            // Clamp input to [-1,1].
            bool nonnegative = (value >= 0.0f);
            float x = fabsf(value);
            float omx = 1.0f - x;
            if (omx < 0.0f)
            {
                omx = 0.0f;
            }
            float root = sqrtf(omx);

            // 3-degree minimax approximation
            float result = ((-0.0187293f * x + 0.0742610f) * x - 0.2121144f) * x + 1.5707288f;
            result *= root;  // acos(|x|)

            // acos(x) = pi - acos(-x) when x < 0, asin(x) = pi/2 - acos(x)
            return (nonnegative ? _PIOVER2 - result : result - _PIOVER2);
        }

        FORCE_INLINE float ScalarACos(float value) noexcept
        {
            // Clamp input to [-1,1].
            bool nonnegative = (value >= 0.0f);
            float x = fabsf(value);
            float omx = 1.0f - x;
            if (omx < 0.0f)
            {
                omx = 0.0f;
            }
            float root = sqrtf(omx);

            // 7-degree minimax approximation
            float result = ((((((-0.0012624911f * x + 0.0066700901f) * x - 0.0170881256f) * x + 0.0308918810f) * x - 0.0501743046f) * x + 0.0889789874f) * x - 0.2145988016f) * x + 1.5707963050f;
            result *= root;

            // acos(x) = pi - acos(-x) when x < 0
            return (nonnegative ? result : _PI - result);
        }

        FORCE_INLINE float ScalarACosEst(float value) noexcept
        {
            // Clamp input to [-1,1].
            bool nonnegative = (value >= 0.0f);
            float x = fabsf(value);
            float omx = 1.0f - x;
            if (omx < 0.0f)
            {
                omx = 0.0f;
            }
            float root = sqrtf(omx);

            // 3-degree minimax approximation
            float result = ((-0.0187293f * x + 0.0742610f) * x - 0.2121144f) * x + 1.5707288f;
            result *= root;

            // acos(x) = pi - acos(-x) when x < 0
            return (nonnegative ? result : _PI - result);
        }
    } // namespace Miscellaneous
    
} // namespace UltReality::Math

#endif // !ULTREALITY_MATH_MISCELLANEOUS_INL