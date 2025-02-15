#ifndef ULTREALITY_MATH_UTILITY_INL
#define ULTREALITY_MATH_UTILITY_INL

#if defined(__GNUC__) or defined(__clang__)
#define FORCE_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline
#endif

#include <cmath>
#include <assert.h>

#include <Constants.h>

namespace UltReality::Math
{
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

    // High-accuracy tangent polynomial approximation (11-degree minimax)
    // Approximates tan(x) for x in [-pi/4, pi/4]
    static inline float tanPolynomial(float x) noexcept
    {
        // Precompute powers (tan is an odd function)
        float x2 = x * x;
        float x3 = x * x2;
        float x5 = x3 * x2;
        float x7 = x5 * x2;
        float x9 = x7 * x2;
        float x11 = x9 * x2;
        // Using coefficients similar to the Taylor expansion:
        // tan(x) ~ x + x^3/3 + (2x^5)/15 + (17x^7)/315 + (62x^9)/2835 + (1382x^11)/155925
        return x
            + 0.3333333333f * x3
            + 0.1333333333f * x5
            + 0.053968254f * x7
            + 0.021854f * x9
            + 0.008863f * x11;
    }

    // Estimated tangent polynomial approximation (7-degree minimax)
    static inline float tanPolynomialEst(float x) noexcept
    {
        float x2 = x * x;
        float x3 = x * x2;
        float x5 = x3 * x2;
        float x7 = x5 * x2;
        return x
            + 0.3333333333f * x3
            + 0.1333333333f * x5
            + 0.053968254f * x7;
    }

    // Helper: Reduce an angle modulo pi so that the result is in [-pi/2, pi/2]
    // (tan has period pi)
    static inline float reduceTanAngle(float x) noexcept
    {
        // Round to nearest integer k, then subtract k*pi.
        int k = static_cast<int>(std::round(x / _PI));
        return x - k * _PI;
    }

    // High-accuracy ScalarTan implementation
    FORCE_INLINE float ScalarTan(float x) noexcept
    {
        // Step 1: Reduce x to [-pi/2, pi/2]
        float xr = reduceTanAngle(x);
        // Step 2: Further restrict to [-pi/4, pi/4] for the polynomial.
        // For |xr| > pi/4, use the identity:
        //    tan(x) = sign(x) * (1 / tan(pi/2 - |x|))
        if (std::fabs(xr) > (_PIOVER2 * 0.5f)) // i.e. |xr| > pi/4
        {
            float sign = (xr >= 0.0f) ? 1.0f : -1.0f;
            float y = _PIOVER2 - std::fabs(xr);
            // Avoid division by zero near the asymptote.
            if (std::fabs(y) < 1e-6f)
            {
                return sign * 1e30f; // a large number to indicate asymptotic behavior
            }
            return sign / tanPolynomial(y);
        }
        else
        {
            return tanPolynomial(xr);
        }
    }

    // Estimated version of ScalarTan (using the 7-degree polynomial)
    FORCE_INLINE float ScalarTanEst(float x) noexcept
    {
        float xr = reduceTanAngle(x);
        if (std::fabs(xr) > (_PIOVER2 * 0.5f))
        {
            float sign = (xr >= 0.0f) ? 1.0f : -1.0f;
            float y = _PIOVER2 - std::fabs(xr);
            if (std::fabs(y) < 1e-6f)
            {
                return sign * 1e30f;
            }
            return sign / tanPolynomialEst(y);
        }
        else
        {
            return tanPolynomialEst(xr);
        }
    }

    //--------------------------------------------------
    // ATAN IMPLEMENTATION
    //--------------------------------------------------

    // High-accuracy arctan polynomial approximation (11-degree minimax)
    // Approximates atan(x) for x in [0, 1]
    static inline float atanPolynomial(float x) noexcept
    {
        float x2 = x * x;
        float x3 = x * x2;
        float x5 = x3 * x2;
        float x7 = x5 * x2;
        float x9 = x7 * x2;
        float x11 = x9 * x2;
        // Using coefficients from the Taylor series for arctan:
        // atan(x) ~ x - x^3/3 + x^5/5 - x^7/7 + x^9/9 - x^11/11
        return x
            - (x3) / 3.0f
            + (x5) / 5.0f
            - (x7) / 7.0f
            + (x9) / 9.0f
            - (x11) / 11.0f;
    }

    // Estimated arctan polynomial approximation (7-degree minimax)
    static inline float atanPolynomialEst(float x) noexcept
    {
        float x2 = x * x;
        float x3 = x * x2;
        float x5 = x3 * x2;
        float x7 = x5 * x2;
        return x
            - (x3) / 3.0f
            + (x5) / 5.0f
            - (x7) / 7.0f;
    }

    // High-accuracy ScalarATan implementation
    FORCE_INLINE float ScalarATan(float x) noexcept
    {
        // arctan is an odd function
        if (x < 0.0f)
            return -ScalarATan(-x);
        // For x > 1, use the identity:
        //    atan(x) = pi/2 - atan(1/x)
        if (x > 1.0f)
            return _PIOVER2 - ScalarATan(1.0f / x);
        // Now x is in [0, 1]
        return atanPolynomial(x);
    }

    // Estimated version of ScalarATan (using the 7-degree polynomial)
    FORCE_INLINE float ScalarATanEst(float x) noexcept
    {
        if (x < 0.0f)
            return -ScalarATanEst(-x);
        if (x > 1.0f)
            return _PIOVER2 - ScalarATanEst(1.0f / x);
        return atanPolynomialEst(x);
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
        *pSine = (((((-2.3889859e-08f * y2 + 2.7525562e-06f) * y2 - 0.00019840874f) * y2 + 0.0083333310f) * y2 - 0.16666667f) * y2 + 1.0f) * y;

        // 10-degree minimax approximation
        float p = ((((-2.6051615e-07f * y2 + 2.4760495e-05f) * y2 - 0.0013888378f) * y2 + 0.041666638f) * y2 - 0.5f) * y2 + 1.0f;
        *pCos = sign * p;
    }

    _Use_decl_annotations_
    FORCE_INLINE void ScalarSineCosEst(float* pSine, float* pCos, float value) noexcept
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
        *pSine = (((-0.00018524670f * y2 + 0.0083139502f) * y2 - 0.16665852f) * y2 + 1.0f) * y;

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
}

#endif // !ULTREALITY_MATH_UTILITY_INL
