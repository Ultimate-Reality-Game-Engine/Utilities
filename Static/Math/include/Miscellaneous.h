#ifndef ULTREALITY_MATH_MISCELLANEOUS_H
#define ULTREALITY_MATH_MISCELLANEOUS_H

#include <SSE2VectorConfig.h>
#include <VectorMathConstants.h>
#include <Float4.h>

namespace UltReality::Math
{
    namespace Miscellaneous
    {
        bool VerifyCPUSupport() noexcept;

        VECTOR VEC_CALLCONV FresnelTerm(A_VECTOR cosIncidentAngle, A_VECTOR refractionIndex) noexcept;

        bool ScalarNearEqual(float S1, float S2, float epsilon) noexcept;
        // Modulo the range of the given angle such that -PI <= angle <= PI
        float ScalarModAngle(float value) noexcept;

        float ScalarSine(float value) noexcept;
        float ScalarSineEst(float value) noexcept;

        float ScalarCos(float value) noexcept;
        float ScalarCosEst(float value) noexcept;

        void ScalarSineCos(_Out_ float* pSine, _Out_ float* pCos, float value) noexcept;
        void ScalarSineCosEst(_Out_ float* pSine, _Out_ float* pCos, float value) noexcept;

        float ScalarASine(float value) noexcept;
        float ScalarASineEst(float value) noexcept;

        float ScalarACos(float value) noexcept;
        float ScalarACosEst(float value) noexcept;
    } // namespace Miscellaneous
    
} // namespace UltReality::Math

#endif // !ULTREALITY_MATH_MISCELLANEOUS_H