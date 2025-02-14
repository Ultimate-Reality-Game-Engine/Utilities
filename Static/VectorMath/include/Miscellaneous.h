#ifndef ULTREALITY_MATH_MISCELLANEOUS_H
#define ULTREALITY_MATH_MISCELLANEOUS_H

#include <SIMDVectorConfig.h>
#include <VectorMathConstants.h>
#include <Float4.h>

namespace UltReality::Math
{
    namespace Miscellaneous
    {
        bool VerifyCPUSupport() noexcept;

        VECTOR VEC_CALLCONV FresnelTerm(A_VECTOR cosIncidentAngle, A_VECTOR refractionIndex) noexcept;
    } // namespace Miscellaneous
    
} // namespace UltReality::Math

#endif // !ULTREALITY_MATH_MISCELLANEOUS_H