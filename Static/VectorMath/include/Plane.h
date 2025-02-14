#ifndef ULTREALITY_MATH_PLANE_H
#define ULTREALITY_MATH_PLANE_H

#include <SIMDVectorConfig.h>

namespace UltReality::Math
{
    struct Float3;
    struct Float4;

    namespace Plane
    {
        bool VEC_CALLCONV Equal(A_VECTOR P1, A_VECTOR P2) noexcept;
        bool VEC_CALLCONV NearEqual(A_VECTOR P1, A_VECTOR P2, A_VECTOR epsilon) noexcept;
        bool VEC_CALLCONV NotEqual(A_VECTOR P1, A_VECTOR P2) noexcept;

        bool VEC_CALLCONV IsNaN(A_VECTOR p) noexcept;
        bool VEC_CALLCONV IsInfinite(A_VECTOR p) noexcept;

        VECTOR VEC_CALLCONV Dot(A_VECTOR p, A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV DotCoord(A_VECTOR p, A_VECTOR v) noexcept;
        VECTOR VEC_CALLCONV DotNormal(A_VECTOR p, A_VECTOR v) noexcept;
        // Uses a reciprocal estimate and return QNaN on zero and infinite vectors
        VECTOR VEC_CALLCONV NormalizeEst(A_VECTOR p) noexcept;
        VECTOR VEC_CALLCONV Normalize(A_VECTOR p) noexcept;
        VECTOR VEC_CALLCONV IntersectLine(A_VECTOR p, A_VECTOR linePoint1, A_VECTOR linePoint2) noexcept;
        void VEC_CALLCONV IntersectPlane(_Out_ VECTOR* pLinePoint1, _Out_ VECTOR* pLinePoint2, _In_ A_VECTOR P1, _In_ A_VECTOR P2) noexcept;

        // Transforms a plane given an inverse transpose matrix
        VECTOR VEC_CALLCONV Transform(A_VECTOR p, A_MATRIX ITM) noexcept;

        // Transforms an array of planes given an inverse transpose matrix
        Float4* VEC_CALLCONV TransformStream(_Out_writes_bytes_(sizeof(Float4) + outputStride * (planeCount - 1)) Float4* pOutputStream, 
            _In_ size_t outputStride, 
            _In_reads_bytes_(sizeof(Float4) + inputStride * (planeCount - 1)) const Float4* pInputStream, 
            _In_ size_t inputStride, _In_ size_t planeCount, _In_ A_MATRIX ITM) noexcept;

        VECTOR VEC_CALLCONV FromPointNormal(A_VECTOR point, A_VECTOR normal) noexcept;
        VECTOR VEC_CALLCONV FromPoints(A_VECTOR point1, A_VECTOR point2, A_VECTOR point3) noexcept;
    } // namespace Plane
    
} // namespace UltReality::Math


#endif // !ULTREALITY_MATH_PLANE_H