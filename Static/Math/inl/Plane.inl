#ifndef ULTREALITY_MATH_PLANE_INL
#define ULTREALITY_MATH_PLANE_INL

#if defined(__GNUC__) or defined(__clang__)
#define FORCE_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline
#endif

namespace UltReality::Math
{
    namespace Plane
    {
        FORCE_INLINE bool VEC_CALLCONV Equal(A_VECTOR P1, A_VECTOR P2) noexcept
        {
            return VEC4::Equal(P1, P2);
        }

        FORCE_INLINE bool VEC_CALLCONV NearEqual(A_VECTOR P1, A_VECTOR P2, A_VECTOR epsilon) noexcept
        {
            VECTOR NP1 = Normalize(P1);
            VECTOR NP2 = Normalize(P2);

            return VEC4::NearEqual(NP1, NP2, epsilon);
        }

        FORCE_INLINE bool VEC_CALLCONV NotEqual(A_VECTOR P1, A_VECTOR P2) noexcept
        {
            return VEC4::NotEqual(P1, P2);
        }

        FORCE_INLINE bool VEC_CALLCONV IsNaN(A_VECTOR p) noexcept
        {
            return VEC4::IsNaN(p);
        }

        FORCE_INLINE bool VEC_CALLCONV IsInfinite(A_VECTOR p) noexcept
        {
            return VEC4::IsInfinite(p);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV Dot(A_VECTOR p, A_VECTOR v) noexcept
        {
            return VEC4::Dot(p, v);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV DotCoord(A_VECTOR p, A_VECTOR v) noexcept
        {
            // Result = p[0] + v[0] + p[1] + v[1] + p[2] + v[2] + p[3]
            VECTOR V3 = VEC::Select(g_One.v, v, g_Select1110.v);
            
            return VEC4::Dot(p, V3);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV DotNormal(A_VECTOR p, A_VECTOR v) noexcept
        {
            return VEC3::Dot(p, v);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV NormalizeEst(A_VECTOR p) noexcept
        {
#if defined(_NO_INTRINSICS_)
            VECTOR Result = VEC3::ReciprocalLengthEst(p);
            
            return VEC::Multiply(p, Result);

#elif defined(_SSE4_INTRINSICS_)
            VECTOR vTemp = _mm_dp_ps(p, p, 0x7f);
            VECTOR vResult = _mm_rsqrt_ps(vTemp);
            
            return _mm_mul_ps(vResult, p);

#elif defined(_SSE2_INTRINSICS_)
            // Perform the dot product
            VECTOR vDot = _mm_mul_ps(p, p);
            // x=Dot.y, y=Dot.z
            VECTOR vTemp = PERMUTE_PS(vDot, _MM_SHUFFLE(2, 1, 2, 1));
            // Result.x = x+y
            vDot = _mm_add_ss(vDot, vTemp);
            // x=Dot.z
            vTemp = PERMUTE_PS(vTemp, _MM_SHUFFLE(1, 1, 1, 1));
            // Result.x = (x+y)+z
            vDot = _mm_add_ss(vDot, vTemp);
            // Splat x
            vDot = PERMUTE_PS(vDot, _MM_SHUFFLE(0, 0, 0, 0));
            // Get the reciprocal
            vDot = _mm_rsqrt_ps(vDot);
            
            // Get the reciprocal
            return _mm_mul_ps(vDot, P);
#endif
        }

        FORCE_INLINE VECTOR VEC_CALLCONV Normalize(A_VECTOR p) noexcept
        {
#if defined(_NO_INTRINSICS_)
            float fLengthSq = sqrtf((p.vector4_f32[0] * p.vector4_f32[0]) + (p.vector4_f32[1] * p.vector4_f32[1]) + (p.vector4_f32[2] * p.vector4_f32[2]));
            // Prevent divide by zero
            if (fLengthSq > 0)
            {
                fLengthSq = 1.0f / fLengthSq;
            }
            VECTOR_F32 vResult = { { {
                    p.vector4_f32[0] * fLengthSq,
                    p.vector4_f32[1] * fLengthSq,
                    p.vector4_f32[2] * fLengthSq,
                    p.vector4_f32[3] * fLengthSq
                } } };
            
            return vResult.v;

#elif defined(_SSE4_INTRINSICS_)
            VECTOR vLengthSq = _mm_dp_ps(p, p, 0x7f);
            // Prepare for the division
            VECTOR vResult = _mm_sqrt_ps(vLengthSq);
            // Failsafe on zero (Or epsilon) length planes
            // If the length is infinity, set the elements to zero
            vLengthSq = _mm_cmpneq_ps(vLengthSq, g_Infinity);
            // Reciprocal mul to perform the normalization
            vResult = _mm_div_ps(p, vResult);
            
            // Any that are infinity, set to zero
            return _mm_and_ps(vResult, vLengthSq);

#elif defined(_SSE2_INTRINSICS_)
            // Perform the dot product on x,y and z only
            VECTOR vLengthSq = _mm_mul_ps(p, p);
            VECTOR vTemp = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(2, 1, 2, 1));
            vLengthSq = _mm_add_ss(vLengthSq, vTemp);
            vTemp = PERMUTE_PS(vTemp, _MM_SHUFFLE(1, 1, 1, 1));
            vLengthSq = _mm_add_ss(vLengthSq, vTemp);
            vLengthSq = PERMUTE_PS(vLengthSq, _MM_SHUFFLE(0, 0, 0, 0));
            // Prepare for the division
            VECTOR vResult = _mm_sqrt_ps(vLengthSq);
            // Failsafe on zero (Or epsilon) length planes
            // If the length is infinity, set the elements to zero
            vLengthSq = _mm_cmpneq_ps(vLengthSq, g_Infinity);
            // Reciprocal mul to perform the normalization
            vResult = _mm_div_ps(p, vResult);
            
            // Any that are infinity, set to zero
            return _mm_and_ps(vResult, vLengthSq);
#endif
        }

        FORCE_INLINE VECTOR VEC_CALLCONV IntersectLine(A_VECTOR p, A_VECTOR linePoint1, A_VECTOR linePoint2) noexcept
        {
            VECTOR V1 = VEC3::Dot(p, linePoint1);
            VECTOR V2 = VEC3::Dot(p, linePoint2);
            VECTOR vD = VEC::Subtract(V1, V2);

            VECTOR vT = DotCoord(p, linePoint1);
            vT = VEC::Divide(vT, vD);

            VECTOR point = VEC::Subtract(linePoint2, linePoint1);
            point = VEC::MultiplyAdd(point, linePoint1);

            const VECTOR zero = VEC::Zero();
            VECTOR control = VEC::NearEqual(vD, zero, g_Epsilon.v);

            return VEC::Select(point, g_QNaN.v, control);
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV IntersectPlane(VECTOR* pLinePoint1, VECTOR* pLinePoint2, A_VECTOR P1, A_VECTOR P2) noexcept
        {
            assert(pLinePoint1);
            assert(pLinePoint2);

            VECTOR V1 = VEC3::Cross(P2, P1);

            VECTOR lengthSq = VEC3::LengthSq(V1);

            VECTOR V2 = VEC3::Cross(P2, V1);

            VECTOR P1W = VEC::SplatW(P1);
            VECTOR point = VEC::Multiply(V2, P1W);

            VECTOR V3 = VEC3::Cross(V1, P1);

            VECTOR P2W = VEC::SplatW(P2);
            point = VEC::MultiplyAdd(V3, P2W, point);

            VECTOR linePoint1 = VEC::Divide(point, lengthSq);

            VECTOR linePoint2 = VEC::Add(linePoint1, V1);

            VECTOR control = VEC::LessOrEqual(lengthSq, g_Epsilon.v);
            *pLinePoint1 = VEC::Select(linePoint1, g_QNaN.v, control);
            *pLinePoint2 = VEC::Select(linePoint2, g_QNaN.v, control);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV Transform(A_VECTOR p, A_MATRIX ITM) noexcept
        {
            VECTOR W = VEC::SplatW(p);
            VECTOR Z = VEC::SplatZ(p);
            VECTOR Y = VEC::SplatY(p);
            VECTOR X = VEC::SplatX(p);

            VECTOR Result = VEC::Multiply(W, ITM.r[3]);
            Result = VEC::MultiplyAdd(Z, ITM.r[2], Result);
            Result = VEC::MultiplyAdd(Y, ITM.r[1], Result);
            
            return VEC::MultiplyAdd(X, ITM.r[0], Result);
        }

        _Use_decl_annotations_
        FORCE_INLINE Float4* VEC_CALLCONV TransformStream(
            Float4* pOutputStream, 
            size_t outputStride, 
            const Float4* pInputStream, 
            size_t inputStride, 
            size_t planeCount, 
            A_MATRIX ITM
        ) noexcept
        {
            return VEC4::TransformStream(
                pOutputStream, outputStride, 
                pInputStream, inputStride, 
                planeCount, 
                ITM);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV FromPointNormal(A_VECTOR point, A_VECTOR normal) noexcept
        {
            VECTOR W = VEC3::Dot(point, normal);
            W = VEC::Negate(W);

            return VEC::Select(W, normal, g_Select1110.v);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV FromPoints(A_VECTOR point1, A_VECTOR point2, A_VECTOR point3) noexcept
        {
            VECTOR V21 = VEC::Subtract(point1, point2);
            VECTOR V31 = VEC::Subtract(point1, point3);

            VECTOR N = VEC3::Cross(V21, V31);
            N = VEC3::Normalize(N);

            VECTOR D = DotNormal(N, point1);
            D = VEC::Negate(D);

            return VEC::Select(D, N, g_Select1110.v);
        }
    } // namespace Plane
    
} // namespace UltReality::Math

#endif // !ULTREALITY_MATH_PLANE_INL