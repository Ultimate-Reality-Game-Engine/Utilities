#ifndef ULTREALITY_MATH_QUATERNION_INL
#define ULTREALITY_MATH_QUATERNION_INL

#if defined(__GNUC__) or defined(__clang__)
#define FORCE_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline
#endif

namespace UltReality::Math
{
    namespace Quaternion
    {
        FORCE_INLINE bool VEC_CALLCONV Equal(A_VECTOR Q1, A_VECTOR Q2) noexcept
        {
            return VEC4::Equal(Q1, Q2);
        }

        FORCE_INLINE bool VEC_CALLCONV NotEqual(A_VECTOR Q1, A_VECTOR Q2) noexcept
        {
            return VEC4::NotEqual(Q1, Q2);
        }

        FORCE_INLINE bool VEC_CALLCONV IsNaN(A_VECTOR Q1, A_VECTOR Q2) noexcept
        {
            return VEC4::IsNaN(Q1, Q2);
        }

        FORCE_INLINE bool VEC_CALLCONV IsInfinite(A_VECTOR q) noexcept
        {
            return VEC4::IsInfinite(q);
        }

        FORCE_INLINE bool VEC_CALLCONV IsIdentity(A_VECTOR q) noexcept
        {
            return VEC4::Equal(q, g_IdentityR3.v);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV Dot(A_VECTOR Q1, A_VECTOR Q2) noexcept
        {
            return VEC4::Dot(Q1, Q2);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV Multiply(A_VECTOR Q1, A_VECTOR Q2) noexcept
        {
            // Returns the product Q2*Q1 (which is the concatenation of a rotation Q1 followed by the rotation Q2)

            // [ (Q2.w * Q1.x) + (Q2.x * Q1.w) + (Q2.y * Q1.z) - (Q2.z * Q1.y),
            //   (Q2.w * Q1.y) - (Q2.x * Q1.z) + (Q2.y * Q1.w) + (Q2.z * Q1.x),
            //   (Q2.w * Q1.z) + (Q2.x * Q1.y) - (Q2.y * Q1.x) + (Q2.z * Q1.w),
            //   (Q2.w * Q1.w) - (Q2.x * Q1.x) - (Q2.y * Q1.y) - (Q2.z * Q1.z) ]

#if defined(_NO_INTRINSICS_)
            VECTOR_F32 Result = { { {
                    (Q2.vector4_f32[3] * Q1.vector4_f32[0]) + (Q2.vector4_f32[0] * Q1.vector4_f32[3]) + (Q2.vector4_f32[1] * Q1.vector4_f32[2]) - (Q2.vector4_f32[2] * Q1.vector4_f32[1]),
                    (Q2.vector4_f32[3] * Q1.vector4_f32[1]) - (Q2.vector4_f32[0] * Q1.vector4_f32[2]) + (Q2.vector4_f32[1] * Q1.vector4_f32[3]) + (Q2.vector4_f32[2] * Q1.vector4_f32[0]),
                    (Q2.vector4_f32[3] * Q1.vector4_f32[2]) + (Q2.vector4_f32[0] * Q1.vector4_f32[1]) - (Q2.vector4_f32[1] * Q1.vector4_f32[0]) + (Q2.vector4_f32[2] * Q1.vector4_f32[3]),
                    (Q2.vector4_f32[3] * Q1.vector4_f32[3]) - (Q2.vector4_f32[0] * Q1.vector4_f32[0]) - (Q2.vector4_f32[1] * Q1.vector4_f32[1]) - (Q2.vector4_f32[2] * Q1.vector4_f32[2])
                } } };
            
            return Result.v;

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 ControlWZYX = { { { 1.0f, -1.0f, 1.0f, -1.0f } } };
            static const VECTOR_F32 ControlZWXY = { { { 1.0f, 1.0f, -1.0f, -1.0f } } };
            static const VECTOR_F32 ControlYXWZ = { { { -1.0f, 1.0f, 1.0f, -1.0f } } };
            // Copy to SSE registers and use as few as possible for x86
            VECTOR Q2X = Q2;
            VECTOR Q2Y = Q2;
            VECTOR Q2Z = Q2;
            VECTOR vResult = Q2;
            // Splat with one instruction
            vResult = PERMUTE_PS(vResult, _MM_SHUFFLE(3, 3, 3, 3));
            Q2X = PERMUTE_PS(Q2X, _MM_SHUFFLE(0, 0, 0, 0));
            Q2Y = PERMUTE_PS(Q2Y, _MM_SHUFFLE(1, 1, 1, 1));
            Q2Z = PERMUTE_PS(Q2Z, _MM_SHUFFLE(2, 2, 2, 2));
            // Retire Q1 and perform Q1*Q2W
            vResult = _mm_mul_ps(vResult, Q1);
            VECTOR Q1Shuffle = Q1;
            // Shuffle the copies of Q1
            Q1Shuffle = PERMUTE_PS(Q1Shuffle, _MM_SHUFFLE(0, 1, 2, 3));
            // Mul by Q1WZYX
            Q2X = _mm_mul_ps(Q2X, Q1Shuffle);
            Q1Shuffle = PERMUTE_PS(Q1Shuffle, _MM_SHUFFLE(2, 3, 0, 1));
            // Flip the signs on y and z
            vResult = FMADD_PS(Q2X, ControlWZYX, vResult);
            // Mul by Q1ZWXY
            Q2Y = _mm_mul_ps(Q2Y, Q1Shuffle);
            Q1Shuffle = PERMUTE_PS(Q1Shuffle, _MM_SHUFFLE(0, 1, 2, 3));
            // Flip the signs on z and w
            Q2Y = _mm_mul_ps(Q2Y, ControlZWXY);
            // Mul by Q1YXWZ
            Q2Z = _mm_mul_ps(Q2Z, Q1Shuffle);
            // Flip the signs on x and w
            Q2Y = FMADD_PS(Q2Z, ControlYXWZ, Q2Y);
            vResult = _mm_add_ps(vResult, Q2Y);
            
            return vResult;
#endif
        }

        FORCE_INLINE VECTOR VEC_CALLCONV LengthSq(A_VECTOR q) noexcept
        {
            return LengthSq(q);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV ReciprocalLength(A_VECTOR q) noexcept
        {
            return VEC4::ReciprocalLength(q);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV Length(A_VECTOR q) noexcept
        {
            return VEC4::Length(q);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV NormalizeEst(A_VECTOR q) noexcept
        {
            return VEC4::NormalizeEst(q);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV Normalize(A_VECTOR q) noexcept
        {
            return VEC4::Normalize(q);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV Conjugate(A_VECTOR q) noexcept
        {
#if defined(_NO_INTRINSICS_)
            VECTOR_F32 Result = { { {
                    -q.vector4_f32[0],
                    -q.vector4_f32[1],
                    -q.vector4_f32[2],
                    q.vector4_f32[3]
                } } };
            
            return Result.v;

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 negativeOne3 = { { { -1.0f, -1.0f, -1.0f, 1.0f } } };
            
            return _mm_mul_ps(q, negativeOne3);
#endif
        }

        FORCE_INLINE VECTOR VEC_CALLCONV Inverse(A_VECTOR q) noexcept
        {
            VECTOR l = VEC4::LengthSq(q);
            VECTOR conjugate = Conjugate(q);

            VECTOR control = VEC::LessOrEqual(l, g_Epsilon.v);

            VECTOR Result = VEC::Divide(conjugate, l);

            Result = VEC::Select(Result, g_Zero, control);

            return Result;
        }

        FORCE_INLINE VECTOR VEC_CALLCONV Ln(A_VECTOR q) noexcept
        {
            static const VECTOR_F32 oneMinusEpsilon = { { { 1.0f - 0.00001f, 1.0f - 0.00001f, 1.0f - 0.00001f, 1.0f - 0.00001f } } };

            VECTOR QW = VEC::SplatW(Q);
            VECTOR Q0 = VEC::Select(g_Select1110.v, q, g_Select1110.v);

            VECTOR controlW = VEC::InBounds(QW, oneMinusEpsilon.v);

            VECTOR theta = VEC::ACos(QW);
            VECTOR sinTheta = VEC::Sine(theta);

            VECTOR S = VEC::Divide(theta, sinTheta);

            VECTOR Result = VEC::Multiply(Q0, S);
            Result = VEC::Select(Q0, Result, controlW);

            return Result;
        }

        FORCE_INLINE VECTOR VEC_CALLCONV Exp(A_VECTOR q) noexcept
        {
            VECTOR theta = VEC3::Length(q);

            VECTOR sineTheta, cosTheta;
            VEC::SineCos(&sineTheta, &cosTheta, theta);

            VETOR s = VEC::Divide(sineTheta, theta);

            VECTOR Result = VEC::Multiply(q, s);

            cost VECTOR zero = VEC::Zero();
            VECTOR control = VEC::NearEqual(theta, zero, g_Epsilon.v);
            Result = VEC::Select(Result, q, control);

            return VEC::Select(cosTheta, Result, g_Select1110.v);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV Slerp(A_VECTOR Q1, A_VECTOR Q2, float t) noexcept
        {
            VECTOR T = VEC::Replicate(t);

            return SlerpV(Q1, Q2, T);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV SlerpV(A_VECTOR Q1, A_VECTOR Q2, A_VECTOR T) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert((VEC::GetY(T) == VEC::GetX(T)) && (VEC::GetZ(T) == VEC::GetX(T)) && (VEC::GetW(T) == VEC::GetX(T)));
#endif

            // Result = Q0 * sine((1.0 - t) * Omega) / sine(Omega) + Q1 * sine(t * Omega) / sine(Omega)
#if defined(_NO_INTRINSICS_)
            const VECTOR_F32 oneMinusEpsilon = { { { 1.0f - 0.00001f, 1.0f - 0.00001f, 1.0f - 0.00001f, 1.0f - 0.00001f } } };

            VECTOR cosOmega = Dot(Q1, Q2);

            const VECTOR zero = VEC::Zero();
            VECTOR control = VEC::Less(cosOmega, zero);
            VECTOR sign = VEC::Select(g_One.v, g_NegativeOne.v, control);

            cosOmega = VEC::Multiply(cosOmega, sign);

            control = VEC::Less(cosOmega, oneMinusEpsilon);

            VECTOR sineOmega = VEC::NegativeMultiplySubtract(cosOmega, cosOmega, g_One.v);
            sineOmega = VEC::Sqrt(sineOmega);

            VECTOR omega = VEC::ATan2(sineOmega, cosOmega);

            VECTOR signMask = VEC::SplatSignMask();
            VECTOR V01 = VEC::ShiftLeft(T, zero, 2);
            signMask = VEC::ShiftLeft(signMask, zero, 3);
            V01 = VEC::XorInt(V01, signMask);
            V01 = VEC::Add(g_IdentityR0.v, V01);

            VECTOR invSinOmega = VEC::Reciprocal(sineOmega);

            VECTOR S0 = VEC::Multiply(V01, omega);
            S0 = VEC::Sin(S0);
            S0 = VEC::Multiply(S0, invSinOmega);

            S0 = VEC::Select(V01, S0, control);

            VECTOR S1 = VEC::SplatY(S0);
            S0 = VEC::SplatX(S0);

            S1 = VEC::Multiply(S1, sign);

            VECTOR Result = VEC::Multiply(Q0, S0);
            
            return VEC::MultiplyAdd(Q1, S1, Result);

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 oneMinusEpsilon = { { { 1.0f - 0.00001f, 1.0f - 0.00001f, 1.0f - 0.00001f, 1.0f - 0.00001f } } };
            static const VECTOR_U32 signMask2 = { { { 0x80000000, 0x00000000, 0x00000000, 0x00000000 } } };

            VECTOR cosOmega = Dot(Q1, Q2);

            const VECTOR zero = VEC::Zero();
            VECTOR control = VEC::Less(cosOmega, zero);
            VECTOR sign = VEC::Select(g_One, g_NegativeOne, control);

            cosOmega = _mm_mul_ps(cosOmega, sign);

            control = VEC::Less(cosOmega, oneMinusEpsilon);

            VECTOR sineOmega = _mm_mul_ps(cosOmega, cosOmega);
            sineOmega = _mm_sub_ps(g_One, sineOmega);
            sineOmega = _mm_sqrt_ps(sineOmega);

            VECTOR omega = VEC::ATan2(sineOmega, cosOmega);

            VECTOR V01 = PERMUTE_PS(T, _MM_SHUFFLE(2, 3, 0, 1));
            V01 = _mm_and_ps(V01, g_MaskXY);
            V01 = _mm_xor_ps(V01, signMask2);
            V01 = _mm_add_ps(g_IdentityR0, V01);

            VECTOR S0 = _mm_mul_ps(V01, omega);
            S0 = VEC::Sin(S0);
            S0 = _mm_div_ps(S0, sineOmega);

            S0 = VEC::Select(V01, S0, control);

            VECTOR S1 = VEC::SplatY(S0);
            S0 = VEC::SplatX(S0);

            S1 = _mm_mul_ps(S1, sign);
            VECTOR Result = _mm_mul_ps(Q1, S0);
            S1 = _mm_mul_ps(S1, Q2);
            
            return _mm_add_ps(Result, S1);
#endif
        }

        FORCE_INLINE VECTOR VEC_CALLCONV Squad(A_VECTOR Q1, A_VECTOR Q2, A_VECTOR Q3, B_VECTOR Q4, float t) noexcept
        {
            VECTOR T = VEC::Replicate(t);

            return SquadV(Q1, Q2, Q3, Q4, T);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV SquadV(A_VECTOR Q1, A_VECTOR Q2, A_VECTOR Q3, B_VECTOR Q4, C_VECTOR T) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert((VEC::GetY(T) == VEC::GetX(T)) && (VEC::GetZ(T) == VEC::GetX(T)) && (VEC::GetW(T) == VEC::GetX(T)));
#endif

            VECTOR TP = T;
            const VECTOR two = VEC::SplatConstant(2, 0);

            VECTOR Q03 = SlerpV(Q1, Q4, T);
            VECTOR Q12 = SlerpV(Q2, Q3, T);

            TP = VEC::NegativeMultiplySubtract(TP, TP, TP);
            TP = VEC::Multiply(TP, two);

            return SlerpV(Q03, Q12, TP);
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV SquadSetup(VECTOR* pA, VECTOR* pB, VECTOR* pC, 
            A_VECTOR Q1, A_VECTOR Q2, A_VECTOR Q3, B_VECTOR Q4
        ) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pA != nullptr);
            assert(pB != nullptr);
            assert(pC != nullptr);
#endif

            VECTOR LS23 = LengthSq(VEC::Add(Q2, Q3));
            VECTOR LD23 = LengthSq(VEC::Subtract(Q2, Q3));
            VECTOR SQ3 = VEC::Negate(Q3);

            VECTOR control1 = VEC::Less(LS23, LD23);
            SQ3 = VEC::Select(Q3, SQ3, control1);

            VECTOR LS12 = LengthSq(VEC::Add(Q1, Q2));
            VECTOR LD12 = LengthSq(VEC::Subtract(Q1, Q2));
            VECTOR SQ1 = VEC::Negate(Q1);

            VECTOR LS34 = LengthSq(VEC::Add(SQ3, Q4));
            VECTOR LD34 = LengthSq(VEC::Subtract(SQ3, Q4));
            VECTOR SQ4 = VEC::Negate(Q4);

            VECTOR control0 = VEC::Less(LS12, LD12);
            VECTOR control2 = VEC::Less(LS34, LD34);

            SQ1 = VEC::Select(Q1, SQ1, control0);
            SQ4 = VEC::Select(Q4, SQ4, control2);

            VECTOR invQ2 = Inverse(Q2);
            VECTOR invQ3 = Inverse(SQ3);

            VECTOR LnQ1 = Ln(Multiply(invQ2, SQ1));
            VECTOR LnQ3 = Ln(Multiply(invQ2, SQ3));
            VECTOR LnQ2 = Ln(Multiply(invQ3, Q2));
            VECTOR LnQ4 = Ln(Multiply(invQ3, SQ4));

            const VECTOR negativeOneQuarter = VEC::SplatConstant(-1, 2);

            VECTOR expQ13 = VEC::Multiply(VEC::Add(LnQ1, LnQ3), negativeOneQuarter);
            VECTOR expQ24 = VEC::Multiply(VEC::Add(LnQ2, LnQ4), negativeOneQuarter);
            expQ13 = Exp(expQ13);
            expQ24 = Exp(expQ24);

            *pA = Multiply(Q2, expQ13);
            *pB = Multiply(SQ3, expQ24);
            *pC = SQ3;
        }

        FORCE_INLINE VECTOR VEC_CALLCONV BaryCentric(A_VECTOR Q1, A_VECTOR Q2, A_VECTOR Q3, float f, float g) noexcept
        {
            float s = f + g;

            VECTOR Result;
            if((s < 0.00001f) && (s > -0.00001f))
            {
                Result = Q1;
            }
            else
            {
                VECTOR Q12 = Slerp(Q1, Q2, s);
                VECTOR Q13 = Slerp(Q1, Q3, s);

                Result = Slerp(Q12, Q13, g / s);
            }

            return Result;
        }

        FORCE_INLINE VECTOR VEC_CALLCONV BaryCentricV(A_VECTOR Q1, A_VECTOR Q2, A_VECTOR Q3, B_VECTOR F, C_VECTOR G) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert((VEC::GetY(F) == VEC::GetX(F)) && (VEC::GetZ(F) == VEC::GetX(F)) && (VEC::GetW(F) == VEC::GetX(F)));
            assert((VEC::GetY(G) == VEC::GetX(G)) && (VEC::GetZ(G) == VEC::GetX(G)) && (VEC::GetW(G) == VEC::GetX(G)));
#endif

            const VECTOR epsilon = VEC::SplatConstant(1, 16);

            VECTOR s = VEC::Add(F, G);

            VECTOR Result;
            if(VEC4::InBounds(s, epsilon))
            {
                Result = Q1;
            }
            else
            {
                VECTOR Q12 = SlerpV(Q1, Q2, s);
                VECTOR Q13 = SlerpV(Q1, Q3, s);
                VECTOR GS = VEC::Reciprocal(s);
                GS = VEC::Multiply(G, GS);

                Result = SlerpV(Q12, Q13, GS);
            }

            return Result;
        }

        FORCE_INLINE VECTOR VEC_CALLCONV Identity() noexcept
        {
            return g_IdentityR3.v;
        }

        FORCE_INLINE VECTOR VEC_CALLCONV RotationPitchYawRoll(float pitch, float yaw, float roll) noexcept
        {
#if defined(_NO_INTRINSICS_)
            const float halfPitch = pitch * 0.5f;
            float cp = cosf(halfPitch);
            float sp = sinf(halfPitch);

            const float halfYaw = yaw * 0.5f;
            float cy = cosf(halfYaw);
            float sy = sinf(halfYaw);

            const float halfRoll = roll * 0.5f;
            float cr = cosf(halfRoll);
            float sr = sinf(halfRoll);

            VECTOR_F32 vResult = { { {
                    cr * sp * cy + sr * cp * sy,
                    cr * cp * sy - sr * sp * cy,
                    sr * cp * cy - cr * sp * sy,
                    cr * cp * cy + sr * sp * sy
                } } };

            return vResult;

#else
            VECTOR angles = VEC::Set(pitch, yaw, roll, 0.0f);

            return RotationPitchYawRollFromVector(angles);
#endif
        }

        FORCE_INLINE VECTOR VEC_CALLCONV RotationPitchYawRollFromVector(A_VECTOR angles) noexcept
        {
#if defined(_NO_INTRINSICS_)
            const float halfPitch = angles.vector4_f32[0] * 0.5f;
            float cp = cosf(halfPitch);
            float sp = sinf(halfPitch);

            const float halfYaw = angles.vector4_f32[1] * 0.5f;
            float cy = cosf(halfYaw);
            float sy = sinf(halfYaw);

            const float halfRoll = angles.vector4_f32[2] * 0.5f;
            float cr = cosf(halfRoll);
            float sr = sinf(halfRoll);

            VECTOR_F32 vResult = { { {
                    cr * sp * cy + sr * cp * sy,
                    cr * cp * sy - sr * sp * cy,
                    sr * cp * cy - cr * sp * sy,
                    cr * cp * cy + sr * sp * sy
                } } };

            return vResult;
#else
            static const VECTOR_F32  sign = { { { 1.0f, -1.0f, -1.0f, 1.0f } } };

            VECTOR halfAngles = VEC::Multiply(angles, g_OneHalf.v);

            VECTOR sineAngles, cosAngles;
            VEC::SinCos(&sineAngles, &cosAngles, halfAngles);

            VECTOR P0 = VEC::Permute<PERMUTE_0X, PERMUTE_1X, PERMUTE_1X, PERMUTE_1X>(sineAngles, cosAngles);
            VECTOR Y0 = VEC::Permute<PERMUTE_1Y, PERMUTE_0Y, PERMUTE_1Y, PERMUTE_1Y>(sineAngles, cosAngles);
            VECTOR R0 = VEC::Permute<PERMUTE_1Z, PERMUTE_1Z, PERMUTE_0Z, PERMUTE_1Z>(sineAngles, cosAngles);
            VECTOR P1 = VEC::Permute<PERMUTE_0X, PERMUTE_1X, PERMUTE_1X, PERMUTE_1X>(cosAngles, sineAngles);
            VECTOR Y1 = VEC::Permute<PERMUTE_1Y, PERMUTE_0Y, PERMUTE_1Y, PERMUTE_1Y>(cosAngles, sineAngles);
            VECTOR R1 = VEC::Permute<PERMUTE_1Z, PERMUTE_1Z, PERMUTE_0Z, PERMUTE_1Z>(cosAngles, sineAngles);

            VECTOR Q1 = VEC::Multiply(P1, sign.v);
            VECTOR Q0 = VEC::Multiply(P0, Y0);
            Q1 = VEC::Multiply(Q1, Y1);
            Q0 = VEC::Multiply(Q0, R0);
            
            return VEC::MultiplyAdd(Q1, R1, Q0);
#endif
        }

        FORCE_INLINE VECTOR VEC_CALLCONV RotationNormal(A_VECTOR normalAxis, float angle) noexcept
        {
#if defined(_NO_INTRINSICS_)
            VECTOR N = VEC::Select(g_One.v, normalAxis, g_Select1110.v);

            float sineV, cosV;
            VEC::SinCos(&sineV, &cosV, 0.5f * angle);

            VECTOR scale = VEC::Set(sineV, sineV, sineV, cosV);
            
            return Multiply(N, scale);

#elif defined(_SSE2_INTRINSICS_)
            VECTOR N = _mm_and_ps(normalAxis, g_Mask3);
            N = _mm_or_ps(N, g_IdentityR3);
            VECTOR scale = _mm_set_ps1(0.5f * angle);
            VECTOR vSine;
            VECTOR vCosine;
            VEC::SinCos(&vSine, &vCosine, scale);
            scale = _mm_and_ps(vSine, g_Mask3);
            vCosine = _mm_and_ps(vCosine, g_MaskW);
            scale = _mm_or_ps(scale, vCosine);
            
            return _mm_mul_ps(N, scale);
#endif
        }

        FORCE_INLINE VECTOR VEC_CALLCONV RotationAxis(A_VECTOR axis, float angle) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(!VEC3::Equal(axis, VEC::Zero()));
            assert(!VEC3::IsInfinite(axis));
#endif

            VECTOR normal = VEC3::Normalize(axis);
            
            return RotationNormal(normal, angle);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV RotationMatrix(A_MATRIX m) noexcept
        {
#if defined(_NO_INTRINSICS_)
            VECTOR_F32 q;
            float r22 = m.m[2][2];
            if(r22 <= 0.0f) // x^2 + y^2 >= z^2 + w^2
            {
                float dif10 = m.m[1][1] - m.m[0][0];
                float omr22 = 1.0f - r22;
                if (dif10 <= 0.0f)  // x^2 >= y^2
                {
                    float fourXSqr = omr22 - dif10;
                    float inv4x = 0.5f / sqrtf(fourXSqr);
                    q.f[0] = fourXSqr * inv4x;
                    q.f[1] = (m.m[0][1] + m.m[1][0]) * inv4x;
                    q.f[2] = (m.m[0][2] + m.m[2][0]) * inv4x;
                    q.f[3] = (m.m[1][2] - m.m[2][1]) * inv4x;
                }
                else  // y^2 >= x^2
                {
                    float fourYSqr = omr22 + dif10;
                    float inv4y = 0.5f / sqrtf(fourYSqr);
                    q.f[0] = (m.m[0][1] + m.m[1][0]) * inv4y;
                    q.f[1] = fourYSqr * inv4y;
                    q.f[2] = (m.m[1][2] + m.m[2][1]) * inv4y;
                    q.f[3] = (m.m[2][0] - m.m[0][2]) * inv4y;
                }
            }
            else  // z^2 + w^2 >= x^2 + y^2
            {
                float sum10 = m.m[1][1] + m.m[0][0];
                float opr22 = 1.0f + r22;
                if (sum10 <= 0.0f)  // z^2 >= w^2
                {
                    float fourZSqr = opr22 - sum10;
                    float inv4z = 0.5f / sqrtf(fourZSqr);
                    q.f[0] = (m.m[0][2] + m.m[2][0]) * inv4z;
                    q.f[1] = (m.m[1][2] + m.m[2][1]) * inv4z;
                    q.f[2] = fourZSqr * inv4z;
                    q.f[3] = (m.m[0][1] - m.m[1][0]) * inv4z;
                }
                else  // w^2 >= z^2
                {
                    float fourWSqr = opr22 + sum10;
                    float inv4w = 0.5f / sqrtf(fourWSqr);
                    q.f[0] = (m.m[1][2] - m.m[2][1]) * inv4w;
                    q.f[1] = (m.m[2][0] - m.m[0][2]) * inv4w;
                    q.f[2] = (m.m[0][1] - m.m[1][0]) * inv4w;
                    q.f[3] = fourWSqr * inv4w;
                }
            }

            return q.v;

#elif defined(_SSE2_INTRINSICS_)
            static const VECTOR_F32 XMPMMP = { { { +1.0f, -1.0f, -1.0f, +1.0f } } };
            static const VECTOR_F32 XMMPMP = { { { -1.0f, +1.0f, -1.0f, +1.0f } } };
            static const VECTOR_F32 XMMMPP = { { { -1.0f, -1.0f, +1.0f, +1.0f } } };

            VECTOR r0 = m.r[0];  // (r00, r01, r02, 0)
            VECTOR r1 = m.r[1];  // (r10, r11, r12, 0)
            VECTOR r2 = m.r[2];  // (r20, r21, r22, 0)

            // (r00, r00, r00, r00)
            VECTOR r00 = PERMUTE_PS(r0, _MM_SHUFFLE(0, 0, 0, 0));
            // (r11, r11, r11, r11)
            VECTOR r11 = PERMUTE_PS(r1, _MM_SHUFFLE(1, 1, 1, 1));
            // (r22, r22, r22, r22)
            VECTOR r22 = PERMUTE_PS(r2, _MM_SHUFFLE(2, 2, 2, 2));

            // x^2 >= y^2 equivalent to r11 - r00 <= 0
            // (r11 - r00, r11 - r00, r11 - r00, r11 - r00)
            VECTOR r11mr00 = _mm_sub_ps(r11, r00);
            VECTOR x2gey2 = _mm_cmple_ps(r11mr00, g_Zero);

            // z^2 >= w^2 equivalent to r11 + r00 <= 0
            // (r11 + r00, r11 + r00, r11 + r00, r11 + r00)
            VECTOR r11pr00 = _mm_add_ps(r11, r00);
            VECTOR z2gew2 = _mm_cmple_ps(r11pr00, g_Zero);

            // x^2 + y^2 >= z^2 + w^2 equivalent to r22 <= 0
            VECTOR x2py2gez2pw2 = _mm_cmple_ps(r22, g_Zero);

            // (4*x^2, 4*y^2, 4*z^2, 4*w^2)
            VECTOR t0 = FMADD_PS(XMPMMP, r00, g_One);
            VECTOR t1 = _mm_mul_ps(XMMPMP, r11);
            VECTOR t2 = FMADD_PS(XMMMPP, r22, t0);
            VECTOR x2y2z2w2 = _mm_add_ps(t1, t2);

            // (r01, r02, r12, r11)
            t0 = _mm_shuffle_ps(r0, r1, _MM_SHUFFLE(1, 2, 2, 1));
            // (r10, r10, r20, r21)
            t1 = _mm_shuffle_ps(r1, r2, _MM_SHUFFLE(1, 0, 0, 0));
            // (r10, r20, r21, r10)
            t1 = PERMUTE_PS(t1, _MM_SHUFFLE(1, 3, 2, 0));
            // (4*x*y, 4*x*z, 4*y*z, unused)
            VECTOR xyxzyz = _mm_add_ps(t0, t1);

            // (r21, r20, r10, r10)
            t0 = _mm_shuffle_ps(r2, r1, _MM_SHUFFLE(0, 0, 0, 1));
            // (r12, r12, r02, r01)
            t1 = _mm_shuffle_ps(r1, r0, _MM_SHUFFLE(1, 2, 2, 2));
            // (r12, r02, r01, r12)
            t1 = PERMUTE_PS(t1, _MM_SHUFFLE(1, 3, 2, 0));
            // (4*x*w, 4*y*w, 4*z*w, unused)
            VECTOR xwywzw = _mm_sub_ps(t0, t1);
            xwywzw = _mm_mul_ps(XMMPMP, xwywzw);

            // (4*x^2, 4*y^2, 4*x*y, unused)
            t0 = _mm_shuffle_ps(x2y2z2w2, xyxzyz, _MM_SHUFFLE(0, 0, 1, 0));
            // (4*z^2, 4*w^2, 4*z*w, unused)
            t1 = _mm_shuffle_ps(x2y2z2w2, xwywzw, _MM_SHUFFLE(0, 2, 3, 2));
            // (4*x*z, 4*y*z, 4*x*w, 4*y*w)
            t2 = _mm_shuffle_ps(xyxzyz, xwywzw, _MM_SHUFFLE(1, 0, 2, 1));

            // (4*x*x, 4*x*y, 4*x*z, 4*x*w)
            VECTOR tensor0 = _mm_shuffle_ps(t0, t2, _MM_SHUFFLE(2, 0, 2, 0));
            // (4*y*x, 4*y*y, 4*y*z, 4*y*w)
            VECTOR tensor1 = _mm_shuffle_ps(t0, t2, _MM_SHUFFLE(3, 1, 1, 2));
            // (4*z*x, 4*z*y, 4*z*z, 4*z*w)
            VECTOR tensor2 = _mm_shuffle_ps(t2, t1, _MM_SHUFFLE(2, 0, 1, 0));
            // (4*w*x, 4*w*y, 4*w*z, 4*w*w)
            VECTOR tensor3 = _mm_shuffle_ps(t2, t1, _MM_SHUFFLE(1, 2, 3, 2));

            // Select the row of the tensor-product matrix that has the largest
            // magnitude.
            t0 = _mm_and_ps(x2gey2, tensor0);
            t1 = _mm_andnot_ps(x2gey2, tensor1);
            t0 = _mm_or_ps(t0, t1);
            t1 = _mm_and_ps(z2gew2, tensor2);
            t2 = _mm_andnot_ps(z2gew2, tensor3);
            t1 = _mm_or_ps(t1, t2);
            t0 = _mm_and_ps(x2py2gez2pw2, t0);
            t1 = _mm_andnot_ps(x2py2gez2pw2, t1);
            t2 = _mm_or_ps(t0, t1);

            // Normalize the row.  No division by zero is possible because the
            // quaternion is unit-length (and the row is a nonzero multiple of
            // the quaternion).
            t0 = VEC4::Length(t2);
            
            return _mm_div_ps(t2, t0);
#endif
        }

        _Use_decl_annotations_
        FORCE_INLINE void VEC_CALLCONV ToAxisAngle(VECTOR* pAxis, float* pAngle, A_VECTOR q) noexcept
        {
#if defined(DEBUG) || defined(_DEBUG)
            assert(pAxis != nullptr);
            assert(pAngle != nullptr);
#endif

            *pAxis = q;
            *pAngle = 2.0f * ScalarACos(VEC::GetW(q));
        }
    }
}

#endif // !ULTREALITY_MATH_QUATERNION_INL