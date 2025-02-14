#ifndef ULTREALITY_MATH_COLOR_INL
#define ULTREALITY_MATH_COLOR_INL

#if defined(__GNUC__) or defined(__clang__)
#define FORCE_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline
#endif

namespace UltReality::Math
{
    namespace Color
    {
        FORCE_INLINE bool VEC_CALLCONV Equal(A_VECTOR C1, A_VECTOR C2) noexcept
        {
            return Vector4::Equal(C1, C2);
        }

        FORCE_INLINE bool VEC_CALLCONV NotEqual(A_VECTOR C1, A_VECTOR C2) noexcept
        {
            return Vector4::NotEqual(C1, C2);
        }

        FORCE_INLINE bool VEC_CALLCONV Greater(A_VECTOR C1, A_VECTOR C2) noexcept
        {
            return Vector4::Greater(C1, C2);
        }

        FORCE_INLINE bool VEC_CALLCONV GreaterOrEqual(A_VECTOR C1, A_VECTOR C2) noexcept
        {
            return Vector4::GreaterOrEqual(C1, C2);
        }

        FORCE_INLINE bool VEC_CALLCONV Less(A_VECTOR C1, A_VECTOR C2) noexcept
        {
            return Vector4::Less(C1, C2);
        }

        FORCE_INLINE bool VEC_CALLCONV LessOrEqual(A_VECTOR C1, A_VECTOR C2) noexcept
        {
            return Vector4::LessOrEqual(C1, C2);
        }

        FORCE_INLINE bool VEC_CALLCONV IsNaN(A_VECTOR c) noexcept
        {
            return Vector4::IsNaN(c);
        }

        FORCE_INLINE bool VEC_CALLCONV IsInfinite(A_VECTOR c) noexcept
        {
            return Vector4::IsInfinite(c);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV Negative(A_VECTOR c) noexcept
        {
#if defined(_NO_INTRINSICS_)
            VECTOR_F32 Result = { { {
                1.0f - c.vector4_f32[0], 
                1.0f - c.vector4_f32[1], 
                1.0f - c.vector4_f32[2], 
                c.vector4_f32[3]
            } } };

            return Result.v;

#elif defined(_SSE2_INTRINSICS_)
            // Negate only x, y, and z
            VECTOR vTemp = _mm_xor_ps(c, g_Negate3);

            // Add 1, 1, 1, 0 to -x, -y, -z, w
            return _mm_add_ps(vTemp, g_One3);
#endif
        }

        FORCE_INLINE VECTOR VEC_CALLCONV Modulate(A_VECTOR C1, A_VECTOR C2) noexcept
        {
            return Vector::Multiply(C1, C2);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV AdjustSaturation(A_VECTOR c, float saturation) noexcept
        {
            // Luminance = 0.2125f * C[0] + 0.7154f * C[1] + 0.0721f * C[2];
            // Result = (C - Luminance) * Saturation + Luminance;

            const VECTOR_F32 gvLuminance = { { { 0.2125f, 0.7154f, 0.0721f, 0.0f } } };
#if defined(_NO_INTRINSICS_)
            float fLuminance = (c.vector4_f32[0] * gvLuminance.f[0]) + (c.vector4_f32[1] * gvLuminance.f[1]) + (c.vector4_f32[2] * gvLuminance.f[2]);
            VECTOR vResult;
            vResult.vector4_f32[0] = ((c.vector4_f32[0] - fLuminance) * fSaturation) + fLuminance;
            vResult.vector4_f32[1] = ((c.vector4_f32[1] - fLuminance) * fSaturation) + fLuminance;
            vResult.vector4_f32[2] = ((c.vector4_f32[2] - fLuminance) * fSaturation) + fLuminance;
            vResult.vector4_f32[3] = c.vector4_f32[3];
            
            return vResult;

#elif defined(_SSE2_INTRINSICS_)
            VECTOR vLuminance = Vector3::Dot(c, gvLuminance);
            // Splat saturation
            VECTOR vSaturation = _mm_set_ps1(saturation);
            // vResult = ((vColor-vLuminance)*vSaturation)+vLuminance;
            VECTOR vResult = _mm_sub_ps(c, vLuminance);
            vResult = FMADD_PS(vResult, vSaturation, vLuminance);
            // Retain w from the source color
            vLuminance = _mm_shuffle_ps(vResult, c, _MM_SHUFFLE(3, 2, 2, 2));   // x = vResult.z,y = vResult.z,z = vColor.z,w=vColor.w
            
            return _mm_shuffle_ps(vResult, vLuminance, _MM_SHUFFLE(3, 0, 1, 0));  // x = vResult.x,y = vResult.y,z = vResult.z,w=vColor.w
#endif
        }

        FORCE_INLINE VECTOR VEC_CALLCONV AdjustContrast(A_VECTOR c, float contrast) noexcept
        {
            // Result = (vColor - 0.5f) * fContrast + 0.5f;

#if defined(_NO_INTRINSICS_)
            VECTOR_F32 Result = { { {
                    ((c.vector4_f32[0] - 0.5f) * fContrast) + 0.5f,
                    ((c.vector4_f32[1] - 0.5f) * fContrast) + 0.5f,
                    ((c.vector4_f32[2] - 0.5f) * fContrast) + 0.5f,
                    c.vector4_f32[3]        // Leave W untouched
                } } };
            
            return Result.v;

#elif defined(_SSE2_INTRINSICS_)
            VECTOR vScale = _mm_set_ps1(contrast); // Splat the scale
            VECTOR vResult = _mm_sub_ps(c, g_OneHalf); // Subtract 0.5f from the source (Saving source)
            vResult = FMADD_PS(vResult, vScale, g_OneHalf);
            // Retain w from the source color
            vScale = _mm_shuffle_ps(vResult, c, _MM_SHUFFLE(3, 2, 2, 2)); // x = vResult.z,y = vResult.z,z = c.z,w=vColor.w
            
            return _mm_shuffle_ps(vResult, vScale, _MM_SHUFFLE(3, 0, 1, 0)); // x = vResult.x,y = vResult.y,z = vResult.z,w=c.w
#endif
        }

        FORCE_INLINE VECTOR VEC_CALLCONV RGBToHSL(A_VECTOR rgb) noexcept
        {
            VECTOR r = Vector::SplatX(rgb);
            VECTOR g = Vector::SplatY(rgb);
            VECTOR b = Vector::SplatZ(rgb);

            VECTOR min = Vector::Min(r, Vector::Min(g, b));
            VECTOR max = Vector::Max(r, Vector::Max(g, b));

            VECTOR l = Vector::Multiply(Vector::Add(min, max), g_OneHalf);

            VECTOR d = Vector::Subtract(max, min);

            VECTOR la = Vector::Select(rgb, l, g_Select1110);

            if (Vector3::Less(d, g_Epsilon))
            {
                // Achromatic, assume H and S of 0
                return Vector::Select(la, g_Zero, g_Select1100);
            }
            else
            {
                VECTOR s, h;

                VECTOR d2 = Vector::Add(min, max);

                if (Vector3::Greater(l, g_OneHalf))
                {
                    // d / (2-max-min)
                    s = Vector::Divide(d, Vector::Subtract(g_Two, d2));
                }
                else
                {
                    // d / (max+min)
                    s = Vector::Divide(d, d2);
                }

                if (Vector3::Equal(r, max))
                {
                    // Red is max
                    h = Vector::Divide(Vector::Subtract(g, b), d);
                }
                else if (Vector3::Equal(g, max))
                {
                    // Green is max
                    h = Vector::Divide(Vector::Subtract(b, r), d);
                    h = Vector::Add(h, g_Two);
                }
                else
                {
                    // Blue is max
                    h = Vector::Divide(Vector::Subtract(r, g), d);
                    h = Vector::Add(h, g_Four);
                }

                h = Vector::Divide(h, g_Six);

                if (Vector3::Less(h, g_Zero))
                    h = Vector::Add(h, g_One);

                VECTOR lha = Vector::Select(la, h, g_Select1100);
                
                return Vector::Select(s, lha, g_Select1011);
            }
        }

        namespace
        {
            FORCE_INLINE VECTOR VEC_CALLCONV Hue2Clr(A_VECTOR p, A_VECTOR q, A_VECTOR h) noexcept
            {
                static const VECTOR_F32 oneSixth = { { { 1.0f / 6.0f, 1.0f / 6.0f, 1.0f / 6.0f, 1.0f / 6.0f } } };
                static const VECTOR_F32 twoThirds = { { { 2.0f / 3.0f, 2.0f / 3.0f, 2.0f / 3.0f, 2.0f / 3.0f } } };

                VECTOR t = h;

                if (Vector3::Less(t, g_Zero))
                    t = Vector::Add(t, g_One);

                if (Vector3::Greater(t, g_One))
                    t = Vector::Subtract(t, g_One);

                if (Vector3::Less(t, oneSixth))
                {
                    // p + (q - p) * 6 * t
                    VECTOR t1 = Vector::Subtract(q, p);
                    VECTOR t2 = Vector::Multiply(g_Six, t);
                    return Vector::MultiplyAdd(t1, t2, p);
                }

                if (Vector3::Less(t, g_OneHalf))
                    return q;

                if (Vector3::Less(t, twoThirds))
                {
                    // p + (q - p) * 6 * (2/3 - t)
                    VECTOR t1 = Vector::Subtract(q, p);
                    VECTOR t2 = Vector::Multiply(g_Six, Vector::Subtract(twoThirds, t));
                    return Vector::MultiplyAdd(t1, t2, p);
                }

                return p;
            }
        }

        FORCE_INLINE VECTOR VEC_CALLCONV HSLToRGB(A_VECTOR hsl) noexcept
        {
            static const VECTOR_F32 oneThird = { { { 1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 3.0f } } };

            VECTOR s = Vector::SplatY(hsl);
            VECTOR l = Vector::SplatZ(hsl);

            if (Vector3::NearEqual(s, g_Zero, g_Epsilon))
            {
                // Achromatic
                return Vector::Select(hsl, l, g_Select1110);
            }
            else
            {
                VECTOR h = Vector::SplatX(hsl);

                VECTOR q;
                if (Vector3::Less(l, g_OneHalf))
                {
                    q = Vector::Multiply(l, Vector::Add(g_One, s));
                }
                else
                {
                    q = Vector::Subtract(Vector::Add(l, s), Vector::Multiply(l, s));
                }

                VECTOR p = Vector::Subtract(Vector::Multiply(g_Two, l), q);

                VECTOR r = Hue2Clr(p, q, Vector::Add(h, oneThird));
                VECTOR g = Hue2Clr(p, q, h);
                VECTOR b = Hue2Clr(p, q, Vector::Subtract(h, oneThird));

                VECTOR rg = Vector::Select(g, r, g_Select1000);
                VECTOR ba = Vector::Select(hsl, b, g_Select1110);

                return Vector::Select(ba, rg, g_Select1100);
            }
        }

        FORCE_INLINE VECTOR VEC_CALLCONV RGBToHSV(A_VECTOR rgb) noexcept
        {
            VECTOR r = Vector::SplatX(rgb);
            VECTOR g = Vector::SplatY(rgb);
            VECTOR b = Vector::SplatZ(rgb);

            VECTOR min = Vector::Min(r, Vector::Min(g, b));
            VECTOR v = Vector::Max(r, Vector::Max(g, b));

            VECTOR d = Vector::Subtract(v, min);

            VECTOR s = (Vector3::NearEqual(v, g_Zero, g_Epsilon)) ? g_Zero : Vector::Divide(d, v);

            if (Vector3::Less(d, g_Epsilon))
            {
                // Achromatic, assume H of 0
                VECTOR hv = Vector::Select(v, g_Zero, g_Select1000);
                VECTOR hva = Vector::Select(rgb, hv, g_Select1110);
                return Vector::Select(s, hva, g_Select1011);
            }
            else
            {
                VECTOR h;

                if (Vector3::Equal(r, v))
                {
                    // Red is max
                    h = Vector::Divide(Vector::Subtract(g, b), d);

                    if (Vector3::Less(g, b))
                        h = Vector::Add(h, g_Six);
                }
                else if (Vector3::Equal(g, v))
                {
                    // Green is max
                    h = Vector::Divide(Vector::Subtract(b, r), d);
                    h = Vector::Add(h, g_Two);
                }
                else
                {
                    // Blue is max
                    h = Vector::Divide(Vector::Subtract(r, g), d);
                    h = Vector::Add(h, g_Four);
                }

                h = Vector::Divide(h, g_Six);

                VECTOR hv = Vector::Select(v, h, g_Select1000);
                VECTOR hva = Vector::Select(rgb, hv, g_Select1110);
                return Vector::Select(s, hva, g_Select1011);
            }
        }

        FORCE_INLINE VECTOR VEC_CALLCONV HSVToRGB(A_VECTOR hsv) noexcept
        {
            VECTOR h = Vector::SplatX(hsv);
            VECTOR s = Vector::SplatY(hsv);
            VECTOR v = Vector::SplatZ(hsv);

            VECTOR h6 = Vector::Multiply(h, g_Six);

            VECTOR i = Vector::Floor(h6);
            VECTOR f = Vector::Subtract(h6, i);

            // p = v* (1-s)
            VECTOR p = Vector::Multiply(v, Vector::Subtract(g_One, s));

            // q = v*(1-f*s)
            VECTOR q = Vector::Multiply(v, Vector::Subtract(g_One, Vector::Multiply(f, s)));

            // t = v*(1 - (1-f)*s)
            VECTOR t = Vector::Multiply(v, Vector::Subtract(g_One, Vector::Multiply(Vector::Subtract(g_One, f), s)));

            auto ii = static_cast<int>(Vector::GetX(Vector::Mod(i, g_Six)));

            VECTOR _rgb;

            switch (ii)
            {
            case 0: // rgb = vtp
            {
                VECTOR vt = Vector::Select(t, v, g_Select1000);
                _rgb = Vector::Select(p, vt, g_Select1100);
                break;
            }
            case 1: // rgb = qvp
            {
                VECTOR qv = Vector::Select(v, q, g_Select1000);
                _rgb = Vector::Select(p, qv, g_Select1100);
                break;
            }
            case 2: // rgb = pvt
            {
                VECTOR pv = Vector::Select(v, p, g_Select1000);
                _rgb = Vector::Select(t, pv, g_Select1100);
                break;
            }
            case 3: // rgb = pqv
            {
                VECTOR pq = Vector::Select(q, p, g_Select1000);
                _rgb = Vector::Select(v, pq, g_Select1100);
                break;
            }
            case 4: // rgb = tpv
            {
                VECTOR tp = Vector::Select(p, t, g_Select1000);
                _rgb = Vector::Select(v, tp, g_Select1100);
                break;
            }
            default: // rgb = vpq
            {
                VECTOR vp = Vector::Select(p, v, g_Select1000);
                _rgb = Vector::Select(q, vp, g_Select1100);
                break;
            }
            }

            return Vector::Select(hsv, _rgb, g_Select1110);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV RGBToYUV(A_VECTOR rgb) noexcept
        {
            static const VECTOR_F32 scale0 = { { { 0.299f, -0.147f, 0.615f, 0.0f } } };
            static const VECTOR_F32 scale1 = { { { 0.587f, -0.289f, -0.515f, 0.0f } } };
            static const VECTOR_F32 scale2 = { { { 0.114f, 0.436f, -0.100f, 0.0f } } };

            MATRIX M(scale0, scale1, scale2, g_Zero);
            VECTOR clr = Vector3::Transform(rgb, M);

            return Vector::Select(rgb, clr, g_Select1110);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV YUVToRGB(A_VECTOR yuv) noexcept
        {
            static const VECTOR_F32 scale1 = { { { 0.0f, -0.395f, 2.032f, 0.0f } } };
            static const VECTOR_F32 scale2 = { { { 1.140f, -0.581f, 0.0f, 0.0f } } };

            MATRIX M(g_One, scale1, scale2, g_Zero);
            VECTOR clr = Vector3::Transform(yuv, M);

            return Vector::Select(yuv, clr, g_Select1110);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV RGBToYUV_HD(A_VECTOR rgb) noexcept
        {
            static const VECTOR_F32 scale0 = { { { 0.2126f, -0.0997f, 0.6150f, 0.0f } } };
            static const VECTOR_F32 scale1 = { { { 0.7152f, -0.3354f, -0.5586f, 0.0f } } };
            static const VECTOR_F32 scale2 = { { { 0.0722f, 0.4351f, -0.0564f, 0.0f } } };

            MATRIX M(scale0, scale1, scale2, g_Zero);
            VECTOR clr = Vector3::Transform(rgb, M);

            return Vector::Select(rgb, clr, g_Select1110);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV YUVToRGB_HD(A_VECTOR yuv) noexcept
        {
            static const VECTOR_F32 scale1 = { { { 0.0f, -0.2153f, 2.1324f, 0.0f } } };
            static const VECTOR_F32 scale2 = { { { 1.2803f, -0.3806f, 0.0f, 0.0f } } };

            MATRIX M(g_One, scale1, scale2, g_Zero);
            VECTOR clr = Vector3::Transform(yuv, M);

            return Vector::Select(yuv, clr, g_Select1110);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV RGBToYUV_UHD(A_VECTOR rgb) noexcept
        {
            static const VECTOR_F32 scale0 = { { { 0.2627f, -0.1215f,  0.6150f, 0.0f } } };
            static const VECTOR_F32 scale1 = { { { 0.6780f, -0.3136f, -0.5655f, 0.0f } } };
            static const VECTOR_F32 scale2 = { { { 0.0593f,  0.4351f, -0.0495f, 0.0f } } };

            MATRIX M(scale0, scale1, scale2, g_Zero);
            VECTOR clr = Vector3::Transform(rgb, M);

            return Vector::Select(rgb, clr, g_Select1110);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV YUVToRGB_UHD(A_VECTOR yuv) noexcept
        {
            static const VECTOR_F32 scale1 = { { {    0.0f, -0.1891f, 2.1620f, 0.0f } } };
            static const VECTOR_F32 scale2 = { { { 1.1989f, -0.4645f,    0.0f, 0.0f } } };

            MATRIX M(g_One, scale1, scale2, g_Zero);
            VECTOR clr = Vector3::Transform(yuv, M);

            return Vector::Select(yuv, clr, g_Select1110);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV RGBToXYZ(A_VECTOR rgb) noexcept
        {
            static const VECTOR_F32 scale0 = { { { 0.4887180f, 0.1762044f, 0.0000000f, 0.0f } } };
            static const VECTOR_F32 scale1 = { { { 0.3106803f, 0.8129847f, 0.0102048f, 0.0f } } };
            static const VECTOR_F32 scale2 = { { { 0.2006017f, 0.0108109f, 0.9897952f, 0.0f } } };
            static const VECTOR_F32 scale = { { { 1.f / 0.17697f, 1.f / 0.17697f, 1.f / 0.17697f, 0.0f } } };

            MATRIX M(scale0, scale1, scale2, g_Zero);
            VECTOR clr = Vector::Multiply(Vector3::Transform(rgb, M), scale);

            return Vector::Select(rgb, clr, g_Select1110);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV XYZToRGB(A_VECTOR xyz) noexcept
        {
            static const VECTOR_F32 scale0 = { { { 2.3706743f, -0.5138850f, 0.0052982f, 0.0f } } };
            static const VECTOR_F32 scale1 = { { { -0.9000405f, 1.4253036f, -0.0146949f, 0.0f } } };
            static const VECTOR_F32 scale2 = { { { -0.4706338f, 0.0885814f, 1.0093968f, 0.0f } } };
            static const VECTOR_F32 scale = { { { 0.17697f, 0.17697f, 0.17697f, 0.0f } } };

            MATRIX M(scale0, scale1, scale2, g_Zero);
            VECTOR clr = Vector3::Transform(Vector::Multiply(xyz, scale), M);

            return Vector::Select(xyz, clr, g_Select1110);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV XYZToSRGB(A_VECTOR xyz) noexcept
        {
            static const VECTOR_F32 scale0 = { { { 3.2406f, -0.9689f, 0.0557f, 0.0f } } };
            static const VECTOR_F32 scale1 = { { { -1.5372f, 1.8758f, -0.2040f, 0.0f } } };
            static const VECTOR_F32 scale2 = { { { -0.4986f, 0.0415f, 1.0570f, 0.0f } } };
            static const VECTOR_F32 cutoff = { { { 0.0031308f, 0.0031308f, 0.0031308f, 0.0f } } };
            static const VECTOR_F32 exp = { { { 1.0f / 2.4f, 1.0f / 2.4f, 1.0f / 2.4f, 1.0f } } };

            MATRIX M(scale0, scale1, scale2, g_Zero);
            VECTOR lclr = Vector3::Transform(xyz, M);

            VECTOR sel = Vector::Greater(lclr, cutoff);

            // clr = 12.92 * lclr for lclr <= 0.0031308f
            VECTOR smallC = Vector::Multiply(lclr, g_srgbScale);

            // clr = (1+a)*pow(lclr, 1/2.4) - a for lclr > 0.0031308 (where a = 0.055)
            VECTOR largeC = Vector::Subtract(Vector::Multiply(g_srgbA1, Vector::Pow(lclr, exp)), g_srgbA);

            VECTOR clr = Vector::Select(smallC, largeC, sel);

            return Vector::Select(xyz, clr, g_Select1110);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV SRGBToXYZ(A_VECTOR srgb) noexcept
        {
            static const VECTOR_F32 scale0 = { { { 0.4124f, 0.2126f, 0.0193f, 0.0f } } };
            static const VECTOR_F32 scale1 = { { { 0.3576f, 0.7152f, 0.1192f, 0.0f } } };
            static const VECTOR_F32 scale2 = { { { 0.1805f, 0.0722f, 0.9505f, 0.0f } } };
            static const VECTOR_F32 cutoff = { { { 0.04045f, 0.04045f, 0.04045f, 0.0f } } };
            static const VECTOR_F32 exp = { { { 2.4f, 2.4f, 2.4f, 1.0f } } };

            VECTOR sel = Vector::Greater(srgb, cutoff);

            // lclr = clr / 12.92
            VECTOR smallC = Vector::Divide(srgb, g_srgbScale);

            // lclr = pow( (clr + a) / (1+a), 2.4 )
            VECTOR largeC = Vector::Pow(Vector::Divide(Vector::Add(srgb, g_srgbA), g_srgbA1), exp);

            VECTOR lclr = Vector::Select(smallC, largeC, sel);

            MATRIX M(scale0, scale1, scale2, g_Zero);
            VECTOR clr = Vector3::Transform(lclr, M);

            return Vector::Select(srgb, clr, g_Select1110); 
        }

        FORCE_INLINE VECTOR VEC_CALLCONV RGBToSRGB(A_VECTOR rgb) noexcept
        {
            static const VECTOR_F32 cutoff = { { { 0.0031308f, 0.0031308f, 0.0031308f, 1.f } } };
            static const VECTOR_F32 linear = { { { 12.92f, 12.92f, 12.92f, 1.f } } };
            static const VECTOR_F32 scale = { { { 1.055f, 1.055f, 1.055f, 1.f } } };
            static const VECTOR_F32 bias = { { { 0.055f, 0.055f, 0.055f, 0.f } } };
            static const VECTOR_F32 invGamma = { { { 1.0f / 2.4f, 1.0f / 2.4f, 1.0f / 2.4f, 1.f } } };

            VECTOR V = Vector::Saturate(rgb);
            VECTOR V0 = Vector::Multiply(V, linear);
            VECTOR V1 = Vector::Subtract(Vector::Multiply(scale, Vector::Pow(V, invGamma)), bias);
            VECTOR select = Vector::Less(V, cutoff);
            V = Vector::Select(V1, V0, select);
            
            return Vector::Select(rgb, V, g_Select1110);
        }

        FORCE_INLINE VECTOR VEC_CALLCONV SRGBToRGB(A_VECTOR srgb) noexcept
        {
            static const VECTOR_F32 cutoff = { { { 0.04045f, 0.04045f, 0.04045f, 1.f } } };
            static const VECTOR_F32 iLinear = { { { 1.f / 12.92f, 1.f / 12.92f, 1.f / 12.92f, 1.f } } };
            static const VECTOR_F32 scale = { { { 1.f / 1.055f, 1.f / 1.055f, 1.f / 1.055f, 1.f } } };
            static const VECTOR_F32 bias = { { { 0.055f, 0.055f, 0.055f, 0.f } } };
            static const VECTOR_F32 gamma = { { { 2.4f, 2.4f, 2.4f, 1.f } } };

            VECTOR V = Vector::Saturate(srgb);
            VECTOR V0 = Vector::Multiply(V, iLinear);
            VECTOR V1 = Vector::Pow(Vector::Multiply(Vector::Add(V, bias), scale), gamma);
            VECTOR select = Vector::Greater(V, cutoff);
            V = Vector::Select(V0, V1, select);
            
            return Vector::Select(srgb, V, g_Select1110);
        }
    } // namespace Color
    
} // namespace UltReality::Math


#endif // !ULTREALITY_MATH_COLOR_INL