#ifndef ULTREALITY_MATH_COLOR_H
#define ULTREALITY_MATH_COLOR_H

#include <SSE2VectorConfig.h>
#include <Float3.h>
#include <Float4.h>

namespace UltReality::Math
{
    namespace Color
    {
        bool VEC_CALLCONV Equal(A_VECTOR C1, A_VECTOR C2) noexcept;
        bool VEC_CALLCONV NotEqual(A_VECTOR C1, A_VECTOR C2) noexcept;
        bool VEC_CALLCONV Greater(A_VECTOR C1, A_VECTOR C2) noexcept;
        bool VEC_CALLCONV GreaterOrEqual(A_VECTOR C1, A_VECTOR C2) noexcept;
        bool VEC_CALLCONV Less(A_VECTOR C1, A_VECTOR C2) noexcept;
        bool VEC_CALLCONV LessOrEqual(A_VECTOR C1, A_VECTOR C2) noexcept;

        bool VEC_CALLCONV IsNaN(A_VECTOR c) noexcept;
        bool VEC_CALLCONV IsInfinite(A_VECTOR c) noexcept;

        VECTOR VEC_CALLCONV Negative(A_VECTOR c) noexcept;
        VECTOR VEC_CALLCONV Modulate(A_VECTOR C1, A_VECTOR C2) noexcept;
        VECTOR VEC_CALLCONV AdjustSaturation(A_VECTOR c, float saturation) noexcept;
        VECTOR VEC_CALLCONV AdjustContrast(A_VECTOR c, float contrast) noexcept;

        VECTOR VEC_CALLCONV RGBToHSL(A_VECTOR rgb) noexcept;
        VECTOR VEC_CALLCONV HSVToRGB(A_VECTOR hsl) noexcept;

        VECTOR VEC_CALLCONV RGBToHSV(A_VECTOR rgb) noexcept;
        VECTOR VEC_CALLCONV HSVToRGB(A_VECTOR hsv) noexcept;

        VECTOR VEC_CALLCONV RGBToYUV(A_VECTOR rgb) noexcept;
        VECTOR VEC_CALLCONV YUVToRGB(A_VECTOR yuv) noexcept;

        VECTOR VEC_CALLCONV RGBToYUV_HD(A_VECTOR rgb) noexcept;
        VECTOR VEC_CALLCONV YUVToRGB_HD(A_VECTOR yuv) noexcept;

        VECTOR VEC_CALLCONV RGBToYUV_UHD(A_VECTOR rgb) noexcept;
        VECTOR VEC_CALLCONV YUVToRGB_UHD(A_VECTOR yuv) noexcept;

        VECTOR VEC_CALLCONV RGBToXYZ(A_VECTOR rgb) noexcept;
        VECTOR VEC_CALLCONV XYZToRGB(A_VECTOR xyz) noexcept;

        VECTOR VEC_CALLCONV XYZToSRGB(A_VECTOR xyz) noexcept;
        VECTOR VEC_CALLCONV SRGBToXYZ(A_VECTOR srgb) noexcept;

        VECTOR VEC_CALLCONV RGBToSRGB(A_VECTOR rgb) noexcept;
        VECTOR VEC_CALLCONV SRGBToRGB(A_VECTOR srgb) noexcept;
    } // namespace Color

    namespace Colors
    {
        // Standard colors (Red/Green/Blue/Alpha) in sRGB colorspace
        
        VEC_GLOBCONST VECTOR_F32 AliceBlue = { { { 0.941176534f, 0.972549081f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 AntiqueWhite = { { { 0.980392218f, 0.921568692f, 0.843137324f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Aqua = { { { 0.f, 1.f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Aquamarine = { { { 0.498039246f, 1.f, 0.831372619f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Azure = { { { 0.941176534f, 1.f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Beige = { { { 0.960784376f, 0.960784376f, 0.862745166f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Bisque = { { { 1.f, 0.894117713f, 0.768627524f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Black = { { { 0.f, 0.f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 BlanchedAlmond = { { { 1.f, 0.921568692f, 0.803921640f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Blue = { { { 0.f, 0.f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 BlueViolet = { { { 0.541176498f, 0.168627456f, 0.886274576f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Brown = { { { 0.647058845f, 0.164705887f, 0.164705887f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 BurlyWood = { { { 0.870588303f, 0.721568644f, 0.529411793f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 CadetBlue = { { { 0.372549027f, 0.619607866f, 0.627451003f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Chartreuse = { { { 0.498039246f, 1.f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Chocolate = { { { 0.823529482f, 0.411764741f, 0.117647067f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Coral = { { { 1.f, 0.498039246f, 0.313725501f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 CornflowerBlue = { { { 0.392156899f, 0.584313750f, 0.929411829f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Cornsilk = { { { 1.f, 0.972549081f, 0.862745166f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Crimson = { { { 0.862745166f, 0.078431375f, 0.235294133f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Cyan = { { { 0.f, 1.f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkBlue = { { { 0.f, 0.f, 0.545098066f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkCyan = { { { 0.f, 0.545098066f, 0.545098066f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkGoldenrod = { { { 0.721568644f, 0.525490224f, 0.043137256f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkGray = { { { 0.662745118f, 0.662745118f, 0.662745118f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkGreen = { { { 0.f, 0.392156899f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkKhaki = { { { 0.741176486f, 0.717647076f, 0.419607878f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkMagenta = { { { 0.545098066f, 0.f, 0.545098066f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkOliveGreen = { { { 0.333333343f, 0.419607878f, 0.184313729f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkOrange = { { { 1.f, 0.549019635f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkOrchid = { { { 0.600000024f, 0.196078449f, 0.800000072f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkRed = { { { 0.545098066f, 0.f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkSalmon = { { { 0.913725555f, 0.588235319f, 0.478431404f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkSeaGreen = { { { 0.560784340f, 0.737254918f, 0.545098066f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkSlateBlue = { { { 0.282352954f, 0.239215702f, 0.545098066f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkSlateGray = { { { 0.184313729f, 0.309803933f, 0.309803933f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkTurquoise = { { { 0.f, 0.807843208f, 0.819607913f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkViolet = { { { 0.580392182f, 0.f, 0.827451050f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DeepPink = { { { 1.f, 0.078431375f, 0.576470613f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DeepSkyBlue = { { { 0.f, 0.749019623f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DimGray = { { { 0.411764741f, 0.411764741f, 0.411764741f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DodgerBlue = { { { 0.117647067f, 0.564705908f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Firebrick = { { { 0.698039234f, 0.133333340f, 0.133333340f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 FloralWhite = { { { 1.f, 0.980392218f, 0.941176534f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 ForestGreen = { { { 0.133333340f, 0.545098066f, 0.133333340f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Fuchsia = { { { 1.f, 0.f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Gainsboro = { { { 0.862745166f, 0.862745166f, 0.862745166f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 GhostWhite = { { { 0.972549081f, 0.972549081f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Gold = { { { 1.f, 0.843137324f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Goldenrod = { { { 0.854902029f, 0.647058845f, 0.125490203f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Gray = { { { 0.501960814f, 0.501960814f, 0.501960814f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Green = { { { 0.f, 0.501960814f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 GreenYellow = { { { 0.678431392f, 1.f, 0.184313729f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Honeydew = { { { 0.941176534f, 1.f, 0.941176534f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 HotPink = { { { 1.f, 0.411764741f, 0.705882370f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 IndianRed = { { { 0.803921640f, 0.360784322f, 0.360784322f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Indigo = { { { 0.294117659f, 0.f, 0.509803951f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Ivory = { { { 1.f, 1.f, 0.941176534f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Khaki = { { { 0.941176534f, 0.901960850f, 0.549019635f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Lavender = { { { 0.901960850f, 0.901960850f, 0.980392218f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LavenderBlush = { { { 1.f, 0.941176534f, 0.960784376f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LawnGreen = { { { 0.486274540f, 0.988235354f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LemonChiffon = { { { 1.f, 0.980392218f, 0.803921640f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightBlue = { { { 0.678431392f, 0.847058892f, 0.901960850f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightCoral = { { { 0.941176534f, 0.501960814f, 0.501960814f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightCyan = { { { 0.878431439f, 1.f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightGoldenrodYellow = { { { 0.980392218f, 0.980392218f, 0.823529482f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightGray = { { { 0.827451050f, 0.827451050f, 0.827451050f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightGreen = { { { 0.564705908f, 0.933333397f, 0.564705908f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightPink = { { { 1.f, 0.713725507f, 0.756862819f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightSalmon = { { { 1.f, 0.627451003f, 0.478431404f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightSeaGreen = { { { 0.125490203f, 0.698039234f, 0.666666687f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightSkyBlue = { { { 0.529411793f, 0.807843208f, 0.980392218f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightSlateGray = { { { 0.466666698f, 0.533333361f, 0.600000024f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightSteelBlue = { { { 0.690196097f, 0.768627524f, 0.870588303f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightYellow = { { { 1.f, 1.f, 0.878431439f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Lime = { { { 0.f, 1.f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LimeGreen = { { { 0.196078449f, 0.803921640f, 0.196078449f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Linen = { { { 0.980392218f, 0.941176534f, 0.901960850f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Magenta = { { { 1.f, 0.f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Maroon = { { { 0.501960814f, 0.f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MediumAquamarine = { { { 0.400000036f, 0.803921640f, 0.666666687f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MediumBlue = { { { 0.f, 0.f, 0.803921640f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MediumOrchid = { { { 0.729411781f, 0.333333343f, 0.827451050f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MediumPurple = { { { 0.576470613f, 0.439215720f, 0.858823597f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MediumSeaGreen = { { { 0.235294133f, 0.701960802f, 0.443137288f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MediumSlateBlue = { { { 0.482352972f, 0.407843173f, 0.933333397f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MediumSpringGreen = { { { 0.f, 0.980392218f, 0.603921592f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MediumTurquoise = { { { 0.282352954f, 0.819607913f, 0.800000072f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MediumVioletRed = { { { 0.780392230f, 0.082352944f, 0.521568656f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MidnightBlue = { { { 0.098039225f, 0.098039225f, 0.439215720f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MintCream = { { { 0.960784376f, 1.f, 0.980392218f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MistyRose = { { { 1.f, 0.894117713f, 0.882353008f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Moccasin = { { { 1.f, 0.894117713f, 0.709803939f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 NavajoWhite = { { { 1.f, 0.870588303f, 0.678431392f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Navy = { { { 0.f, 0.f, 0.501960814f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 OldLace = { { { 0.992156923f, 0.960784376f, 0.901960850f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Olive = { { { 0.501960814f, 0.501960814f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 OliveDrab = { { { 0.419607878f, 0.556862772f, 0.137254909f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Orange = { { { 1.f, 0.647058845f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 OrangeRed = { { { 1.f, 0.270588249f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Orchid = { { { 0.854902029f, 0.439215720f, 0.839215755f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 PaleGoldenrod = { { { 0.933333397f, 0.909803987f, 0.666666687f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 PaleGreen = { { { 0.596078455f, 0.984313786f, 0.596078455f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 PaleTurquoise = { { { 0.686274529f, 0.933333397f, 0.933333397f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 PaleVioletRed = { { { 0.858823597f, 0.439215720f, 0.576470613f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 PapayaWhip = { { { 1.f, 0.937254965f, 0.835294187f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 PeachPuff = { { { 1.f, 0.854902029f, 0.725490212f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Peru = { { { 0.803921640f, 0.521568656f, 0.247058839f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Pink = { { { 1.f, 0.752941251f, 0.796078503f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Plum = { { { 0.866666734f, 0.627451003f, 0.866666734f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 PowderBlue = { { { 0.690196097f, 0.878431439f, 0.901960850f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Purple = { { { 0.501960814f, 0.f, 0.501960814f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Red = { { { 1.f, 0.f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 RosyBrown = { { { 0.737254918f, 0.560784340f, 0.560784340f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 RoyalBlue = { { { 0.254901975f, 0.411764741f, 0.882353008f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 SaddleBrown = { { { 0.545098066f, 0.270588249f, 0.074509807f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Salmon = { { { 0.980392218f, 0.501960814f, 0.447058856f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 SandyBrown = { { { 0.956862807f, 0.643137276f, 0.376470625f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 SeaGreen = { { { 0.180392161f, 0.545098066f, 0.341176480f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 SeaShell = { { { 1.f, 0.960784376f, 0.933333397f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Sienna = { { { 0.627451003f, 0.321568638f, 0.176470593f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Silver = { { { 0.752941251f, 0.752941251f, 0.752941251f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 SkyBlue = { { { 0.529411793f, 0.807843208f, 0.921568692f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 SlateBlue = { { { 0.415686309f, 0.352941185f, 0.803921640f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 SlateGray = { { { 0.439215720f, 0.501960814f, 0.564705908f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Snow = { { { 1.f, 0.980392218f, 0.980392218f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 SpringGreen = { { { 0.f, 1.f, 0.498039246f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 SteelBlue = { { { 0.274509817f, 0.509803951f, 0.705882370f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Tan = { { { 0.823529482f, 0.705882370f, 0.549019635f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Teal = { { { 0.f, 0.501960814f, 0.501960814f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Thistle = { { { 0.847058892f, 0.749019623f, 0.847058892f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Tomato = { { { 1.f, 0.388235331f, 0.278431386f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Transparent = { { { 0.f, 0.f, 0.f, 0.f } } };
        VEC_GLOBCONST VECTOR_F32 Turquoise = { { { 0.250980407f, 0.878431439f, 0.815686345f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Violet = { { { 0.933333397f, 0.509803951f, 0.933333397f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Wheat = { { { 0.960784376f, 0.870588303f, 0.701960802f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 White = { { { 1.f, 1.f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 WhiteSmoke = { { { 0.960784376f, 0.960784376f, 0.960784376f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Yellow = { { { 1.f, 1.f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 YellowGreen = { { { 0.603921592f, 0.803921640f, 0.196078449f, 1.f } } };
    } // namespace Colors

    namespace ColorsLinear
    {
        // Standard colors (Red/Green/Blue/Alpha) in linear colorspace
        
        VEC_GLOBCONST VECTOR_F32 AliceBlue = { { { 0.871367335f, 0.938685894f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 AntiqueWhite = { { { 0.955973506f, 0.830770075f, 0.679542601f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Aqua = { { { 0.f, 1.f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Aquamarine = { { { 0.212230787f, 1.f, 0.658374965f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Azure = { { { 0.871367335f, 1.f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Beige = { { { 0.913098991f, 0.913098991f, 0.715693772f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Bisque = { { { 1.f, 0.775822461f, 0.552011609f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Black = { { { 0.f, 0.f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 BlanchedAlmond = { { { 1.f, 0.830770075f, 0.610495746f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Blue = { { { 0.f, 0.f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 BlueViolet = { { { 0.254152179f, 0.024157630f, 0.760524750f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Brown = { { { 0.376262218f, 0.023153365f, 0.023153365f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 BurlyWood = { { { 0.730461001f, 0.479320228f, 0.242281199f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 CadetBlue = { { { 0.114435382f, 0.341914445f, 0.351532698f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Chartreuse = { { { 0.212230787f, 1.f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Chocolate = { { { 0.644479871f, 0.141263321f, 0.012983031f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Coral = { { { 1.f, 0.212230787f, 0.080219828f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 CornflowerBlue = { { { 0.127437726f, 0.300543845f, 0.846873462f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Cornsilk = { { { 1.f, 0.938685894f, 0.715693772f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Crimson = { { { 0.715693772f, 0.006995410f, 0.045186214f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Cyan = { { { 0.f, 1.f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkBlue = { { { 0.f, 0.f, 0.258182913f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkCyan = { { { 0.f, 0.258182913f, 0.258182913f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkGoldenrod = { { { 0.479320228f, 0.238397658f, 0.003346536f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkGray = { { { 0.396755308f, 0.396755308f, 0.396755308f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkGreen = { { { 0.f, 0.127437726f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkKhaki = { { { 0.508881450f, 0.473531544f, 0.147027299f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkMagenta = { { { 0.258182913f, 0.f, 0.258182913f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkOliveGreen = { { { 0.090841733f, 0.147027299f, 0.028426038f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkOrange = { { { 1.f, 0.262250721f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkOrchid = { { { 0.318546832f, 0.031896040f, 0.603827536f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkRed = { { { 0.258182913f, 0.f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkSalmon = { { { 0.814846814f, 0.304987371f, 0.194617867f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkSeaGreen = { { { 0.274677366f, 0.502886593f, 0.258182913f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkSlateBlue = { { { 0.064803280f, 0.046665095f, 0.258182913f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkSlateGray = { { { 0.028426038f, 0.078187428f, 0.078187428f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkTurquoise = { { { 0.f, 0.617206752f, 0.637597024f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DarkViolet = { { { 0.296138316f, 0.f, 0.651405811f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DeepPink = { { { 1.f, 0.006995410f, 0.291770697f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DeepSkyBlue = { { { 0.f, 0.520995677f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DimGray = { { { 0.141263321f, 0.141263321f, 0.141263321f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 DodgerBlue = { { { 0.012983031f, 0.278894335f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Firebrick = { { { 0.445201248f, 0.015996292f, 0.015996292f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 FloralWhite = { { { 1.f, 0.955973506f, 0.871367335f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 ForestGreen = { { { 0.015996292f, 0.258182913f, 0.015996292f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Fuchsia = { { { 1.f, 0.f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Gainsboro = { { { 0.715693772f, 0.715693772f, 0.715693772f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 GhostWhite = { { { 0.938685894f, 0.938685894f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Gold = { { { 1.f, 0.679542601f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Goldenrod = { { { 0.701102138f, 0.376262218f, 0.014443844f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Gray = { { { 0.215860531f, 0.215860531f, 0.215860531f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Green = { { { 0.f, 0.215860531f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 GreenYellow = { { { 0.417885154f, 1.f, 0.028426038f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Honeydew = { { { 0.871367335f, 1.f, 0.871367335f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 HotPink = { { { 1.f, 0.141263321f, 0.456411064f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 IndianRed = { { { 0.610495746f, 0.107023112f, 0.107023112f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Indigo = { { { 0.070360109f, 0.f, 0.223227978f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Ivory = { { { 1.f, 1.f, 0.871367335f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Khaki = { { { 0.871367335f, 0.791298151f, 0.262250721f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Lavender = { { { 0.791298151f, 0.791298151f, 0.955973506f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LavenderBlush = { { { 1.f, 0.871367335f, 0.913098991f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LawnGreen = { { { 0.201556295f, 0.973445475f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LemonChiffon = { { { 1.f, 0.955973506f, 0.610495746f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightBlue = { { { 0.417885154f, 0.686685443f, 0.791298151f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightCoral = { { { 0.871367335f, 0.215860531f, 0.215860531f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightCyan = { { { 0.745404482f, 1.f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightGoldenrodYellow = { { { 0.955973506f, 0.955973506f, 0.644479871f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightGray = { { { 0.651405811f, 0.651405811f, 0.651405811f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightGreen = { { { 0.278894335f, 0.854992807f, 0.278894335f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightPink = { { { 1.f, 0.467783839f, 0.533276618f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightSalmon = { { { 1.f, 0.351532698f, 0.194617867f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightSeaGreen = { { { 0.014443844f, 0.445201248f, 0.401977867f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightSkyBlue = { { { 0.242281199f, 0.617206752f, 0.955973506f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightSlateGray = { { { 0.184475034f, 0.246201396f, 0.318546832f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightSteelBlue = { { { 0.434153706f, 0.552011609f, 0.730461001f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LightYellow = { { { 1.f, 1.f, 0.745404482f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Lime = { { { 0.f, 1.f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 LimeGreen = { { { 0.031896040f, 0.610495746f, 0.031896040f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Linen = { { { 0.955973506f, 0.871367335f, 0.791298151f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Magenta = { { { 1.f, 0.f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Maroon = { { { 0.215860531f, 0.f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MediumAquamarine = { { { 0.132868364f, 0.610495746f, 0.401977867f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MediumBlue = { { { 0.f, 0.f, 0.610495746f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MediumOrchid = { { { 0.491020888f, 0.090841733f, 0.651405811f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MediumPurple = { { { 0.291770697f, 0.162029430f, 0.708376050f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MediumSeaGreen = { { { 0.045186214f, 0.450785846f, 0.165132239f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MediumSlateBlue = { { { 0.198069349f, 0.138431653f, 0.854992807f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MediumSpringGreen = { { { 0.f, 0.955973506f, 0.323143244f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MediumTurquoise = { { { 0.064803280f, 0.637597024f, 0.603827536f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MediumVioletRed = { { { 0.571125031f, 0.007499032f, 0.234550655f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MidnightBlue = { { { 0.009721218f, 0.009721218f, 0.162029430f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MintCream = { { { 0.913098991f, 1.f, 0.955973506f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 MistyRose = { { { 1.f, 0.775822461f, 0.752942443f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Moccasin = { { { 1.f, 0.775822461f, 0.462077051f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 NavajoWhite = { { { 1.f, 0.730461001f, 0.417885154f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Navy = { { { 0.f, 0.f, 0.215860531f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 OldLace = { { { 0.982250869f, 0.913098991f, 0.791298151f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Olive = { { { 0.215860531f, 0.215860531f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 OliveDrab = { { { 0.147027299f, 0.270497859f, 0.016807375f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Orange = { { { 1.f, 0.376262218f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 OrangeRed = { { { 1.f, 0.059511241f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Orchid = { { { 0.701102138f, 0.162029430f, 0.672443330f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 PaleGoldenrod = { { { 0.854992807f, 0.806952477f, 0.401977867f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 PaleGreen = { { { 0.313988745f, 0.964686573f, 0.313988745f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 PaleTurquoise = { { { 0.428690553f, 0.854992807f, 0.854992807f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 PaleVioletRed = { { { 0.708376050f, 0.162029430f, 0.291770697f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 PapayaWhip = { { { 1.f, 0.863157392f, 0.665387452f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 PeachPuff = { { { 1.f, 0.701102138f, 0.485149980f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Peru = { { { 0.610495746f, 0.234550655f, 0.049706575f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Pink = { { { 1.f, 0.527115345f, 0.597202003f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Plum = { { { 0.723055363f, 0.351532698f, 0.723055363f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 PowderBlue = { { { 0.434153706f, 0.745404482f, 0.791298151f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Purple = { { { 0.215860531f, 0.f, 0.215860531f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Red = { { { 1.f, 0.f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 RosyBrown = { { { 0.502886593f, 0.274677366f, 0.274677366f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 RoyalBlue = { { { 0.052860655f, 0.141263321f, 0.752942443f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 SaddleBrown = { { { 0.258182913f, 0.059511241f, 0.006512091f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Salmon = { { { 0.955973506f, 0.215860531f, 0.168269455f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 SandyBrown = { { { 0.904661357f, 0.371237785f, 0.116970696f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 SeaGreen = { { { 0.027320892f, 0.258182913f, 0.095307484f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 SeaShell = { { { 1.f, 0.913098991f, 0.854992807f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Sienna = { { { 0.351532698f, 0.084376216f, 0.026241222f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Silver = { { { 0.527115345f, 0.527115345f, 0.527115345f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 SkyBlue = { { { 0.242281199f, 0.617206752f, 0.830770075f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 SlateBlue = { { { 0.144128501f, 0.102241747f, 0.610495746f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 SlateGray = { { { 0.162029430f, 0.215860531f, 0.278894335f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Snow = { { { 1.f, 0.955973506f, 0.955973506f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 SpringGreen = { { { 0.f, 1.f, 0.212230787f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 SteelBlue = { { { 0.061246071f, 0.223227978f, 0.456411064f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Tan = { { { 0.644479871f, 0.456411064f, 0.262250721f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Teal = { { { 0.f, 0.215860531f, 0.215860531f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Thistle = { { { 0.686685443f, 0.520995677f, 0.686685443f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Tomato = { { { 1.f, 0.124771863f, 0.063010029f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Transparent = { { { 0.f, 0.f, 0.f, 0.f } } };
        VEC_GLOBCONST VECTOR_F32 Turquoise = { { { 0.051269468f, 0.745404482f, 0.630757332f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Violet = { { { 0.854992807f, 0.223227978f, 0.854992807f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Wheat = { { { 0.913098991f, 0.730461001f, 0.450785846f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 White = { { { 1.f, 1.f, 1.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 WhiteSmoke = { { { 0.913098991f, 0.913098991f, 0.913098991f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 Yellow = { { { 1.f, 1.f, 0.f, 1.f } } };
        VEC_GLOBCONST VECTOR_F32 YellowGreen = { { { 0.323143244f, 0.610495746f, 0.031896040f, 1.f } } };
    } // namespace ColorsLinear
} // namespace UltReality::Math

#endif // !ULTREALITY_MATH_COLOR_H
