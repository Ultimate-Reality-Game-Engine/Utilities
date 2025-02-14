#ifndef ULTREALITY_MATH_PACKED_VECTOR_H
#define ULTREALITY_MATH_PACKED_VECTOR_H

#include <SIMDVectorConfig.h>

namespace UltReality::Math
{
    namespace PackedVector
    {
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4201 4365 4324 4996)
        // C4201: nonstandard extension used
        // C4365: Off by default noise
        // C4324: alignment padding warnings
        // C4996: deprecation warnings
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wgnu-anonymous-struct"
#pragma clang diagnostic ignored "-Wnested-anon-types"
#endif

        //------------------------------------------------------------------------------
        // ARGB Color; 8-8-8-8 bit unsigned normalized integer components packed into
        // a 32 bit integer.  The normalized color is packed into 32 bits using 8 bit
        // unsigned, normalized integers for the alpha, red, green, and blue components.
        // The alpha component is stored in the most significant bits and the blue
        // component in the least significant bits (A8R8G8B8):
        // [32] aaaaaaaa rrrrrrrr gggggggg bbbbbbbb [0]
        struct COLOR
        {
            union
            {
                struct
                {
                    uint8_t b; // Blue:     0/255 to 255/255
                    uint8_t g; // Green:    0/255 to 255/255
                    uint8_t r; // Red:      0/255 to 255/255
                    uint8_t a; // Alpha:    0/255 to 255/255
                };
                uint32_t c;
            };

            COLOR() = default;

            COLOR(const COLOR&) = default;
            COLOR& operator=(const COLOR&) = default;

            COLOR(COLOR&&) = default;
            COLOR& operator=(COLOR&&) = default;

            constexpr COLOR(uint32_t color) noexcept: c(color) {}
            COLOR(float _r, float _g, float _b, float _a) noexcept;
            explicit COLOR(_In_reads_(4) const float* pArray) noexcept;

            operator uint32_t() const noexcept { return c; }

            COLOR& operator=(const uint32_t color) noexcept { c = color; return *this; }
        };

        //------------------------------------------------------------------------------
        // 16 bit floating point number consisting of a sign bit, a 5 bit biased
        // exponent, and a 10 bit mantissa
        using HALF = uint16_t;

        //------------------------------------------------------------------------------
        // 2D Vector; 16 bit floating point components
        struct HALF2
        {
            union
            {
                struct
                {
                    HALF x;
                    HALF y;
                };
                uint32_t v;
            };
            
            HALF2() = default;

            HALF2(const HALF2&) = default;
            HALF2& operator=(const HALF2&) = default;

            HALF2(HALF2&&) = default;
            HALF2& operator=(HALF2&&) = default;

            explicit constexpr HALF2(uint32_t packed) noexcept : v(packed) {}
            constexpr HALF2(HALF _x, HALF _y) noexcept : x(_x), y(_y) {}
            explicit HALF2(_In_reads_(2) const HALF* pArray) noexcept : x(pArray[0]), y(pArray[1]) {}
            HALF2(float _x, float _y) noexcept;
            explicit HALF2(_In_reads_(2) const float* pArray) noexcept;

            HALF2& operator=(uint32_t packed) noexcept { v = packed; return *this; }
        };

        //------------------------------------------------------------------------------
        // 2D Vector; 16 bit signed normalized integer components
        struct SHORTN2
        {
            union
            {
                struct
                {
                    int16_t x;
                    int16_t y;
                };
                uint32_t v;
            };
            
            SHORTN2() = default;

            SHORTN2(const SHORTN2&) = default;
            SHORTN2& operator=(const SHORTN2&) = default;

            SHORTN2(SHORTN2&&) = default;
            SHORTN2& operator=(SHORTN2&&) = default;

            explicit constexpr SHORTN2(uint32_t packed) noexcept : v(packed) {}
            constexpr SHORTN2(int16_t _x, int16_t _y) noexcept : x(_x), y(_y) {}
            explicit SHORTN2(_In_reads_(2) const uint16_t* pArray) noexcept : x(pArray[0]), y(pArray[1]) {}
            SHORTN2(float _x, float _y) noexcept;
            explicit SHORTN2(_In_reads_(2) const float* pArray) noexcept;

            SHORTN2 operator=(uint32_t packed) noexcept { v = packed; return *this; }
        };

        // 2D Vector; 16 bit signed integer components
        struct SHORT2
        {
            union
            {
                struct
                {
                    int16_t x;
                    int16_t y;
                };
                uint32_t v;
            };
            
            SHORT2() = default;

            SHORT2(const SHORT2&) = default;
            SHORT2& operator=(const SHORT2&) = default;

            SHORT2(SHORT2&&) = default;
            SHORT2& operator=(SHORT2&&) = default;

            explicit constexpr SHORT2(uint32_t packed) noexcept : v(packed) {}
            constexpr SHORT2(int16_t _x, int16_t _y) noexcept : x(_x), y(_y) {}
            explicit SHORT2(_In_reads_(2) const int16_t* pArray) noexcept : x(pArray[0]), y(pArray[1]) {}
            SHORT2(float _x, float _y) noexcept;
            explicit SHORT2(_In_reads_(2) const float* pArray) noexcept;

            SHORT2& operator=(uint32_t packed) noexcept { v = packed; return *this; }
        };

        // 2D Vector; 16 bit unsinged normalized integer components
        struct USHORTN2
        {
            union
            {
                struct
                {
                    uint16_t x;
                    uint16_t y;
                };
                uint32_t v;
            };
            
            USHORTN2() = default;

            USHORTN2(const USHORTN2&) = default;
            USHORTN2& operator=(const USHORTN2&) = default;

            USHORTN2(USHORTN2&&) = default;
            USHORTN2& operator=(USHORTN2&&) = default;

            explicit constexpr USHORTN2(uint32_t packed) noexcept : v(packed) {}
            constexpr USHORTN2(uint16_t _x, uint16_t _y) noexcept : x(_x), y(_y) {}
            explicit USHORTN2(_In_reads_(2) const uint16_t* pArray) noexcept : x(pArray[0]), y(pArray[1]) {}
            USHORTN2(float _x, float _y) noexcept;
            explicit USHORTN2(_In_reads_(2) const float* pArray) noexcept;

            USHORTN2& operator=(uint32_t packed) noexcept { v = packed; return *this; }
        };

        // 2D Vector; 16 bit unsigned integer components
        struct USHORT2
        {
            union
            {
                struct
                {
                    uint16_t x;
                    uint16_t y;
                };
                uint32_t v;
            };
            
            USHORT2() = default;

            USHORT2(const USHORT2&) = default;
            USHORT2& operator=(const USHORT2&) = default;

            USHORT2(USHORT2&&) = default;
            USHORT2& operator=(USHORT2&&) = default;

            explicit constexpr USHORT2(uint32_t packed) noexcept : v(packed) {}
            constexpr USHORT2(uint16_t _x, uint16_t _y) noexcept : x(_x), y(_y) {}
            explicit USHORT2(_In_reads_(2) const uint16_t* pArray) noexcept : x(pArray[0]), y(pArray[1]) {}
            USHORT2(float _x, float _y) noexcept;
            explicit USHORT2(_In_reads_(2) const float* pArray) noexcept;

            USHORT2& operator=(uint32_t packed) noexcept { v = packed; return *this; }
        };

        //------------------------------------------------------------------------------
        // 2D Vector; 8 bit signed normalized integer components
        struct BYTEN2
        {
            union
            {
                struct
                {
                    int8_t x;
                    int8_t y;
                };
                uint16_t v;
            };
            
            BYTEN2() = default;

            BYTEN2(const BYTEN2&) = default;
            BYTEN2& operator=(const BYTEN2&) = default;

            BYTEN2(BYTEN2&&) = default;
            BYTEN2& operator=(BYTEN2&&) = default;

            explicit constexpr BYTEN2(uint16_t packed) noexcept : v(packed) {}
            constexpr BYTEN2(int8_t _x, int8_t _y) noexcept : x(_x), y(_y) {}
            explicit BYTEN2(_In_reads_(2) const int8_t* pArray) noexcept : x(pArray[0]), y(pArray[1]) {}
            BYTEN2(float _x, float _y) noexcept;
            explicit BYTEN2(_In_reads_(2) const float * pArray) noexcept;

            BYTEN2& operator=(uint16_t packed) noexcept { v = packed; return *this; }
        };

        // 2D Vector; 8 bit signed integer components
        struct BYTE2
        {
            union
            {
                struct
                {
                    int8_t x;
                    int8_t y;
                };
                uint16_t v;
            };
            
            BYTE2() = default;

            BYTE2(const BYTE2&) = default;
            BYTE2& operator=(const BYTE2&) = default;

            BYTE2(BYTE2&&) = default;
            BYTE2& operator=(BYTE2&&) = default;

            explicit constexpr BYTE2(uint16_t packed) noexcept : v(packed) {}
            constexpr BYTE2(int8_t _x, int8_t _y) noexcept : x(_x), y(_y) {}
            explicit BYTE2(_In_reads_(2) const int8_t* pArray) noexcept : x(pArray[0]), y(pArray[1]) {}
            BYTE2(float _x, float _y) noexcept;
            explicit BYTE2(_In_reads_(2) const float* pArray) noexcept;

            BYTE2& operator=(uint16_t packed) noexcept { v = packed; return *this; }
        };
        
        // 2D Vector; 8 bit unsigned normalized integer components
        struct UBYTEN2
        {
            union
            {
                struct
                {
                    uint8_t x;
                    uint8_t y;
                };
                uint16_t v;
            };
            
            UBYTEN2() = default;

            UBYTEN2(const UBYTEN2&) = default;
            UBYTEN2& operator=(const UBYTEN2&) = default;

            UBYTEN2(UBYTEN2&&) = default;
            UBYTEN2& operator=(UBYTEN2&&) = default;

            explicit constexpr UBYTEN2(uint16_t packed) noexcept : v(packed) {}
            constexpr UBYTEN2(uint8_t _x, uint8_t _y) noexcept : x(_x), y(_y) {}
            explicit UBYTEN2(_In_reads_(2) const uint8_t* pArray) noexcept : x(pArray[0]), y(pArray[1]) {}
            UBYTEN2(float _x, float _y) noexcept;
            explicit UBYTEN2(_In_reads_(2) const float* pArray) noexcept;

            UBYTEN2& operator=(uint16_t packed) noexcept { v = packed; return *this; }
        };
        
        // 2D Vector; 8 bit unsinged integer components
        struct UBYTE2
        {
            union
            {
                struct
                {
                    uint8_t x;
                    uint8_t y;
                };
                uint16_t v;
            };
            
            UBYTE2() = default;

            UBYTE2(UBYTE2&) = default;
            UBYTE2& operator=(const UBYTE2&) = default;

            UBYTE2(UBYTE2&&) = default;
            UBYTE2& operator=(UBYTE2&&) = default;

            explicit constexpr UBYTE2(uint16_t packed) noexcept : v(packed) {}
            constexpr UBYTE2(uint8_t _x, uint8_t _y) noexcept : x(_x), y(_y) {}
            explicit UBYTE2(_In_reads_(2) const uint8_t* pArray) noexcept : x(pArray[0]), y(pArray[1]) {}
            UBYTE2(float _x, float _y) noexcept;
            explicit UBYTE2(_In_reads_(2) const float* pArray) noexcept;

            UBYTE2& operator=(uint16_t packed) noexcept { v = packed; return *this; }
        };
        
        //------------------------------------------------------------------------------
        // 3D vector: 5/6/5 unsigned integer components
        struct U565
        {
            union
            {
                struct
                {
                    uint16_t x : 5; // 0 to 31
                    uint16_t y : 6; // 0 to 63
                    uint16_t z : 5; // 0 to 31
                };
                uint16_t v;
            };
            
            U565() = default;

            U565(const U565&) = default;
            U565& operator=(const U565&) = default;

            U565(U565&&) = default;
            U565& operator=(U565&&) = default;

            explicit constexpr U565(uint16_t packed) noexcept : v(packed) {}
            constexpr U565(uint8_t _x, uint8_t _y, uint8_t _z) noexcept : x(_x), y(_y), z(_z) {}
            explicit U565(_In_reads_(3) const uint8_t* pArray) noexcept : x(pArray[0]), y(pArray[1]), z(pArray[2]) {}
            U565(float _x, float _y, float _z) noexcept;
            explicit U565(_In_reads_(3) const float* pArray) noexcept;

            operator uint16_t() const noexcept { return v; }

            U565& operator=(uint16_t packed) noexcept { v = packed; return *this; }
        };

        //------------------------------------------------------------------------------
        // 3D vector: 11/11/10 floating-point components
        // The 3D vector is packed into 32 bits as follows: a 5-bit biased exponent
        // and 6-bit mantissa for x component, a 5-bit biased exponent and
        // 6-bit mantissa for y component, a 5-bit biased exponent and a 5-bit
        // mantissa for z. The z component is stored in the most significant bits
        // and the x component in the least significant bits. No sign bits so
        // all partial-precision numbers are positive.
        // (Z10Y11X11): [32] ZZZZZzzz zzzYYYYY yyyyyyXX XXXxxxxx [0]
        struct FLOAT3PK
        {
            union
            {
                struct
                {
                    uint32_t xm : 6; // x-mantissa
                    uint32_t xe : 5; // x-exponent
                    uint32_t ym : 6; // y-mantissa
                    uint32_t ye : 5; // y-exponent
                    uint32_t zm : 5; // z-mantissa
                    uint32_t ze : 5; // z-exponent
                };
                uint32_t v;
            };
            
            FLOAT3PK() = default;

            FLOAT3PK(const FLOAT3PK&) = default;
            FLOAT3PK& operator=(const FLOAT3PK&) = default;

            FLOAT3PK(FLOAT3PK&&) = default;
            FLOAT3PK& operator=(FLOAT3PK&&) = default;

            explicit constexpr FLOAT3PK(uint32_t packed) noexcept : v(packed) {}
            FLOAT3PK(float _x, float _y, float _z) noexcept;
            explicit FLOAT3PK(_In_reads_(3) const float* pArray) noexcept;

            operator uint32_t () const noexcept { return v; }

            FLOAT3PK& operator=(uint32_t packed) noexcept { v = packed; return *this; }
        };
        
        //------------------------------------------------------------------------------
        // 3D vector: 9/9/9 floating-point components with shared 5-bit exponent
        // The 3D vector is packed into 32 bits as follows: a 5-bit biased exponent
        // with 9-bit mantissa for the x, y, and z component. The shared exponent
        // is stored in the most significant bits and the x component mantissa is in
        // the least significant bits. No sign bits so all partial-precision numbers
        // are positive.
        // (E5Z9Y9X9): [32] EEEEEzzz zzzzzzyy yyyyyyyx xxxxxxxx [0]
        struct FLOAT3SE
        {
            union
            {
                struct
                {
                    uint32_t xm : 9; // x-mantissa
                    uint32_t ym : 9; // y-mantissa
                    uint32_t zm : 9; // z-mantissa
                    uint32_t e : 5; // shared exponent
                };
                uint32_t v;
            };

            FLOAT3SE() = default;

            FLOAT3SE(const FLOAT3SE&) = default;
            FLOAT3SE& operator=(const FLOAT3SE&) = default;

            FLOAT3SE(FLOAT3SE&&) = default;
            FLOAT3SE& operator=(FLOAT3SE&&) = default;

            explicit constexpr FLOAT3SE(uint32_t packed) noexcept : v(packed) {}
            FLOAT3SE(float _x, float _y, float _z) noexcept;
            explicit FLOAT3SE(_In_reads_(3) const float* pArray) noexcept;

            operator uint32_t () const noexcept { return v; }

            FLOAT3SE& operator=(uint32_t packed) noexcept { v = packed; return *this; }
        };

        //------------------------------------------------------------------------------
        // 4D Vector; 16 bit floating point components
        struct HALF4
        {
            union
            {
                struct
                {
                    HALF x;
                    HALF y;
                    HALF z;
                    HALF w;
                };
                uint64_t v;
            };

            HALF4() = default;

            HALF4(const HALF4&) = default;
            HALF4& operator=(const HALF4&) = default;

            HALF4(HALF4&&) = default;
            HALF4& operator=(HALF4&&) = default;

            explicit constexpr HALF4(uint64_t packed) noexcept : v(packed) {}
            constexpr HALF4(HALF _x, HALF _y, HALF _z, HALF _w) noexcept : x(_x), y(_y), z(_z), w(_w) {}
            explicit HALF4(_In_reads_(4) const HALF* pArray) noexcept : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}
            HALF4(float _x, float _y, float _z, float _w) noexcept;
            explicit HALF4(_In_reads_(4) const float* pArray) noexcept;

            HALF4& operator=(uint64_t packed) noexcept { v = packed; return *this; }
        };

        //------------------------------------------------------------------------------
        // 4D Vector; 16 bit signed normalized integer components
        struct SHORTN4
        {
            union
            {
                struct
                {
                    int16_t x;
                    int16_t y;
                    int16_t z;
                    int16_t w;
                };
                uint64_t v;
            };

            SHORTN4() = default;

            SHORTN4(const SHORTN4&) = default;
            SHORTN4& operator=(const SHORTN4&) = default;

            SHORTN4(SHORTN4&&) = default;
            SHORTN4& operator=(SHORTN4&&) = default;

            explicit constexpr SHORTN4(uint64_t packed) noexcept : v(packed) {}
            constexpr SHORTN4(int16_t _x, int16_t _y, int16_t _z, int16_t _w) noexcept : x(_x), y(_y), z(_z), w(_w) {}
            explicit SHORTN4(_In_reads_(4) const int16_t* pArray) noexcept : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}
            SHORTN4(float _x, float _y, float _z, float _w) noexcept;
            explicit SHORTN4(_In_reads_(4) const float* pArray) noexcept;

            SHORTN4& operator=(uint64_t packed) noexcept { v = packed; return *this; }
        };

        // 4D Vector; 16 bit signed integer components
        struct SHORT4
        {
            union
            {
                struct
                {
                    int16_t x;
                    int16_t y;
                    int16_t z;
                    int16_t w;
                };
                uint64_t v;
            };

            SHORT4() = default;

            SHORT4(const SHORT4&) = default;
            SHORT4& operator=(const SHORT4&) = default;

            SHORT4(SHORT4&&) = default;
            SHORT4& operator=(SHORT4&&) = default;

            explicit constexpr SHORT4(uint64_t packed) noexcept : v(packed) {}
            constexpr SHORT4(int16_t _x, int16_t _y, int16_t _z, int16_t _w) noexcept : x(_x), y(_y), z(_z), w(_w) {}
            explicit SHORT4(_In_reads_(4) const int16_t* pArray) noexcept : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}
            SHORT4(float _x, float _y, float _z, float _w) noexcept;
            explicit SHORT4(_In_reads_(4) const float* pArray) noexcept;

            SHORT4& operator=(uint64_t packed) noexcept { v = packed; return *this; }
        };

        // 4D Vector; 16 bit unsigned normalized integer components
        struct USHORTN4
        {
            union
            {
                struct
                {
                    uint16_t x;
                    uint16_t y;
                    uint16_t z;
                    uint16_t w;
                };
                uint64_t v;
            };

            USHORTN4() = default;

            USHORTN4(const USHORTN4&) = default;
            USHORTN4& operator=(const USHORTN4&) = default;

            USHORTN4(USHORTN4&&) = default;
            USHORTN4& operator=(USHORTN4&&) = default;

            explicit constexpr USHORTN4(uint64_t packed) noexcept : v(packed) {}
            constexpr USHORTN4(uint16_t _x, uint16_t _y, uint16_t _z, uint16_t _w) noexcept : x(_x), y(_y), z(_z), w(_w) {}
            explicit USHORTN4(_In_reads_(4) const uint16_t* pArray) noexcept : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}
            USHORTN4(float _x, float _y, float _z, float _w) noexcept;
            explicit USHORTN4(_In_reads_(4) const float* pArray) noexcept;

            USHORTN4& operator=(uint64_t packed) noexcept { v = packed; return *this; }
        };

        // 4D Vector; 16 bit unsigned integer components
        struct USHORT4
        {
            union
            {
                struct
                {
                    uint16_t x;
                    uint16_t y;
                    uint16_t z;
                    uint16_t w;
                };
                uint64_t v;
            };

            USHORT4() = default;

            USHORT4(const USHORT4&) = default;
            USHORT4& operator=(const USHORT4&) = default;

            USHORT4(USHORT4&&) = default;
            USHORT4& operator=(USHORT4&&) = default;

            explicit constexpr USHORT4(uint64_t packed) noexcept : v(packed) {}
            constexpr USHORT4(uint16_t _x, uint16_t _y, uint16_t _z, uint16_t _w) noexcept : x(_x), y(_y), z(_z), w(_w) {}
            explicit USHORT4(_In_reads_(4) const uint16_t* pArray) noexcept : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}
            USHORT4(float _x, float _y, float _z, float _w) noexcept;
            explicit USHORT4(_In_reads_(4) const float* pArray) noexcept;

            USHORT4& operator=(uint32_t packed) noexcept { v = packed; return *this; }
        };

        //------------------------------------------------------------------------------
        // 4D Vector; 10-10-10-2 bit normalized components packed into a 32 bit integer
        // The normalized 4D Vector is packed into 32 bits as follows: a 2 bit unsigned,
        // normalized integer for the w component and 10 bit signed, normalized
        // integers for the z, y, and x components.  The w component is stored in the
        // most significant bits and the x component in the least significant bits
        // (W2Z10Y10X10): [32] wwzzzzzz zzzzyyyy yyyyyyxx xxxxxxxx [0]
        struct XDECN4
        {
            union
            {
                struct
                {
                    int32_t x : 10;    // -511/511 to 511/511
                    int32_t y : 10;    // -511/511 to 511/511
                    int32_t z : 10;    // -511/511 to 511/511
                    uint32_t w : 2;    //      0/3 to     3/3
                };
                uint32_t v;
            };

            XDECN4() = default;

            XDECN4(const XDECN4&) = default;
            XDECN4& operator=(const XDECN4&) = default;

            XDECN4(XDECN4&&) = default;
            XDECN4& operator=(XDECN4&&) = default;

            explicit constexpr XDECN4(uint32_t packed) : v(packed) {}
            XDECN4(float _x, float _y, float _z, float _w) noexcept;
            explicit XDECN4(_In_reads_(4) const float* pArray) noexcept;

            operator uint32_t () const noexcept { return v; }

            XDECN4& operator=(uint32_t packed) noexcept { v = packed; return *this; }
        };

        // 4D Vector; 10-10-10-2 bit components packed into a 32 bit integer
        // The normalized 4D Vector is packed into 32 bits as follows: a 2 bit unsigned
        // integer for the w component and 10 bit signed integers for the
        // z, y, and x components.  The w component is stored in the
        // most significant bits and the x component in the least significant bits
        // (W2Z10Y10X10): [32] wwzzzzzz zzzzyyyy yyyyyyxx xxxxxxxx [0]
        struct DEPRECATED XDEC4
        {
            union
            {
                struct
                {
                    int32_t x : 10;    // -511 to 511
                    int32_t y : 10;    // -511 to 511
                    int32_t z : 10;    // -511 to 511
                    uint32_t w : 2;    // 0 to 3
                };
                uint32_t v;
            };

            XDEC4() = default;

            XDEC4(const XDEC4&) = default;
            XDEC4& operator=(const XDEC4&) = default;

            XDEC4(XDEC4&&) = default;
            XDEC4& operator=(XDEC4&&) = default;

            explicit constexpr XDEC4(uint32_t packed) noexcept : v(packed) {}
            XDEC4(float _x, float _y, float _z, float _w) noexcept;
            explicit XDEC4(_In_reads_(4) const float* pArray) noexcept;

            operator uint32_t () const noexcept { return v; }

            XDEC4& operator=(uint32_t packed) noexcept { v = packed; return *this; }
        };

        // 4D Vector; 10-10-10-2 bit normalized components packed into a 32 bit integer
        // The normalized 4D Vector is packed into 32 bits as follows: a 2 bit signed,
        // normalized integer for the w component and 10 bit signed, normalized
        // integers for the z, y, and x components.  The w component is stored in the
        // most significant bits and the x component in the least significant bits
        // (W2Z10Y10X10): [32] wwzzzzzz zzzzyyyy yyyyyyxx xxxxxxxx [0]
        struct DEPRECATED DECN4
        {
            union
            {
                struct
                {
                    int32_t x : 10;    // -511/511 to 511/511
                    int32_t y : 10;    // -511/511 to 511/511
                    int32_t z : 10;    // -511/511 to 511/511
                    int32_t w : 2;     //     -1/1 to     1/1
                };
                uint32_t v;
            };

            DECN4() = default;

            DECN4(const DECN4&) = default;
            DECN4& operator=(const DECN4&) = default;

            DECN4(DECN4&&) = default;
            DECN4& operator=(DECN4&&) = default;

            explicit constexpr DECN4(uint32_t packed) noexcept : v(packed) {}
            DECN4(float _x, float _y, float _z, float _w) noexcept;
            explicit DECN4(_In_reads_(4) const float* pArray) noexcept;

            operator uint32_t () const noexcept { return v; }

            DECN4& operator=(uint32_t packed) noexcept { v = packed; return *this; }
        };

        // 4D Vector; 10-10-10-2 bit components packed into a 32 bit integer
        // The 4D Vector is packed into 32 bits as follows: a 2 bit signed,
        // integer for the w component and 10 bit signed integers for the
        // z, y, and x components.  The w component is stored in the
        // most significant bits and the x component in the least significant bits
        // (W2Z10Y10X10): [32] wwzzzzzz zzzzyyyy yyyyyyxx xxxxxxxx [0]
        struct DEPRECATED DEC4
        {
            union
            {
                struct
                {
                    int32_t  x : 10;    // -511 to 511
                    int32_t  y : 10;    // -511 to 511
                    int32_t  z : 10;    // -511 to 511
                    int32_t  w : 2;     //   -1 to   1
                };
                uint32_t v;
            };

            DEC4() = default;

            DEC4(const DEC4&) = default;
            DEC4& operator=(const DEC4&) = default;

            DEC4(DEC4&&) = default;
            DEC4& operator=(DEC4&&) = default;

            explicit constexpr DEC4(uint32_t packed) noexcept : v(packed) {}
            DEC4(float _x, float _y, float _z, float _w) noexcept;
            explicit DEC4(_In_reads_(4) const float* pArray) noexcept;

            operator uint32_t () const noexcept { return v; }

            DEC4& operator=(uint32_t packed) noexcept { v = packed; return *this; }
        };

        // 4D Vector; 10-10-10-2 bit normalized components packed into a 32 bit integer
        // The normalized 4D Vector is packed into 32 bits as follows: a 2 bit unsigned,
        // normalized integer for the w component and 10 bit unsigned, normalized
        // integers for the z, y, and x components.  The w component is stored in the
        // most significant bits and the x component in the least significant bits
        // (W2Z10Y10X10): [32] wwzzzzzz zzzzyyyy yyyyyyxx xxxxxxxx [0]
        struct UDECN4
        {
            union
            {
                struct
                {
                    uint32_t x : 10;    // 0/1023 to 1023/1023
                    uint32_t y : 10;    // 0/1023 to 1023/1023
                    uint32_t z : 10;    // 0/1023 to 1023/1023
                    uint32_t w : 2;     //    0/3 to       3/3
                };
                uint32_t v;
            };

            UDECN4() = default;

            UDECN4(const UDECN4&) = default;
            UDECN4& operator=(const UDECN4&) = default;

            UDECN4(UDECN4&&) = default;
            UDECN4& operator=(UDECN4&&) = default;

            explicit constexpr UDECN4(uint32_t packed) noexcept : v(packed) {}
            UDECN4(float _x, float _y, float _z, float _w) noexcept;
            explicit UDECN4(_In_reads_(4) const float* pArray) noexcept;

            operator uint32_t () const noexcept { return v; }

            UDECN4& operator=(uint32_t packed) noexcept { v = packed; return *this; }
        };

        // 4D Vector; 10-10-10-2 bit components packed into a 32 bit integer
        // The 4D Vector is packed into 32 bits as follows: a 2 bit unsigned,
        // integer for the w component and 10 bit unsigned integers
        // for the z, y, and x components.  The w component is stored in the
        // most significant bits and the x component in the least significant bits
        // (W2Z10Y10X10): [32] wwzzzzzz zzzzyyyy yyyyyyxx xxxxxxxx [0]
        struct UDEC4
        {
            union
            {
                struct
                {
                    uint32_t x : 10;    // 0 to 1023
                    uint32_t y : 10;    // 0 to 1023
                    uint32_t z : 10;    // 0 to 1023
                    uint32_t w : 2;     // 0 to    3
                };
                uint32_t v;
            };

            UDEC4() = default;

            UDEC4(const UDEC4&) = default;
            UDEC4& operator=(const UDEC4&) = default;

            UDEC4(UDEC4&&) = default;
            UDEC4& operator=(UDEC4&&) = default;

            explicit constexpr UDEC4(uint32_t packed) noexcept : v(packed) {}
            UDEC4(float _x, float _y, float _z, float _w) noexcept;
            explicit UDEC4(_In_reads_(4) const float* pArray) noexcept;

            operator uint32_t () const noexcept { return v; }

            UDEC4& operator=(uint32_t packed) noexcept { v = packed; return *this; }
        };

        //------------------------------------------------------------------------------
        // 4D Vector; 8 bit signed normalized integer components
        struct BYTEN4
        {
            union
            {
                struct
                {
                    int8_t x;
                    int8_t y;
                    int8_t z;
                    int8_t w;
                };
                uint32_t v;
            };

            BYTEN4() = default;

            BYTEN4(const BYTEN4&) = default;
            BYTEN4& operator=(const BYTEN4&) = default;

            BYTEN4(BYTEN4&&) = default;
            BYTEN4& operator=(BYTEN4&&) = default;

            constexpr BYTEN4(int8_t _x, int8_t _y, int8_t _z, int8_t _w) noexcept : x(_x), y(_y), z(_z), w(_w) {}
            explicit constexpr BYTEN4(uint32_t packed) noexcept : v(packed) {}
            explicit BYTEN4(_In_reads_(4) const int8_t* pArray) noexcept : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}
            BYTEN4(float _x, float _y, float _z, float _w) noexcept;
            explicit BYTEN4(_In_reads_(4) const float* pArray) noexcept;

            BYTEN4& operator=(uint32_t packed) noexcept { v = packed; return *this; }
        };

        // 4D Vector; 8 bit signed integer components
        struct BYTE4
        {
            union
            {
                struct
                {
                    int8_t x;
                    int8_t y;
                    int8_t z;
                    int8_t w;
                };
                uint32_t v;
            };

            BYTE4() = default;

            BYTE4(const BYTE4&) = default;
            BYTE4& operator=(const BYTE4&) = default;

            BYTE4(BYTE4&&) = default;
            BYTE4& operator=(BYTE4&&) = default;

            constexpr BYTE4(int8_t _x, int8_t _y, int8_t _z, int8_t _w) noexcept : x(_x), y(_y), z(_z), w(_w) {}
            explicit constexpr BYTE4(uint32_t packed) noexcept : v(packed) {}
            explicit BYTE4(_In_reads_(4) const int8_t* pArray) noexcept : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}
            BYTE4(float _x, float _y, float _z, float _w) noexcept;
            explicit BYTE4(_In_reads_(4) const float* pArray) noexcept;

            BYTE4& operator=(uint32_t packed) noexcept { v = packed; return *this; }
        };

        // 4D Vector; 8 bit unsigned normalized integer components
        struct UBYTEN4
        {
            union
            {
                struct
                {
                    uint8_t x;
                    uint8_t y;
                    uint8_t z;
                    uint8_t w;
                };
                uint32_t v;
            };

            UBYTEN4() = default;

            UBYTEN4(const UBYTEN4&) = default;
            UBYTEN4& operator=(const UBYTEN4&) = default;

            UBYTEN4(UBYTEN4&&) = default;
            UBYTEN4& operator=(UBYTEN4&&) = default;

            constexpr UBYTEN4(uint8_t _x, uint8_t _y, uint8_t _z, uint8_t _w) noexcept : x(_x), y(_y), z(_z), w(_w) {}
            explicit constexpr UBYTEN4(uint32_t packed) noexcept : v(packed) {}
            explicit UBYTEN4(_In_reads_(4) const uint8_t* pArray) noexcept : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}
            UBYTEN4(float _x, float _y, float _z, float _w) noexcept;
            explicit UBYTEN4(_In_reads_(4) const float* pArray) noexcept;

            UBYTEN4& operator=(uint32_t packed) noexcept { v = packed; return *this; }
        };

        // 4D Vector; 8 bit unsigned integer components
        struct UBYTE4
        {
            union
            {
                struct
                {
                    uint8_t x;
                    uint8_t y;
                    uint8_t z;
                    uint8_t w;
                };
                uint32_t v;
            };

            UBYTE4() = default;

            UBYTE4(const UBYTE4&) = default;
            UBYTE4& operator=(const UBYTE4&) = default;

            UBYTE4(UBYTE4&&) = default;
            UBYTE4& operator=(UBYTE4&&) = default;

            constexpr UBYTE4(uint8_t _x, uint8_t _y, uint8_t _z, uint8_t _w) noexcept : x(_x), y(_y), z(_z), w(_w) {}
            explicit constexpr UBYTE4(uint32_t packed)  noexcept : v(packed) {}
            explicit UBYTE4(_In_reads_(4) const uint8_t* pArray)  noexcept : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}
            UBYTE4(float _x, float _y, float _z, float _w) noexcept;
            explicit UBYTE4(_In_reads_(4) const float* pArray) noexcept;

            UBYTE4& operator=(uint32_t packed)  noexcept { v = packed; return *this; }
        };

        //------------------------------------------------------------------------------
        // 4D vector; 4 bit unsigned integer components
        struct UNIBBLE4
        {
            union
            {
                struct
                {
                    uint16_t x : 4;    // 0 to 15
                    uint16_t y : 4;    // 0 to 15
                    uint16_t z : 4;    // 0 to 15
                    uint16_t w : 4;    // 0 to 15
                };
                uint16_t v;
            };

            UNIBBLE4() = default;

            UNIBBLE4(const UNIBBLE4&) = default;
            UNIBBLE4& operator=(const UNIBBLE4&) = default;

            UNIBBLE4(UNIBBLE4&&) = default;
            UNIBBLE4& operator=(UNIBBLE4&&) = default;

            explicit constexpr UNIBBLE4(uint16_t packed)  noexcept : v(packed) {}
            constexpr UNIBBLE4(uint8_t _x, uint8_t _y, uint8_t _z, uint8_t _w)  noexcept : x(_x), y(_y), z(_z), w(_w) {}
            explicit UNIBBLE4(_In_reads_(4) const uint8_t* pArray)  noexcept : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}
            UNIBBLE4(float _x, float _y, float _z, float _w) noexcept;
            explicit UNIBBLE4(_In_reads_(4) const float* pArray) noexcept;

            operator uint16_t() const  noexcept { return v; }

            UNIBBLE4& operator=(uint16_t packed) noexcept { v = packed; return *this; }
        };

        //------------------------------------------------------------------------------
        // 4D vector: 5/5/5/1 unsigned integer components
        struct U555
        {
            union
            {
                struct
                {
                    uint16_t x : 5;    // 0 to 31
                    uint16_t y : 5;    // 0 to 31
                    uint16_t z : 5;    // 0 to 31
                    uint16_t w : 1;    // 0 or 1
                };
                uint16_t v;
            };

            U555() = default;

            U555(const U555&) = default;
            U555& operator=(const U555&) = default;

            U555(U555&&) = default;
            U555& operator=(U555&&) = default;

            explicit constexpr U555(uint16_t packed) noexcept : v(packed) {}
            constexpr U555(uint8_t _x, uint8_t _y, uint8_t _z, bool _w) noexcept : x(_x), y(_y), z(_z), w(_w ? 0x1 : 0) {}
            U555(_In_reads_(3) const uint8_t* pArray, _In_ bool _w) noexcept : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(_w ? 0x1 : 0) {}
            U555(float _x, float _y, float _z, bool _w) noexcept;
            U555(_In_reads_(3) const float* pArray, _In_ bool _w) noexcept;

            operator uint16_t() const noexcept { return v; }

            U555& operator=(uint16_t packed) noexcept { v = packed; return *this; }
        };

#ifdef __clang__
#pragma clang diagnostic pop
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif

        float ConvertHalfToFloat(HALF Value) noexcept;
        float* ConvertHalfToFloatStream(_Out_writes_bytes_(sizeof(float) + outputStride * (halfCount - 1)) float* pOutputStream,
            _In_ size_t outputStride,
            _In_reads_bytes_(sizeof(HALF) + inputStride * (halfCount - 1)) const HALF* pInputStream,
            _In_ size_t inputStride, _In_ size_t halfCount) noexcept;
        HALF ConvertFloatToHalf(float value) noexcept;
        HALF* ConvertFloatToHalfStream(_Out_writes_bytes_(sizeof(HALF) + outputStride * (floatCount - 1)) HALF* pOutputStream,
            _In_ size_t outputStride,
            _In_reads_bytes_(sizeof(float) + inputStride * (floatCount - 1)) const float* pInputStream,
            _In_ size_t inputStride, _In_ size_t floatCount) noexcept;

        VECTOR VEC_CALLCONV LoadColor(_In_ const COLOR* pSource) noexcept;

        VECTOR VEC_CALLCONV LoadHalf2(_In_ const HALF2* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadShortN2(_In_ const SHORTN2* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadShort2(_In_ const SHORT2* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadUShortN2(_In_ const USHORTN2* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadUShort2(_In_ const USHORT2* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadByteN2(_In_ const BYTEN2* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadByte2(_In_ const BYTE2* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadUByteN2(_In_ const UBYTEN2* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadUByte2(_In_ const UBYTE2* pSource) noexcept;

        VECTOR VEC_CALLCONV LoadU565(_In_ const U565* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadFloat3PK(_In_ const FLOAT3PK* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadFloat3SE(_In_ const FLOAT3SE* pSource) noexcept;

        VECTOR VEC_CALLCONV LoadHalf4(_In_ const HALF4* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadShortN4(_In_ const SHORTN4* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadShort4(_In_ const SHORT4* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadUShortN4(_In_ const USHORTN4* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadUShort4(_In_ const USHORT4* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadXDecN4(_In_ const XDECN4* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadUDecN4(_In_ const UDECN4* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadUDecN4_XR(_In_ const UDECN4* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadUDec4(_In_ const UDEC4* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadByteN4(_In_ const BYTEN4* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadByte4(_In_ const BYTE4* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadUByteN4(_In_ const UBYTEN4* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadUByte4(_In_ const UBYTE4* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadUNibble4(_In_ const UNIBBLE4* pSource) noexcept;
        VECTOR VEC_CALLCONV LoadU555(_In_ const U555* pSource) noexcept;

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4996)
        // C4996: ignore deprecation warning
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

        DEPRECATED
        VECTOR VEC_CALLCONV LoadDecN4(_In_ const DECN4* pSource) noexcept;

        DEPRECATED
        VECTOR VEC_CALLCONV LoadDec4(_In_ const DEC4* pSource) noexcept;

        DEPRECATED
        VECTOR VEC_CALLCONV LoadXDec4(_In_ const XDEC4* pSource) noexcept;

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif

        void VEC_CALLCONV StoreColor(_Out_ COLOR* pDestination, _In_ A_VECTOR V) noexcept;

        void VEC_CALLCONV StoreHalf2(_Out_ HALF2* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreShortN2(_Out_ SHORTN2* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreShort2(_Out_ SHORT2* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreUShortN2(_Out_ USHORTN2* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreUShort2(_Out_ USHORT2* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreByteN2(_Out_ BYTEN2* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreByte2(_Out_ BYTE2* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreUByteN2(_Out_ UBYTEN2* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreUByte2(_Out_ UBYTE2* pDestination, _In_ A_VECTOR V) noexcept;

        void VEC_CALLCONV StoreU565(_Out_ U565* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreFloat3PK(_Out_ FLOAT3PK* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreFloat3SE(_Out_ FLOAT3SE* pDestination, _In_ A_VECTOR V) noexcept;

        void VEC_CALLCONV StoreHalf4(_Out_ HALF4* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreShortN4(_Out_ SHORTN4* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreShort4(_Out_ SHORT4* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreUShortN4(_Out_ USHORTN4* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreUShort4(_Out_ USHORT4* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreXDecN4(_Out_ XDECN4* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreUDecN4(_Out_ UDECN4* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreUDecN4_XR(_Out_ UDECN4* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreUDec4(_Out_ UDEC4* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreByteN4(_Out_ BYTEN4* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreByte4(_Out_ BYTE4* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreUByteN4(_Out_ UBYTEN4* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreUByte4(_Out_ UBYTE4* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreUNibble4(_Out_ UNIBBLE4* pDestination, _In_ A_VECTOR V) noexcept;
        void VEC_CALLCONV StoreU555(_Out_ U555* pDestination, _In_ A_VECTOR V) noexcept;

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4996)
        // C4996: ignore deprecation warning
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

        DEPRECATED
        void VEC_CALLCONV StoreDecN4(_Out_ DECN4* pDestination, _In_ A_VECTOR V) noexcept;

        DEPRECATED
        void VEC_CALLCONV StoreDec4(_Out_ DEC4* pDestination, _In_ A_VECTOR V) noexcept;

        DEPRECATED
        void VEC_CALLCONV StoreXDec4(_Out_ XDEC4* pDestination, _In_ A_VECTOR V) noexcept;

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4068 4214 4204 4365 4616 6001 6101)
         // C4068/4616: ignore unknown pragmas
         // C4214/4204: nonstandard extension used
         // C4365: Off by default noise
         // C6001/6101: False positives
#endif

#ifdef _PREFAST_
#pragma prefast(push)
#pragma prefast(disable : 25000, "A_VECTOR is 16 bytes")
#pragma prefast(disable : 26495, "Union initialization confuses /analyze")
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif

#include "PackedVector.inl"

#ifdef __clang__
#pragma clang diagnostic pop
#endif
#ifdef _PREFAST_
#pragma prefast(pop)
#endif
#ifdef _MSC_VER
#pragma warning(pop)
#endif
    } // namespace PackedVector
    
} // namespace UltReality::Math


#endif // !ULTREALITY_MATH_PACKED_VECTOR_H