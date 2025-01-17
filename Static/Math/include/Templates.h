#ifndef ULTREALITY_MATH_SIMD_TEMPLATES_H
#define ULTREALITY_MATH_SIMD_TEMPLATES_H

#include <SSE2VectorConfig.h>

#if defined(__GNUC__) or defined(__clang__)
#define FORCE_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline
#endif

namespace UltReality::Math
{
#if defined(__XNAMATH_H__) && defined(XMMin)
#undef XMMin
#undef XMMax
#endif

    template<class T> FORCE_INLINE T XMMin(T a, T b) noexcept { return (a < b) ? a : b; }
    template<class T> FORCE_INLINE T XMMax(T a, T b) noexcept { return (a > b) ? a : b; }

#if defined(_SSE2_INTRINSICS_) && !defined(_NO_INTRINSICS_)
    // PermuteHelper internal template (SSE only)
    namespace
    {
        // Slow path fallback for permutes that do not map to a single SSE shuffle opcode.
        template<uint32_t Shuffle, bool WhichX, bool WhichY, bool WhichZ, bool WhichW> struct PermuteHelper
        {
            static VECTOR VEC_CALLCONV Permute(A_VECTOR v1, A_VECTOR v2) noexcept
            {
                static const VECTOR_U32 selectMask =
                { { {
                        WhichX ? 0xFFFFFFFF : 0,
                        WhichY ? 0xFFFFFFFF : 0,
                        WhichZ ? 0xFFFFFFFF : 0,
                        WhichW ? 0xFFFFFFFF : 0,
                } } };

                VECTOR shuffled1 = PERMUTE_PS(v1, Shuffle);
                VECTOR shuffled2 = PERMUTE_PS(v2, Shuffle);

                VECTOR masked1 = _mm_andnot_ps(selectMask, shuffled1);
                VECTOR masked2 = _mm_and_ps(selectMask, shuffled2);

                return _mm_or_ps(masked1, masked2);
            }
        };

        // Fast path for permutes that only read from the first vector.
        template<uint32_t Shuffle> struct PermuteHelper<Shuffle, false, false, false, false>
        {
            static VECTOR VEC_CALLCONV Permute(A_VECTOR v1, A_VECTOR) noexcept { return PERMUTE_PS(v1, Shuffle); }
        };

        // Fast path for permutes that only read from the second vector.
        template<uint32_t Shuffle> struct PermuteHelper<Shuffle, true, true, true, true>
        {
            static VECTOR VEC_CALLCONV Permute(A_VECTOR, A_VECTOR v2) noexcept { return PERMUTE_PS(v2, Shuffle); }
        };

        // Fast path for permutes that read XY from the first vector, ZW from the second.
        template<uint32_t Shuffle> struct PermuteHelper<Shuffle, false, false, true, true>
        {
            static VECTOR VEC_CALLCONV Permute(A_VECTOR v1, A_VECTOR v2) noexcept { return _mm_shuffle_ps(v1, v2, Shuffle); }
        };
    }

#endif // _SSE2_INTRINSICS_ && !_NO_INTRINSICS_

    namespace VEC
    {
        // General permute template
        template<uint32_t PermuteX, uint32_t PermuteY, uint32_t PermuteZ, uint32_t PermuteW>
        inline VECTOR VEC_CALLCONV Permute(A_VECTOR V1, A_VECTOR V2) noexcept
        {
            static_assert(PermuteX <= 7, "PermuteX template parameter out of range");
            static_assert(PermuteY <= 7, "PermuteY template parameter out of range");
            static_assert(PermuteZ <= 7, "PermuteZ template parameter out of range");
            static_assert(PermuteW <= 7, "PermuteW template parameter out of range");

#if defined(_SSE2_INTRINSICS_) && !defined(_NO_INTRINSICS_)
            constexpr uint32_t shuffle = _MM_SHUFFLE(PermuteW & 3, PermuteZ & 3, PermuteY & 3, PermuteX & 3);

            constexpr bool WhichX = PermuteX > 3;
            constexpr bool WhichY = PermuteY > 3;
            constexpr bool WhichZ = PermuteZ > 3;
            constexpr bool WhichW = PermuteW > 3;

            return PermuteHelper<shuffle, WhichX, WhichY, WhichZ, WhichW>::Permute(V1, V2);

#else
            return Permute(V1, V2, PermuteX, PermuteY, PermuteZ, PermuteW);

#endif
        }

        // Special-case permute templates
        template<> constexpr VECTOR VEC_CALLCONV Permute<0, 1, 2, 3>(A_VECTOR V1, A_VECTOR) noexcept { return V1; }
        template<> constexpr VECTOR VEC_CALLCONV Permute<4, 5, 6, 7>(A_VECTOR, A_VECTOR V2) noexcept { return V2; }

#if defined(_SSE2_INTRINSICS_) && !defined(_NO_INTRINSICS_)
        template<> inline VECTOR VEC_CALLCONV Permute<0, 1, 4, 5>(A_VECTOR V1, A_VECTOR V2) noexcept { return _mm_movelh_ps(V1, V2); }
        template<> inline VECTOR VEC_CALLCONV Permute<6, 7, 2, 3>(A_VECTOR V1, A_VECTOR V2) noexcept { return _mm_movehl_ps(V1, V2); }
        template<> inline VECTOR VEC_CALLCONV Permute<0, 4, 1, 5>(A_VECTOR V1, A_VECTOR V2) noexcept { return _mm_unpacklo_ps(V1, V2); }
        template<> inline VECTOR VEC_CALLCONV Permute<2, 6, 3, 7>(A_VECTOR V1, A_VECTOR V2) noexcept { return _mm_unpackhi_ps(V1, V2); }
        template<> inline VECTOR VEC_CALLCONV Permute<2, 3, 6, 7>(A_VECTOR V1, A_VECTOR V2) noexcept { return _mm_castpd_ps(_mm_unpackhi_pd(_mm_castps_pd(V1), _mm_castps_pd(V2))); }
#endif

#if defined(_SSE4_INTRINSICS_) && !defined(_NO_INTRINSICS_)
        template<> inline VECTOR VEC_CALLCONV Permute<4, 1, 2, 3>(A_VECTOR V1, A_VECTOR V2) noexcept { return _mm_blend_ps(V1, V2, 0x1); }
        template<> inline VECTOR VEC_CALLCONV Permute<0, 5, 2, 3>(A_VECTOR V1, A_VECTOR V2) noexcept { return _mm_blend_ps(V1, V2, 0x2); }
        template<> inline VECTOR VEC_CALLCONV Permute<4, 5, 2, 3>(A_VECTOR V1, A_VECTOR V2) noexcept { return _mm_blend_ps(V1, V2, 0x3); }
        template<> inline VECTOR VEC_CALLCONV Permute<0, 1, 6, 3>(A_VECTOR V1, A_VECTOR V2) noexcept { return _mm_blend_ps(V1, V2, 0x4); }
        template<> inline VECTOR VEC_CALLCONV Permute<4, 1, 6, 3>(A_VECTOR V1, A_VECTOR V2) noexcept { return _mm_blend_ps(V1, V2, 0x5); }
        template<> inline VECTOR VEC_CALLCONV Permute<0, 5, 6, 3>(A_VECTOR V1, A_VECTOR V2) noexcept { return _mm_blend_ps(V1, V2, 0x6); }
        template<> inline VECTOR VEC_CALLCONV Permute<4, 5, 6, 3>(A_VECTOR V1, A_VECTOR V2) noexcept { return _mm_blend_ps(V1, V2, 0x7); }
        template<> inline VECTOR VEC_CALLCONV Permute<0, 1, 2, 7>(A_VECTOR V1, A_VECTOR V2) noexcept { return _mm_blend_ps(V1, V2, 0x8); }
        template<> inline VECTOR VEC_CALLCONV Permute<4, 1, 2, 7>(A_VECTOR V1, A_VECTOR V2) noexcept { return _mm_blend_ps(V1, V2, 0x9); }
        template<> inline VECTOR VEC_CALLCONV Permute<0, 5, 2, 7>(A_VECTOR V1, A_VECTOR V2) noexcept { return _mm_blend_ps(V1, V2, 0xA); }
        template<> inline VECTOR VEC_CALLCONV Permute<4, 5, 2, 7>(A_VECTOR V1, A_VECTOR V2) noexcept { return _mm_blend_ps(V1, V2, 0xB); }
        template<> inline VECTOR VEC_CALLCONV Permute<0, 1, 6, 7>(A_VECTOR V1, A_VECTOR V2) noexcept { return _mm_blend_ps(V1, V2, 0xC); }
        template<> inline VECTOR VEC_CALLCONV Permute<4, 1, 6, 7>(A_VECTOR V1, A_VECTOR V2) noexcept { return _mm_blend_ps(V1, V2, 0xD); }
        template<> inline VECTOR VEC_CALLCONV Permute<0, 5, 6, 7>(A_VECTOR V1, A_VECTOR V2) noexcept { return _mm_blend_ps(V1, V2, 0xE); }
#endif

         // General swizzle template
        template<uint32_t SwizzleX, uint32_t SwizzleY, uint32_t SwizzleZ, uint32_t SwizzleW>
        inline VECTOR VEC_CALLCONV Swizzle(A_VECTOR V) noexcept
        {
            static_assert(SwizzleX <= 3, "SwizzleX template parameter out of range");
            static_assert(SwizzleY <= 3, "SwizzleY template parameter out of range");
            static_assert(SwizzleZ <= 3, "SwizzleZ template parameter out of range");
            static_assert(SwizzleW <= 3, "SwizzleW template parameter out of range");

#if defined(_SSE2_INTRINSICS_) && !defined(_NO_INTRINSICS_)
            return PERMUTE_PS(V, _MM_SHUFFLE(SwizzleW, SwizzleZ, SwizzleY, SwizzleX));

#else
            return Swizzle(V, SwizzleX, SwizzleY, SwizzleZ, SwizzleW);

#endif
        }

        // Specialized swizzles
        template<> constexpr VECTOR VEC_CALLCONV Swizzle<0, 1, 2, 3>(A_VECTOR V) noexcept { return V; }

#if defined(_SSE2_INTRINSICS_) && !defined(_NO_INTRINSICS_)
        template<> inline VECTOR VEC_CALLCONV Swizzle<0, 1, 0, 1>(A_VECTOR V) noexcept { return _mm_movelh_ps(V, V); }
        template<> inline VECTOR VEC_CALLCONV Swizzle<2, 3, 2, 3>(A_VECTOR V) noexcept { return _mm_movehl_ps(V, V); }
        template<> inline VECTOR VEC_CALLCONV Swizzle<0, 0, 1, 1>(A_VECTOR V) noexcept { return _mm_unpacklo_ps(V, V); }
        template<> inline VECTOR VEC_CALLCONV Swizzle<2, 2, 3, 3>(A_VECTOR V) noexcept { return _mm_unpackhi_ps(V, V); }
#endif

#if defined(_SSE3_INTRINSICS_) && !defined(_NO_INTRINSICS_)
        template<> inline VECTOR VEC_CALLCONV Swizzle<0, 0, 2, 2>(A_VECTOR V) noexcept { return _mm_moveldup_ps(V); }
        template<> inline VECTOR VEC_CALLCONV Swizzle<1, 1, 3, 3>(A_VECTOR V) noexcept { return _mm_movehdup_ps(V); }
#endif

#if defined(_AVX2_INTRINSICS_) && !defined(_NO_INTRINSICS_) && defined(_FAVOR_INTEL_)
        template<> inline VECTOR VEC_CALLCONV Swizzle<0, 0, 0, 0>(A_VECTOR V) noexcept { return _mm_broadcastss_ps(V); }
#endif

        template<uint32_t Elements>
        inline VECTOR VEC_CALLCONV ShiftLeft(A_VECTOR V1, A_VECTOR V2) noexcept
        {
            static_assert(Elements < 4, "Elements template parameter out of range");
            return Permute<Elements, (Elements + 1), (Elements + 2), (Elements + 3)>(V1, V2);
        }

        template<uint32_t Elements>
        inline VECTOR VEC_CALLCONV RotateLeft(A_VECTOR V) noexcept
        {
            static_assert(Elements < 4, "Elements template parameter out of range");
            return Swizzle<Elements & 3, (Elements + 1) & 3, (Elements + 2) & 3, (Elements + 3) & 3>(V);
        }

        template<uint32_t Elements>
        inline VECTOR VEC_CALLCONV RotateRight(A_VECTOR V) noexcept
        {
            static_assert(Elements < 4, "Elements template parameter out of range");
            return Swizzle<(4 - Elements) & 3, (5 - Elements) & 3, (6 - Elements) & 3, (7 - Elements) & 3>(V);
        }

        template<uint32_t VSLeftRotateElements, uint32_t Select0, uint32_t Select1, uint32_t Select2, uint32_t Select3>
        inline VECTOR VEC_CALLCONV Insert(A_VECTOR VD, A_VECTOR VS) noexcept
        {
            VECTOR control = SelectControl(Select0 & 1, Select1 & 1, Select2 & 1, Select3 & 1);
            return Select(VD, RotateLeft<VSLeftRotateElements>(VS), control);
        }
    } // namespace VEC
}

#endif // !ULTREALITY_MATH_SIMD_TEMPLATES_H