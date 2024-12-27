#ifndef ULTREALITY_MATH_SIMD_LEVEL_H
#define ULTREALITY_MATH_SIMD_LEVEL_H

// Determine calling convention
#if defined(_MSC_VER) && !defined(_M_ARM) && !defined(_M_ARM64) && !defined(_M_HYBRID_X86_ARM64) && !defined(_M_ARM64EC) && (!_MANAGED) && (!_M_CEE) && (!defined(_M_IX86_FP) || (_M_IX86_FP > 1)) && !defined(_NO_INTRINSICS_) && !defined(_VECTORCALL_)
#define _VECTORCALL_ 1
#endif

#if defined(_VECTORCALL_)
#define VEC_CALLCONV __vectorcall
#elif defined(__GNUC__)
#define VEC_CALLCONV 
#else
#define VEC_CALLCONV __fastcall
#endif

// Set deprecated macro value
#if !defined(_DEPRECATED)
#if (__cplusplus >= 201402L)
#define _DEPRECATED [[deprecated]]
#elif defined(__GNUC__) || defined(__clang__)
#define _DEPRECATED __attribute__ ((deprecated))
#else
#define _DEPRECATED __declspec(deprecated("This is deprecated and will be removed in a future version"))
#endif
#endif

// Set AVX2 Intrinsic level
#if !defined(_AVX2_INTRINSICS_) && defined(__AVX2__) && !defined(_NO_INTRINSICS_)
#define _AVX2_INTRINSICS_
#endif

// Set FMA3 Intrinsic level
#if !defined(_FMA3_INTRINSICS_) && defined(_AVX2_INTRINSICS_) && !defined(_NO_INTRINSICS_)
#define _FMA3_INTRINSICS_
#endif

// Set F16C Intrinsic level
#if !defined(_F16C_INTRINSICS_) && defined(_AVX2_INTRINSICS_) && !defined(_NO_INTRINSICS_)
#define _F16C_INTRINSICS_
#endif
#if !defined(_F16C_INTRINSICS_) && defined(__F16C__) && !defined(_NO_INTRINSICS_)
#define _F16C_INTRINSICS_
#endif

// Set AVX Intrinsic level
#if defined(_FMA3_INTRINSICS_) && !defined(_AVX_INTRINSICS_)
#define _AVX_INTRINSICS_
#endif
#if defined(_F16C_INTRINSICS_) && !defined(_AVX_INTRINSICS_)
#define _AVX_INTRINSICS_
#endif
#if !defined(_AVX_INTRINSICS_) && defined(__AVX__) && !defined(_NO_INTRINSICS_)
#define _AVX_INTRINSICS_
#endif

// Set SSE4 Intrinsic level
#if defined(_AVX_INTRINSICS_) && !defined(_SSE4_INTRINSICS_)
#define _SSE4_INTRINSICS_
#endif

// Set SSE3 Intrinsic level
#if defined(_SSE4_INTRINSICS_) && !defined(_SSE3_INTRINSICS_)
#define _SSE3_INTRINSICS_
#endif

// Set SSE2 Intrinsic level
#if defined(_SSE3_INTRINSICS_) && !defined(_SSE2_INTRINSICS_)
#define _SSE2_INTRINSICS_
#endif
#if !defined(_SSE2_INTRINSICS_) && !defined(_NO_INTRINSICS_)
#if (defined(_M_IX86) || defined(_M_X64) || __i386__ || __x86_64__) && !defined(_M_HYBRID_X86_ARM64) && !defined(_M_ARM64EC)
#define _SSE2_INTRINSICS_
#elif !defined(_NO_INTRINSICS_)
#error SIMD not supported on this target. Set _NO_INTRINSICS_ to explicitly compile without SIMD support
#endif
#endif

#if defined(_SSE2_INTRINSICS_) && defined(_MSC_VER) && (_MSC_VER >= 1920) && !defined(__clang__) && !defined(_SVML_INTRINSICS_) && !defined(_DISABLE_SVML_)
#define _SVML_INTRINSICS_
#endif

#endif // !ULTREALITY_MATH_SIMD_LEVEL_H
