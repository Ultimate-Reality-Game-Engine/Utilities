#ifndef ULTREALITY_MATH_VECTOR_CONFIG_H
#define ULTREALITY_MATH_VECTOR_CONFIG_H

#ifndef __cplusplus
#error This library requires C++
#endif

// Library version
#define VECTOR_MATH_LIB_VERSION 100 // 1.0.0

// Determine calling convention
#if defined(_MSC_VER) && !defined(_M_ARM) && !defined(_M_ARM64) && !defined(_M_HYBRID_X86_ARM64) && !defined(_M_ARM64EC) && (!_MANAGED) && (!_M_CEE) && (!defined(_M_IX86_FP) || (_M_IX86_FP > 1)) && !defined(_NO_INTRINSICS) && !defined(_VECTORCALL_)
#define _VECTORCALL_ 1
#endif

#if defined(_VECTORCALL_)
#define VECMATH_CALLCONV __vectorcall
#elif defined(__GNUC__) or defined(__clang__)
#define VECMATH_CALLCONV 
#else
#define VECMATH_CALLCONV __fastcall
#endif

#if !defined(_DEPRECATED)
#if (__cplusplus >= 201402L)
#define _DEPRECATED [[deprecated]]
#elif defined(__GNUC__) or defined(__clang__)
#define _DEPRECATED __attribute__ ((deprecated))
#else
#define _DEPRECATED __declspec(deprecated("This is deprecated and will be removed in a future version"))
#endif
#endif




#endif // !ULTREALITY_MATH_VECTOR_CONFIG_H
