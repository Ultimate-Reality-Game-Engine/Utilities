#ifndef ULTREALITY_UTILITIES_LIBRARY_EXPORT_H
#define ULTREALITY_UTILITIES_LIBRARY_EXPORT_H

#if defined(_WIN_TARGET)
#if defined(LIBRARY_EXPORTS)
#define LIBRARY_ABI __declspec(dllexport)
#else
#define LIBRARY_ABI __declspec(dllimport)
#endif

#define LIBRARY_CALL __cdecl
#elif defined(_LINUX_TARGET) or defined(_MAC_TARGET)
#if defined(LIBRARY_EXPORTS)
#define LIBRARY_ABI __attribute__((visibility("default")))
#else
#define LIBRARY_ABI __attribute__((visibility("default")))
#endif

#define LIBRARY_CALL __cdecl
#endif

#endif // !ULTREALITY_UTILITIES_LIBRARY_EXPORT_H
