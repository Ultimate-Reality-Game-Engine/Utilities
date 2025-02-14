#ifndef ULTREALITY_MATH_VECTOR_CONSTANTS_H
#define ULTREALITY_MATH_VECTOR_CONSTANTS_H

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4005 4668)
// C4005/4668: Old header issue
#endif
#include <stdint.h>
#ifdef _MSC_VER
#pragma warning(pop)
#endif

namespace UltReality::Math
{
	constexpr uint32_t SELECT_NONE = 0x00000000;
	constexpr uint32_t SELECT_ALL = 0xFFFFFFFF;

	constexpr uint32_t PERMUTE_0X = 0;
	constexpr uint32_t PERMUTE_0Y = 1;
	constexpr uint32_t PERMUTE_0Z = 2;
	constexpr uint32_t PERMUTE_0W = 3;

	constexpr uint32_t PERMUTE_1X = 4;
	constexpr uint32_t PERMUTE_1Y = 5;
	constexpr uint32_t PERMUTE_1Z = 6;
	constexpr uint32_t PERMUTE_1W = 7;

	constexpr uint32_t SWIZZLE_X = 0;
	constexpr uint32_t SWIZZLE_Y = 1;
	constexpr uint32_t SWIZZLE_Z = 2;
	constexpr uint32_t SWIZZLE_W = 3;

	constexpr uint32_t CRMASK_CR6 = 0x000000F0;
	constexpr uint32_t CRMASK_CR6TRUE = 0x00000080;
	constexpr uint32_t CRMASK_CR6FALSE = 0x00000020;
	constexpr uint32_t CRMASK_CR6BOUNDS = CRMASK_CR6FALSE;

	constexpr size_t CACHE_LINE_SIZE = 64;

	constexpr bool CompareAllTrue(uint32_t cr) noexcept
	{
		return (cr & CRMASK_CR6TRUE) == CRMASK_CR6TRUE;
	}

	constexpr bool CompareAnyTrue(uint32_t cr) noexcept
	{
		return (cr & CRMASK_CR6FALSE) != CRMASK_CR6FALSE;
	}

	constexpr bool CompareAllFalse(uint32_t cr) noexcept
	{
		return (cr & CRMASK_CR6FALSE) == CRMASK_CR6FALSE;
	}

	constexpr bool CompareAnyFalse(uint32_t cr) noexcept
	{
		return (cr & CRMASK_CR6TRUE) != CRMASK_CR6TRUE;
	}

	constexpr bool CompareMixed(uint32_t cr) noexcept
	{
		return (cr & CRMASK_CR6) == 0;
	}

	constexpr bool CompareAllInBounds(uint32_t cr) noexcept
	{
		return (cr & CRMASK_CR6BOUNDS) == CRMASK_CR6BOUNDS;
	}

	constexpr bool CompareAnyInBounds(uint32_t cr) noexcept
	{
		return (cr & CRMASK_CR6BOUNDS) != CRMASK_CR6BOUNDS;
	}
}

#endif // !ULTREALITY_MATH_VECTOR_CONSTANTS_H
