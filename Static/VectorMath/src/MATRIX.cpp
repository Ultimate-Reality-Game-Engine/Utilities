#include <VectorTypes.h>

namespace UltReality::Math
{
	const MATRIX MATRIX::Zero(
		Vector::Zero(),
		Vector::Zero(),
		Vector::Zero(),
		Vector::Zero()
	);

	const MATRIX MATRIX::Identity(
		Vector::Set(1.0f, 0.0f, 0.0f, 0.0f),
		Vector::Set(0.0f, 1.0f, 0.0f, 0.0f),
		Vector::Set(0.0f, 0.0f, 1.0f, 0.0f),
		Vector::Set(0.0f, 0.0f, 0.0f, 1.0f)
	);
}