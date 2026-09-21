#include <cmath>
#include <iostream>
#ifndef PICOLV2
#include "Sinh.h"
#include "ASinh.h"
#include "Cosh.h"
#endif
#include "HyperbolicTables.h"


using namespace std;

#ifdef PICOLV2
static float fast_exp(float x)
{
	/* Range reduction to 2^n * 2^f followed by a fourth-order expansion.
	 * This keeps the distortion port compact and is sufficiently accurate for
	 * its diode model without the original 71 KiB lookup tables. */
	x = x < -87.0f ? -87.0f : (x > 87.0f ? 87.0f : x);
	const float scaled = x * 1.4426950408889634f;
	int exponent = (int)scaled;
	if ((float)exponent > scaled) --exponent;
	const float fraction = scaled - (float)exponent;
	const float mantissa = 1.0f + fraction *
		(0.6931471806f + fraction *
		(0.2402265070f + fraction *
		(0.0555041087f + fraction * 0.0096181291f)));
	union { unsigned int bits; float value; } scale;
	scale.bits = (unsigned int)(exponent + 127) << 23;
	return scale.value * mantissa;
}

void SINHCOSH(float x, float *sinh_value, float *cosh_value)
{
	const float positive = fast_exp(x);
	const float negative = fast_exp(-x);
	*sinh_value = 0.5f * (positive - negative);
	*cosh_value = 0.5f * (positive + negative);
}

float SINH(float x)
{
	float sinh_value, cosh_value;
	SINHCOSH(x, &sinh_value, &cosh_value);
	return sinh_value;
}

float COSH(float x)
{
	float sinh_value, cosh_value;
	SINHCOSH(x, &sinh_value, &cosh_value);
	return cosh_value;
}

float ASINH(float x)
{
	return logf(x + sqrtf(x*x + 1.0f));
}
#else
float SINH( float x)
{
	int flag = 1;
	
	if (x < 0)
	{
		flag = -1;
		x = -x;
	}
	
	float SinH;
	
	if( x > static_cast<float>(SINH_fim))
	{
		SinH = (Sinh[SINH_N-1])*flag;
	}
	else
	{
		float naux = x*SINH_Idx;
		int n = round(naux);
		SinH = Sinh[n]*flag;
	}
	
	return SinH;
}

float COSH( float x)
{
	if (x < 0)
	{
		x = -x;
	}
	
	float CosH;
	
	if( x > static_cast<float>(COSH_fim))
	{
		CosH = Cosh[COSH_N-1];
	}
	else
	{
		float naux = x*COSH_Idx;
		int n = round(naux);
		CosH = Cosh[n];
	}
	
	return CosH;
}

float ASINH( float x)
{
	int flag = 1;
	
	if (x < 0)
	{
		flag = -1;
		x = -x;
	}
	
	float ASinH;
	
	if( x > static_cast<float>(ASINH_fim))
	{
		ASinH = (ASinh[ASINH_N-1])*flag;
	}
	else
	{
		float naux = x*ASINH_Idx;
		int n = round(naux);
		ASinH = ASinh[n]*flag;
	}
	
	return ASinH;
}
#endif
