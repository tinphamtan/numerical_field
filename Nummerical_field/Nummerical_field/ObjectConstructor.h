#ifndef OBJECTCONSTRUCTOR_H
#define OBJECTCONSTRUCTOR_H
#pragma once

#ifndef SMALLSQUARE_H
#include "SmallSquare.h"
#endif

#ifndef RVEKTOR_H
#include "RVektor.h"
#endif

class ObjectConstructor
{
private:
	int N; // Total number of SmallSquares
	int nX; // the Number of element in x coordinate
	int nY; // the Number of element in y coordinate
	int nZ;
	double unitLength;
	Vektor<SmallSquare> Plate;
public:
	ObjectConstructor(void);
	~ObjectConstructor(void);
	void MakeRectangle();
};
#endif