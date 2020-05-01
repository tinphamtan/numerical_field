#ifndef RECTANGULAR_H
#define RECTANGULAR_H
#pragma once

#ifndef RVEKTOR_H
#include "RVektor.h"
#endif


#ifndef RMATRIX_H
#include "RMatrix.h"
#endif


#ifndef SMALLSQUARE_H
#include "SmallSquare.h"
#endif

class RectangularPlate
{
private:
	double plateLengh;// length of the plate
	double plateWidth;//width of the plate
	int wRes; //width resolution
	int lRes;// length resolution
	double oX;//x coordinate of the plate
	double oY;//y coordinate of the plate
	double oZ;//z coordinate of the plate
	double oA;//angle
	Vektor<SmallSquare> vPlate;
	typedef struct {double x, y, z;} point;
	point coordinate[];
	

public:
	RectangularPlate(void);
	~RectangularPlate(void);
	void CreatePlate(const double plateWidth, const int wRes, const int lRes, const double oA, const double oX, const double oY, const double oZ);
};
#endif;