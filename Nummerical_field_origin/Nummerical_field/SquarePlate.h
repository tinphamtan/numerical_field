#ifndef SQUAREPLATE_H
#define SQUAREPLATE_H
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

class SquarePlate
{
private:
	int nX; // Must be not even
	int nY; // How many SmallSquares are used to approximate the plate
	int nN; // Total number of SmallSquares
	double side; //length of the side in meter
	Vektor<SmallSquare> vPlate; //Plate composed of SmallSquares
	RVektor gmVec;
	RVektor gmVec_pos;
	RVektor chargeVec;
	RMatrix myMatrix;
	double potPlate;
	double totalCharge;

public:
	SquarePlate(void);
	~SquarePlate(void);
	void SetSide(const double aSide){side = aSide;}; //Set the side length of the plate
	void CreatePlate(); //Fillup Vektor to describe the plate
	void CreateMatrix();
	void CreateGmVector();
	void FindChargeDistribution();
	void CalculateTotalCharge();

	
};
#endif

