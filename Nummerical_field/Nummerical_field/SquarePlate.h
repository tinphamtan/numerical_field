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
	int nX; // Must be not even because there will be no center smallsquare then
	int nY; // How many SmallSquares are used to approximate the plate
	int nZ; // define quantity of smallsquare of the z dimention
	int nN; // Total number of SmallSquares
	int slide;
	double side;
	double width;
	double height;
	double swidth;
	double length; 
	Vektor<SmallSquare> vPlate;
	Vektor<SmallSquare> pPlate;//Plate composed of SmallSquares
	RVektor gmVec;
	RVektor chargeVec;

	
	RMatrix myMatrix;
	
	double potPlate;
	double totalCharge;

	//--------------
	int potSide;
	RVektor potVec; //contains coordinate of the potential plate
	RVektor potentialVec;
	float potSideLength;
	RMatrix chargeMatrix;
	RMatrix potMatrix;

public:
	SquarePlate(void);
	~SquarePlate(void);
	void SetSide(const double aSide){side = aSide;}; //Set the side length of the plate
	void Create_U_Plate();
	void CreatePlate();  
	void CreateMatrix();
	void CreateGmVector();
	void FindChargeDistribution();
	void CalculateTotalCharge();
	//--------------------------
	void CreatePotMatrix();
	void CreateAVector();
	void CalculatePotential();

	
};
#endif

