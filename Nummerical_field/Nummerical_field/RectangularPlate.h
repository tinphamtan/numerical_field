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
	int N;
	double oX;//x coordinate of the plate
	double oY;//y coordinate of the plate
	double oZ;//z coordinate of the plate
	double oA;//angle
	Vektor<SmallSquare> vObject; //temporary vector space, will be deleted when better solution found
	Vektor<SmallSquare> cObject; // the combined object 
	Vektor<SmallSquare> pPlate; // contains potential plate
	RMatrix myMatrix;
	RVektor gmVec;
	RVektor defVec;
	RVektor potentialVec;

	RVektor chargeVec;
	RMatrix chargeMatrix_top;
	RMatrix potMatrix;

	double totalCharge;
	double totalCharge2;

public:
	RectangularPlate(void);
	~RectangularPlate(void);
	void CreatePlate(const double plateWidth, const int wRes, const int lRes, 
		const double oA, const double oX, const double oY, const double oZ);
	void ConstructObject(void);

	void AssignVectors(const int dim,const double res); //turn the recPlate.m file which contains bunch of coordinates of the object to a vector field
	void CreateMatrix(const int size);
	void CreateGmVector(const int size, const int sizeP1, const int sizeP2, const int sizeP3);
	void FindChargeDistribution(const int chargeVecSize);
	void FindTotalCharge(const double sideOfSmallSquare, const int totalDim);
	void FindPotential(const int sizeOfObject);
	void CreateYZPlate(const int yRes,const int zRes, const double yLength,
		const double oA, const double oX, const double oY, const double oZ);
	void CreatePotMatrix(const int yRes, const int zRes, const int size);
	void CreateAVector(const int size);
	void CalculatePotential(const int yRes, const int zRes);
	void TurnThemToMatrix(const int dimOfObject,const int dimP1w, const int dimP1l,const int dimP2w,const int dimP2l,const int dimP3w, const int dimP3l);
};

#endif;