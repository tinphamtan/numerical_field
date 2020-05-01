#include "stdafx.h"
#include <iostream>
#include <fstream>

using namespace std;
#ifndef RVEKTOR_H
#include "RVektor.h"
#endif

#ifndef SMALLSQUARE_H
#include "SmallSquare.h"
#endif

#ifndef VEKTOR_H
#include "Vector.h"
#endif

#ifndef SQUAREPLATE_H
#include "SquarePlate.h"
#endif

#ifndef RMATRIX_H
#include "RMatrix.h"
#endif

void (*Vektor<SmallSquare>::ErrHdlr)(char*) = Vektor<SmallSquare>::DefaultErrHdlr;
char Vektor<SmallSquare>::ErrFlag = 0;

//-----------------------------------
SquarePlate::SquarePlate(void)
{
	nX   = 5;
	nY   = nX;
	nN   = nX*nY;
	side = 1;
	potPlate = 1.0;
	totalCharge = 0.;

}
//-----------------------------------
SquarePlate::~SquarePlate(void)
{
}
//---------------------------------------
void SquarePlate::CreatePlate()
{

	
	vPlate.SetDim(nN);
	double delX = side/nX;
	double delY = delX; 
	int iN      = 0;
	for (int iX = -nX/2; iX <= nX/2; iX++)
	{
		for(int iY = -nY/2; iY <= nY/2; iY++)
		{
			double x = delX*iX;
			double y = delY*iY;
			double z = 0;
			vPlate[iN].MoveTo(x,y,z);
			vPlate[iN].SetSide(delX);
			std::cout << iN << "" << vPlate[iN];
			iN++;
		}
	}

}
void SquarePlate::CreateMatrix()
{
	std::cout << "\nCreate Matrix\n";
	
	myMatrix.SetDimension(nN,nN);
	for (int i = 0; i < nN; i++)
	{
		for (int j = 0; j < nN; j++)
		{
			myMatrix[i][j] = vPlate[i].CouplingCoefficient(vPlate[j]);
		}
	}
	std::cout << myMatrix;
}

void SquarePlate::CreateGmVector()
{
	std::cout << "\n\ngmVector\n\n";
	gmVec.SetDim(nN);
	for (int i = 0; i < nN; i++)
	{
		gmVec[i] = potPlate;
	}
	
	
	std::cout << gmVec;
}
void SquarePlate::FindChargeDistribution()
{
	ofstream file;
	file.open("c:/software/files/ChargeData_2.m");
	chargeVec.SetDim(nN);
	chargeVec = RealGaussElimination(myMatrix, gmVec);
	std::cout << "\n\n Charge Distribution \n\n";
	int iN = 0;
	for (int i = 0; i < nX; i++)
	{
		file << "\n"<<"  ";
		std::cout <<"\n"<< i<<"  ";
		for (int j = 0; j < nY; j++)
		{
			file << chargeVec[iN]<<"  ";
			std::cout << chargeVec[iN] << "  ";
			iN++;
		}
	}
	file.close();
}
void SquarePlate::CalculateTotalCharge()
{
	std::cout << "\n\nTotalCharge: ";
	int iN = 0;
	
	for (int i = 0; i < nX; i++){
		for (int j = 0; j < nY; j++)
		{

			totalCharge += chargeVec[iN]*(side/nX)*(side/nY);

			iN++;
		}
	}
	
	std::cout << totalCharge<<"\n";
}