// Nummerical_field.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Vector.h"
#include <iostream>
#include <fstream>

#ifndef CMPLX_H
#include "cmplx.h"
#endif



#ifndef RVEKTOR_H
#include "RVektor.h"
#endif

#ifndef RMATRIX_H
#include "RMatrix.h"
#endif

#ifndef SQUAREPLATE_H
#include "SquarePlate.h"
#endif

#ifndef SMALLSQUARE_H
#include "SmallSquare.h"
#endif

#ifndef RECTANGULARPLATE_H
#include "RectangularPlate.h"
#endif

using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	printf("\n\nHello World\n\n");
	/*
	Cmplx z(1, 2);
	std::cout << "z=" << z << "\n\n";
	//////////////////////////////////////
	int N = 5;
	int f1 = ((N +1)+(N + 1));
	int f2 = (-(N + 1));
	
	RMatrix lMatrix(N,N);
	RVektor gVector(N);
	for (int c = 0; c < N; c++)
	{
		
		for (int r=0; r < N; r++)
		{
			
			if (c == r)
			{
				lMatrix[c][r] = f1;
			}
			if (abs(c - r) == 1)
			{
				lMatrix[c][r] = f2;

			}
			
		}	
	}

	for (int r = 1; r < N+1; r++)
	{
		//gVector[r] = (1.+(4*r*r+1./3)/((N+1)*(N+1)))/(N+1);
		gVector[r] = (1. / (N + 1)) * (1. + ((4 * (r*r) + (1. / 3)) / ((N + 1.)*(N + 1.))));
	}
	
	RVektor rMatrix;
	rMatrix = RealGaussElimination(lMatrix, gVector);

	cout << "fMatrix\n"<<lMatrix << "\n\n";
	cout << "gMatrix\n" << gVector << "\n\n";
	cout << "rMatrix\n" << rMatrix << "\n\n";
	
	
	



	//////////////////////////////////////
	RVektor aVec(1,2,3);
	RVektor bVec( 1,2,3);
	RVektor cVec( 7,8,9);
	
	
	RMatrix aMa(3, 3);
	aMa[0] = aVec;
	aMa[1] = bVec;
	aMa[2] = cVec;
	
	RMatrix xMa(5, 5);





	cout << "aVec=" << aVec << "\n\n";
	cout << aMa << "\n\n";*/

	//SmallSquare mySquare;
	//SquarePlate myPlate;
	RectangularPlate recPlate;

	
	
	//-----TOP---------------------

	double tpWidth = 0.5; //length of top plate is 2 "unit"?????????????
	int tpwRes = 10;
	int tplRes = 20;

	//-----BOTTOM------------------
	double bpWidth = 1;
	int bpwRes = 20;
	int bplRes = 20;
	
	//-----------------------------
	int sizeP1 = tpwRes*tplRes;
	int sizeP3 = bplRes*bpwRes;
	int dim = 2*(tpwRes*tplRes)+bpwRes*bplRes;
	double res = tpWidth/tpwRes;
	//-----------------------------
	//recPlate.TurnThemToMatrix(dim,tpwRes,tplRes,tpwRes,tplRes,bpwRes,bplRes);
	//-----------------------------

	//Construct object with plates, an text file of all coordinates of object will be created.
	recPlate.CreatePlate(tpWidth, tpwRes, tplRes, 5, 0., tpWidth, 0.0125);
	recPlate.CreatePlate(tpWidth, tpwRes, tplRes, 175, 0., tpWidth, 0.0125);
	recPlate.CreatePlate(bpWidth, bpwRes, bplRes, 0, 0, 0, 0);
	//----------------------------------
	recPlate.AssignVectors(dim, res);

	recPlate.CreateMatrix(dim);
	recPlate.CreateGmVector(dim, sizeP1, sizeP1, sizeP3);
	recPlate.FindChargeDistribution(dim);
	recPlate.FindTotalCharge(res,dim);
	recPlate.FindPotential(dim);
	
	/*
	myPlate.Create_U_Plate();
	myPlate.CreateMatrix();
	myPlate.CreateGmVector();
	myPlate.FindChargeDistribution();
	myPlate.CalculateTotalCharge();
	*/
	
	//myPlate.CreatePlate();
	/*
	myPlate.CreatePotMatrix();
	myPlate.CreateAVector();
	myPlate.CalculatePotential();
	*/
	std::cin.get();
	return 0;
}
