// Nummerical_field.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Vector.h"
#include <iostream>

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

#ifndef RECTANGULARPLATE_H
#include "RectangularPlate.h"
#endif

#ifndef SMALLSQUARE_H
#include "SmallSquare.h"
#endif

using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	printf("\n\nProgram is started and running...\n\n");
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
	
	//myPlate.CreatePlate();
	recPlate.CreatePlate(1,5,10,90,0,0,0);
	//myPlate.CreateMatrix();
	//myPlate.CreateGmVector();
	
	//myPlate.FindChargeDistribution();
	//myPlate.CalculateTotalCharge();	

	std::cin.get();
	return 0;
}
