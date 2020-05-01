#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <algorithm>
//#include "math.h""
#include "physicon.h"


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

#ifndef RECTANGULARPLATE_H
#include "RectangularPlate.h"
#endif

#ifndef RMATRIX_H
#include "RMatrix.h"
#endif



void(*Vektor<SmallSquare>::ErrHdlr)(char*) = Vektor<SmallSquare>::DefaultErrHdlr;
char Vektor<SmallSquare>::ErrFlag = 0;


RectangularPlate::RectangularPlate(void)
{

	plateLengh = 2;
	plateWidth = 1;
	wRes = 10;
	lRes = 20;
	N = 0;
	oA = 90;
	oX = 0;
	oY = 0;
	oZ = 0;
}

RectangularPlate::~RectangularPlate(void){}
void RectangularPlate::CreatePlate(const double plateWidth, const int wRes, const int lRes, const double oA, const double oX, const double oY, const double oZ)
{
	vObject.SetDim(wRes*lRes);
	double plateLength = plateWidth*lRes / wRes;
	double res = plateWidth / wRes;
	//double ress = plateLength / lRes;
	double delX = res;
	double delY = (res * cos(oA*PI/180));
	double delZ = (res * sin(oA*PI/180));
	ofstream file;
	
	//remove("c:/software/files/recPlate.m");
	file.open("c:/software/files/recPlate.m", std::ios::app);

	int iN = 0;
	for (int iX = 0; iX < lRes; iX++)
	{
		for (int iY = 1; iY <= wRes; iY++)
		{
			int iZ = iY;

			double x = delX*iX + 0.5*res + oX;
			double y = delY*iY + oY - 0.5*delY;
			double z = delZ*iZ + oZ - 0.5*delZ;
			vObject[iN].MoveTo(x, y, z);
			vObject[iN].SetSide(res);
			std::cout << iN << " " << vObject[iN];
			file << vObject[iN];
			iN++;


		}
	}
	file.close();
}
//construct object with three rectangular plates
void RectangularPlate::ConstructObject()
{
	int dim = 10;
	double res = 1.0;
	CreatePlate(2, 2, 2, 45, 0, 2, 2);
	CreatePlate(2, 2, 2, 135, 0, 2, 2);
	CreatePlate(4, 4, 2, 0, 0, 0, 0);
	
	AssignVectors(dim, res);
}
void RectangularPlate::AssignVectors(int dim, double res)
{
	ofstream file;

	file.open("c:/software/files/recPlate_copy.m", std::ios::app);
	cObject.SetDim(dim);
	//double res = 1;
	std::fstream myfile("c:/software/files/recPlate.m", std::ios_base::in);
	double a;
	double b;
	double c;
	int iN = 0;
	while (myfile >> a >> b >> c)
	{
		/*
		printf("a: %f ", a);
		printf("b: %f ", b);
		printf("c: %f ", c);
		printf("\n");
		*/
		//std::cout << cObject[iN];
		cObject[iN].MoveTo(a, b, c);
		cObject[iN].SetSide(res);
		file << cObject[iN];
		//std::cout << cObject[iN];
		iN++;
	}
	file.close();

}
void RectangularPlate::CreateMatrix(int size)
{
	std::cout << "\nCreate Matrix\n";
	N = size;
	myMatrix.SetDimension(N, N);
	
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			myMatrix[i][j] = cObject[i].CouplingCoefficient(cObject[j]);
			//cout << myMatrix[i][j]<<"\n";
		}
	}
	/*
	for (int i = 0; i < 0.5*N; i++)
	{
		for (int j = 0; j < 0.5*N; j++)
		{
			// L11 no Plus
			myMatrix[i][j] = cObject[i].CouplingCoefficient(cObject[j]);
		}

		for (int j = 0.5*N; j <N; j++)
		{
			// L12
			myMatrix[i][j] = cObject[i].CouplingCoefficient_plus(cObject[j]);
		}
	}
	for (int i = 0.5*N; i <N; i++)
	{
		for (int j = 0; j < 0.5*N; j++)
		{
			// L21
			myMatrix[i][j] = cObject[i].CouplingCoefficient_plus(cObject[j]);
		}
	
		for (int j = 0.5*N; j <N; j++)
		{
			// L122 no Plus
			myMatrix[i][j] = cObject[i].CouplingCoefficient(cObject[j]);
		}
	}
	*/
	//std::cout << myMatrix;
	
}

void RectangularPlate::CreateGmVector(int size,int sizeP1, int sizeP2, int sizeP3)
{
	N = size;
	std::cout << "\n\ngmVector\n\n";
	gmVec.SetZero();
	gmVec.SetDim(N);
	int q = 0;
	for (int i = 0; i < sizeP1; i++)
	{
		gmVec[q] = 0.5;
		q++;
	}
	for (int i = 0; i < sizeP2; i++)
	{
		gmVec[q] = 0.5;
		q++;
	}
	for (int i = 0; i < sizeP3; i++)
	{
		gmVec[q] = -0.5;
		q++;
	}
	std::cout << gmVec;
}
void RectangularPlate::FindChargeDistribution(int chargeVecSize)
{
	N = chargeVecSize;
	ofstream file;
	file.open("c:/software/files/U_Charge.m");
	chargeVec.SetDim(N);
	chargeVec.SetZero();

	//chargeMatrix_top.SetDimension(5, 10);						//matrix to store the charge density
	chargeVec = RealGaussElimination(myMatrix, gmVec);
	std::cout << "\n\n Charge Distribution \n\n"; 
	int iN = 0;
	for (int i = 0; i < N; i++)
	{
		file << chargeVec[iN] << "\n";
		iN++;
	}
	file.close();
}

void RectangularPlate::FindTotalCharge(double sideOfSmallSquare, int dim)
{
	std::cout << "\n\nTotalCharge: \n";
	totalCharge = 0.;

	int iN = 0;
	double side = sideOfSmallSquare;
	for (int i = 0; i < dim/2; i++){
		totalCharge += chargeVec[iN]*(side)*(side);
		iN++;
		cout << i << "P2: " << totalCharge << "\n";
		
	}
	std::cout << "P1: " << totalCharge << "\n";
	totalCharge = 0.;
	for (int i = dim / 2; i < dim; i++){
		totalCharge += chargeVec[iN] * (side)*(side);
		iN++;
		cout << i << "\n";
	}
	std::cout <<"P2: "<< totalCharge << "\n";
}

void RectangularPlate::CreateYZPlate(int yRes, int zRes, double yLength, double oA , double oX, double oY, double oZ)
{
	ofstream file;							//helps create and write to file

	file.open("c:/software/files/potPlate.m");

	pPlate.SetDim(yRes*zRes);

	double res = yLength / yRes;
	double delX = (res * cos(oA*PI / 180));
	double delY = res;
	double delZ = (res * sin(oA*PI / 180));
	int iN = 0;
	
	for (int iY = 0; iY < yRes; iY++)
	{
		for (int iZ = 0; iZ < zRes; iZ++)
		{
			int iX = iZ;
			double x = delX*iX + oX - 0.5*delX;
			double y = delY*iY + oY - 0.5*delY;
			double z = delZ*iZ + oZ - 0.5*delZ;
			pPlate[iN].MoveTo(x, y, z);
			pPlate[iN].SetSide(res);
			std::cout << iN << " " << pPlate[iN];
			file << pPlate[iN];
			iN++;
			//std::cout << "\n" << potVec;
		}

	}

	file.close();

}

void RectangularPlate::CreatePotMatrix(int yRes, int zRes, int sizeOfObject)
{
	std::cout << "\n\nCreating Potential Matrix\n\n";
	N = sizeOfObject;
	potMatrix.SetDimension(yRes*zRes, N);
	for (int i = 0; i < yRes*zRes; i++)
	{
		for (int j = 0; j < N; j++)
		{
			potMatrix[i][j] = pPlate[i].GetPotential(cObject[j], chargeVec[j]);
		}
	}
	std::cout << potMatrix;
}
void RectangularPlate::CreateAVector(int size)
{
	std::cout << "\n\ngmVector\n\n";
	N = size;
	defVec.SetDim(N);
	int q = 0;
	for (int i = 0; i < N; i++)
	{
		defVec[i] = 1;
	}
	std::cout << defVec;
}

void RectangularPlate::CalculatePotential(int yRes, int zRes)
{
	ofstream file;
	file.open("c:/software/files/Potential.m");
	potentialVec.SetDim(yRes* zRes);
	potentialVec = potMatrix * defVec;
	std::cout << "\n\n Charge Distribution \n\n";
	int iN = 0;
	for (int i = 0; i < zRes; i++)
	{
		file << "\n" << "			";
		std::cout << "\n" << "			";
		for (int j = 0; j < yRes; j++)
		{

			std::cout << potentialVec[iN] << "	";
			file << potentialVec[iN] << "			";
			iN++;

		}
		
	}

	file.close();
}

void RectangularPlate::TurnThemToMatrix(int dimOfObject, int dimP1w, int dimP1l, int dimP2w, int dimP2l, int dimP3w, int dimP3l)
{
	/*
	int dimP1w = 1;
	int dimP1l = 1;
	int dimP2w = 1;
	int dimP2l = 1;

	int dimP3w = 2;
	int dimP3l = 2;
	*/
	int iN = 0;
	//check for missmatch
	if ((dimP1w * dimP1l + dimP2w * dimP2l + dimP3w * dimP3l) != dimOfObject)
	{
		cout << "Dimention not match";
		throw std::exception();
	}
	double tempvar;
	ofstream fileP1;
	ofstream fileP2;
	ofstream fileP3;
	fileP1.open("c:/software/files/P1.m");
	fileP2.open("c:/software/files/P2.m");
	fileP3.open("c:/software/files/P3.m");

	chargeVec.SetZero();
	chargeVec.SetDim(dimOfObject);
	std::fstream myfile("c:/software/files/U_Charge.m", std::ios_base::in);
	while (myfile >> tempvar)
	{
		
		chargeVec[iN] = tempvar;
		iN++;
	}
	iN = 0;
	for (int i = 0; i < dimP1w+dimP2w; i++)
	{
		fileP1 << "\n" << "		";
		for (int j = 0 ; j < dimP1l; j++)
		{
			//myfile >> tempvar;
			fileP1 << chargeVec[iN] << "		";
			iN++;
		}
		
	}
	/*	
	for (int i = 0; i < dimP2w; i++)
	{
		fileP2 << "\n" << "		";
		for (int j = 0; j < dimP2l; j++)
		{
			//myfile >> tempvar;
			fileP2 << chargeVec[iN] << "		";
			iN++;
		}
	}*/

	for (int i = 0; i <dimP3w; i++)
	{
		fileP3 << "\n" << "		";
		for (int j = 0; j < dimP3l; j++)
		{
			//myfile >> tempvar;
			fileP3 << chargeVec[iN] << "		";
			iN++;
		}
	}
	fileP1.close();
	fileP2.close();
	fileP3.close();
}

void RectangularPlate::FindPotential(int sizeOfObject)
{
	int yRes = 10;
	int zRes = 10;
	double yLength = 3.;
	double oA = 90.;
	double oX = 0.6;
	double oY = -0.75;
	double oZ = -0.75;
	CreateYZPlate(yRes, zRes, yLength, oA, oX, oY, oZ);
	CreatePotMatrix(yRes, zRes, sizeOfObject);
	CreateAVector(sizeOfObject);
	CalculatePotential(yRes, zRes);
}