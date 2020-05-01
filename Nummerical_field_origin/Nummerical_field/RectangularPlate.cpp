#include "stdafx.h"
#include <iostream>
#include <fstream>
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
	oA = 90;
	oX = 0;
	oY = 0;
	oZ = 0;
	//point coordinate[6];

}

RectangularPlate::~RectangularPlate(void){}
void RectangularPlate::CreatePlate(const double plateWidth, const int wRes, const int lRes, const double oA, const double oX, const double oY, const double oZ)
{

	vPlate.SetDim(wRes*lRes);
	double plateLength = plateWidth*lRes / wRes;
	double res = plateWidth / wRes;
	//double ress = plateLength / lRes;
	double delX = res;
	double delY = (res * cos(oA*PI / 180));
	double delZ = (res * sin(oA*PI / 180));
	ofstream file;
	file.open("c:/software/files/recPlate.m");

	int iN = 0;
	for (int iX = 0; iX < lRes; iX++)
	{
		for (int iY = 1; iY <= wRes; iY++)
		{
			int iZ = iY;

			double x = delX*iX + 0.5*res + oX;
			double y = delY*iY + oY - 0.5*delY;
			double z = delZ*iZ + oZ - 0.5*delZ;
			vPlate[iN].MoveTo(x, y, z);
			vPlate[iN].SetSide(res);
			std::cout << iN << " " << vPlate[iN];
			file << vPlate[iN];
			iN++;

		}
	}
	file.close();
}