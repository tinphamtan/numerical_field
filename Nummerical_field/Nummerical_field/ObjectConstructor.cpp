#include "stdafx.h"
#include <iostream>
#include <fstream> //ability to create and write text file

#ifndef OBJECTCONSTRUCTOR_H
#include "ObjectConstructor.h"
#endif

using namespace std;

void(*Vektor<SmallSquare>::ErrHdlr)(char*) = Vektor<SmallSquare>::DefaultErrHdlr;
char Vektor<SmallSquare>::ErrFlag = 0;

ObjectConstructor::ObjectConstructor(void)
{
	N = 30;
	nX = 5;
	nY = 6;
	unitLength = 1.;

}

ObjectConstructor::~ObjectConstructor(void){

}

void ObjectConstructor::MakeRectangle()
{

	ofstream file;											//helps create and write to file
	file.open("c:/software/files/U_Structure.m");			//create file with the name
	Plate.SetDim(N);
	int iN = 0;
	int iY = 0;
	int iZ = 0;
	for (int iX = 0; iX <= nX; iX++)
	{
		for (int iY = 0; iY <= nY; iY++)
		{
			//x, y are coordinate of the center of every small square element
			double x = iX*unitLength + unitLength / 2;
			double y = iY*unitLength + unitLength / 2;
			double z = 0;
			Plate[iN].MoveTo(x, y, z);
			Plate[iN].SetSide(unitLength);
			std::cout << Plate[iN];
			file << Plate[iN];				//write data to the file
			iN++;
		}
	}
}