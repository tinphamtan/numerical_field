#include "stdafx.h"
#include <iostream>
#include <fstream> //ability to create and write text file

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

using namespace std;

void (*Vektor<SmallSquare>::ErrHdlr)(char*) = Vektor<SmallSquare>::DefaultErrHdlr;
char Vektor<SmallSquare>::ErrFlag = 0;

//-----------------------------------
SquarePlate::SquarePlate(void)
{
	nX   = 55*3;  // Length-----dimension of the plate
	nY   = 55;	// width
	nZ   = 5;	//height
	slide = 70; // total number of square plate in one slide
	nN   = slide*nX;  // what is this?
	side = 1;
	height = 1;
	swidth = 1;
	width = swidth*11;
	length = (swidth / 5) * nX;
	potPlate = 1.0;
	totalCharge = 0.;
	potSide = 60;
	potSideLength = 20;
}
//-----------------------------------

SquarePlate::~SquarePlate(void)
{
}
//---------------------------------------
void SquarePlate::Create_U_Plate( )
{

	
	vPlate.SetDim(nN);
	
	double delY = side / nZ;
	double delZ = delY;
	double delX = delY;
	int iN = 0;
	int iY = 0;
	int iZ = 0;
	ofstream file;							//helps create and write to file
	
	file.open("c:/software/files/U_Structure.m");			//create file with the name
	
	for (int iX = 1; iX <= nX; iX++)
	{

		for(iY =1; iY <= nY; iY++)    //z=0, bottom y side
		{
			iZ = 0;
			double x = delX*iX;
			double y = delY*iY;
			double z = delZ*iZ;
			vPlate[iN].MoveTo(x,y,z);
			vPlate[iN].SetSide(delX);
			std::cout << vPlate[iN];
			file << vPlate[iN];				//write data to the file
			iN++;
			
		}
		
		for (iZ =1; iZ <= nZ; iZ++)   //y=nY, left z side
		{
			iY = 56;
			double x = delX*iX;
			double y = delY*iY-delY/2;
			double z = delZ*iZ-delZ/2;
			vPlate[iN].MoveTo(x, y, z);
			vPlate[iN].SetSide(delX);
			std::cout << vPlate[iN];
			file << vPlate[iN];				
			iN++;			
		}
		
		for (iY = iY - 26; iY >= 26; iY--)	//z=nZ, top y side
		{
			
			iZ = 6;
			double x = delX*iX;
			double y = delY*iY;
			double z = delZ*iZ - delZ;
			vPlate[iN].MoveTo(x, y, z);
			vPlate[iN].SetSide(delX);
			std::cout << vPlate[iN];
			file << vPlate[iN];				
			iN++;			
		}
		
		

		for (iZ = iZ-1; iZ >=1; iZ--)	//y=0, right z side
		{
			iY = 0;
			double x = delX*iX;
			double y = delY*iY + delY / 2;
			double z = delZ*iZ - delZ / 2;
			vPlate[iN].MoveTo(x, y, z);
			vPlate[iN].SetSide(delX);
			std::cout << vPlate[iN];
			file << vPlate[iN];				
			iN++;
			

		}
		

	}
	file.close();							//close the file

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
	gmVec.SetZero();
	gmVec.SetDim(nN);
	int q = 0;
	for (int i = 0; i < nX; i++)
	{
		for (int j = 0; j < slide; j++)
		{
			if (j == 61 || j == 62 || j == 63 || j == 64 || j == 65 )
			{
				gmVec[q] = potPlate;
				q++;
			}
			
			else
			{
				gmVec[q] = 0;
				q++;
			}
			
		}
	}
	std::cout << gmVec;
}
void SquarePlate::FindChargeDistribution()
{
	ofstream file;
	file.open("c:/software/files/U_Charge.m");
	chargeVec.SetDim(nN);
	
	chargeMatrix.SetDimension(nX,slide);						//matrix to store the charge density
	chargeVec = RealGaussElimination(myMatrix, gmVec);
	std::cout << "\n\n Charge Distribution \n\n";
	int iN = 0;
	for (int i = 0; i < nX; i++)
	{
		//std::cout <<"\n"<< i<<"			";
		//file << "\n" << "		";
		for (int j = 0; j < slide; j++)
		{
				chargeMatrix[i][j] = chargeVec[iN];
				//std::cout << chargeVec[iN] << "			";
				//file << chargeVec[iN] << "			";
				iN++;
				
		}
		
	}
	file<<"\n"<< chargeMatrix;
	std::cout << chargeMatrix;
	file.close();
}


void SquarePlate::CalculateTotalCharge()
{
	std::cout << "\n\nTotalCharge: ";
	int iN = 0;
	
	for (int i = 0; i < nX; i++){
		for (int j = 0; j < nY; j++)
		{

			totalCharge += chargeVec[iN]*((side/nX)*(side/nX));

			iN++;
		}
	}
	
	std::cout << totalCharge<<"\n";
}

void SquarePlate::CreatePlate()
{
	ofstream file;							//helps create and write to file

	file.open("c:/software/files/potPlate.m");

	pPlate.SetDim(potSide*potSide);

	double delX = potSideLength / potSide;
	double delY = delX;
	double delZ = delX;
	int iN = 0;
	float iX = 16.5;
	int iY = 0;
	int iZ = 0;
	for (iY = 0; iY < potSide; iY++)
	{
		for (iZ = 0; iZ < potSide; iZ++)
		{
			double x = 0.2*iX;
			double y = delY*iY-1;
			double z = delZ*iZ-1;
			pPlate[iN].MoveTo(x, y, z);
			pPlate[iN].SetSide(delX);
			std::cout << iN << " " << pPlate[iN];
			file << pPlate[iN];
			iN++;
			//std::cout << "\n" << potVec;
		}

	}

	file.close();

}
void SquarePlate::CreatePotMatrix()
{
	std::cout << "\n\nCreating Potential Matrix\n\n";
	potMatrix.SetDimension(potSide*potSide, nN);
	for (int i = 0; i < potSide*potSide; i++)
	{
		for (int j = 0; j < nN; j++)
		{
			potMatrix[i][j] = pPlate[i].GetPotential(vPlate[j],chargeVec[j]);
		}
	}
	std::cout << potMatrix;
}

void SquarePlate::CreateAVector()
{
	std::cout << "\n\ngmVector\n\n";
	gmVec.SetDim(nN);
	int q = 0;
	for (int i = 0; i < nN; i++)
	{
		gmVec[i] = potPlate;
	}
	std::cout << gmVec;
}
void SquarePlate::CalculatePotential()
{
	ofstream file;
	file.open("c:/software/files/Potential.m");
	potentialVec.SetDim(potSide*potSide);
	potentialVec = potMatrix * gmVec;
	std::cout << "\n\n Charge Distribution \n\n";
	int iN = 0;
	for (int i = 0; i < potSide; i++)
	{
		std::cout << "\n" << "		";
		file << "\n" << "		";
		for (int j = 0; j < potSide; j++)
		{
			
			std::cout << potentialVec[iN] << "	";
			file << potentialVec[iN] << "		";
			iN++;

		}

	}
	
	file.close();
}