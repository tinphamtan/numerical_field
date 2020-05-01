#ifndef SMALLSQUARE_H
#define SMALLSQUARE_H
#pragma once
class SmallSquare
{
private:
	RVektor rCenter; //Center of the square in meter
	double  side;    //Lenght of a side in meter


public:
	SmallSquare(void);
	~SmallSquare(void);

	void MoveTo(const RVektor& aVec){ rCenter = aVec; }; // Moves Center to aVec
	void MoveTo(const double aX, const double aY, const double aZ);     // Moves in the x-y-z plane
	void SetSide(const double aSide){ side = aSide; };   // Set the side
	double GetSide() const { return side; };			//Get the side of the small square
	double CouplingCoefficient(const SmallSquare& aSq);//calculate the coupling coefficient
	bool IsSelfCoupling(const SmallSquare & aSq);
	double SelfCoupling();
	double DistanceCoupling(const SmallSquare & aSq);
	RVektor GetCenter(const SmallSquare& a);
	////////////////////
	RVektor GetCenter()const { return rCenter; };
	double GetPotential(const SmallSquare& aSq, double chargeVal);
	double DistancePotential(const SmallSquare & aSq, double chargeVal);
};

std::ostream& operator <<(std::ostream& os, const SmallSquare& a);
#endif

