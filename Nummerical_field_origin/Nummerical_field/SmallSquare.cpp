#include "stdafx.h"
#include <iostream>
#include "Physicon.h"
#ifndef RVEKTOR_H
#include "RVektor.h"
#endif

#ifndef SMALLSQUARE_H
#include "SmallSquare.h"
#endif
#

//---------------------------------------
SmallSquare::SmallSquare(void)
{
	rCenter = RVektor(0., 0., 0.);
	side = 0.;
}
//---------------------------------------
SmallSquare::~SmallSquare(void)
{
}
//--------------------------------------------
void SmallSquare::MoveTo(const double aX, const double aY, const double aZ)
{
	rCenter.SetX(aX);
	rCenter.SetY(aY);
	rCenter.SetZ(aZ);
}
//--------------------------------------------

double SmallSquare::CouplingCoefficient(const SmallSquare& aSq)
{
	if (IsSelfCoupling(aSq)) return SelfCoupling();
	else return DistanceCoupling(aSq);
}

bool SmallSquare::IsSelfCoupling(const SmallSquare& aSq)
{
	return (rCenter == aSq.rCenter);
}

double SmallSquare::SelfCoupling()
{
	double retVal = (side) / (PI*EPS_0);
	retVal *= log(sqrt(2.) + 1);
	return retVal;
}
double SmallSquare::DistanceCoupling(const SmallSquare& aSq)
{
	RVektor disVec = rCenter - aSq.GetCenter();
	double dis = disVec.Norm();
	double retVal = side*side;
	retVal /= 4.*PI*EPS_0*dis;
	return retVal;
}



//--------------------------------------------
std::ostream& operator <<(std::ostream& os, const SmallSquare& aSq)
{
	//os << '(';
	//os << aSq.GetSide() << ",";     // not printing out the side and brackets
	os << aSq.GetCenter() << "\n";
	return os;
}



/////////////////////////--------------------------------------------
RVektor SmallSquare::GetCenter(const SmallSquare& aSq)
{

	return  aSq.GetCenter();
}
//--------------------------------------------
double SmallSquare::GetPotential(const SmallSquare& aSq, double chargeVal)
{
	if (IsSelfCoupling(aSq)) return 0;
	else return DistancePotential(aSq, chargeVal);
}

double SmallSquare::DistancePotential(const SmallSquare& aSq, double chargeVal)
{
	RVektor disVec = rCenter - aSq.GetCenter();
	double dis = disVec.Norm();

	chargeVal /= 4.*PI*EPS_0*dis;
	return chargeVal;
}