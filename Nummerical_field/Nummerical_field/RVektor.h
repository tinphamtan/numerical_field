#ifndef RVEKTOR_H
#define RVEKTOR_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifndef  CMPLX_H
#include "cmplx.h"
#endif

#ifndef VEKTOR_H
#include "Vector.h"
#endif

class RVektor : public Vektor<REAL>
{

public:
	 RVektor();
     RVektor(int n);
     RVektor(const RVektor &aVec):Vektor<REAL>(aVec){};
     RVektor(REAL a, REAL b):Vektor<REAL>(a,b){};
     RVektor(REAL a, REAL b, REAL c):Vektor<REAL>(a,b,c){};

	 void SetZero();//Vektor auf Null setzen
	 REAL Norm();   //Betrag des Vektors
	 void SetX(const double aX){a[0] = aX;};
	 void SetY(const double aY){a[1] = aY;};
	 void SetZ(const double aZ){a[2] = aZ;};

#ifndef MY_KONSOL_APP
	 void Serialize(CArchive& a);
#endif

	 RVektor& operator/=(const double aDouble);
	 RVektor& operator-=(const RVektor& aVec);
     RVektor& operator*=(const double&);
     RVektor& operator+=(const RVektor&);    // Addition eines Vektors: v += w


	 //Operatoren die nicht verändern
	 RVektor operator-(const RVektor& aVec);
     RVektor operator*(const double&) const;


      // Vektor+Vektor: v+w
    RVektor  operator+(const RVektor& w) const;

      // Vektor-Vektor: v-w
    RVektor  operator-(const RVektor& w) const;

     //RVektor  operator=(const Vektor<REAL>& aVec);
     //RVektor  operator=(const RVektor& aVec);
};
//-----------------------------------------------------
// Stream Ausgabe

std::ostream& operator<<(std::ostream& os, const RVektor& m);

// Stream Eingabe

std::istream& operator>>(std::istream& is, RVektor& m);

#endif