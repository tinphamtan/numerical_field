#ifndef CVEKTOR_H
#define CVEKTOR_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#ifndef  CMPLX_H
#include "cmplx.h"
#endif

#ifndef VEKTOR_H
#include "Vector.h"
#endif

//-------------------------------------------------------------
class CVektor : public Vektor<Cmplx>
{

public:
	 CVektor(): Vektor<Cmplx>(1){};
     CVektor(int n): Vektor<Cmplx>(n){};
     CVektor(const RVektor &aVec);
     CVektor(const CVektor &aVec):Vektor<Cmplx>(aVec){};
     CVektor(const Vektor<Cmplx> &aVec):Vektor<Cmplx>(aVec){};
     CVektor(Cmplx a, Cmplx b):Vektor<Cmplx>(a,b){};
     CVektor(Cmplx a, Cmplx b, Cmplx c):Vektor<Cmplx>(a,b,c){};

     CVektor  operator=(const Vektor<Cmplx>& aVec);
     CVektor  operator=(const CVektor& aVec);
#ifndef MY_KONSOL_APP
	 void FFTRealFnkt(const RVektor &aVec);//FFT-Transformierte reeller Abtastwerte
	 void FFTRealFnkt1(const RVektor & aVec);//FFT mit komplexer Transformation
	 void FFTTestFnkt1(const RVektor & aVec);//FFT mit komplexer Transformation
	 RVektor InversFFT1();//Rücktransformation
	 void DFTRealFnkt(const RVektor & aVec); //Diskrete Fourieranalyse einer reellen Funktion
	 void Serialize(CArchive& ar);//Ins Archive
#endif
	 void SetZero();//Vektor auf Null setzen
};
//-----------------------------------------------------
// Stream Ausgabe

std::ostream& operator<<(std::ostream& os, const CVektor& m);

// Stream Eingabe

std::istream& operator>>(std::istream& is, CVektor& m);

#endif
