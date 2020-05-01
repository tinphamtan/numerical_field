#include "stdafx.h"

#ifndef RVEKTOR_H
#include "RVektor.h"
#endif

#ifndef SMALLSQUARE_H
#include "SmallSquare.h"
#endif

#ifndef VEKTOR_H
#include "Vector.h"
#endif

#ifndef RMATRIX_H
#include "RMatrix.h"
#endif

#ifndef UPLATE_H
#include "UPlate.h"
#endif

using namespace std;

void(*Vektor<SmallSquare>::ErrHdlr)(char*) = Vektor<SmallSquare>::DefaultErrHdlr;
char Vektor<SmallSquare>::ErrFlag = 0;

//-----------------------------------------------
UPlate::UPlate(void)
{

}

