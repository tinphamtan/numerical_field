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

#ifndef POTENTIAL_H
#include "Potential.h"
#endif

#ifndef RMATRIX_H
#include "RMatrix.h"
#endif

using namespace std;

void(*Vektor<SmallSquare>::ErrHdlr)(char*) = Vektor<SmallSquare>::DefaultErrHdlr;
char Vektor<SmallSquare>::ErrFlag = 0;

//-----------------------------------