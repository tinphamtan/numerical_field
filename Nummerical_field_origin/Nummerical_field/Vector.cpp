// Vector.cpp
// Definition der Elementfunktionen von Vector

#include "stdafx.h"
#include <iostream>

typedef double Type;

#ifndef VEKTOR_H
#include "Vector.h"
#endif


// ErrHdlr ist eine statische Elementfunktion
// zu einer void Funktion mit char* Argument
void (*Vektor<Type>::ErrHdlr)(char*) = Vektor<Type>::DefaultErrHdlr;
char Vektor<Type>::ErrFlag = 0;
//<Type> ZERO;
//------------------------------------------------------------
Vektor<Type>::Vektor(void)
{
  Size = 1;

  a = new Type[1];

  a[0] = Type();

}
//------------------------------------------------------------
Vektor<Type>::Vektor( int n )
{
  if ( n <= 0 )
  {
    Error( "Vektor muss größer als Null sein" );
    n = 0;
  }
  Size = n;
  a = new Type[ Size ];
  for ( int i=0; i<Size; i++ )
                          a[i] = 0;
}
//------------------------------------------------------
Vektor<Type>::Vektor(const Vektor<Type>& v)
{
  Size = v.Size;
  a = new Type[Size];
  for ( int i=0; i<Size; i++ )
    a[i] = v.a[i];
}
//--------------------------------------------------------
Vektor<Type>& Vektor<Type>::operator=(const Vektor& v)
{
  if ( Size < v.Size )
  {
    delete  a;
    a = new Type[v.Size];
  }
  Size = v.Size;
  for ( int i=0; i<Size; i++ )
    a[i] = v.a[i];
  return *this;
}
//-----------------------------------------------------
Vektor<Type>& Vektor<Type>::operator/=(const Type& r)
{
  for ( int i=0; i<Size; i++ )
    a[i] /= r;
  return *this;
}
//-----------------------------------------------------
Vektor<Type>& Vektor<Type>::operator/=(const float& r)
{
  for ( int i=0; i<Size; i++ )
    a[i] /= r;
  return *this;
}
//-----------------------------------------------------
Vektor<Type>& Vektor<Type>::operator*=(const Type& r)
{
  for ( int i=0; i<Size; i++ )
    a[i] *= r;
  return *this;
}
//---------------------------------------------------------
Vektor<Type>& Vektor<Type>::operator+=(const Vektor<Type>& v)
{
  if ( v.Size != Size )
    Error( "Vektoren verschiedener Größen können nicht addiert werden" );
  else
    for ( int i=0; i<Size; i++ )
      a[i] += v.a[i];
  return *this;
}
//----------------------------------------------------------
Vektor<Type>& Vektor<Type>::operator-=(const Vektor<Type>& v)
{
  if ( v.Size != Size )
    Error( "Vektoren verschiedener Größen können nicht subtrahiert werden" );
  else
    for ( int i=0; i<Size; i++ )
      a[i] -= v.a[i];
  return *this;
}
//----------------------------------------------------------------
Vektor<Type> Vektor<Type>::operator-() const
{
  Vektor Result( Size );

  for ( int i=0; i<Size; i++ )
    Result.a[i] = -a[i];
  return Result;
}
//------------------------------------------------------------
/*Type Vektor<Type>::operator*(const Vektor<Type>& v) const
{
  Type Result = 0;
  if ( Size != v.Size )
    Error( "Vektoren verschiedener Größe können nicht multipliziert werden" );
  else
    for ( int i=0; i<Size; i++ )
      Result += a[i]*v.a[i];
  return Result;
}
*/
//---------------------------------------------------------------
Vektor<Type> Vektor<Type>::operator*(const Type& t) const
{
  Vektor<Type> Result = *this;
  Result *= t;
  return Result;
}
//---------------------------------------------------------------
Vektor<Type> operator*(const Type& t, const Vektor<Type>& v)
{ 
	return v * t; 
}
//---------------------------------------------------------------
Vektor<Type> Vektor<Type>::operator/(const Type& t) const
{
  Vektor<Type> Result = *this;
  Result /= t;
  return Result;
}
//---------------------------------------------------------------
Vektor<Type> Vektor<Type>::operator+(const Vektor<Type>& v) const
{
  Vektor Result = *this;
  Result += v;
  return Result;
}
//-------------------------------------------------------------
Vektor<Type> Vektor<Type>::operator-(const Vektor<Type>& v) const
{
  Vektor Result = *this;
  Result -= v;
  return Result;
}
//---------------------------------------------------------------
int Vektor<Type>::operator==(const Vektor<Type>& v) const
{
  int Result = 1;
  if ( Size != v.Size )
  {
    Error( "Vektoren verschiedener Größe können nicht verglichen werden" );
    Result = 0;
  }
  else
    for ( int i=0; i<Size && Result; i++ )
      if ( a[i] != v.a[i] )
        Result = 0;
  return Result;
}
/*
//--------------------------------------------------------------------
Type Vektor<Type>::Min() const
{
  Type Result = a[0];
  for ( int i=1; i<Size; i++ )
    if ( a[i] < Result )
      Result = a[i];
  return Result;
}

Type Vektor::Max() const
{
  Type Result = a[0];
  for ( int i=1; i<Size; i++ )
    if ( a[i] > Result )
      Result = a[i];
  return Result;
}

float Vektor::Mean() const
{
  float Result = 0;
  for ( int i=0; i<Size; i++ )
    Result += a[i];
  return Result / float(Size);
}
*/
//-----------------------------------------------
Vektor<Type> Vektor<Type>::operator/(const float& t) const
{
  Vektor<Type> Result = *this;
  Result /= t;
  return Result;
}
//------------------------------------------------------------------
/*
Vektor<Type> Vektor<Type>::operator+(const Vektor<Type>& v) const
{
  Vektor<Type> Result = *this;
  Result += v;
  return Result;
}
//-------------------------------------------------------------------
Vektor<Type> Vektor<Type>::operator-(const Vektor<Type>& v) const
{
  Vektor<Type> Result = *this;
  Result -= v;
  return Result;
}
//---------------------------------------------------------------
int Vektor<Type>::operator==(const Vektor<Type>& v) const
{
  int Result = 1;
  if ( Size != v.Size )
  {
    Result = 0;
  }
  else
    for ( int i=0; i<Size && Result; i++ )
      if ( a[i] != v.a[i] )
        Result = 0;
  return Result;
}*/

/*
Type Vektor::Min() const
{
  Type Result = a[0];
  for ( int i=1; i<Size; i++ )
    if ( a[i] < Result )
      Result = a[i];
  return Result;
}

Type Vektor::Max() const
{
  Type Result = a[0];
  for ( int i=1; i<Size; i++ )
    if ( a[i] > Result )
      Result = a[i];
  return Result;
}

float Vektor::Mean() const
{
  float Result = 0;
  for ( int i=0; i<Size; i++ )
    Result += a[i];
  return Result / float(Size);
}
*/
// Standard-Fehlerbehandlung
//void Vektor<Type>::DefaultErrHdlr( char *Msg )
//------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const Vektor<Type>& v)
{
  os << '(';
  int n = v.nElems() - 1;
  int i;
  for (i=0; i<n; i++ )
          os << v.Elem(i) << ' ';
  os << v.Elem(i) << ')';
  return os;
}
//------------------------------------------------------------
std::istream& operator>>(std::istream& is, Vektor<Type>& v)
{
  char c;
  int n = v.nElems();
  is >> c;  // lese '('
  for (int i=0; i<n; i++ )
    is >> v[i];
  is >> c;  // lese ')'
  return is;
}
//------------------------------------------------------------
/*
Vektor<Type>& Vektor<Type>::operator*=(const Type& r)
{
  for ( int i=0; i<Size; i++ )
    a[i] *= r;
  return Vektor<type> *const this;
}
*/