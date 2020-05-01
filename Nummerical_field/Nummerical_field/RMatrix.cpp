
#include "stdafx.h"
#include <iostream>
#include "Physicon.h"
using namespace std;
#ifndef RMATRIX_H
#include "RMatrix.h"
#endif
extern int g_N;
extern double g_b;
extern double g_d;



//-------------------------------------------------------------
// Implementierung der Klasse RMatrix
//-------------------------------------------------------------
// Statische Elementfunktionen

void (*RMatrix::ErrHdlr)(char*) = RMatrix::DefaultErrHdlr;

char RMatrix::ErrFlag = 0;

//-----------------------------------------------------
// Konstruktor:
//
// Default-Konstrutor
RMatrix::RMatrix()
{
   Size = 1;
   v    = new RVektorPtr[Size];
   v[0] = new RVektor();
}
//--------------------------------------------
// RMatrix anlegen mit M Zeilen und N Spalten
RMatrix::RMatrix( int M, int N )
{
  //ASSERT(M > 0);
  Size = M;
  v = new RVektorPtr[ Size ];
  for (int i=0; i<Size; i++)
               v[i] = new RVektor( N );
}
//---------------------------------------------------
// Konstruktor: RMatrix kopieren

RMatrix::RMatrix( const RMatrix& m )
{
  Size = m.Size;
  v = new RVektorPtr[ Size ];
  for ( int i=0; i<Size; i++ )
    *(v+i) = new RVektor( m.Elem(i) );
}
//----------------------------------------------------
void RMatrix::SetDimension(int z, int s)
{
	//Vektoren löschen
	for(int i = 0; i < Size; i++)
	{
		   if(v[i] != NULL) 
			         delete v[i];
	}

	//Neukonstruktion
  Size = z;
  v = new RVektorPtr[ Size ];
  for ( int i=0; i<Size; i++)
               v[i] = new RVektor( s );
}
//----------------------------------------------------
// Destruktor

RMatrix::~RMatrix()
{
  for (int i=0; i<Size; i++)
    delete *(v+i);  // Vektoren löschen
  delete v;         // Zeiger löschen
}
//--------------------------------------------------------------
void RMatrix::AddToElement(REAL aVal,int z, int s)
{
	//ASSERT(z >= 0);
	//ASSERT(s >= 0);
	//ASSERT(z < Size);
	//ASSERT(s < v[z]->Size);
	v[z]->a[s] += aVal;
}
//--------------------------------------------------------------
void RMatrix::InsertMatrix(const RMatrix& aMa, int z, int s)
{
	//ASSERT( z == 0);//nur z= 0 und s = 0 implementiert
	//ASSERT( s == 0);

	//ASSERT(Rows() >= aMa.Rows());
	//ASSERT(Cols() >= aMa.Cols());

	for(int i = 0; i < aMa.Rows(); i++)
	{
		for(int k = 0; k < aMa.Cols(); k++)
		{
          v[i]->a[k] = aMa.v[i]->a[k];
		}
	}
}
//--------------------------------------------------------------
void RMatrix::InsertSpalte(const RVektor& aVec, int indx)
{
	int nCol = Cols();
	//ASSERT( indx < nCol);
	//ASSERT(Rows() == aVec.nElems());

	for(int i = 0; i < aVec.nElems(); i++)
	{
	   v[i]->a[indx] = aVec.a[i];
	}
}
//--------------------------------------------------------------
RVektor RealGaussElimination(const RMatrix& aMa,const RVektor& bVec)
{
	//ASSERT(aMa.Rows() == aMa.Cols());    //quadratische Matrix
	//ASSERT(aMa.Cols() == bVec.nElems()); //richtige Dimension des Vektors

	int nDim = aMa.Rows();
	int n    = nDim;
	int   i, k;
	int   kmax, j, m = n + 1;

	RVektor aResult(nDim);
	        aResult.SetZero();

	RMatrix nMa(nDim,nDim + 1);

	nMa.InsertMatrix(aMa,0,0);
	nMa.InsertSpalte(bVec, nDim);
//---------------------------------
/*    Funktion realisiert den Gauss-Algorithmus
		ma - Matrix der Dimension n x n+1
		n  - Dimension der Matrix
*/
	REAL  smax;
	REAL  pil, zz;
	REAL  det = 1.0;

	for( i=0 ; i<n ; i++ ) /* Groesstes Element der Spalte bestimmen */
	{
	    kmax = 0; smax = 0.;
	    for( k=i ; k<n ; k++ )
		{
		 if( fabs(nMa[k][i]) > smax)
			{
			   kmax = k; smax = fabs(nMa[k][i]);
			}
		 }
	  if( smax == 0)
		 {
  		    return aResult;
		 }
	  else
		 {
		 if(kmax != i)
			{            /*   Zeilentausch  */
			for( j=i ; j<m ; j++ )
			  {
			    zz           = nMa[i][j];
			    nMa[i][j]    = nMa[kmax][j];
			    nMa[kmax][j] = zz;
			  }
			  det = det*(-1);
			}
	  }

	/* Dividieren der Pilotzeile */
	pil = nMa[i][i];
	det = pil*det;
	for( j=i+1 ; j<m ; j++ )
		nMa[i][j] = nMa[i][j]/pil;

	/* Eliminationsschritt       */
	for( k=i+1 ; k<n ; k++ )
		{
		  zz = nMa[k][i];
		  for( j=i+1 ; j<n+1 ; j++ )
		  {
		    nMa[k][j] = nMa[k][j] - zz*nMa[i][j];
		  }
		}
	}
	/*  Rueckwaertseinsetzen     */
	for( i=n-1 ; i>=0 ; i-- )
	  {
	    for( k=i-1 ; k>=0 ; k-- )
		 {
		   nMa[k][n] = nMa[k][n] - nMa[k][i]*nMa[i][n];
		 }
	  }
//---------------------------------
    for(i = 0; i < n; i++)
	{
		aResult[i] = nMa[i][n];
	}
	return aResult;
}
//--------------------------------------------------------------
/*RVektor RMatrix::GaussAlgo(RVektor& inVec)
{
  if( Rows() != Cols())
    {
      Error("Fehler GaussAlgo: Keine quadratische  Matrix");
      return inVec;
    }
  if( Rows() != inVec.nElems())
    {
      Error("Fehler GaussAlgo: Matrixdimension != Vektordimension");
      return inVec;
    }

  int dim = Rows();
  RVektor outVec(dim);

  int iZeil, nSpalt;
  Cmplx *ma = NULL;
  ma = new Cmplx[dim*(dim + 1)];
  //Cmplx ma[6*9];
  if(ma != NULL)
   { // Uebernehmen der neuen Matrix
     for( iZeil = 0; iZeil < dim; iZeil++)
       {
	     for( nSpalt = 0; nSpalt < dim; nSpalt++)
	        {
	          //ma[iZeil*(dim + 1) + nSpalt] = in_ma[iZeil*maxDim + nSpalt];
	          ma[iZeil*(dim + 1) + nSpalt] = Elem(iZeil).Elem(nSpalt);
	        }
       }
     for( iZeil = 0; iZeil < dim; iZeil++)
       { //Uebernahme der Erregung
	        ma[iZeil*(dim + 1) + dim] = inVec[iZeil];
       }
     if(gauss(ma,dim))
      {
	     for( iZeil = 0; iZeil < dim; iZeil++)
	        { //Uebernahme des Ergebnisses
	           outVec[iZeil] = ma[iZeil*(dim + 1) + dim];
	        }
        delete ma;
	     return outVec;
      }
     delete ma;
     Error("Fehler GausAlgo: Determinante verschwindet");
     return inVec;
   }
  Error("Fehler GausAlgo: Nicht genEend Speicherplatz");
  return inVec;
}
//--------------------------------------------------------------
RVektor RMatrix::operator/(const RVektor& aVec)
{
   return this->GaussAlgo(aVec);
}*/
//--------------------------------------------------------------
RMatrix& RMatrix::operator=( const RMatrix& m2 )
{
  for ( int i=0; i<Size; i++ )
    delete *(v+i);
  delete v;

  Size = m2.Size;
  v = new RVektorPtr[ Size ];
  for (int  i=0; i<Size; i++ )
             v[i] = new RVektor( m2.Elem(i) );
  return *this;
}
//----------------------------------------------------
// I/O stream Operatoren
std::ostream& operator<<(std::ostream& os, const RMatrix& m)
{
  for ( int i=0; i<m.Rows(); i++ )
                    os << m.Elem(i) << "\n";
  return os;
}
//-----------------------------------------------------
std::istream& operator>>(std::istream& is, RMatrix& m)
{
  for ( int i=0; i<m.Rows(); i++ )
    is >> m[i];
  return is;
}
//-----------------------------------------------------
// mathematische Operatoren
RMatrix& RMatrix::operator*=(const REAL& r)
{
  for ( int i=0; i<Size; i++ )
    *(v[i]) *= r;
  return *this;
}
//----------------------------------------------------
/*RMatrix& RMatrix::operator*=(const Cmplx& r)
{
  for ( int i=0; i<Size; i++ )
    *(v[i]) *= r;
  return *this;
}*/
//--------------------------------------------------------
RMatrix& RMatrix::operator/=(const REAL& r)
{
  for ( int i=0; i<Size; i++ )
    *(v[i]) /= r;
  return *this;
}
//----------------------------------------------------------
/*RMatrix& RMatrix::operator/=(const Cmplx& r)
{
  for ( int i=0; i<Size; i++ )
    *(v[i]) /= r;
  return *this;
}*/
//---------------------------------------------------------

RMatrix& RMatrix::operator+=(const RMatrix& m2)
{
  //ASSERT(Rows() >= m2.Rows() );
  //ASSERT(Cols() >= m2.Cols() );
  for ( int i=0; i< m2.Rows(); i++ )
      (*this)[i] += m2.Elem(i);
  return *this;
}
//-----------------------------------------------------------

RMatrix& RMatrix::operator-=(const RMatrix& m2)
{
  if ( !EqualSize(*this,m2) )
    Error( "verschiedene RMatrixgrößen in operator-" );
  else
    for ( int i=0; i<Rows(); i++ )
      (*this)[i] -= m2.Elem(i);
  return *this;
}
//------------------------------------------------------------
RMatrix& RMatrix::operator*=(const RMatrix& m2)
{
  if ( Cols() != m2.Rows() )
    Error( "Fehler bei (A * B): A.Cols != B.Rows" );
  else
  {
    RMatrix Result( Rows(), m2.Cols() );
    for ( int i=0; i<Rows(); i++ )
      for ( int j=0; j<m2.Cols(); j++ )
      {
        REAL Sum = 0;
        for ( int k=0; k<Cols(); k++ )
        Sum += (*v[i])[k] * m2.Elem(k).Elem(j);
        Result[i][j] = Sum;
      }
    *this = Result;
  }
  return *this;
}
//--------------------------------------------------------
// Konstanten
//--------------------------------------------------------
RVektor RMatrix::operator*(const RVektor& aVec) const
{
  RVektor Result(Rows() );
  if ( Cols() !=  aVec.nElems() )
    Error( "Fehler bei (A * B): A.Cols != B.Rows" );
  else
  {
    for ( int i=0; i<Rows(); i++ )
      {
        REAL Sum = 0;
        for ( int k=0; k<Cols(); k++ )
        Sum += (*v[i])[k] * aVec.Elem(k);
        Result[i] = Sum;
      }
    return Result;
  }
  return Result;
}
//---------------------------------------------------
/*
CVektor RMatrix::operator*(const CVektor& aVec) const
{
  CVektor Result(Rows() );
  if ( Cols() !=  aVec.nElems() )
    Error( "Fehler bei (A * B): A.Cols != B.Rows" );
  else
  {
    for ( int i=0; i<Rows(); i++ )
      {
        Cmplx Sum = 0.;
        for ( int k=0; k<Cols(); k++ )
          {
             Cmplx aC = aVec.Elem(k);
             Sum +=  aC*((*v[i])[k]);
          }
        Result[i] = Sum;
      }
    return Result;
  }
  return Result;
}
*/
//---------------------------------------------------
RMatrix RMatrix::operator*(const REAL& r) const
{
  RMatrix Result = *this;
  Result *= r;
  return Result;
}
//------------------------------------------------------

RMatrix RMatrix::operator/(const REAL& r) const
{
  RMatrix Result = *this;
  Result /= r;
  return Result;
}
//----------------------------------------------------

RMatrix RMatrix::operator+(const RMatrix& m2) const
{
  RMatrix Result = *this;
  Result += m2;
  return Result;
}
//-----------------------------------------------------

RMatrix RMatrix::operator-(const RMatrix& m2) const
{
  RMatrix Result = *this;
  Result -= m2;
  return Result;
}
//-----------------------------------------------------

RMatrix RMatrix::operator*(const RMatrix& m2) const
{
  RMatrix Result = *this;
  Result *= m2;
  return Result;
}
//--------------------------------------------------------

RMatrix RMatrix::operator-() const
{
  RMatrix Result = *this;
  for ( int i=0; i<Rows(); i++ )
    for ( int j=0; j<Cols(); j++ )
      Result[i][j] = -((*v[i])[j]);
  return Result;
}
//-----------------------------------------------------
// Verschiedene Funktionen und Operatoren

void RMatrix::DefaultErrHdlr( char *Msg )
{
	std::cout << "RMatrix Fehler: " << Msg << "\n";
}
//-----------------------------------------------------

int RMatrix::EqualSize(const RMatrix& m1,
		      const RMatrix& m2)
{
  return m1.Cols() == m2.Cols()
      && m1.Rows() == m2.Rows();
}
/*
Type RMatrix::Min() const
{
  Type Result = v[0]->Min();
  for ( int i=1; i<Size; i++ )
  {
    Type Temp = v[i]->Min();
    if ( Temp < Result )
      Result = Temp;
  }
  return Result;
}

Type RMatrix::Max() const
{
  Type Result = v[0]->Max();
  for ( int i=1; i<Size; i++ )
  {
    Type Temp = v[i]->Max();
    if ( Temp > Result )
      Result = Temp;
  }
  return Result;
}

Type RMatrix::Mean() const
{
  Type Result = 0;
  for ( int i=0; i<Size; i++ )
    Result += v[i]->Mean();
  Result /= Type(nElems());
  return Result;
}
*/
//-----------------------------------------------------
RMatrix RMatrix::Identity( int n )
{
  RMatrix Result( n, n );
  for ( int i=0; i<n; i++ )
    Result[i][i] = 1;
  return Result;
}
