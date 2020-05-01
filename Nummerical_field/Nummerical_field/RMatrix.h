#ifndef RMATRIX_H
#define RMATRIX_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifndef RVEKTOR_H
#include "RVektor.h"
#endif

typedef RVektor *RVektorPtr;

class RMatrix
{
  private:
    RVektorPtr *v;
    int Size;


      // Statische (Klassen-) Variablen
      // ------------------------------

    static char ErrFlag;
    static void (*ErrHdlr)(char*);
    static void DefaultErrHdlr(char*);
    static int EqualSize(const RMatrix& m1,
                         const RMatrix& m2);


  public:
	RMatrix();
      // Anlegen einer Matrix mit M Zeilen (Rows) und N Spalten (Cols)
	RMatrix( int M, int N );
      // Copy Konstruktor
    RMatrix( const RMatrix& );
      // Destruktor
    ~RMatrix();

    friend class CMatrix;

    //Lösung eines Gleichungssystems
    //CVektor GaussAlgo(CVektor& inVec);

      // Zuweisungs-Operator
    RMatrix& operator=(const RMatrix& m2);

    // Index Operator, für Schreib- und Lesezugriffe
    RVektor& operator[](int row) const
    {
       if ( row >= 0 && row < Size )
	                     return *v[row];
      else
      {
	       Error( "ungültiger Zeilen-Index" );
	       return *v[0];
      }
    }
    // Nur-Lese Index Operator
    const RVektor& Elem(int row) const
    {
      if ( row >= 0 && row < Size )
	                        return *v[row];
      else
      {
	              Error( "ungültiger Zeilen-Index" );
	              return *v[0];
      }
    }

    // Operatoren, welche die aktuelle RMatrix verändern
    // ------------------------------------------------
	// z=1, s=1 -> 0,0
	void AddToElement(REAL aReal,int z, int s);
    // Multiplikation mit Skalar: m *= r
    RMatrix& operator*=(const REAL& r);
      // Division durch Skalar: m /= r
    RMatrix& operator/=(const REAL& r);
      // RMatrix addieren: m1 += m2
    RMatrix& operator+=(const RMatrix& m2);
      // RMatrix subtrahieren: m1 -= m2
    RMatrix& operator-=(const RMatrix& m1);
      // Multiplikation mit RMatrix: m1 *= m2
    RMatrix& operator*=(const RMatrix& m2);
      // Multiplikation mit Vektor: m *= v
    //RVektor& operator*=(const RVektor& v);
    //CVektor& operator*=(const CVektor& v);

      // konstante Operatoren
      // --------------------
      // Multiplikation mit Skalar: m * r
    RMatrix operator*(const REAL&) const;
      // Multiplikation mit Skalar: r * m
    friend RMatrix operator*(const REAL& r,
                            const RMatrix& m)
                      { return m * r; }
      // Division durch Skalar: m / r
    RMatrix operator/(const REAL&) const;
      // RMatrix addieren: m1 + m2
    RMatrix operator+(const RMatrix& m2) const;
      // RMatrix subtrahieren: m1 - m2
    RMatrix operator-(const RMatrix& m2) const;
      // Multiplikation mit RMatrix: m1 * m2
    RMatrix operator*(const RMatrix& m2) const;
      // Un„res minus: -m
    RMatrix operator-() const;
    //Multiplikation mit einem Vektor
    RVektor operator*(const RVektor& aVec) const;

    //Division durch einen Vektor (GausAlgo)
    RVektor operator/(const RVektor& aVec);

    void InsertSpalte(const RVektor& aVec, int indx);
    void InsertMatrix(const RMatrix& aMa, int z, int s);
	void SetDimension(int z, int s);
    int Rows() const { return Size; }           //Zeilen
    int Cols() const { return v[0]->nElems(); } //Spalten

    static int Error()
    {
      if ( ErrFlag )
      {
	     ErrFlag = 0;
	     return 1;
      }
      else
	       return RVektor::Error();
    }
    static void Error( char *Msg )
    {
      ErrFlag = 1;
      if ( ErrHdlr && Msg )
	       ErrHdlr( Msg );
    }
    static void SetErrHdlr( void (*Hdlr)(char*) =0 )
    {
      ErrHdlr = Hdlr;
    }

    //Type Min() const;
    //Type Max() const;
    //Type Mean() const;
    int nElems() const { return Cols() * Rows(); }

      // Einheits-RMatrix erzeugen
      // statische Funktion: Aufruf mit 'RMatrix::Identity(n)'
    static RMatrix Identity(int n);

};

RVektor RealGaussElimination(const RMatrix& aMa,const RVektor& bVec);
#endif

//-----------------------------------------------------
// Stream Ausgabe

std::ostream& operator<<(std::ostream& os, const RMatrix& m);

// Stream Eingabe

std::istream& operator>>(std::istream& is, RMatrix& m);
