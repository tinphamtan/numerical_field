// Vector.h
// Vector Klassendefinition

#ifndef VEKTOR_H
#define VEKTOR_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000





//extern Type Zero;
template <class Type> class Vektor
{
  protected:
    Type *a;
    int Size;

      // Statische (Klassen-) Variablen
      // ------------------------------

      // Zeiger auf Funktion zur Fehlerbehandlung
    static void (*ErrHdlr)(char*);

      // Default-Fehlerbehandlung
    static void DefaultErrHdlr(char* Msg)
	{
       std::cout << "Vektor Fehler: " << Msg << "\n";
    };

      // Fehler Flag: 1 = letzte Operation fehlerhaft
    static char ErrFlag;

  public:
      // Erzeugen eines Vektors mit 'n' Elementen
//    friend class CVektor;
	friend class RVektor;
	friend class RMatrix;
//	friend class CMatrix;
//	friend class RealSignal;
//	friend class CmplxSignal;
//-----------------------------------------------
    Vektor()
    {
      Size = 1;
      a = new Type[1];
      a[0] = Type();
	};
//-----------------------------------------------
    Vektor( int n )
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
   };
//----------------------------------------------------------------
      // Konstruktor kopieren
    Vektor(const Vektor& v)
    {
       Size = v.Size;
       a = new Type[Size];
       for ( int i=0; i<Size; i++ )
                           a[i] = v.a[i];
    };
//------------------------------------------------------------------
   // Erzeugen eines Vektors mit 2 Elementen
    Vektor( Type x1, Type x2  )
    {
      Size = 2;
      a = new Type[ Size ];
      a[0] = x1; a[1] = x2;
    };
//-----------------------------------------------------------------
      // Erzeugen eines Vektors mit 3 Elementen
    Vektor( Type x1, Type x2, Type x3 )
    {
      Size = 3;
      a = new Type[ Size ];
      a[0] = x1; a[1] = x2; a[2] = x3;
    }
//----------------------------------------------------------------
      // Destruktor
    ~Vektor()
    {
      delete [Size]a;
    }
//----------------------------------------------------------------
      // Zuweisungs-Operator
    Vektor& operator=(const Vektor& v)
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
    };
//-----------------------------------------------------------------
	void SetDim(int nDim)//Löschen und neue Dimension
	{
		if(Size)
		      delete [Size]a;
	    Size = nDim;
		a    = new Type[Size];
    };
//-----------------------------------------------------------------
      // Index Operator
    Type& operator[]( int n ) const
    {
      if ( n>=0 && n<Size )
        return *(a+n);
      else
      {
        Error( "ungültiger Index" );
        return *(a);
      }
    }
//-----------------------------------------------------------------
      // Nur-Lese Index Operator, fr konstante Objekte (const)
    const Type& Elem( int n ) const
    {
      if ( n>=0 && n<Size )
        return a[n];
      else
      {
        Error( "ungültiger Index" );
        return a[0];
      }
    }
//-----------------------------------------------------------------
      // Operatoren, die das Objekt modifizieren
      // ---------------------------------------
      // Multiplikation mit Skalar: v *= r
//    Vektor& operator*=(const float&);
//    Vektor& operator*=(const Type&);

      // Division durch Skalar: v /= r
//    Vektor& operator/=(const float&);
//    Vektor& operator/=(const Type&);

 
      // Subtraktion eines Vektors: v -= w
 //   Vektor& operator-=(const Vektor&);


      // Konstante Operatoren
      // --------------------

      // Un„res Minus: -v
 //   Vektor  operator-() const;

      // Vektor*Skalar: v * r
    Vektor  operator*(const float&) const;
    Vektor  operator*(const Type&) const;

      // Vektor/Skalar: v / r
    Vektor  operator/(const float&) const;
    Vektor  operator/(const Type&) const;


      // Vergleich zweiter Vektoren: v==w
    int operator==(const Vektor& v) const
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
   };
   int operator!=(const Vektor& w) const
      { return ! (*this==w); }

      // Prüfen ob die letzte Operation zu einem Fehler führte
    static int Error()
    {
      int Result = ErrFlag;
      ErrFlag = 0;
      return Result;
    }
      // Fehlermeldung ausgeben
    static void Error( char *Msg )
    {
      ErrFlag = 1;
      if ( ErrHdlr && Msg )
                  ErrHdlr( Msg );
    };
      // Routine zur Fehlerbehandlung setzen
    static void SetErrHdlr( void (*Hdlr)(char*) = 0 ) //const
    {
      ErrHdlr = Hdlr;
    };


/*      // Rückgabe des Minimums aller Elemente
    Type Min() const;

      // Rückgabe des Maximus aller Elemente
    Type Max() const;

      // Rückgabe des Mittelwertes aller Elemente
    float Mean() const;  */

      // Rückgabe der Anzahl der Elemente
    int nElems() const { return Size; }

};

//typedef Vektor<Type> *VektorPtr;
// Ausgabe zu einem ostream
/*
std::ostream& operator<<(std::ostream& os,  const Vektor<Type>& v)
{
  os << '(';
  int n = v.nElems() - 1;
  int i;
  for (i=0; i<n; i++ )
          os << v.Elem(i) << ' ';
  os << v.Elem(i) << ')';
  return os;
};
*/
// Eingabe von einem istream
//std::istream& operator>>(std::istream& is,  Vektor<Type>& v );


//-------------------------------------------------------------------
#endif
// eof Vektor.h
