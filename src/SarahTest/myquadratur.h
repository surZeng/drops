// Meine Quadraturklassen für Quadratur auf Raum-Zeit-Interface


#ifndef DROPS_MYQUADRATUR_H
#define DROPS_MYQUADRATUR_H

#include "misc/container.h"
#include "SarahTest/mygeometry.h"
#include "num/discretize.h"

namespace DROPS
{

class Quad5SpaceTimeInterfaceDataCL; //Quadraturklasse mit 7 Stützstellen und Exaktheitsgrad 5
//class Quad5SpaceTimeInterfaceCL;
class Quad9Space1DDataCL; //Quadraturklasse für Linienstücke mit 5 Stützstellen und Exaktheitsgrad 9 (Gauß-Quadratur)



//**************************************************************************
// Class:   MyLocalP1CL                                                    *
// Template Parameter:                                                     *
//          T - The result-type of the finite-element-function             *
// Purpose: Evaluate a P1-function on a triangle and calculate with such   *
//          functions. As MyLocalP1CL is derived from valarray, arithmetic *
//          operations are carried out efficiently.                        *
//          The valarray holds the values in the 3 degrees of freedom,     *
//          vertex_0,..., vertex_2.                                        *
//**************************************************************************
template<class T= double>
class MyLocalP1CL: public GridFunctionCL<T>
{
  public:
    typedef GridFunctionCL<T> base_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::instat_fun_ptr instat_fun_ptr;
    
  protected:
    typedef MyLocalP1CL<T> self_;

  public:
    MyLocalP1CL() : base_type( value_type(), 3) {}  // 3 DoFs
    MyLocalP1CL(const value_type& t): base_type( t, 3) {}
    // Initialize from a given function
    MyLocalP1CL(const MyTriangleCL&, instat_fun_ptr , double= 0.0);
	// Initialize from a given Vector with 3 entries (for instance Hat-Function)
	MyLocalP1CL(const MyTriangleCL&, T a[3]);
    

DROPS_DEFINE_VALARRAY_DERIVATIVE(MyLocalP1CL, T, base_type)

    // These "assignment-operators" correspond to the constructors
    // with multiple arguments
    inline self_&
    assign(const MyTriangleCL&, instat_fun_ptr, double= 0.0); // instat_fun_ptr sollte lineare Funktion auf MyTriangle sein
	inline self_&
	assign(const MyTriangleCL&, T a[3]);
    

    // pointwise evaluation in barycentric coordinates
    inline value_type operator()(const Point3DCL& p) const;
};



//**************************************************************************
// Class:   MyLocalP2CL                                                    *
// Template Parameter:                                                     *
//          T - The result-type of the finite-element-function             *
// Purpose: Evaluate a P2-function on a triangle and calculate with such   *
//          functions. As LocalP2CL is derived from valarray, arithmetic   *
//          operations are carried out efficiently.                        *
//          The valarray holds the values in the 6 degrees of freedom,     *
//          vertex_0,..., vertex_2, edge_0,..., edge2.                     *
//**************************************************************************
template<class T= double>
class MyLocalP2CL: public GridFunctionCL<T>
{
  public:
    typedef GridFunctionCL<T> base_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::instat_fun_ptr instat_fun_ptr;
    

  protected:
    typedef MyLocalP2CL<T> self_;

  public:
    MyLocalP2CL() : base_type( value_type(), 6) {}
    MyLocalP2CL(const value_type& t): base_type( t, 6) {}
    // Initialize from a given function
    MyLocalP2CL(const MyTriangleCL&, instat_fun_ptr , double= 0.0);
	// Initialize from a given Vector with 6 entries (for instance P2-Basis-Function)
	MyLocalP2CL(const MyTriangleCL&, T a[6]);
	// Initialize from MyLocalP1CL (braucht man wenn man 2 MyLocalP1CL miteinander multipliziert und dann MyLocalP2CL erhält)
	MyLocalP2CL(const MyLocalP1CL<T>&);
    

DROPS_DEFINE_VALARRAY_DERIVATIVE(MyLocalP2CL, T, base_type)

    // These "assignment-operators" correspond to the constructors
    // with multiple arguments
    inline self_&
    assign(const MyTriangleCL&, instat_fun_ptr, double= 0.0);
	inline self_&
	assign(const MyTriangleCL&, T a[6]);
	inline self_&
    assign(const MyLocalP1CL<T>&);

	// P2-Basisfunktionen auf Referenzdreieck in baryzentrischen Koordinaten
	// WICHTIG: es gelten folgenden baryzentrischen Koordinaten für T_hut:
	//    Knoten0=(0,0) entspricht baryzentrischen Koordinaten (1,0,0)
	//    Knoten1=(1,0) entspricht baryzentrischen Koordinaten (0,1,0)
	//    Knoten2=(0,1) entspricht baryzentrischen Koordinaten (0,0,1)
	//    Knoten3=(1/2,0) entspricht (1/2,1/2,0)
	//    Knoten4=(1/2,1/2)
	//    Knoten5=(0,1/2)
    static double H0(const Point3DCL& p)	{ const double sum= p[1] + p[2]; return 1. +sum*(2.*sum -3.); }
    static double H1(const Point3DCL& p)    { return p[1]*(2.*p[1] -1.); }
    static double H2(const Point3DCL& p)    { return p[2]*(2.*p[2] -1.); }
    static double H3(const Point3DCL& p)    { return 4.*p[1]*( 1. -(p[1] + p[2]) ); }
    static double H4(const Point3DCL& p)    { return 4.*p[1]*p[2]; }
    static double H5(const Point3DCL& p)    { return 4.*p[2]*( 1. -(p[1] + p[2]) ); }

    // pointwise evaluation in barycentric coordinates
    inline value_type operator()(const Point3DCL&) const;
};

//*************************************************************
//
//   Klasse für die Quadratur-Daten (2d-Quadratur)
//
//*************************************************************
class Quad5SpaceTimeInterfaceDataCL
{
  public:
    Quad5SpaceTimeInterfaceDataCL ();

    enum { NumNodes= 7 };

	// statische Variablen, damit nur einmal initialisiert
    static Point3DCL  Node[NumNodes];    //Stützstellen in baryzentrischen Koordinaten (lambda1, lambda2, lambda3) für dreieck
    static const double Weight[NumNodes];  //Gewichte

    ///  M enthält die Koordinaten der 3 Knoten des Dreiecks
    ///  p: darin werden Quadraturstützstellen für dieses Dreieck gespeichert
    ///          If p == 0, the array is new[]-allocated
    /// \return  adress of the array of quadrature points
    Point3DCL* TransformNodes (const SArrayCL<Point3DCL,3>& M, Point3DCL* p= 0);
};


//*************************************************************
//
//   Klasse für die Quadratur-Daten (1d-Quadratur)
//
//*************************************************************
class Quad9Space1DDataCL
{
  public:
    Quad9Space1DDataCL();

    enum { NumNodes= 5 };

	// statische Variablen, damit nur einmal initialisiert
    static Point2DCL  Node[NumNodes];    //Stützstellen in baryzentrischen Koordinaten (lambda1, lambda2) für intervall
	static const double Weight[NumNodes];  //Gewichte

    ///  M enthält die Koordinaten der 2 Knoten des Liniensegments
    ///  p: darin werden Quadraturstützstellen für dieses Liniensegment gespeichert
    ///          If p == 0, the array is new[]-allocated
    /// \return  adress of the array of quadrature points
    Point3DCL* TransformNodes (const SArrayCL<Point3DCL,2>& M, Point3DCL* p= 0);
};


//*******************************************************************
//
//  Quadraturklasse 2D (T ist Template-Parameter für die Funktionswerte)
//
//*******************************************************************



template<class T=double>
class Quad5SpaceTimeInterfaceCL: public GridFunctionCL<T>
{
  public:
    typedef GridFunctionCL<T> base_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::instat_fun_ptr instat_fun_ptr;
    
    Quad5SpaceTimeInterfaceDataCL QuadraturData; // Daten der Quadraturregel


  protected:
    typedef Quad5SpaceTimeInterfaceCL<T> self_;

  public:
    Quad5SpaceTimeInterfaceCL(): base_type( value_type(), Quad5SpaceTimeInterfaceDataCL::NumNodes) {}
    Quad5SpaceTimeInterfaceCL(const value_type& t): base_type( t, Quad5SpaceTimeInterfaceDataCL::NumNodes) {}

    Quad5SpaceTimeInterfaceCL(const MyTriangleCL&, instat_fun_ptr, double= 0.0);
	Quad5SpaceTimeInterfaceCL(const MyLocalP1CL<value_type>&);
    Quad5SpaceTimeInterfaceCL(const MyLocalP2CL<value_type>&);
    

DROPS_DEFINE_VALARRAY_DERIVATIVE(Quad5SpaceTimeInterfaceCL, T, base_type)// ??????????

    inline self_&
    assign(const MyTriangleCL&, instat_fun_ptr, double= 0.0);
    inline self_&
    assign(const MyLocalP1CL<value_type>&);
    inline self_&
    assign(const MyLocalP2CL<value_type>&);
    
    

    // Integration:
    // absdet wird als Parameter uebergeben, damit dieser Faktor bei der
    // Diskretisierung nicht vergessen wird (beliebter folgenschwerer Fehler :-)
    T quad (double absdet) const
    {
      double sum=0.0;
	  for(int i=0; i<Quad5SpaceTimeInterfaceDataCL::NumNodes; i++)
	  {
		  sum+= Quad5SpaceTimeInterfaceDataCL::Weight[i]*(*this)[i];
	  }
	  return sum*absdet;
    }

    
};



//*******************************************************************
//
//  Quadraturklasse 1D (T ist Template-Parameter für die Funktionswerte)
//
//*******************************************************************



template<class T=double>
class Quad9Space1DCL: public GridFunctionCL<T>
{
  public:
    typedef GridFunctionCL<T> base_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::instat_fun_ptr instat_fun_ptr;
    
    Quad9Space1DDataCL QuadraturData; // Daten der Quadraturregel


  protected:
    typedef Quad9Space1DCL<T> self_;

  public:
    Quad9Space1DCL(): base_type( value_type(), Quad9Space1DDataCL::NumNodes) {}
    Quad9Space1DCL(const value_type& t): base_type( t, Quad9Space1DDataCL::NumNodes) {}

    Quad9Space1DCL(const MyLineCL&, instat_fun_ptr, double= 0.0);
	//übergebe die 2 FE-Funktionswerte an den Knoten des Liniensegments, lineare Interpolation dazwischen
	Quad9Space1DCL(value_type values[2]);
   
    

DROPS_DEFINE_VALARRAY_DERIVATIVE(Quad9Space1DCL, T, base_type)// ??????????

    inline self_&
    assign(const MyLineCL&, instat_fun_ptr, double= 0.0);
    inline self_&
    assign(value_type values[2]);
    
    
    

    // Integration:
    // absdet wird als Parameter uebergeben, damit dieser Faktor bei der
    // Diskretisierung nicht vergessen wird (beliebter folgenschwerer Fehler :-)
    T quad (double absdet) const
    {
      double sum=0.0;
	  for(int i=0; i<Quad9Space1DDataCL::NumNodes; i++)
	  {
		  sum+= Quad9Space1DDataCL::Weight[i]*(*this)[i];
	  }
	  return sum*absdet;
    }

    
};





}// end of namespace DROPS

#include "SarahTest/myquadratur.tpp"

#endif
