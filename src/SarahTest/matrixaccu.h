// Klassen für die Aufstellung der Matrizen

#ifndef DROPS_MATRIXACCU_H
#define DROPS_MATRIXACCU_H

#include "SarahTest/mygeometry.h"
#include "SarahTest/myfespace.h"
#include "SarahTest/myquadratur.h"
#include "misc/container.h"
#include <map>

namespace DROPS
{

class MassAccumulatorP1CL;
class RHSAccumulatorP1CL;
class StiffAccumulatorP1CL;
class DerivAccumulatorP1CL;



//****************************************************
//
//	Klasse: MassAccumulatorP1CL
//
//****************************************************
class MassAccumulatorP1CL
{
public:
	// typedef für Koeffizient-Funktionen (entweder skalar, oder vektoriell für Geschwindigkeitsfeld)
	typedef double (*scalar_fct_ptr) (const Point3DCL& P, double t); // P[2] ist Zeit, t wird nicht gebraucht
	typedef Point2DCL (*vec2D_fct_ptr)(const Point3DCL& P, double t);// P[2] ist Zeit, t wird nicht gebraucht
	typedef SMatrixCL<2,2> (*mat2_fct_ptr) (const Point3DCL& P, double t); // gibt 2x2-Matrix zurück (für grad(w))

private:
	const MySpaceTimeGridCL& stg_;
    //const BndDataCL<> * BndData_;
    MatrixCL* Mmat_;
    //VecDescCL* b_;
    SpaceTimeInterfaceP1CL fe_;
	scalar_fct_ptr f_; // Koeffizient
	vec2D_fct_ptr w_;  // Geschwindigkeitsfeld
	mat2_fct_ptr gradw_; //Gradient des Geschwindigkeitsfeld 
						 // Eintrag (0,0): d/dx w1, Eintrag (1,0): d/dy w1
						 // Eintrag (0,1): d/dx w2, Eintrag (1,1): d/dy w2
	scalar_fct_ptr divw_; //Divergenz des Geschwindigkeitsfeld div(w)

    MatrixBuilderCL * M_;

    
    double absdet;
    IdxT Numb[3];
    MyLocalP1CL<double> phi[3];
    //Quad5SpaceTimeInterfaceCL<> phiQuad[3];
    
    double coup[3][3];
    
    static const Uint idx=0;
	
	bool surfactant_;

    const double t;

  public:
	void begin_accumulation ();
	void finalize_accumulation();
	void visit ( MySpaceTimeGridCL&, const MyTriangleCL&);
	void local_setup( MySpaceTimeGridCL&, const MyTriangleCL&);
    void update_global_matrix(MySpaceTimeGridCL& , const MyTriangleCL&);
	
   //  void update_coupling(const MyTriangleCL& triang);

 
    MassAccumulatorP1CL (const MySpaceTimeGridCL & stg, /*const BndDataCL<> * BndData,*/ MatrixCL* Mmat, /*VecDescCL* b,*/
            SpaceTimeInterfaceP1CL& fe, scalar_fct_ptr f, vec2D_fct_ptr w, mat2_fct_ptr gradw, scalar_fct_ptr divw, 
			bool surfactant, const double t_=0);

    
    
};

//****************************************************
//
//	RHSAccumulatorP1CL
//
//****************************************************
class compare { // lexikographische ordnung für Point3DCL
   public:
      bool operator()(const DROPS::Point3DCL& P, const DROPS::Point3DCL& Q) // ist P<Q?
	  {
		  if(P[0]>Q[0])
		  {
			  return false;
		  }else{
			  if(P[0]==Q[0])
			  {
				  if(P[1]>Q[1])
				  {
					  return false;
				  }else{
					  if(P[1]==Q[1])
					  {
						  if(P[2]>=Q[2]) return false;
					  }
				  }
			  }
		  }
		  return true;
		  
	  } 
};

// typdef für Map für die Abbildung von Randknoten zu Fktswerten
typedef std::map<Point3DCL, double, compare> boundarymap;

class RHSAccumulatorP1CL
{
public:
	// typedef für skalare Randdaten-Funktionen und vektorielles Geschwindigkeitsfeld 
	
	typedef double (*scalar_fct_ptr) (const Point3DCL& P, double t); // P[2] ist Zeit, t wird nicht gebraucht
	typedef Point2DCL (*vec2D_fct_ptr)(const Point3DCL& P, double t);// P[2] ist Zeit, t wird nicht gebraucht

private:
	const MySpaceTimeGridCL& stg_;
	boundarymap BndMap_;	// Map für die Abbildung von Randknoten zu Funktionswerten von u_h für Kopplung von Time Slabs
    scalar_fct_ptr BndData_; // Randdaten
	scalar_fct_ptr f_;		 // Quellterm der rechten Seite
	vec2D_fct_ptr  w_;       // Geschwindigkeitsfeld
    VecDescCL* b_;
    SpaceTimeInterfaceP1CL fe_;
    
	// zeiger, welche matrizen aufgebaut werden können, damit kann man anzahlMatrizen setzen
	MassAccumulatorP1CL* mass_;
	// stiff und konv werden nicht mass-accumulator sonder andere accumulator
	StiffAccumulatorP1CL* stiff_;
	//double alpha_; //ist diffusionskoeffizient vor der steifigkeitsmatrix
	DerivAccumulatorP1CL* deriv_;
    
    double absdet;
    IdxT Numb[3];
    MyLocalP1CL<double> phi[3];

	double anzahlMatrizen; //gibt anzahl der matrizen != 0, die aufgestellt werden sollen (brauche ich bei local_setup), beachtet auch Diffusionskoeffizient alpha
   
    
    double coup[3];
    
    static const Uint idx=0;
	
	
    const double t;

  public:
	void begin_accumulation ();
	void visit ( MySpaceTimeGridCL&, const MyTriangleCL&);
	void local_setup( MySpaceTimeGridCL&, const MyTriangleCL&);
    void update_global_rhs(MySpaceTimeGridCL& , const MyTriangleCL&);
	
   

 
   RHSAccumulatorP1CL (const MySpaceTimeGridCL & stg, scalar_fct_ptr BndData, boundarymap BndMap, scalar_fct_ptr f, VecDescCL* b,
            SpaceTimeInterfaceP1CL& fe,  vec2D_fct_ptr w, 
			MassAccumulatorP1CL* mass, StiffAccumulatorP1CL* stiff, /*double alpha,*/ DerivAccumulatorP1CL* deriv, const double t_=0);

    
    
};



//****************************************************
//
//	Klasse: StiffAccumulatorP1CL
//
//****************************************************
class StiffAccumulatorP1CL
{
public:
	// typedef für Koeffizient-Funktionen (entweder skalar, oder vektoriell für Geschwindigkeitsfeld)
	typedef double (*scalar_fct_ptr) (const Point3DCL& P, double t); // P[2] ist Zeit, t wird nicht gebraucht
	typedef Point2DCL (*vec2D_fct_ptr)(const Point3DCL& P, double t);// P[2] ist Zeit, t wird nicht gebraucht
	

private:
	const MySpaceTimeGridCL& stg_;
    MatrixCL* Smat_;
    SpaceTimeInterfaceP1CL fe_;
	scalar_fct_ptr f_; // Koeffizient
	vec2D_fct_ptr w_;  // Geschwindigkeitsfeld
    MatrixBuilderCL * S_;

    
    double absdet;
    IdxT Numb[3];
    MyLocalP1CL<double> phi[3];
	static const Point2DCL GradRef[3]; // Gradienten der 3 Hutfunktionen auf Referenzdreieck
    
    
    double coup[3][3];
    
    static const Uint idx=0;
	
	bool surfactant_;

    const double t;

  public:
	void begin_accumulation ();
	void finalize_accumulation();
	void visit ( MySpaceTimeGridCL&, const MyTriangleCL&);
	void local_setup( MySpaceTimeGridCL&, const MyTriangleCL&);
    void update_global_matrix(MySpaceTimeGridCL& , const MyTriangleCL&);
	
   
 
    StiffAccumulatorP1CL (const MySpaceTimeGridCL & stg,  MatrixCL* Smat, 
            SpaceTimeInterfaceP1CL& fe, scalar_fct_ptr f, vec2D_fct_ptr w, 
			bool surfactant, const double t_=0);

    
    
};




//****************************************************
//
//	Klasse: DerivAccumulatorP1CL // stellt int_\Gamma_stern dot(u)v auf
//
//****************************************************
class DerivAccumulatorP1CL
{
public:
	// typedef für Koeffizient-Funktionen (entweder skalar, oder vektoriell für Geschwindigkeitsfeld)
	typedef double (*scalar_fct_ptr) (const Point3DCL& P, double t); // P[2] ist Zeit, t wird nicht gebraucht
	typedef Point2DCL (*vec2D_fct_ptr)(const Point3DCL& P, double t);// P[2] ist Zeit, t wird nicht gebraucht
	

private:
	const MySpaceTimeGridCL& stg_;
    MatrixCL* Umat_;
    SpaceTimeInterfaceP1CL fe_;
	scalar_fct_ptr f_; // Koeffizient
	vec2D_fct_ptr w_;  // Geschwindigkeitsfeld
    MatrixBuilderCL * U_;

    
    double absdet;
    IdxT Numb[3];
    MyLocalP1CL<double> phi[3];
	static const Point2DCL GradRef[3]; // Gradienten der 3 Hutfunktionen auf Referenzdreieck
    
    
    double coup[3][3];
    
    static const Uint idx=0;
	
	bool surfactant_;

    const double t;

  public:
	void begin_accumulation ();
	void finalize_accumulation();
	void visit ( MySpaceTimeGridCL&, const MyTriangleCL&);
	void local_setup( MySpaceTimeGridCL&, const MyTriangleCL&);
    void update_global_matrix(MySpaceTimeGridCL&, const MyTriangleCL&);
	
   
 
    DerivAccumulatorP1CL (const MySpaceTimeGridCL & stg,  MatrixCL* Umat, 
            SpaceTimeInterfaceP1CL& fe, scalar_fct_ptr f, vec2D_fct_ptr w, 
			bool surfactant, const double t_=0);

    
    
};


// stellt matrizen und rechte seite auf, indem einmal über gitter gelaufen wird
void accumulate_matrix( MySpaceTimeGridCL& stg, MassAccumulatorP1CL * mass, StiffAccumulatorP1CL * stiff,
						DerivAccumulatorP1CL * deriv, RHSAccumulatorP1CL * rhs);

// prüft, ob Knoten auf Dirichlet-Rand (unterer Rand) liegt
bool IsOnBoundary(MySpaceTimeGridCL& stg, const MyVertexCL& vert);

// prüft, ob Knoten auf dem oberen Rand liegt
bool IsOnUpperBoundary(MySpaceTimeGridCL& stg, const MyVertexCL& vert);

// baut boundarymap auf für kopplung zwischen zwei time slabs
boundarymap BuildBoundaryMap(MySpaceTimeGridCL& stg, VectorCL& solution);

// berechnet Normale n_gamma an \Gamma_h(t)  und berechnet vh und vh*n_gamma
void GetN_gamma(const MyTriangleCL& triang, Point2DCL& ngamma, double& vh, Point2DCL& vhngamma);

// berechnet w_h für das Dreieck triang und speichert es in wh
void GetW_h(const MyTriangleCL& triang, Point2DCL& wh);

// gibt Determinante auf Dereieck zurück
double GetAbsDet(const MyTriangleCL& triang);

// gibt Trafo-Matrix zurück
SMatrixCL<3,2> GetGradTrafoMatrix ( const MyTriangleCL& triang);

// gibt Determinante auf Linienstück zurück (Det= norm(P-Q))
double GetAbsDet1D(const MyLineCL & line);

}



#endif