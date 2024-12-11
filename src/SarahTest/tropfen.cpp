// Experiment zu verschmelzenden Tropfen

#include "SarahTest/matrixaccu.h"
#include "SarahTest/vtkoutput.h"
#include "levelset/levelset.h"
#include "levelset/surfacetension.h"
#include "levelset/adaptriang.h"
#include "out/vtkOut.h"
// für Direkte Löser
#include "num/directsolver.h"
#include "num/precond.h"
#include "num/krylovsolver.h"
#include <fstream>

double sigma= 0.022; // sigma aus film.json
double sigmaf (const DROPS::Point3DCL&, double) { return sigma; } // sigma-funktion (pointer-to-function)

// Level-Set-Funktion ...
//verschmelzende Tropfen
//phi(x,y,t)= (x-0.05)^2*(x+0.05)^2 + 0.5*y^2 - t
double lsetf (const DROPS::Point3DCL& P, double t) { return (P[0]-0.5)*(P[0]-0.5)*(P[0]+0.5)*(P[0]+0.5) + 0.5*P[1]*P[1] - 0.1*(P[2]+0.2); }

// Geschwindigkeitsfeld w
//verschmelzende Tropfen
DROPS::Point2DCL wf (const DROPS::Point3DCL& P, double t) {
	if((P[0]==0.5 && P[1]==0.0) || (P[0]==-0.5 && P[1]==0.0) || (P[0]==0.0 && P[1]==0.0)){
		return DROPS::MakePoint2D(0.0, 0.0);
	}else{
		double nenner= 10.*(4.*P[0]*P[0]*P[0]- P[0])*(4.*P[0]*P[0]*P[0]- P[0]) + 10.*P[1]*P[1];
		DROPS::Point2DCL vec= DROPS::MakePoint2D(4.*P[0]*P[0]*P[0]- P[0], P[1]);
		return vec/nenner;
	}
}


// div(w)
// verschmelzende Tropfen
double divwf (const DROPS::Point3DCL& P, double t) { 
	if((P[0]==0.5 && P[1]==0.0) || (P[0]==-0.5 && P[1]==0.0) || (P[0]==0.0 && P[1]==0.0)){
		return 0.0;
	}else{
		double x=P[0];
		double y=P[1];
		double zaehler= -(96.*x*x*x*x*x*x*x*x - 64.*x*x*x*x*x*x + 14.*x*x*x*x - 6.*x*x*y*y - x*x + y*y);
		double nenner= 5.*(16.*x*x*x*x*x*x - 8.*x*x*x*x + x*x + y*y)*(16.*x*x*x*x*x*x - 8.*x*x*x*x + x*x + y*y);
		return zaehler/nenner;	
	}
}

// grad(w)
//verschmelzende Tropfen
DROPS::SMatrixCL<2,2> gradwf (const DROPS::Point3DCL& P, double t) 
{	DROPS::SMatrixCL<2,2> M;
	if((P[0]==0.5 && P[1]==0.0) || (P[0]==-0.5 && P[1]==0.0) || (P[0]==0.0 && P[1]==0.0)){
		M(0,0)=0.0;
		M(1,0)=0.0;
		M(0,1)=0.0;
		M(1,1)=0.0;
	}else{
		double x=P[0];
		double y=P[1];
		double nenner= (16.*x*x*x*x*x*x - 8.*x*x*x*x + x*x + y*y)*(16.*x*x*x*x*x*x - 8.*x*x*x*x + x*x + y*y);
		// erste Spalte von M ist grad(w_1)
		M(0,0)= -(12.*x*x - 1.)*(16.*x*x*x*x*x*x - 8.*x*x*x*x + x*x - y*y)/(10.*nenner); 
		M(1,0)= -x*y*(4.*x*x - 1.)/(5.*nenner); 
		// zweite Spalte von M ist grad(w_2)
		M(0,1)= -x*y*(4.*x*x-1.)*(12.*x*x-1.)/(5.*nenner); 
		M(1,1)= (16.*x*x*x*x*x*x - 8.*x*x*x*x + x*x - y*y)/(10.*nenner); 

	}
	return M;
}


// Quellterm f für rechte Seite
double quellf (const DROPS::Point3DCL& P, double t) { return 0.0; }

// Rand-Daten
// bndDataf= 1
//double bndDataf (const DROPS::Point3DCL& P, double t) {return 1.0;}

//bndDataf= cos(alpha)^2
double bndDataf (const DROPS::Point3DCL& P, double t) {
	if(P[0]>=0){
		double sin= (P[1]-0.0)/std::sqrt((P[0]-0.5)*(P[0]-0.5) + (P[1]-0.0)*(P[1]-0.0));
		double alpha= std::asin(sin);
		return std::cos(alpha)*std::cos(alpha);
	}else{
		double sin= (P[1]-0.0)/std::sqrt((P[0]+0.5)*(P[0]+0.5) + (P[1]-0.0)*(P[1]-0.0));
		double alpha= std::asin(sin);
		return std::cos(alpha)*std::cos(alpha);
		
	}
}

//bndDataf= sin(alpha)^2
/*
double bndDataf (const DROPS::Point3DCL& P, double t) {
	if(P[0]>=0){
		double sin= (P[1]-0.0)*(P[1]-0.0)/((P[0]-0.5)*(P[0]-0.5) + (P[1]-0.0)*(P[1]-0.0));
		return sin;
	}else{
		double sin= (P[1]-0.0)*(P[1]-0.0)/((P[0]+0.5)*(P[0]+0.5) + (P[1]-0.0)*(P[1]-0.0));
		return sin;
		
	}
	
}
*/

// Koeffizient alpha(s,t)
//double alphaf (const DROPS::Point3DCL& P, double t) {return 1.0;}
double alphaf (const DROPS::Point3DCL& P, double t) {return 0.001;}

/*
double bndDataf (const DROPS::Point3DCL& P, double t){
	return std::sin(norm(P));
}
*/



int main(){
			DROPS::MySpaceTimeGridCL newstg;

			// Punkte für den BrickBuilder		
			DROPS::Point3DCL origin = DROPS::MakePoint3D (-1.0, -1.0, 0.0);   
			DROPS::Point3DCL e1 = DROPS::MakePoint3D (2.0, 0.0, 0.0);     
			DROPS::Point3DCL e2 = DROPS::MakePoint3D (0.0, 2.0, 0.0);     
			DROPS::Point3DCL e3 = DROPS::MakePoint3D (0.0, 0.0, 1.0);
		
			// baue MultiGrid
			DROPS::BrickBuilderCL builder( origin, e1, e2, e3, 20, 20, 10 ); 
			DROPS::MultiGridCL* mg;
			mg = new DROPS::MultiGridCL(builder);


			//adaptive Verfeinerung	
			double width= 1e-25;
			int clevel= 0; 
			int flevel1= 1; 
			DROPS::AdapTriangCL adap(*mg, width, clevel, flevel1);
			adap.MakeInitialTriang(lsetf);


			// baue Levelset-Funktion
			DROPS::BndCondT bc_ls[6]= { DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC };
			DROPS::SurfaceTensionCL sf( sigmaf, 0 );
			DROPS::LevelsetP2CL lset( *mg, DROPS::LsetBndDataCL( 6, bc_ls ), sf); // default: SD=0, CurvDiff=-1
			DROPS::IdxDescCL* lidx= &lset.idx; // FE-Raum für levelset funktion
			lset.CreateNumbering( mg->GetLastLevel(), lidx ); // oder .CreateNumbering(*mg.GetLastLevel(),lidx)?
			lset.Phi.SetIdx( lidx); // wozu brauche ich das? -> länge des vektors wird gesetzt
			lset.Init(lsetf); // levelset-function wird initialisiert

			//baue Gitter
			std::cout << "Baue Gitter...\n";
			newstg.BuildMyGeometry( *mg, lset.Phi);

			// lege FE-Raum an
			DROPS::SpaceTimeInterfaceP1CL fe;
			fe.NumberingVertices(newstg);

			//lege Matrizen an
			DROPS::MatrixCL* Mmatsurf2;
			Mmatsurf2= new DROPS::MatrixCL;

			DROPS::MatrixCL* Smatsurf2;
			Smatsurf2= new DROPS::MatrixCL;

			DROPS::MatrixCL* Umatsurf2;
			Umatsurf2= new DROPS::MatrixCL;

			DROPS::VecDescCL* rhskoefftrafo;
			rhskoefftrafo= new DROPS::VecDescCL;


			//massematrix 
			DROPS::MassAccumulatorP1CL masssurf2(newstg, Mmatsurf2, fe, 0, wf, gradwf, divwf, true);


			//stiffmatrix 
			DROPS::StiffAccumulatorP1CL stiffsurf2(newstg, Smatsurf2, fe, alphaf, wf, true);

			//deriv-Matrix
			DROPS::DerivAccumulatorP1CL derivsurf2(newstg, Umatsurf2, fe, 0, wf, true);

			//rechte Seite
			DROPS::boundarymap bnd;
			DROPS::RHSAccumulatorP1CL rhs(newstg, bndDataf, bnd, quellf, rhskoefftrafo, fe, wf, &masssurf2, &stiffsurf2, &derivsurf2);
				
			// accumulate:
			std::cout << "Accumulate...\n";
			DROPS::accumulate_matrix(newstg, &masssurf2, &stiffsurf2, &derivsurf2, &rhs);


			// Löse LGS (M+S+U)x=rhs
			DROPS::MatrixCL MSU;
			MSU.LinComb(1.0, (*Mmatsurf2), 1.0, (*Smatsurf2), 1.0, (*Umatsurf2)); // M+S+U
		

			DROPS::VectorCL x( fe.GetNumUnknowns());

			//Direkte Löser
			DROPS::DirectNonSymmSolverCL MSUSolver(MSU);
			std::cout << "Loese LGS...\n";
			MSUSolver.Solve(MSU, x, rhskoefftrafo->Data);




			//Visualisierung
			DROPS::MyVTKOutCL* vtk = NULL;
			vtk = new DROPS::MyVTKOutCL (newstg, x, x, std::string("vtk"), std::string("TropfenB1") );
			vtk->PutGeometryAndSolution();

			


			delete vtk;
			

	return 0;
}