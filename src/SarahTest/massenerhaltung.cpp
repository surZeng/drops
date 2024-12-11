// hier wird experimentell untersucht, ob massenerhaltung auch im diskreten gilt

#include "SarahTest/matrixaccu.h"
#include "SarahTest/vtkoutput.h"
#include "levelset/levelset.h"
#include "levelset/surfacetension.h"
#include "levelset/adaptriang.h"
// für Direkte Löser
#include "num/directsolver.h"
#include "num/precond.h"
#include "num/krylovsolver.h"
#include <fstream>

double sigma= 0.022; // sigma aus film.json
double sigmaf (const DROPS::Point3DCL&, double) { return sigma; } // sigma-funktion (pointer-to-function)

// Level-Set-Funktion ...
//schiefe Möhre
double lsetf (const DROPS::Point3DCL& P, double t) { return 1./((1.0e-3- 1./12.*P[2])*(1.0e-3- 1./12.*P[2]))*( (P[0]-3.5e-3-1./3.*P[2])*(P[0]-3.5e-3-1./3.*P[2]) + (P[1]-3.5e-3)*(P[1]-3.5e-3) ) - 1.0 ; }

// Geschwindigkeitsfeld w
//schiefe Möhre
DROPS::Point2DCL wf (const DROPS::Point3DCL& P, double t) { 
	DROPS::Point2DCL res= -1./(12.* (1.0e-3 - 1./12.*P[2]))* DROPS::MakePoint2D(P[0] - 3.5e-3 - 1./3.*P[2], P[1] - 3.5e-3) + DROPS::MakePoint2D(1./3., 0.0);
	return res;
}

// div(w)
// gerade Möhre, schiefe Möhre
double divwf (const DROPS::Point3DCL& P, double t) { 
	if(1.0e-3 - 1./12.*P[2] < 1.0e-14) {std::cout << "division durch 0 in divwf\n"; }
	return -1./(6.* (1.0e-3 - 1./12.*P[2])); }

// grad(w)
//gerade Möhre, schiefe Möhre
DROPS::SMatrixCL<2,2> gradwf (const DROPS::Point3DCL& P, double t) 
{	DROPS::SMatrixCL<2,2> M;
	double fac= -1./(12.* (1.0e-3 - 1./12.*P[2]));
	if(1.0e-3- 1./12.*P[2] < 1.0e-14) {std::cout << "division durch 0 in gradwf\n"; }
	// erste Spalte von M ist grad(w_1)
	M(0,0)= fac; 
	M(1,0)= 0.0; 
	// zweite Spalte von M ist grad(w_2)
	M(0,1)= 0.0; 
	M(1,1)= fac; 

	return M;
}

// Quellterm f für rechte Seite
double quellf (const DROPS::Point3DCL& P, double t) { return 0.0; }

// Rand-Daten
// randdaten für u_{h,0}=sin(alpha(x,y,t))^2
double bndDataf (const DROPS::Point3DCL& P, double t) {return ((P[1]-3.5e-3)*(P[1]-3.5e-3))/((P[0]-3.5e-3)*(P[0]-3.5e-3) + (P[1]-3.5e-3)*(P[1]-3.5e-3));}

int main(){

	std::fstream file;
	file.open("Massenerhaltung5.dat", std::ios::out);

	int N= 10; // Anzahl der TimeSlabs
	int initial=8;


	for(int i=5; i<=5; i++){ // Verfeinerungen
		DROPS::boundarymap bnd;
		DROPS::MySpaceGridCL sgold;
		DROPS::MySpaceGridCL sgnew;
		double masseold =0.0;
		double massenew =0.0;
		std::cout << "Schleifendurchlauf " << i << "\n\n";

		for(int k=0; k<N; k++){ // gehe N TimeSlabs durch		
			std::cout << "Time Slab " << k+1 << "\n";
			DROPS::MySpaceTimeGridCL newstg;

			// Punkte für den BrickBuilder		
			DROPS::Point3DCL origin = DROPS::MakePoint3D (0.0, 0.0, (k/(double)(N))*6.0e-3);   
			DROPS::Point3DCL e1 = DROPS::MakePoint3D (6.0e-3, 0.0, 0.0);     
			DROPS::Point3DCL e2 = DROPS::MakePoint3D (0.0, 6.0e-3, 0.0);     
			DROPS::Point3DCL e3 = DROPS::MakePoint3D (0.0, 0.0, (6.0e-3)/(double)(N));
		
			// baue MultiGrid
			DROPS::BrickBuilderCL builder( origin, e1, e2, e3, initial, initial, 1 ); 
			DROPS::MultiGridCL* mg;
			if(k==0) std::cout << "Baue MultiGrid...\n";
			mg = new DROPS::MultiGridCL(builder);

			//adaptive Verfeinerung
			//braucht man für anisotrope Gitter nicht		
			double width= 1e-25;
			int clevel= 0; // gröbstes gitter hat level 1
			int flevel1= i; // feinstes gitter hat level 2 (da width= hinitial/2)

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

			if(k==0) std::cout << "Baue Gitter...\n";
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




			//massematrix mit surfactant ohne koeff
			DROPS::MassAccumulatorP1CL masssurf2(newstg, Mmatsurf2, fe, 0, wf, gradwf, divwf, true);


			//stiffmatrix mit surfactant und koeff
			DROPS::StiffAccumulatorP1CL stiffsurf2(newstg, Smatsurf2, fe, 0, wf, true);

			//deriv-Matrix mit surfactant ohne koeff
			DROPS::DerivAccumulatorP1CL derivsurf2(newstg, Umatsurf2, fe, 0, wf, true);

		
			// rechte Seite + accumulate M, S, U, rhs
		
			if(k==0){ // für ersten Time Slab Randdaten aus Funktionenpointer
				DROPS::RHSAccumulatorP1CL rhs(newstg, bndDataf, bnd, quellf, rhskoefftrafo, fe, wf, &masssurf2, &stiffsurf2, &derivsurf2);
				// accumulate:
				std::cout << "accumulate...\n";
				DROPS::accumulate_matrix(newstg, &masssurf2, &stiffsurf2, &derivsurf2, &rhs);
			}else{ // für nächste Time Slabs Randdaten aus Boundarymap
				DROPS::RHSAccumulatorP1CL rhs(newstg, 0,  bnd, quellf, rhskoefftrafo, fe, wf, &masssurf2, &stiffsurf2, &derivsurf2);
				// accumulate:
				DROPS::accumulate_matrix(newstg, &masssurf2, &stiffsurf2, &derivsurf2, &rhs);
			}

			// Löse LGS (M+S+U)x=rhs
			DROPS::MatrixCL MSU;
			MSU.LinComb(1.0, (*Mmatsurf2), 1.0, (*Smatsurf2), 1.0, (*Umatsurf2)); // M+S+U
		

			DROPS::VectorCL x( fe.GetNumUnknowns());

			//Direkte Löser
			DROPS::DirectNonSymmSolverCL MSUSolver(MSU);
			if(k==0) std::cout << "löse lgs \n";
			MSUSolver.Solve(MSU, x, rhskoefftrafo->Data);
		

			// baue BoundaryMap
			if(k==0) std::cout << "baue boundarymap...\n";
			bnd.clear();
			bnd= DROPS::BuildBoundaryMap(newstg, x);


			//berechne int_{\Gamma_h(0)}u_h(\cdot, 0)
			if(k==0){
				//Gitter bei t=0
				sgold.BuildSpaceGrid(newstg, false);
				file << "t_old= " << sgold.Get_t() << "\n";
				DROPS::MySpaceGridCL::MyLineIterator it=sgold.GetLineBegin();
				for(it; it!=sgold.GetLineEnd(); it++){
					double absdet= DROPS::GetAbsDet1D(*it);
					DROPS::Quad9Space1DCL<> u_h_0(*it, bndDataf);
					masseold+= u_h_0.quad(absdet);
				}
			}

			//berechne int_{\Gamma_h(T)}u_h(\cdot, T)
			if(k==N-1){
				//Gitter bei t=T
				sgnew.BuildSpaceGrid(newstg, true);
				file << "t_new= " << sgnew.Get_t() << "\n";
				DROPS::MySpaceGridCL::MyLineIterator ot=sgnew.GetLineBegin();
				for(ot; ot!=sgnew.GetLineEnd(); ot++){
					double absdet= DROPS::GetAbsDet1D(*ot);
					double values[2]= {x[(ot->getVertex(0))->Unknowns(0)], x[(ot->getVertex(1))->Unknowns(0)]};
					DROPS::Quad9Space1DCL<> u_h_T(values);
					massenew+= u_h_T.quad(absdet);
				}
			}


		}// Ende: gehe Time Slabs durch

		double h= 6.0e-3/(initial*std::pow(2.0, i));
		double differenz= std::abs(masseold - massenew);

		std::cout << "h= " << h << "\tMasse0= " << masseold << "\tMasseT= " << massenew << "\tDifferenz= " << differenz << "\n";
		file << h << "\t" << masseold << "\t" << massenew << "\t" << differenz << "\n";

	}// Ende: Verfeinerungen
	file.close();
	return 0;


}