// L2-Fehler auf \Gamma_h(T) für schiefe Möhre

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
//alternativ: alte Version
//double lsetf (const DROPS::Point3DCL& P, double t) { return (1.0e-3- 1./12.*P[2])*(1.0e-3- 1./12.*P[2]) - (P[0]-3.0e-3-1./3.*P[2])*(P[0]-3.0e-3-1./3.*P[2]) - (P[1]-3.0e-3)*(P[1]-3.0e-3); }
//gerade Möhre
/*
double lsetf (const DROPS::Point3DCL& P, double t) {
	if(1.0e-3 - 1./12.*P[2] < 1.0e-14) {std::cout << "division durch 0 in phi\n"; }
	return ( (P[0]-3.5e-3)*(P[0]-3.5e-3) + (P[1]-3.5e-3)*(P[1]-3.5e-3) ) / ((1.0e-3 - 1./12.*P[2])*(1.0e-3 - 1./12.*P[2])) - 1.0 ; }
*/
//alternativ: alte Version
//double lsetf (const DROPS::Point3DCL& P, double t) { return (1.0e-3- 1./12.*P[2])*(1.0e-3- 1./12.*P[2]) - (P[0]-3.0e-3)*(P[0]-3.0e-3) - (P[1]-3.0e-3)*(P[1]-3.0e-3); }
//schiefer Zylinder
//double lsetf (const DROPS::Point3DCL& P, double t)  { return ((1.0e-6)-(P[0]-P[2]-1.2e-3)*(P[0]-P[2]-1.2e-3)-(P[1]-3.0e-3)*(P[1]-3.0e-3)); }

// Koeffizient alpha(s,t)
double alphaf (const DROPS::Point3DCL& P, double t) {return 0.0001;}

// für Methode 1
// exakte Lösung u
//double uexakt1(const DROPS::Point3DCL& P, double t) {return 2.;}
//double uexakt1(const DROPS::Point3DCL& P, double t) {return 1.0;} 
//double uexakt(const DROPS::Point3DCL& P, double t) {return 1.0/(1.0-83.333333333333333333 * P[2]);}

//uexakt= 1/(1 -1000/12*t)*sin(alpha(x,y,t))
double uexakt(const DROPS::Point3DCL& P, double t){return 1.0/(1.0-83.333333333333333333 * P[2])*(P[1]-3.5e-3)/std::sqrt((P[0]-3.5e-3-1./3.*P[2])*(P[0]-3.5e-3-1./3.*P[2]) + (P[1]-3.5e-3)*(P[1]-3.5e-3));}

// Tangentialgradient an \Gamma(T) der exakten Lösung (skalar, da auf Geradenstück)
//double uexaktgrad1(const DROPS::Point3DCL& P, double t) {return 0.0;}

//DROPS::Point2DCL uexaktgrad(const DROPS::Point3DCL& P, double t) {return DROPS::MakePoint2D(0.0, 0.0);}

//uexaktgrad von uexakt=1/(1 -1000/12*t)*sin(alpha(x,y,t))

DROPS::Point2DCL uexaktgrad(const DROPS::Point3DCL& P, double t) {
	double ux= -27000.*(2000.*P[1]-7.)*(-6000.*P[0]+21.+2000.*P[2])/std::sqrt(36000000.*P[0]*P[0]-24000000.*P[0]*P[2]+36000000.*P[1]*P[1]+4000000.*P[2]*P[2]-252000.*P[0]-252000.*P[1]+84000.*P[2]+882.)/(-3.+250.*P[2])/(18000000.*P[0]*P[0]-12000000.*P[0]*P[2]+18000000.*P[1]*P[1]+2000000.*P[2]*P[2]-126000.*P[0]-126000.*P[1]+42000.*P[2]+441.);
	double uy= -9000./std::sqrt(36000000.*P[0]*P[0]-24000000.*P[0]*P[2]+36000000.*P[1]*P[1]+4000000.*P[2]*P[2]-252000.*P[0]-252000.*P[1]+84000.*P[2]+882.)*(36000000.*P[0]*P[0]-24000000.*P[0]*P[2]+4000000.*P[2]*P[2]-252000.*P[0]+84000.*P[2]+441.)/(-3.+250.*P[2])/(18000000.*P[0]*P[0]-12000000.*P[0]*P[2]+18000000.*P[1]*P[1]+2000000.*P[2]*P[2]-126000.*P[0]-126000.*P[1]+42000.*P[2]+441.);
	return DROPS::MakePoint2D(ux, uy);
}

DROPS::Point3DCL uexaktgrad3(const DROPS::Point3DCL& P, double t) {
	double ux= -27000.*(2000.*P[1]-7.)*(-6000.*P[0]+21.+2000.*P[2])/std::sqrt(36000000.*P[0]*P[0]-24000000.*P[0]*P[2]+36000000.*P[1]*P[1]+4000000.*P[2]*P[2]-252000.*P[0]-252000.*P[1]+84000.*P[2]+882.)/(-3.+250.*P[2])/(18000000.*P[0]*P[0]-12000000.*P[0]*P[2]+18000000.*P[1]*P[1]+2000000.*P[2]*P[2]-126000.*P[0]-126000.*P[1]+42000.*P[2]+441.);
	double uy= -9000./std::sqrt(36000000.*P[0]*P[0]-24000000.*P[0]*P[2]+36000000.*P[1]*P[1]+4000000.*P[2]*P[2]-252000.*P[0]-252000.*P[1]+84000.*P[2]+882.)*(36000000.*P[0]*P[0]-24000000.*P[0]*P[2]+4000000.*P[2]*P[2]-252000.*P[0]+84000.*P[2]+441.)/(-3.+250.*P[2])/(18000000.*P[0]*P[0]-12000000.*P[0]*P[2]+18000000.*P[1]*P[1]+2000000.*P[2]*P[2]-126000.*P[0]-126000.*P[1]+42000.*P[2]+441.);
	double ut= 2250.*(2000.*P[1]-7.)*(18000000.*P[0]*P[0]-18000000.*P[0]*P[2]+18000000.*P[1]*P[1]+4000000.*P[2]*P[2]-54000.*P[0]-126000.*P[1]+39000.*P[2]+189.)/((18000000.*P[0]*P[0]-12000000.*P[0]*P[2]+18000000.*P[1]*P[1]+2000000.*P[2]*P[2]-126000.*P[0]-126000.*P[1]+42000.*P[2]+441.)*std::sqrt(36000000.*P[0]*P[0]-24000000.*P[0]*P[2]+36000000.*P[1]*P[1]+4000000.*P[2]*P[2]-252000.*P[0]-252000.*P[1]+84000.*P[2]+882.)*(-3.+250.*P[2])*(-3.+250.*P[2]));
	return DROPS::MakePoint3D(ux, uy, ut);
}

//für Methode 2



// ????????? das ist jetzt der örtliche Gradient von u(., T), oder soll da grad(u) (ort + zeit ) hin?
//DROPS::Point2DCL uexaktgrad2(const DROPS::Point3DCL& P, double t) {return DROPS::MakePoint2D(0.0, 0.0);}



// Geschwindigkeitsfeld w
//schiefe Möhre

DROPS::Point2DCL wf (const DROPS::Point3DCL& P, double t) { 
	DROPS::Point2DCL res= -1./(12.* (1.0e-3 - 1./12.*P[2]))* DROPS::MakePoint2D(P[0] - 3.5e-3 - 1./3.*P[2], P[1] - 3.5e-3) + DROPS::MakePoint2D(1./3., 0.0);
	return res;
}





// Test wf (konstant, um zu überprüfen, ob grad(u_h) konstant auf linienstück)
//DROPS::Point2DCL wf (const DROPS::Point3DCL & P, double t) {return DROPS::MakePoint2D(1., 1.);}

// gerade Möhre
/*
DROPS::Point2DCL wf (const DROPS::Point3DCL& P, double t) { 
	if(1.0e-3- 1./12.*P[2] < 1.0e-14) {std::cout << "division durch 0 in wf\n"; }
	DROPS::Point2DCL res=  DROPS::MakePoint2D(P[0]-3.5e-3, P[1]-3.5e-3) / (-12. * (1.0e-3 - 1./12.*P[2]))  ;
	return res;
}
*/



//schiefer Zylinder
//DROPS::Point2DCL wf(const DROPS::Point3DCL& P, double t) {return DROPS::MakePoint2D(1.0, 0.0);}




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



//schiefer Zylinder
/*
DROPS::SMatrixCL<2,2> gradwf (const DROPS::Point3DCL& P, double t) 
{	DROPS::SMatrixCL<2,2> M;
	// erste Spalte von M ist grad(w_1)
	M(0,0)= 0.0; 
	M(1,0)= 0.0; 
	// zweite Spalte von M ist grad(w_2)
	M(0,1)= 0.0; 
	M(1,1)= 0.0;

	return M;
}
*/

// Rand-Daten
//double bndDataf (const DROPS::Point3DCL& P, double t) { return 1.0; }
// randdaten für uexakt=1/( 1-1000/12*t)*sin(alpha(x,y,t))
double bndDataf (const DROPS::Point3DCL& P, double t) {return (P[1]-3.5e-3)/std::sqrt((P[0]-3.5e-3)*(P[0]-3.5e-3) + (P[1]-3.5e-3)*(P[1]-3.5e-3));}



// Quellterm f für rechte Seite
//double quellf (const DROPS::Point3DCL& P, double t) { return 0.0; }
// Erweiterung fe für uexakt=1/(1 -1000/12*t)*sin(alpha(x,y,t))

double quellf (const DROPS::Point3DCL& Q, double t){
	// für projektion p= q - dist*ngamma
	DROPS::Point3DCL mittelpunkt= DROPS::MakePoint3D(3.5e-3 + 1./3.*Q[2], 3.5e-3, Q[2]);
	double radius= 1.0e-3 - 1./12.*Q[2];
	double dist= DROPS::norm(Q - mittelpunkt) - radius;
	DROPS::Point3DCL ngamma= DROPS::MakePoint3D(Q[0]-3.5e-3 - 1./3.*Q[2], Q[1]-3.5e-3, 0.0)/ std::sqrt((Q[0]-3.5e-3 - 1./3.*Q[2])*(Q[0]-3.5e-3 - 1./3.*Q[2]) + ( Q[1]-3.5e-3)*( Q[1]-3.5e-3));
	//projektion P
	DROPS::Point3DCL P= Q - dist*ngamma;


	return -162000000.*(2000.*P[1]-7.)/std::sqrt(36000000.*P[0]*P[0]-24000000.*P[0]*P[2]+36000000.*P[1]*P[1]+4000000.*P[2]*P[2]-252000.*P[0]-252000.*P[1]+84000.*P[2]+882.)/(-3.+250.*P[2])/(18000000.*P[0]*P[0]-12000000.*P[0]*P[2]+18000000.*P[1]*P[1]+2000000.*P[2]*P[2]-126000.*P[0]-126000.*P[1]+42000.*P[2]+441.);

}


DROPS::Point3DCL GetNhut(DROPS::MyTriangleCL triang){

	DROPS::Point3DCL p= triang.getVertex(0)->GetCoordinates();
	DROPS::Point3DCL q= triang.getVertex(1)->GetCoordinates();
	DROPS::Point3DCL r= triang.getVertex(2)->GetCoordinates();

	DROPS::Point3DCL a= q-p;
	DROPS::Point3DCL b= r-p;

	// Kreuzprodukt ist Normale an triang
	DROPS::Point3DCL nhut_unnormiert= DROPS::MakePoint3D(a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]);
	DROPS::Point3DCL nhut= nhut_unnormiert/ norm(nhut_unnormiert);

	return nhut;

}


int main()
{

std::fstream file;
file.open("DiskGanzSchiefeMöhreU.dat", std::ios::out);

int N= 1; // Anzahl der TimeSlabs
int initial=8;


for(int i=0; i<=4; i++){ // Verfeinerungen

	DROPS::boundarymap bnd;
	DROPS::MySpaceGridCL sg;
	DROPS::MySpaceGridCL sg2;
	DROPS::VectorCL xx;
	DROPS::VectorCL yy;
	std::cout << "Schleifendurchlauf " << i << "\n\n";
	
	
	//________________________________________________MEHRERE TIME SLABS___________________________________________________
	std::cout << "Methode mit N= " << N << " Time Slabs\n";
	//file << "Methode mit N= " << N << " Time Slabs\n";
	for(int k=0; k<N; k++){ // gehe N TimeSlabs durch

		
		std::cout << "Time Slab " << k+1 << "\n";
		DROPS::MySpaceTimeGridCL newstg;

		// Punkte für den BrickBuilder
		
		DROPS::Point3DCL origin = DROPS::MakePoint3D (0.0, 0.0, (k/(double)(N))*6.0e-3);   
		DROPS::Point3DCL e1 = DROPS::MakePoint3D (6.0e-3, 0.0, 0.0);     
		DROPS::Point3DCL e2 = DROPS::MakePoint3D (0.0, 6.0e-3, 0.0);     
		DROPS::Point3DCL e3 = DROPS::MakePoint3D (0.0, 0.0, (6.0e-3)/(double)(N));
		
		// baue MultiGrid
		DROPS::BrickBuilderCL builder( origin, e1, e2, e3, initial, initial, 10 ); // (10 unterteilungen in x,y,z-richtung)
		//Builder für anisotrope Gitter
		// hs->0, ht=konst
		//DROPS::BrickBuilderCL builder( origin, e1, e2, e3, initial*std::pow(2.,i), initial*std::pow(2.,i), 1 );
		// hs=konst, ht->0
		//DROPS::BrickBuilderCL builder( origin, e1, e2, e3, initial, initial, 1*std::pow(2.,i) );
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
		DROPS::MySpaceGridCL zaehlergitter;
		zaehlergitter.BuildSpaceGrid(newstg, true);
		int doppelzaehler=0;
		for(DROPS::MySpaceGridCL::MyVertexIterator toto=zaehlergitter.GetVerticesBegin(); toto!=zaehlergitter.GetVerticesEnd(); toto++){
			doppelzaehler++;
		}

		//std::cout << "Anzahl Unbekannte= " << fe.GetNumUnknowns() << "\n";
		if(k==0) {file << "Anzahl Unbekannte= " << fe.GetNumUnknowns() << "\n";}
		else{file << "Anzahl Unbekannte= " << fe.GetNumUnknowns() - doppelzaehler << "\n";}



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
		
		





		
		//std::cout << "M=\n" << (*Mmatsurf2) << "\nS=\n" << (*Smatsurf2) << "\nU=\n" <<(*Umatsurf2) << "\nrhs=\n" << rhskoefftrafo-> Data << "\n";

		// Löse LGS (M+S+U)x=rhs
		DROPS::MatrixCL MSU;
		MSU.LinComb(1.0, (*Mmatsurf2), 1.0, (*Smatsurf2), 1.0, (*Umatsurf2)); // M+S+U
		//std::cout << "M+S+U=\n" << MSU << "\n";


		DROPS::VectorCL x( fe.GetNumUnknowns());

		//Direkte Löser

		DROPS::DirectNonSymmSolverCL MSUSolver(MSU);
		if(k==0) std::cout << "löse lgs \n";
		MSUSolver.Solve(MSU, x, rhskoefftrafo->Data);
		//std::cout << "lösung=\n" << x << "\n";

		// baue BoundaryMap
		if(k==0) std::cout << "baue boundarymap...\n";
		bnd.clear();
		bnd= DROPS::BuildBoundaryMap(newstg, x);


		
// Berechne L2-Fehler auf \Gamma_h(T) und H1-Fehler auf \Gamma_h(T) (Methode 1) (nur im letzten TimeSlab)
		if( k == (N-1) ){//Fehlerberechnung nur für letzten TimeSlab
			
			xx= x; // kopiere x in globale variable xx rein, damit ich das später verwenden kann
			sg.BuildSpaceGrid(newstg, true);
			fe.NumberingVerticesSG(sg, newstg);// Gitter sg wird konsistent zu Gitter newstg nummeriert

			// Methode 1
			std::cout << "berechne fehler\n";
			DROPS::MySpaceGridCL::MyLineIterator line=sg.GetLineBegin();
			double sum=0.0;
			double sum2=0.0;
			for(line; line!=sg.GetLineEnd(); line++){
	
				double absdet= DROPS::GetAbsDet1D(*line);
				//std::cout << "absdet= " << absdet << "\n";

				//______________L2____________________________________________-
	
				// für ue (hier ue=uexakt, da uexakt konstant)
				DROPS::Quad9Space1DCL<> ue(*line, uexakt);
				/*		
				std::cout << "ue= \n";
				for(int i=0; i< ue.size(); i++) std::cout << ue[i] << "\n";
				*/
	
				// für u_h
				// values sind Funktionswerte von u_h an Endknoten der Linie line
				double values[2]= {x[(line->getVertex(0))->Unknowns(0)], x[(line->getVertex(1))->Unknowns(0)]};
				//std::cout << "\n uh=\n" << values[0] << "\n" << values[1] << "\n";

				DROPS::Quad9Space1DCL<> uh(values);
	

				// für (ue-uh)^2
				DROPS::Quad9Space1DCL<> res((ue-uh)*(ue-uh));
	
				sum+= res.quad(absdet);
				//std::cout << "sum= " << sum << "\n";

				//______________H1____________________________________________-

				// Tangentialvektor tang
				DROPS::Point2DCL t= DROPS::MakePoint2D((line->getVertex(1))->GetCoordinates()[0]- (line->getVertex(0))->GetCoordinates()[0],
											(line->getVertex(1))->GetCoordinates()[1]- (line->getVertex(0))->GetCoordinates()[1]);
				DROPS::Point2DCL tang= t/norm(t);
				//std::cout << "Tangentialvektor= " << tang << "\n";
				DROPS::Quad9Space1DCL<DROPS::Point2DCL> tangquad(tang);

				// für tangentialgrad(ue) (hier ue=uexakt)
				DROPS::Quad9Space1DCL<DROPS::Point2DCL> uegrad(*line, uexaktgrad);
				DROPS::Quad9Space1DCL<DROPS::Point2DCL> uetanggrad(dot(tangquad,uegrad)* tangquad);

				// für tangentialgrad(u_h)
				// ableitung= ableitung von u_h auf linie (absdet= norm(p-q))
				double ableitung= (values[1]-values[0])/norm(t);
				//std::cout << "ableitung bei methode1= " << ableitung <<  "\n";
				DROPS::Quad9Space1DCL<DROPS::Point2DCL> uhgrad(ableitung*tang);
	

				// für (grad(ue)-grad(uh))^2
				DROPS::Quad9Space1DCL<DROPS::Point2DCL> minus(uetanggrad - uhgrad);
				DROPS::Quad9Space1DCL<> res2(dot(minus,minus));
	
				sum2+= res2.quad(absdet);
	

			}


		// Berechne H1-Fehler (in Seminorm) auf \Gamma_h(T)
		// Methode 2

			double eps= 0.1e-14;
			double tnew=newstg.Get_t_new();
			double sum22=0.0;
			const DROPS::Point2DCL GradRef[3] = {DROPS::MakePoint2D(-1.0,-1.0), DROPS::MakePoint2D(1.0, 0.0), DROPS::MakePoint2D(0.0, 1.0)};

			for(DROPS::MySpaceTimeGridCL::MyTriangIterator iter= newstg.GetTriangBegin(); iter!=newstg.GetTriangEnd(); iter++){ // gehe alle Dreiecke aus SpaceTimeGrid durch
				for(int ii=0; ii<3; ii++){ // (i=0, j=1); (i=1, j=2); (i=2, j=0)
					int jj= (ii+1)%3;
						int kk= (ii+2)%3; // (i=0, k=2) (i=1, k=0) (i=2, k=1)
						// Knoten i und j liegen auf oberem Rand
						if( std::abs((iter->getVertex(ii))->GetCoordinates()[2] - tnew) < eps && 
								std::abs((iter->getVertex(jj))->GetCoordinates()[2] - tnew) < eps) {
						
				
								const DROPS::MyLineCL line( iter->getVertex(ii), iter->getVertex(jj) );
		
								double absdet= DROPS::GetAbsDet1D(line);
						

								//Hilfszeugs						
								double vh=0.0;
								DROPS::Point2DCL ngamma;
								DROPS::Point2DCL vhngamma;
								DROPS::GetN_gamma(*iter, ngamma, vh, vhngamma);


								// für Grad_Gamma_stern
								DROPS::SMatrixCL<3,2> M= DROPS::GetGradTrafoMatrix(*iter);
								DROPS::Point3DCL grad_gamma_stern_i(M*GradRef[ii]);
								DROPS::Point3DCL grad_gamma_stern_j(M*GradRef[jj]);
								DROPS::Point3DCL grad_gamma_stern_k(M*GradRef[kk]);

								// für grad_\Gamma(Chi_i), grad_\Gamma(Chi_j), grad_\Gamma(Chi_k)
						
								DROPS::Quad9Space1DCL<DROPS::Point2DCL> vektor(vhngamma);

								DROPS::Quad9Space1DCL<DROPS::Point2DCL> grad_gamma_i ((DROPS::MakePoint2D(grad_gamma_stern_i[0], grad_gamma_stern_i[1])
															-grad_gamma_stern_i[2]* vektor));
						

								DROPS::Quad9Space1DCL<DROPS::Point2DCL> grad_gamma_j ((DROPS::MakePoint2D(grad_gamma_stern_j[0], grad_gamma_stern_j[1])
															-grad_gamma_stern_j[2]* vektor));
						

								DROPS::Quad9Space1DCL<DROPS::Point2DCL> grad_gamma_k ((DROPS::MakePoint2D(grad_gamma_stern_k[0], grad_gamma_stern_k[1])
															-grad_gamma_stern_k[2]* vektor));
						


								// Projeziere grad_gamma_i,j,k per hand auf t
								// t= Knoten j - Knoten i (3. Komponente ist 0)
								//DROPS::Point3DCL tfull= (iter->getVertex(jj))->GetCoordinates() - (iter->getVertex(ii))->GetCoordinates();
								//DROPS::Point2DCL t= DROPS::MakePoint2D(tfull[0], tfull[1]);
								//DROPS::Quad9Space1DCL<DROPS::Point2DCL> tang(t/norm(t));


								//grad_gamma_i= dot(grad_gamma_i, tang)* tang;
								//grad_gamma_j= dot(grad_gamma_j, tang)* tang;
								//grad_gamma_k= dot(grad_gamma_k, tang)* tang;

						
								//double value1= x[(iter->getVertex(ii)->Unknowns(0))];
								//double value2= x[(iter->getVertex(jj)->Unknowns(0))];
						

								//DROPS::Point2DCL graduhvergleich= (value2-value1)/norm(t) * t/norm(t);
								//std::cout << "differenz= " << std::abs(value2-value1) << "\n";
								//std::cout << "abstand j - i= " << norm(t) << "\n";
						
								//std::cout << "uhgradvergleich= " << graduhvergleich << "norm= " << norm(graduhvergleich)<< " Richtung= "<< t/norm(t) << "\n";
						

								/*
								std::cout << "Funktionswerte von grad_gamma_i:\n"; 
								for (int r=0; r< grad_gamma_i.size(); r++){
									std::cout << r << " " << grad_gamma_i[r] << " " << norm(grad_gamma_i[r]) << " Richtung= "<< grad_gamma_i[r]/norm(grad_gamma_i[r]) << "\n";
								}
								*/

						
								/*
								std::cout << "Funktionswerte von grad_gamma_j:\n"; 
								for (int p=0; p< grad_gamma_j.size(); p++){
									std::cout << p << " " << grad_gamma_j[p] << " " << norm(grad_gamma_j[p]) << " Richtung= "<< grad_gamma_j[p]/norm(grad_gamma_j[p]) << "\n";
								}

						

								std::cout << "Funktionswerte von grad_gamma_k:\n"; 
								for (int l=0; l< grad_gamma_k.size(); l++){
									std::cout << l << " " << grad_gamma_k[l] << " " << norm(grad_gamma_k[l]) << " Richtung= "<< grad_gamma_k[l]/norm(grad_gamma_k[l]) <<"\n";
								}
								*/
						
						

								// für \Grad_\Gamma_h(T)u_h 
								DROPS::Quad9Space1DCL<DROPS::Point2DCL> gu_h((x[(iter->getVertex(ii))->Unknowns(0)]*grad_gamma_i + x[(iter->getVertex(jj))->Unknowns(0)]*grad_gamma_j
																	+ x[(iter->getVertex(kk))->Unknowns(0)]*grad_gamma_k));
								/*
								std::cout << "Funktionswerte von gu_h:\n"; 
								for (int m=0; m< gu_h.size();m++){
									std::cout << m << " " << gu_h[m] << " " << norm(gu_h[m]) << " Richtung= " << gu_h[m]/norm(gu_h[m]) << "\n";
								}
								*/

						

								// für \Grad_\Gamma_h(T)ue (hier ist ue=u, da u konstant ist)
								DROPS::Quad9Space1DCL<DROPS::Point2DCL> gue(line, uexaktgrad);


								// für (grad(u_h)-grad(ue))^2
								DROPS::Quad9Space1DCL<DROPS::Point2DCL> gueminusuh(gue-gu_h);
								DROPS::Quad9Space1DCL<double> res22( dot(gueminusuh, gueminusuh) );

								sum22+= res22.quad(absdet);


				
				}
			}

		}

		// Berechnung des Diskretisierungsfehler auf \Gamma_{\ast,h}
		double suml2ganz=0.0;
		double sumh1ganz=0.0;
		for(DROPS::MySpaceTimeGridCL::MyTriangIterator toto= newstg.GetTriangBegin(); toto!=newstg.GetTriangEnd(); toto++){

			double absdet= DROPS::GetAbsDet(*toto);
			//Werte der diskreten Lösung an Knoten
			double values[3]= {x[toto->getVertex(0)->Unknowns(0)], x[toto->getVertex(1)->Unknowns(0)], x[toto->getVertex(2)->Unknowns(0)],  };

			//L2-Fehler -----------------------------------------
			// für ue:
			DROPS::Quad5SpaceTimeInterfaceCL<> ue(*toto, uexakt);
			// für u_h
			DROPS::MyLocalP1CL<> u_hhelp(*toto, values);
			DROPS::Quad5SpaceTimeInterfaceCL<> u_h(u_hhelp);
			// für Fehler
			DROPS::Quad5SpaceTimeInterfaceCL<> minusl2((ue - u_h)*(ue - u_h));
			suml2ganz+= minusl2.quad(absdet);

			//H1-Seminorm-Fehler --------------------------------------
			DROPS::SMatrixCL<3,2> M= DROPS::GetGradTrafoMatrix(*toto);
			DROPS::Point3DCL grad_gamma_stern_0(M*GradRef[0]);
			DROPS::Point3DCL grad_gamma_stern_1(M*GradRef[1]);
			DROPS::Point3DCL grad_gamma_stern_2(M*GradRef[2]);
			//für grad_T(u_h)
			DROPS::Point3DCL gradu_h= values[0]*grad_gamma_stern_0 + values[1]*grad_gamma_stern_1 + values[2]*grad_gamma_stern_2;
			DROPS::Quad5SpaceTimeInterfaceCL<DROPS::Point3DCL> grad_Tu_h(gradu_h);
			//für grad_T(ue)
			DROPS::Quad5SpaceTimeInterfaceCL<DROPS::Point3DCL> gradue(*toto, uexaktgrad3);
			DROPS::Point3DCL nhut= GetNhut(*toto);
			DROPS::Quad5SpaceTimeInterfaceCL<DROPS::Point3DCL> grad_Tue( gradue - nhut* dot(nhut, gradue));
			DROPS::Quad5SpaceTimeInterfaceCL<DROPS::Point3DCL> minus(grad_Tu_h - grad_Tue);
			DROPS::Quad5SpaceTimeInterfaceCL<> minush1( dot(minus, minus));
			sumh1ganz+= minush1.quad(absdet);

		}


		//double h=6e-3/std::pow(2.0,i);
		double h= 6.0e-3/(initial* std::pow(2.0,i));
		//double h= 3e-3/(1.0*i);
		//std::cout << "sum= " << sum << "\n";
		//std::cout << "sum2= " << sum2 << "\n";
		//std::cout << "sum22= " << sum22 << "\n";
		double l2error= std::sqrt(sum);
		double h1error= std::sqrt(sum2);
		double h1errormethode2= std::sqrt(sum22);
		double l2errorganz= std::sqrt(suml2ganz);
		double h1errorganz= std::sqrt(sumh1ganz);

		file << h << "\t" << l2error << "\t" << h1error << "\t" << h1errormethode2 << "\t" << l2errorganz << "\t" << h1errorganz << "\n";

		std::cout << "h= " << h << "\tL2-Fehler= " << l2error << "\tH1-Fehler (Methode 1)= "<< h1error << "\tH1-Fehler (Methode2)= " << h1errormethode2 << "\tL2-Fehler ganz= " << l2errorganz << "\tH1-Fehler ganz= " << h1errorganz << "\n";


		}// Ende: Fehlerberechnung für letzten TimeSlab





	

	DROPS::MyVTKOutCL* vtk = NULL;
	vtk = new DROPS::MyVTKOutCL (newstg, x, x, std::string("vtk"), std::string("Diskrete") );
	vtk->PutGeometryAndSolution();

	DROPS::VectorCL interpolwerte(fe.GetNumUnknowns());
	int zahlzahl=0;
	for(DROPS::MySpaceTimeGridCL::MyVertexIterator titi= newstg.GetVerticesBegin(); titi!=newstg.GetVerticesEnd(); titi++){
		interpolwerte[zahlzahl]= uexakt(titi->GetCoordinates(), 0.0);
		zahlzahl++;
	}
	DROPS::MyVTKOutCL* vtk2 = NULL;
	vtk2 = new DROPS::MyVTKOutCL (newstg, interpolwerte, interpolwerte, std::string("vtk"), std::string("Interpolierte") );
	vtk2->PutGeometryAndSolution();


	delete mg;
	delete vtk;
	delete vtk2;
	delete Mmatsurf2;
	delete Smatsurf2;
	delete Umatsurf2;
	delete rhskoefftrafo;


	}// ende gehe TimeSlabs durch
	
	
	/*
	//_______________________________________________NUR EIN TIME SLAB__________________________________________________________

		std::cout << "Methode mit 1 Time Slab\n";
		file << "Methode mit 1 Time Slab\n";
		DROPS::MySpaceTimeGridCL newstg;

		// Punkte für den BrickBuilder		
		DROPS::Point3DCL origin = DROPS::MakePoint3D (0.0, 0.0, 0.0);   
		DROPS::Point3DCL e1 = DROPS::MakePoint3D (6.0e-3, 0.0, 0.0);     
		DROPS::Point3DCL e2 = DROPS::MakePoint3D (0.0, 6.0e-3, 0.0);     
		DROPS::Point3DCL e3 = DROPS::MakePoint3D (0.0, 0.0, 6.0e-3);
		
		// baue MultiGrid
		DROPS::BrickBuilderCL builder( origin, e1, e2, e3, initial, initial, N ); // (10 unterteilungen in x,y,z-richtung)
		DROPS::MultiGridCL* mg;
		mg = new DROPS::MultiGridCL(builder);

		//adaptive Verfeinerung
		double width= 1e-25;
		int clevel= 0; // gröbstes gitter hat level 0
		int flevel1= i; // feinstes gitter hat level i (da width= hinitial/2)

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

		std::cout << "Anzahl Unbekannte= " << fe.GetNumUnknowns() << "\n";
		file << "Anzahl Unbekannte= " << fe.GetNumUnknowns() << "\n";



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
		
		DROPS::RHSAccumulatorP1CL rhs(newstg, bndDataf, bnd, quellf, rhskoefftrafo, fe, wf, &masssurf2, &stiffsurf2, &derivsurf2);

		// accumulate:
		std::cout << "accumulate...\n";
		DROPS::accumulate_matrix(newstg, &masssurf2, &stiffsurf2, &derivsurf2, &rhs);

		// Löse LGS (M+S+U)x=rhs
		DROPS::MatrixCL MSU;
		MSU.LinComb(1.0, (*Mmatsurf2), 1.0, (*Smatsurf2), 1.0, (*Umatsurf2)); // M+S+U
		
		DROPS::VectorCL y( fe.GetNumUnknowns());

		//Direkte Löser
		DROPS::DirectNonSymmSolverCL MSUSolver(MSU);
		std::cout << "löse lgs \n";
		MSUSolver.Solve(MSU, y, rhskoefftrafo->Data);

		//Visualisierung
		DROPS::MyVTKOutCL* vtk = NULL;
		vtk = new DROPS::MyVTKOutCL (newstg, y, y, std::string("vtk"), std::string("Diskrete") );
		vtk->PutGeometryAndSolution();

		//Berechne L2 und H1 Fehler
		sg2.BuildSpaceGrid(newstg, true);
		fe.NumberingVerticesSG(sg2, newstg);// Gitter sg wird konsistent zu Gitter newstg nummeriert

		std::cout << "berechne fehler\n";
		DROPS::MySpaceGridCL::MyLineIterator line=sg2.GetLineBegin();
		double sumneu=0.0;
		double sumneu2=0.0;
			for(line; line!=sg2.GetLineEnd(); line++){
	
				double absdet= DROPS::GetAbsDet1D(*line);
				//std::cout << "absdet= " << absdet << "\n";

				//______________L2____________________________________________-
	
				// für ue (hier ue=uexakt, da uexakt konstant)
				DROPS::Quad9Space1DCL<> ue(*line, uexakt);
				
				// für u_h
				// values sind Funktionswerte von u_h an Endknoten der Linie line
				double values[2]= {y[(line->getVertex(0))->Unknowns(0)], y[(line->getVertex(1))->Unknowns(0)]};
				DROPS::Quad9Space1DCL<> uh(values);
	

				// für (ue-uh)^2
				DROPS::Quad9Space1DCL<> res((ue-uh)*(ue-uh));
	
				sumneu+= res.quad(absdet);
				

				//______________H1____________________________________________-

				// Tangentialvektor tang
				DROPS::Point2DCL t= DROPS::MakePoint2D((line->getVertex(1))->GetCoordinates()[0]- (line->getVertex(0))->GetCoordinates()[0],
											(line->getVertex(1))->GetCoordinates()[1]- (line->getVertex(0))->GetCoordinates()[1]);
				DROPS::Point2DCL tang= t/norm(t);
				//std::cout << "Tangentialvektor= " << tang << "\n";
				DROPS::Quad9Space1DCL<DROPS::Point2DCL> tangquad(tang);

				// für tangentialgrad(ue) (hier ue=uexakt)
				DROPS::Quad9Space1DCL<DROPS::Point2DCL> uegrad(*line, uexaktgrad);
				DROPS::Quad9Space1DCL<DROPS::Point2DCL> uetanggrad(dot(tangquad,uegrad)* tangquad);

				// für tangentialgrad(u_h)
				// ableitung= ableitung von u_h auf linie (absdet= norm(p-q))
				double ableitung= (values[1]-values[0])/norm(t);
				DROPS::Quad9Space1DCL<DROPS::Point2DCL> uhgrad(ableitung*tang);
	

				// für (grad(ue)-grad(uh))^2
				DROPS::Quad9Space1DCL<DROPS::Point2DCL> minus(uetanggrad - uhgrad);
				DROPS::Quad9Space1DCL<> res2(dot(minus,minus));
	
				sumneu2+= res2.quad(absdet);
	

			}

			double h=(6.0e-3)/(initial*std::pow(2.0, i));
			sumneu= std::sqrt(sumneu);
			sumneu2= std::sqrt(sumneu2);
			std::cout << "h= " << h << "\tL2-Fehler= " << sumneu << "\tH1-Fehler= " << sumneu2 << "\n";
			file << h << "\t" << sumneu << "\t" << sumneu2 << "\n";


			delete mg;	
			delete Mmatsurf2;
			delete Smatsurf2;
			delete Umatsurf2;
			delete rhskoefftrafo;

	/*
	//______________________________________VERGLEICH ZWISCHEN U_H1 UND U_H2_________________________________
			std::cout << "Vergleich zwischen uh1 und uh2\n";
			file << "Vergleich zwischen uh1 und uh2\n";
			//Vergleich zwischen sg und sg2 (sollten am besten gleich sein)
			// teste, ob Knoten in gittern gleich sind
			DROPS::MySpaceGridCL::MyVertexIterator knoten1=sg.GetVerticesBegin();
			DROPS::MySpaceGridCL::MyVertexIterator knoten2=sg2.GetVerticesBegin();
			bool gittergleich= true;
			for(knoten1; knoten1!=sg.GetVerticesEnd(); knoten1++){
				bool knotengleich=false;
				for(knoten2; knoten2!=sg2.GetVerticesEnd(); knoten2++){
					if(norm(knoten1->GetCoordinates() - knoten2->GetCoordinates()) < 0.1e-14) {
						knotengleich=true;
						break;
					}
				}
				if(knotengleich==false) {
					gittergleich=false;
					break;
				}
				
			}
			if(gittergleich) {
				std::cout << "Raumgitter sind gleich!\n";
				file << "Raumgitter sind gleich!\n";
			}else{
				std::cout << "Raumgitter sind nicht gleich!\n";
				file << "Raumgitter sind nicht gleich!\n";
			}

			double sumneuneu=0.0;
			double sumneuneu2=0.0;
			DROPS::MySpaceGridCL::MyLineIterator lin=sg.GetLineBegin();
			for(lin; lin!=sg.GetLineEnd(); lin++){

				double absdet= DROPS::GetAbsDet1D(*lin);
				
				//______________L2____________________________________________-
	
				
				// für u_h1
				// values sind Funktionswerte von u_h an Endknoten der Linie line
				//std::cout << "baue values1\n";
				double values1[2]= {xx[(lin->getVertex(0))->Unknowns(0)], xx[(lin->getVertex(1))->Unknowns(0)]};
				DROPS::Quad9Space1DCL<> uh1(values1);

				//für u_h2
				// values sind Funktionswerte von u_h an Endknoten der Linie line
				//std::cout << "baue values2\n";
				double values2[2]= {y[(lin->getVertex(0))->Unknowns(0)], y[(lin->getVertex(1))->Unknowns(0)]};
				DROPS::Quad9Space1DCL<> uh2(values2);
	

				// für (uh1-uh2)^2
				DROPS::Quad9Space1DCL<> resres((uh1-uh2)*(uh1-uh2));
	
				sumneuneu+= resres.quad(absdet);
				

				//______________H1____________________________________________-

					

				// für tangentialgrad(u_h1)
				// ableitung= ableitung von u_h auf linie (absdet= norm(p-q))
				double ableitung1= (values1[1]-values1[0])/absdet;
				//file << "ableitung1= " << ableitung1 << "\n";
				DROPS::Quad9Space1DCL<> uh1grad(ableitung1);

				// für tangentialgrad(u_h2)
				// ableitung= ableitung von u_h auf linie (absdet= norm(p-q))
				double ableitung2= (values2[1]-values2[0])/absdet;
				//file << "ableitung2= " << ableitung2 << "\n";
				DROPS::Quad9Space1DCL<> uh2grad(ableitung2);
	

				// für (grad(uh1)-grad(uh2))^2
				DROPS::Quad9Space1DCL<> minus2(uh1grad- uh2grad);
				DROPS::Quad9Space1DCL<> resres2(minus2*minus2);
	
				sumneuneu2+= resres2.quad(absdet);
			}

			//double h=(6.0e-3)/(initial*std::pow(2.0, i));
			sumneuneu= std::sqrt(sumneuneu);
			sumneuneu2= std::sqrt(sumneuneu2);
			std::cout << "h= " << h << "\tVergleich(L2)= " << sumneuneu << "\tVergleich(H1)= " << sumneuneu2 << "\n";
			file << h << "\t" << sumneuneu << "\t" << sumneuneu2 << "\n\n";
			*/

}// Ende: Verfeinerung (for int i=0...)

file.close();
return 0;

}