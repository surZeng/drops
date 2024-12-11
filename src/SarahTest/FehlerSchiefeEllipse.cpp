// L2-Fehler auf \Gamma_h(T) für schiefe Ellipse

#include "SarahTest/matrixaccu.h"
#include "SarahTest/vtkoutput.h"
#include "levelset/levelset.h"
#include "levelset/surfacetension.h"
#include "SarahTest/myadap.h"
// für Direkte Löser
#include "num/directsolver.h"
#include "num/precond.h"
#include "num/krylovsolver.h"
#include <fstream>

double sigma= 0.022; // sigma aus film.json
double sigmaf (const DROPS::Point3DCL&, double) { return sigma; } // sigma-funktion (pointer-to-function)

// Level-Set-Funktion ...
double lsetf (const DROPS::Point3DCL& P, double t) {return 0.0625*(P[0]-P[2]-1.2e-3)*(P[0]-P[2]-1.2e-3) + (P[1]-3.0e-3)*(P[1]-3.0e-3) - 1.0e-6; }

// Krümmung
double kappaf (const DROPS::Point3DCL& P, double t) {
	if(P[0]== P[2]+ 1.2e-3 && P[1]== 3.0e-3)
		return 0.0;
	else 
		return (1./16. *(P[0]-P[2]-1.2e-3)*(P[0]-P[2]-1.2e-3) + (P[1]-3.0e-3)*(P[1]-3.0e-3))/ (std::pow(1./32. *(P[0]-P[2]-1.2e-3)*(P[0]-P[2]-1.2e-3) + 4.*(P[1]-3.0e-3)*(P[1]-3.0e-3), 1.5));

}

// Koeffizient alpha(s,t)
double alphaf (const DROPS::Point3DCL& P, double t) {return 0.0001;}


// exakte Lösung u
/*
double uexakt(const DROPS::Point3DCL& P, double t) {return std::sin((P[0]-P[2]-1.2e-3)*(P[0]-P[2]-1.2e-3) + (P[1]-3.0e-3)*(P[1]-3.0e-3));}
*/
//für diese uexakt 
//hier ist projektion p vorgeschaltet worden, das heißt uexakt=ue
double uexakt(const DROPS::Point3DCL& Q, double t) {
	// für Projektion p(x,y,t)= (x,y,t) - d(x,y,t)*ngamma
	DROPS::Point3DCL P;
	double xterm= Q[0]-Q[2]-1.2e-3;
	double yterm= Q[1]-3.0e-3;
	DROPS::Point2DCL ngamma= DROPS::MakePoint2D(1./8.*xterm, 2.*yterm)/(std::sqrt(1./64.*xterm*xterm + 4.*yterm*yterm));
	//quadratische Gleichung aus levelsetf(p(x,y,t))=0 -> a*d^2 + b*d + c=0
	double a= (1./1024.*xterm*xterm + 4.*yterm*yterm)/(1./64.*xterm*xterm + 4.*yterm*yterm);
	double b= -std::sqrt(1./64.*xterm*xterm + 4.*yterm*yterm);
	double c= 1./16.*xterm*xterm + yterm*yterm - 1.0e-6;
	if(a==0){
		return 0; 
	}else{
		if(b*b - 4.*a*c < 0){ 
			return 0;
		}else{
			double d1= (-b + std::sqrt(b*b - 4.*a*c))/(2.*a);
			double d2= (-b - std::sqrt(b*b - 4.*a*c))/(2.*a);
			if(std::fabs(d1)<=std::fabs(d2)){
				P=DROPS::MakePoint3D(Q[0]-d1*ngamma[0], Q[1]-d1*ngamma[1], Q[2]);
			}else{
				P=DROPS::MakePoint3D(Q[0]-d2*ngamma[0], Q[1]-d2*ngamma[1], Q[2]);
			}
		}
	}

	// jetzt ue=u(P)
	return (P[1]-3.0e-3)/ std::sqrt((P[0]-P[2]-1.2e-3)*(P[0]-P[2]-1.2e-3) + (P[1]-3.0e-3)*(P[1]-3.0e-3));
}


// Tangentialgradient an \Gamma(T) der exakten Lösung 
/*
DROPS::Point2DCL uexaktgrad(const DROPS::Point3DCL& P, double t) {
	double faktor= std::cos ((P[0]-P[2]-1.2e-3)*(P[0]-P[2]-1.2e-3) + (P[1]-3.0e-3)*(P[1]-3.0e-3));
	return 2.0*faktor * DROPS::MakePoint2D(P[0]-P[2]-1.2e-3, P[1]-3.0e-3);
}
*/

//bei uexaktgrad wird grad(u) berechnet und nicht grad(ue), was es eigentlich sein müsste. das geht aber nicht
//weil ich keine geschlossene form für dist-funnktion d habe, die brauche ich aber für Dp.
DROPS::Point2DCL uexaktgrad(const DROPS::Point3DCL& P, double t){
	double faktor= 1./std::pow((P[0]-P[2]-1.2e-3)*(P[0]-P[2]-1.2e-3) + (P[1]-3.0e-3)*(P[1]-3.0e-3), 1.5);
	return faktor* DROPS::MakePoint2D(-(P[0]-P[2]-1.2e-3)*(P[1]-3.0e-3), (P[0]-P[2]-1.2e-3)*(P[0]-P[2]-1.2e-3));
}


// Geschwindigkeitsfeld w
DROPS::Point2DCL wf(const DROPS::Point3DCL& P, double t){ return DROPS::MakePoint2D(1.0, 0.0); }


// div(w)
double divwf (const DROPS::Point3DCL& P, double t) { 	
	return 0.0; }


// grad(w)
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




// Rand-Daten
//double bndDataf (const DROPS::Point3DCL& P, double t) { return std::sin((P[0]-1.2e-3)*(P[0]-1.2e-3) + (P[1]-3.0e-3)*(P[1]-3.0e-3)); }
//schalte auch bei bndDataf die Projektion p vor
double bndDataf (const DROPS::Point3DCL& Q, double t) {
	// für Projektion p(x,y,t)= (x,y,t) - d(x,y,t)*ngamma
	DROPS::Point3DCL P;
	double xterm= Q[0]-Q[2]-1.2e-3;
	double yterm= Q[1]-3.0e-3;
	DROPS::Point2DCL ngamma= DROPS::MakePoint2D(1./8.*xterm, 2.*yterm)/(std::sqrt(1./64.*xterm*xterm + 4.*yterm*yterm));
	//quadratische Gleichung aus levelsetf(p(x,y,t))=0 -> a*d^2 + b*d + c=0
	double a= (1./1024.*xterm*xterm + 4.*yterm*yterm)/(1./64.*xterm*xterm + 4.*yterm*yterm);
	double b= -std::sqrt(1./64.*xterm*xterm + 4.*yterm*yterm);
	double c= 1./16.*xterm*xterm + yterm*yterm - 1.0e-6;
	if(a==0){
		return 0; 
	}else{
		if(b*b - 4.*a*c < 0){ 
			return 0;
		}else{
			double d1= (-b + std::sqrt(b*b - 4.*a*c))/(2.*a);
			double d2= (-b - std::sqrt(b*b - 4.*a*c))/(2.*a);
			if(std::fabs(d1)<=std::fabs(d2)){
				P=DROPS::MakePoint3D(Q[0]-d1*ngamma[0], Q[1]-d1*ngamma[1], Q[2]);
			}else{
				P=DROPS::MakePoint3D(Q[0]-d2*ngamma[0], Q[1]-d2*ngamma[1], Q[2]);
			}
		}
	}

	//jetzt eigentliche randdaten
	return (P[1]-3.0e-3)/ std::sqrt((P[0]-1.2e-3)*(P[0]-1.2e-3) + (P[1]-3.0e-3)*(P[1]-3.0e-3));
}


// Quellterm f für rechte Seite
//FALSCH, habe hier aus versehen div(grad_gamma(u)) genommen anstatt div_gamma(grad_gamma(u))
/*
double quellf (const DROPS::Point3DCL& P, double t) { 
	double xterm= (P[0]-P[2]-3./2500.);
	double yterm= (P[1]-3./1000.);
	double xterm2= xterm*xterm;
	double yterm2= yterm*yterm;

	double term1= 2.*(std::sin)(xterm2+yterm2)*(2.*P[0]-2.*P[2]-3./1250.)*xterm-4.*(std::cos)(xterm2+yterm2);
	double term2= -((std::sin)(xterm2+yterm2)*(2.*P[0]-2.*P[2]-3./1250.)/(1./8.*xterm*(1./8.*P[0]-1./8.*P[2]-3./20000.)+2.*yterm*(2.*P[1]-3./500.))*(1./32.*xterm*xterm2+1./2.*xterm*yterm2));
	double term3= -((std::cos)(xterm2+yterm2)/(1./8.*xterm*(1./8.*P[0]-1./8.*P[2]-3./20000.)+2.*yterm*(2.*P[1]-3./500.))*(1./8.*xterm*(1./8.*P[0]-1./8.*P[2]-3./20000.)+2.*yterm*(2.*P[1]-3./500.))*(1./32.*xterm*xterm2+1./2.*xterm*yterm2)*(1./32.*P[0]-1./32.*P[2]-3./80000.));
	double term4= (std::cos)(xterm2+yterm2)/(1./8.*xterm*(1./8.*P[0]-1./8.*P[2]-3./20000.)+2.*yterm*(2.*P[1]-3./500.))*(3./32.*xterm2+1./2.*yterm2);
	double term5= 2.*(std::sin)(xterm2+yterm2)*(2.*P[1]-3./500.)*yterm;
	double term6= -((std::sin)(xterm2+yterm2)*(2.*P[1]-3./500.)/(1./8.*xterm*(1./8.*P[0]-1./8.*P[2]-3./20000.)+2.*yterm*(2.*P[1]-3./500.))*(8.*yterm*yterm2+1./2.*xterm2*yterm));
	double term7= -((std::cos)(xterm2+yterm2)/(1./8.*xterm*(1./8.*P[0]-1./8.*P[2]-3./20000.)+2.*yterm*(2.*P[1]-3./500.))*(1./8.*xterm*(1./8.*P[0]-1./8.*P[2]-3./20000.)+2.*yterm*(2.*P[1]-3./500.))*(8.*yterm*yterm2+1./2.*xterm2*yterm)*(8.*P[1]-3./125.));
	double term8= (std::cos)(xterm2+yterm2)/(1./8.*xterm*(1./8.*P[0]-1./8.*P[2]-3./20000.)+2.*yterm*(2.*P[1]-3./500.))*(24.*yterm2+1./2.*xterm2);

	return (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8) ;
}
*/

//hier ist projektion p vorgeschaltet, das heißt quellf=fe
double quellf (const DROPS::Point3DCL& Q, double t){
	// für Projektion p(x,y,t)= (x,y,t) - d(x,y,t)*ngamma
	DROPS::Point3DCL P;
	double xterm= Q[0]-Q[2]-1.2e-3;
	double yterm= Q[1]-3.0e-3;
	DROPS::Point2DCL ngamma= DROPS::MakePoint2D(1./8.*xterm, 2.*yterm)/(std::sqrt(1./64.*xterm*xterm + 4.*yterm*yterm));
	//quadratische Gleichung aus levelsetf(p(x,y,t))=0 -> a*d^2 + b*d + c=0
	double a= (1./1024.*xterm*xterm + 4.*yterm*yterm)/(1./64.*xterm*xterm + 4.*yterm*yterm);
	double b= -std::sqrt(1./64.*xterm*xterm + 4.*yterm*yterm);
	double c= 1./16.*xterm*xterm + yterm*yterm - 1.0e-6;
	if(a==0){
		return 0; 
	}else{
		if(b*b - 4.*a*c < 0){ 
			return 0;
		}else{
			double d1= (-b + std::sqrt(b*b - 4.*a*c))/(2.*a);
			double d2= (-b - std::sqrt(b*b - 4.*a*c))/(2.*a);
			if(std::fabs(d1)<=std::fabs(d2)){
				P=DROPS::MakePoint3D(Q[0]-d1*ngamma[0], Q[1]-d1*ngamma[1], Q[2]);
			}else{
				P=DROPS::MakePoint3D(Q[0]-d2*ngamma[0], Q[1]-d2*ngamma[1], Q[2]);
			}
		}
	}

	//ab hier eigentliche Funktion, angewendet auf Projektion P
	xterm= P[0]-P[2]-1.2e-3;
	yterm= P[1]-3.0e-3;
	double xterm2= xterm*xterm;
	double yterm2= yterm*yterm;
	double xterm3= xterm*xterm2;
	double yterm3= yterm*yterm2;
	double xterm4= xterm*xterm3;
	double yterm4= yterm*yterm3;
	double xterm6= xterm3*xterm3;
	double yterm6= yterm3*yterm3;
	double xterm8= xterm4*xterm4;
	double yterm8= yterm4*yterm4;

	double zaehler=1./625000.*(-38550000000.*xterm*xterm2*xterm3*yterm*yterm2-12000000000.*xterm2*xterm3*yterm*yterm2*P[0]+12000000000.*xterm2*xterm3*yterm*yterm2*P[2]-3360000000.*xterm*xterm4*yterm*yterm2*P[0]+3360000000.*xterm*xterm4*yterm*yterm2*P[2]-480000000.*xterm*xterm4*yterm2*P[0]*P[1]+480000000.*xterm*xterm4*yterm2*P[1]*P[2]-118080000000.*xterm*xterm2*yterm*yterm4*P[0]+118080000000.*xterm*xterm2*yterm*yterm4*P[2]-245760000000.*xterm*xterm2*yterm4*P[0]*P[1]+245760000000.*xterm*xterm2*yterm4*P[1]*P[2]
					+77250000000.*xterm*xterm2*xterm3*yterm2*P[1]+14400000.*xterm2*xterm3*yterm*yterm2-30000000.*xterm*xterm6*yterm*P[0]+30000000.*xterm*xterm6*yterm*P[2]-31457280000000.*xterm*yterm6*P[0]*P[1]+31457280000000.*xterm*yterm6*P[1]*P[2]-150000000.*xterm*xterm4*xterm3*yterm+4032000.*xterm*xterm4*yterm*yterm2+1440000.*xterm*xterm4*yterm2*P[0]+576000.*xterm*xterm4*yterm2*P[1]-1440000.*xterm*xterm4*yterm2*P[2]+141696000.*xterm*xterm2*yterm*yterm4+737280000.*xterm*xterm2*yterm4*P[0]+294912000.*xterm*xterm2*yterm4*P[1]-737280000.*xterm*xterm2*yterm4*P[2]-38400000000.*xterm*xterm3*yterm*yterm4-1848000000000.*xterm3*yterm*yterm4*P[0]
					+1848000000000.*xterm3*yterm*yterm4*P[2]+30000000.*xterm2*xterm4*yterm*P[1]*P[1]+1966080000000.*xterm2*yterm*yterm4*P[1]*P[1]+330000000.*xterm2*xterm4*yterm*yterm2-180000.*xterm2*xterm4*yterm*P[1]-11796480000.*xterm2*yterm*yterm4*P[1]-231750000.*xterm*xterm2*xterm3*yterm2+15360000000.*xterm4*yterm*yterm2*P[1]*P[1]-92160000.*xterm4*yterm*yterm2*P[1]+192000000000.*xterm*xterm3*yterm4*P[1]+138240.*xterm4*yterm*yterm2-576000000.*xterm*xterm3*yterm4+10000000.*xterm2*xterm6*yterm
					+17694720.*xterm2*yterm*yterm4+270.*xterm2*xterm4*yterm+1269760000000.*xterm2*yterm*yterm6-76130000000.*xterm6*yterm2*P[1]+12960000000.*xterm4*yterm4*P[1]-1728.*xterm*xterm4*yterm2+2217600000.*xterm3*yterm*yterm4+10567680000000.*xterm2*yterm6*P[1]+36000.*xterm*xterm6*yterm-94371840000.*xterm*yterm6*P[2]+84300000000.*xterm6*yterm*yterm2-884736.*xterm*xterm2*yterm4+37748736000.*xterm*yterm6*P[1]+94371840000.*xterm*yterm6*P[0]+1354080000000.*xterm4*yterm*yterm4-5625.*xterm8-31457280000.*yterm8+300000000.*xterm8*yterm+1875000.*xterm8*P[1]+10485760000000.*yterm8*P[1]-113246208.*xterm*yterm6+228390000.*xterm6*yterm2-38880000.*xterm4*yterm4-31703040000.*xterm2*yterm6);
	
	double nenner= std::pow((xterm2+yterm2), (5./2.))*std::pow((xterm2+256.*yterm2), 3.);

	return zaehler/nenner;


}


int main()
{

std::fstream file;
file.open("FehlerSchiefeEllipseNeu3_10TS_1t_cgleich.dat", std::ios::out);

int N= 10; // Anzahl der TimeSlabs



for(int i=0; i<=3; i++){ // Verfeinerungen

	DROPS::boundarymap bnd;
	for(int k=0; k<N; k++){ // gehe N TimeSlabs durch

		if(k==0) std::cout << "Schleifendurchlauf " << i << "\n";

		DROPS::MySpaceTimeGridCL newstg;

		// Punkte für den BrickBuilder
		DROPS::Point3DCL origin = DROPS::MakePoint3D (-12.0e-3, 0.0, (k/(double)(N))*6.0e-3);   
		DROPS::Point3DCL e1 = DROPS::MakePoint3D (24.0e-3, 0.0, 0.0);     
		DROPS::Point3DCL e2 = DROPS::MakePoint3D (0.0, 6.0e-3, 0.0);     
		DROPS::Point3DCL e3 = DROPS::MakePoint3D (0.0, 0.0, 6.0e-3/(double)(N));

		// baue MultiGrid
		int initial=4;
		//DROPS::BrickBuilderCL builder( origin, e1, e2, e3, 6*initial, 3*initial, 1*initial ); 
		DROPS::BrickBuilderCL builder( origin, e1, e2, e3, 4*initial, initial, 1 ); 
		DROPS::MultiGridCL* mg;
		mg = new DROPS::MultiGridCL(builder);

		//adaptive Verfeinerung
		double width= 1e-25;
		int clevel= 0; 
		int flevel1= i; 
		int flevel2= i+2;
		//double konst= 1.0/std::pow(2.0, i);
		double konst= 0.9;
		//double konst= 1.6; // (h<= c* 2/kappa -> c ist eigentlich 0.8)
		double initialmeshsize= 6.0e-3/initial;
		DROPS::MyAdapCL adap(*mg, width, clevel, flevel1, flevel2, konst, initialmeshsize);
		adap.MakeInitialTriang(lsetf, kappaf);



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





		// lege Matrizen an
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
			DROPS::MySpaceGridCL sg;
			sg.BuildSpaceGrid(newstg);
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


		



		//double h=6e-3/std::pow(2.0,i);
		double h= 6.0e-3/(initial* std::pow(2.0,i));
		//double h= 3e-3/(1.0*i);
		//std::cout << "sum= " << sum << "\n";
		//std::cout << "sum2= " << sum2 << "\n";
		//std::cout << "sum22= " << sum22 << "\n";
		double l2error= std::sqrt(sum);
		double h1error= std::sqrt(sum2);
		//double h1errormethode2= std::sqrt(sum22);

		file << h << "\t" << l2error << "\t" << h1error /* << "\t" << h1errormethode2/* <<"\t" << std::sqrt(l2error*l2error+ h1error*h1error)<< "\t"<< std::sqrt(l2error*l2error+ h1errormethode2*h1errormethode2)*/ <<"\n";

		std::cout << "h= " << h << "\tL2-Fehler= " << l2error << "\tH1-Fehler (Methode 1)= "<< h1error /*<< "\tH1-Fehler (Methode2)= " << h1errormethode2*/ <<"\n";


		}// Ende: Fehlerberechnung für letzten TimeSlab





	

	DROPS::MyVTKOutCL* vtk = NULL;
	vtk = new DROPS::MyVTKOutCL (newstg, x, x, std::string("vtk"), std::string("FehlerSchiefeEllipse") );
	vtk->PutGeometry();

	delete mg;
	delete Mmatsurf2;
	delete Smatsurf2;
	delete Umatsurf2;
	delete rhskoefftrafo;
	delete vtk;


	}// ende gehe TimeSlabs durch

	

}// Ende: Verfeinerung (for int i=0...)

file.close();
return 0;

}