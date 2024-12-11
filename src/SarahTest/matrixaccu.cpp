// Definitionen aus matrixaccu.h


#include "SarahTest/matrixaccu.h"

namespace DROPS
{

typedef Point2DCL (*vec2D_fct_ptr)(const Point3DCL& P, double t);// P[2] ist Zeit, t wird nicht gebraucht

bool IsOnBoundary(MySpaceTimeGridCL& stg, const MyVertexCL& vert)
{
	double eps= 0.1e-15;

	if(std::abs(vert.GetCoordinates()[2]-stg.Get_t_old()) < eps) return true;

	else return false;

}


bool IsOnUpperBoundary(MySpaceTimeGridCL& stg, const MyVertexCL& vert)
{
	double eps= 0.1e-15;
	if(std::abs(vert.GetCoordinates()[2]-stg.Get_t_new()) < eps ){
		return true;
	}else{ 
		return false;
	}

}

boundarymap BuildBoundaryMap(MySpaceTimeGridCL& stg, VectorCL& solution)
{
	boundarymap bnd;
	MySpaceTimeGridCL::MyVertexIterator it= stg.GetVerticesBegin();
	for(it; it!=stg.GetVerticesEnd(); it++){
		if( IsOnUpperBoundary(stg, *it ) ){
			bnd[it->GetCoordinates()]= solution[it->Unknowns(0)];
		} 
	}

	return bnd;

}

boundarymap::iterator FindInBoundaryMap(boundarymap& bound, const Point3DCL& p )
{
	boundarymap::iterator it= bound.begin();
	for(it; it!= bound.end(); it++){
		if( norm(it->first - p) < 0.1e-15) return it;
	}


}

void accumulate_matrix(MySpaceTimeGridCL& stg, MassAccumulatorP1CL * mass, StiffAccumulatorP1CL * stiff, DerivAccumulatorP1CL * deriv, RHSAccumulatorP1CL * rhs){

	// muss rechte Seite aufbauen, dabei muss ich darauf achten, wieviele Matrizen ich aufbaue
	
	if (mass) mass->begin_accumulation();
	if (stiff) stiff-> begin_accumulation();
	if (deriv) deriv-> begin_accumulation();
	if (rhs) rhs-> begin_accumulation();

	MySpaceTimeGridCL::const_MyTriangIterator it= stg.GetTriangBegin();
	int anzahldreiecke= 0;
	for(DROPS::MySpaceTimeGridCL::MyTriangIterator ot= stg.GetTriangBegin(); ot!=stg.GetTriangEnd(); ot++)
	{
		anzahldreiecke++;
	}

	int dreieck=0;
	for(it; it!=stg.GetTriangEnd(); it++){
		dreieck++;
		//std::cout << "baue matrizen auf dreieck " << dreieck << "von " << anzahldreiecke << "\n";
		if (mass) mass->visit(stg, *it);
		if (stiff) stiff -> visit(stg, *it);
		if (deriv) deriv -> visit(stg, *it);
		if (rhs) rhs ->visit(stg, *it);
	}

	if (mass) mass->finalize_accumulation();
	if (stiff) stiff -> finalize_accumulation();
	if (deriv) deriv -> finalize_accumulation();
	// rhs braucht kein finalize_accumulation()
}

// Funktion, die Determinante zu gegebenem Dreieck T berechnet (det= 2*|T|)
double GetAbsDet(const MyTriangleCL& triang)
{
	Point3DCL p= triang.getVertex(0)->GetCoordinates();
	Point3DCL q= triang.getVertex(1)->GetCoordinates();
	Point3DCL r= triang.getVertex(2)->GetCoordinates();

	double a= DROPS::norm(p-q);
	double b= DROPS::norm(p-r);
	double c= DROPS::norm(q-r);
	double s= 0.5*(a+b+c);

	//std::cout << "absdet(dreieck)= " << 2. * std::sqrt(s*(s-a)*(s-b)*(s-c)) << "\n"; 
	return   2. * std::sqrt(s*(s-a)*(s-b)*(s-c));

}

double GetAbsDet1D (const MyLineCL& line)
{
	return norm(  (line.getVertex(0))->GetCoordinates() - (line.getVertex(1))->GetCoordinates()  );

}

void GetN_gamma(const MyTriangleCL& triang, Point2DCL& ngamma, double& vh, Point2DCL& vhngamma)
{
	Point3DCL p= triang.getVertex(0)->GetCoordinates();
	Point3DCL q= triang.getVertex(1)->GetCoordinates();
	Point3DCL r= triang.getVertex(2)->GetCoordinates();

	Point3DCL a= q-p;
	Point3DCL b= r-p;

	// Kreuzprodukt ist Normale an triang
	Point3DCL nhut_unnormiert= MakePoint3D(a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]);
	Point3DCL nhut= nhut_unnormiert/ norm(nhut_unnormiert);

	

	ngamma= std::sqrt(1.+(nhut[2]*nhut[2])/(1.-nhut[2]*nhut[2])) * MakePoint2D(nhut[0], nhut[1]);
	//std::cout << "n_gamma= " << ngamma << "\n"; 

	vh= -nhut[2]/ (std::sqrt(nhut[0]*nhut[0] + nhut[1]*nhut[1]));
	//std::cout << "vh= " << vh << "\n";

	vhngamma= vh* ngamma;
	//std::cout << "vh*n_gamma= " << vhngamma << "\n";

	
}

//berechnet zu triang das lokale diskrete Geschwindigkeitsfeld w_h
void GetW_h(const MyTriangleCL& triang, Point2DCL& wh)
{
	Point3DCL p= triang.getVertex(0)->GetCoordinates();
	Point3DCL q= triang.getVertex(1)->GetCoordinates();
	Point3DCL r= triang.getVertex(2)->GetCoordinates();

	Point3DCL a= q-p;
	Point3DCL b= r-p;

	// Kreuzprodukt ist Normale an triang
	Point3DCL nhut_unnormiert= MakePoint3D(a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]);
	Point3DCL nhut= nhut_unnormiert/ norm(nhut_unnormiert);

	wh= (-nhut[2])/(nhut[0]*nhut[0] + nhut[1]*nhut[1])* MakePoint2D(nhut[0], nhut[1]);

}

// diese Funktion gibt A*(A^T*A)^[-1] zurück, wobei A die Matrix aus der Transformation T_hut -> T ist
SMatrixCL<3,2> GetGradTrafoMatrix ( const MyTriangleCL& triang)
{
	Point3DCL P= triang.getVertex(0)->GetCoordinates();
	Point3DCL Q= triang.getVertex(1)->GetCoordinates();
	Point3DCL R= triang.getVertex(2)->GetCoordinates();

	SMatrixCL<3,2> A; // A= (Q-P | R-P)
	A(0,0)= Q[0]-P[0];
	A(1,0)= Q[1]-P[1];
	A(2,0)= Q[2]-P[2];
	A(0,1)= R[0]-P[0];
	A(1,1)= R[1]-P[1];
	A(2,1)= R[2]-P[2];

	

	SMatrixCL<2,2> A_tA= GramMatrix(A);

	SMatrixCL<2,2> invA_tA;
	
	double det= (A_tA(0,0)*A_tA(1,1) - A_tA(0,1)*A_tA(1,0));
	/*
	if(std::abs(det) < 1.0e-35) {
		if(det >= 0) det= 1.0e-35;
		else		 det= -1.0e-35;
	}
	*/
	//if(det < 1.0e-16) std::cout << "det= " << det << "\n";

	// herkömmliches invertieren
	
	invA_tA(0,0)=  A_tA(1,1)/ det;
	//if(det < 1.0e-16) std::cout << "1. eintrag= " << invA_tA(0,0) << "\n";
	invA_tA(0,1)=  (-A_tA(0,1))/ det;
	invA_tA(1,0)=  (-A_tA(1,0))/ det;
	invA_tA(1,1)=  A_tA(0,0)/ det;
	

	// invertieren über givens-rotation
	/*
	double r= std::sqrt( A_tA(0,0)*A_tA(0,0) + A_tA(1,0)*A_tA(1,0) );
	double c= A_tA(0,0) / r;
	double s= std::sqrt(1-c*c);
	
	SMatrixCL<2,2> Qmat; // orthogonale Matrix Q
	Qmat(0,0)= c; Qmat(1,0)= -s; Qmat(0,1)= s; Qmat(1,1)=c;

	SMatrixCL<2,2> Rmat= Qmat * A_tA; // obere Dreiecksmatrix R

	// A_tA= Q^T*R -> A_tA*x=e_i äquivalent zu R*x= Q*e_i
	Point2DCL b1= Qmat * MakePoint2D(1.0, 0.0);
	Point2DCL b2= Qmat * MakePoint2D(0.0, 1.0);

	// löse R*x=Q*e_i durch Rückwärtseinsetzen
	std::cout << "Rmat(1,1)= " << std::setprecision(11) << Rmat(1,1) << "\n";
	double x2= b1[1]/ Rmat(1,1);
	double x1= (b1[0]-Rmat(0,1)*x2) / Rmat(0,0);

	double y2= b2[1]/ Rmat(1,1);
	double y1= (b2[0]-Rmat(0,1)*y2) / Rmat(0,0);

	// baue invA_tA auf
	invA_tA(0,0)= x1; 
	invA_tA(1,0)= x2;
	invA_tA(0,1)= y1; 
	invA_tA(1,1)= y2;
	*/


	SMatrixCL<3,2> ret= A*invA_tA;
	//if(det < 1.0e-16) std::cout << "1. eintrag ret= " << ret(0,0) << "\n";

	return ret;

}




//****************************************************
//
//  Klasse: MassAccumulatorP1CL
//
//****************************************************


MassAccumulatorP1CL::MassAccumulatorP1CL(const MySpaceTimeGridCL & stg, /*const BndDataCL<> * BndData,*/ MatrixCL* Mmat, /*VecDescCL* b,*/
        SpaceTimeInterfaceP1CL& fe, scalar_fct_ptr f, vec2D_fct_ptr w, mat2_fct_ptr gradw, scalar_fct_ptr divw, bool surfactant, const double t_):
        stg_(stg), /*BndData_(BndData),*/ Mmat_(Mmat), /*b_(b),*/ fe_(fe), f_(f), w_(w), gradw_(gradw), divw_(divw), surfactant_(surfactant),
        M_(0),
        t(t_)
{
    for(int i=0; i<3; i++)
    {
        phi[i][i]=1.;
        //phiQuad[i].assign(phi[i]);
    }
}

void MassAccumulatorP1CL::begin_accumulation ()
{
   /* if (b_ != 0)
        b_->Clear( t);
   */
    if (Mmat_)
        M_ = new MatrixBuilderCL( Mmat_, fe_.GetNumUnknowns(), fe_.GetNumUnknowns());

}

void MassAccumulatorP1CL::finalize_accumulation ()
{
    if (M_ != 0){
        M_->Build(); // danach ist Mmat_ gesetzt in CRS-Format
        delete M_;
    }
}

void MassAccumulatorP1CL::visit ( MySpaceTimeGridCL& stg, const MyTriangleCL& triang)
{
    local_setup(stg, triang);
    if (M_ != 0)
        update_global_matrix(stg , triang);
    //if (b_ != 0 && BndData_ != 0)
      //  update_coupling(triang);
}

void MassAccumulatorP1CL::local_setup ( MySpaceTimeGridCL& stg, const MyTriangleCL& triang)
{
   
   absdet= GetAbsDet(triang);
   
   //für (1+(V_h)^2)^(-1/2)   
   double vh=0.0; //wird in GetN_gamma richtig gesetzt
   Point2DCL vhngamma; // ist nur Dummy, wird hier nicht gebraucht
   Point2DCL ngamma; // wird in GetN_gamma richtig gesetzt
   GetN_gamma(triang, ngamma, vh, vhngamma);
   Quad5SpaceTimeInterfaceCL<> oberflaechentrafo(std::pow(1.0 + vh*vh, -1./2.));
   
   //ÄNDERUNG IST AB HIER FALSCH !!!!!!!!!!!!!!!!!!!!!!!!!!
   /*
   double dummy1;
   Point2DCL ngammahelp;
   Point2DCL dummy2;
   GetN_gamma(triang, ngammahelp, dummy1, dummy2); // Ngamma wird gesetzt
   Quad5SpaceTimeInterfaceCL<Point2DCL> ngamma(ngammahelp);
   Quad5SpaceTimeInterfaceCL<Point2DCL> w(triang, w_);
   Quad5SpaceTimeInterfaceCL<> oberflaechentrafo( std::pow(1.0 + dot(w,ngamma)*dot(w,ngamma), -1./2.) );
   */
   //ÄNDERUNG ENDE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

   if(!surfactant_) // löse nicht surfactant-Gleichung, sondern Gleichung mit bel. Koeffinzient alpha(s,t)
   {
	   
	   for(int i=0; i<3; i++){
		   for(int j=0; j<=i; j++){ // unterer Dreiecksteil der Kopplungsmatrix (Massematrix ist symmetrisch)

			   
			   
				    // für Chi_i* Chi_j
   					MyLocalP2CL<> phi_i(phi[i]);
					MyLocalP2CL<> phi_j(phi[j]);
					MyLocalP2CL<> phi_ij(phi_i*phi_j);
					Quad5SpaceTimeInterfaceCL<> quadraturij(phi_ij); 

					if(f_){ // massematrix mit Koeffinzient alpha(=f_)
				   
												
							Quad5SpaceTimeInterfaceCL<> quadraturf(triang, f_);
							Quad5SpaceTimeInterfaceCL<> quadraturij_f_trafo(oberflaechentrafo*quadraturf*quadraturij);
							coup[i][j]=quadraturij_f_trafo.quad(absdet); 

					}else{// massematrix ohne Koeffizient
							Quad5SpaceTimeInterfaceCL<> quadraturij_trafo(oberflaechentrafo*quadraturij);
							coup[i][j]=quadraturij_trafo.quad(absdet); 
					}

			   
			   
			}
		}	
   
   }else{ // Löse surfactant-Gleichung: div_\Gamma(w) hinzugefügt
	   for(int i=0; i<3; i++){
		   for(int j=0; j<=i; j++){ // unterer Dreiecksteil (Massematrix ist symmetrisch)

			   
			   
				    // für Chi_i* Chi_j
   					MyLocalP2CL<> phi_i(phi[i]);
					MyLocalP2CL<> phi_j(phi[j]);
					MyLocalP2CL<> phi_ij(phi_i*phi_j);
					Quad5SpaceTimeInterfaceCL<> quadraturij(phi_ij); 

					// für div_\Gamma(w)
					Quad5SpaceTimeInterfaceCL<> divw(triang, divw_);
					Quad5SpaceTimeInterfaceCL<Point2DCL> n_gamma(ngamma);
					Quad5SpaceTimeInterfaceCL< SMatrixCL<2,2> > gradw(triang, gradw_);
					Quad5SpaceTimeInterfaceCL<> div_gamma_w(divw-dot(n_gamma, gradw, n_gamma));

					Quad5SpaceTimeInterfaceCL<> ij_div_gamma_w_trafo(quadraturij*div_gamma_w*oberflaechentrafo);
					coup[i][j]= ij_div_gamma_w_trafo.quad(absdet);


			   

		}
	   }
	  }

   // spiegel unteren Dreiecksteil auf oberen Dreiecksteil
   for(int i=0; i<3; i++){
	   for(int j=0; j<i; j++){
		   coup[j][i]= coup[i][j];
	   }
   }

   //Update von Mapping lokale Nummerierung <-> globale Nummerierung 
   for(int k=0; k<3; k++){
		 Numb[k]= triang.getVertex(k)->Unknowns(idx);
   }
}

void MassAccumulatorP1CL::update_global_matrix(MySpaceTimeGridCL& stg , const MyTriangleCL& triang)
{
    for(int i=0; i<3; ++i){      
		
		if(IsOnBoundary(stg, *(triang.getVertex(i)))){ // Knoten i ist auf Rand -> i-te Zeile: 0,0,...0,1,0,..,0
			for(int j=0; j<3; j++){
				if(i==j) {(*M_)(Numb[i], Numb[j])= 1.0;}
				else{	(*M_)(Numb[i], Numb[j])= 0.0;}
			}

		}else{ // Knoten i ist nicht auf Rand
			for(int j=0; j<3; j++){
				(*M_)(Numb[i], Numb[j])+= coup[i][j];
			}
		}
		
		
    }
        
}


//**********************************************
//
//	Klasse: RHSAccumulatorP1CL
//
//**********************************************

RHSAccumulatorP1CL::RHSAccumulatorP1CL (const MySpaceTimeGridCL & stg, scalar_fct_ptr BndData, boundarymap BndMap, scalar_fct_ptr f, VecDescCL* b,
            SpaceTimeInterfaceP1CL& fe,  vec2D_fct_ptr w, 
			MassAccumulatorP1CL* mass, StiffAccumulatorP1CL* stiff, /*double alpha,*/ DerivAccumulatorP1CL* deriv, const double t_):
			
			stg_(stg), BndData_(BndData), f_(f), b_(b), fe_(fe), w_(w), mass_(mass), stiff_(stiff), /*alpha_(alpha),*/ deriv_(deriv), t(t_)
{
    for(int i=0; i<3; i++)
    {
        phi[i][i]=1.;
    }

	// setze anzahlMatrizen
	anzahlMatrizen= 0.0;
	if(mass_!=0) anzahlMatrizen+= 1.0;
	if(stiff_!=0) anzahlMatrizen+= 1.0;
	if(deriv_!=0) anzahlMatrizen+= 1.0;

	// kopiere boundarymap
	BndMap_.clear();
	boundarymap::iterator it= BndMap.begin();
	for(it; it!=BndMap.end(); it++){
		BndMap_[it->first]=it->second;
		
	}
	/*
	std::cout << "BndMap_ groesse bei konstruktor= " << BndMap_.size() << "\n";
	for(boundarymap::iterator iter= BndMap_.begin(); iter!=BndMap_.end(); iter++){
		std::cout<< "BndMap_first= " << iter->first << " BndMap_second= " << iter->second << "\n";
	}
	*/
}

void RHSAccumulatorP1CL::begin_accumulation()
{
	if(b_!=0){ b_ -> Data.resize(0); b_->Data.resize(fe_.GetNumUnknowns()); }

}

void RHSAccumulatorP1CL::visit(MySpaceTimeGridCL & stg, const MyTriangleCL& triang)
{
	local_setup(stg, triang);
    if (b_ != 0) update_global_rhs(stg, triang);
			
}

void RHSAccumulatorP1CL::update_global_rhs(MySpaceTimeGridCL& stg, const MyTriangleCL& triang)
{
	for(int i=0; i<3; ++i){        
			if( !IsOnBoundary(stg, *(triang.getVertex(i))) ){
				(b_ ->Data)[Numb[i]]+=coup[i];
			}
			else{
				if(BndData_ != 0){ // Randdaten mit Funktionenpointer übergeben
					(b_ ->Data)[Numb[i]]= anzahlMatrizen* BndData_(triang.getVertex(i)->GetCoordinates(), 0.0);
				}else{ // Randdaten mit boundarymap übergeben
					boundarymap::iterator help= FindInBoundaryMap(BndMap_, triang.getVertex(i)->GetCoordinates());
					if( help != BndMap_.end()){ // Vertex ist in Map
						(b_ ->Data)[Numb[i]]= anzahlMatrizen* BndMap_[help->first];
					}
					
				}
			}
    }

}


void RHSAccumulatorP1CL::local_setup(MySpaceTimeGridCL& stg, const MyTriangleCL& triang)
{
   
   absdet= GetAbsDet(triang);
   
   //für (1+(V_h)^2)^(-1/2)
   double vh=0.0; //wird in GetN_gamma richtig gesetzt
   Point2DCL vhngamma; // wird in GetN_gamma richtig gesetzt, ist hier aber nur Dummy
   Point2DCL ngamma; // wird in GetN_gamma richtig gesetzt, ist hier aber nur Dummy
   GetN_gamma(triang, ngamma, vh, vhngamma);
   Quad5SpaceTimeInterfaceCL<> oberflaechentrafo(std::pow(1.0 + vh*vh, -1./2.));
   
   //ÄNDERUNG AB HIER IST FALSCH!!!
   /*
   double dummy1;
   Point2DCL ngammahelp;
   Point2DCL dummy2;
   GetN_gamma(triang, ngammahelp, dummy1, dummy2); // Ngamma wird gesetzt
   Quad5SpaceTimeInterfaceCL<Point2DCL> ngamma(ngammahelp);
   Quad5SpaceTimeInterfaceCL<Point2DCL> w(triang, w_);
   Quad5SpaceTimeInterfaceCL<> oberflaechentrafo( std::pow(1.0 + dot(w,ngamma)*dot(w,ngamma), -1./2.) );
   */
   //ÄNDERUNG ENDE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   for(int i=0; i<3; i++){
	   
		   Quad5SpaceTimeInterfaceCL<> quadraturi(phi[i]);
		   Quad5SpaceTimeInterfaceCL<> quadraturf( triang, f_); 
		   Quad5SpaceTimeInterfaceCL<> f_i_trafo(quadraturi*quadraturf*oberflaechentrafo);

		   coup[i]= f_i_trafo.quad(absdet);

	   
	   

	   Numb[i]= triang.getVertex(i)->Unknowns(idx);
   
   
   }

   

}


//*************************************************************
//
//	Klasse: StiffAccumulatorP1CL
//
//*************************************************************

const Point2DCL StiffAccumulatorP1CL::GradRef[3] = {MakePoint2D(-1.0,-1.0), MakePoint2D(1.0, 0.0), MakePoint2D(0.0, 1.0)};

StiffAccumulatorP1CL::StiffAccumulatorP1CL(const MySpaceTimeGridCL & stg,  MatrixCL* Smat, 
        SpaceTimeInterfaceP1CL& fe, scalar_fct_ptr f, vec2D_fct_ptr w, bool surfactant, const double t_):
        stg_(stg), Smat_(Smat),  fe_(fe), f_(f), w_(w),  surfactant_(surfactant),
        S_(0),
        t(t_)
{
    for(int i=0; i<3; i++)
    {
        phi[i][i]=1.;
       
    }
}

void StiffAccumulatorP1CL::begin_accumulation ()
{
   
    if (Smat_)
        S_ = new MatrixBuilderCL( Smat_, fe_.GetNumUnknowns(), fe_.GetNumUnknowns());

}

void StiffAccumulatorP1CL::finalize_accumulation ()
{
    if (S_ != 0){
        S_->Build(); // danach ist Mmat_ gesetzt in CRS-Format
        delete S_;
    }
}

void StiffAccumulatorP1CL::visit ( MySpaceTimeGridCL& stg, const MyTriangleCL& triang)
{
    local_setup(stg, triang);
    if (S_ != 0)
        update_global_matrix(stg, triang);
    
}

void StiffAccumulatorP1CL::update_global_matrix(MySpaceTimeGridCL& stg, const MyTriangleCL& triang)
{	
	 for(int i=0; i<3; ++i){      
		
		if(IsOnBoundary(stg, *(triang.getVertex(i)))){ // Knoten i ist auf Rand -> i-te Zeile: 0,0,...0,1,0,..,0
			for(int j=0; j<3; j++){
				if(i==j) {(*S_)(Numb[i], Numb[j])= 1.0;}
				else{	(*S_)(Numb[i], Numb[j])= 0.0;}
			}

		}else{ // Knoten i ist nicht auf Rand
			for(int j=0; j<3; j++){
				(*S_)(Numb[i], Numb[j])+= coup[i][j];
			}
		}
		
		
    }
        
        
}


void StiffAccumulatorP1CL::local_setup(MySpaceTimeGridCL& stg, const MyTriangleCL& triang)
{
 
   absdet= GetAbsDet(triang);
   
   //für (1+(V_h)^2)^(-1/2)   
   double vh=0.0; //wird in GetN_gamma richtig gesetzt
   Point2DCL vhngamma; // wird in GetN_gamma richtig gesetzt
   Point2DCL ngamma; // wird in GetN_gamma richtig gesetzt
   GetN_gamma(triang, ngamma, vh, vhngamma);
   Quad5SpaceTimeInterfaceCL<> oberflaechentrafo(std::pow(1.0 + vh*vh, -1./2.));
   
   bool detklein=false;
   double det;

   //ÄNDERUNG IST AB HIER FALSCH !!!!!!!!!!!!!!!!!!!!!!!!!!
   /*
   double dummy1;
   Point2DCL ngammahelp;
   Point2DCL dummy2;
   GetN_gamma(triang, ngammahelp, dummy1, dummy2); // Ngamma wird gesetzt
   Quad5SpaceTimeInterfaceCL<Point2DCL> ngamma(ngammahelp);
   Quad5SpaceTimeInterfaceCL<Point2DCL> w(triang, w_);
   Quad5SpaceTimeInterfaceCL<> oberflaechentrafo( std::pow(1.0 + dot(w,ngamma)*dot(w,ngamma), -1./2.) );
   */
   //ÄNDERUNG ENDE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if(surfactant_){

	   SMatrixCL<3,2> M= GetGradTrafoMatrix(triang);
	   // gucke ob det klein ist
	   /*
	   Point3DCL P= triang.getVertex(0)->GetCoordinates();
	   Point3DCL Q= triang.getVertex(1)->GetCoordinates();
	   Point3DCL R= triang.getVertex(2)->GetCoordinates();

	   SMatrixCL<3,2> A; // A= (Q-P | R-P)
	   A(0,0)= Q[0]-P[0];
	   A(1,0)= Q[1]-P[1];
	   A(2,0)= Q[2]-P[2];
	   A(0,1)= R[0]-P[0];
	   A(1,1)= R[1]-P[1];
	   A(2,1)= R[2]-P[2];

	

	   SMatrixCL<2,2> A_tA= GramMatrix(A);

	   SMatrixCL<2,2> invA_tA;
	
	   det= (A_tA(0,0)*A_tA(1,1) - A_tA(0,1)*A_tA(1,0));
	   
	   
	   if(std::abs(det) < 1.0e-18) detklein= true;
	   if(detklein) {
		   //std::cout << "det= " << det << "\n";
		   //std::cout<< "dreieck= " << std::setprecision(20) << "(" << P << ") (" << Q << ") (" << R << ")\n";
	   }
	   // ende: gucke, ob det klein ist
	   */

	   //if(detklein) std::cout << "coup für dieses dreieck\n";
	   for(int i=0; i<3; i++){
		   for(int j=0; j<=i; j++){ // baue unteres Dreieck der coup-Matrix auf (Steifigkeit ist symmetrisch)

			   		  

			   

				   // für grad_\Gamma_stern(Chi_i), grad_\Gamma_stern(Chi_j)
				   
				   Point3DCL grad_gamma_stern_i(M*GradRef[i]);
				   //std::cout << "grad_gamma_stern_i= " << grad_gamma_stern_i << "\n";
				   Point3DCL grad_gamma_stern_j(M*GradRef[j]);

				   // für grad_\Gamma(Chi_i), grad_\Gamma(Chi_j)
				  
				   //RICHTIG!!!!!!!!!!!!!!!
				   Quad5SpaceTimeInterfaceCL<Point2DCL> vektor(vhngamma);
				   //AB HIER IST ÄNDERUNG FALSCH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				   //Quad5SpaceTimeInterfaceCL<Point2DCL> vektor(dot(w,ngamma)*ngamma);
				   //ENDE ÄNDERUNG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				   Quad5SpaceTimeInterfaceCL<Point2DCL> grad_gamma_i ((MakePoint2D(grad_gamma_stern_i[0], grad_gamma_stern_i[1])
															-grad_gamma_stern_i[2]* vektor));
				   

				   Quad5SpaceTimeInterfaceCL<Point2DCL> grad_gamma_j ((MakePoint2D(grad_gamma_stern_j[0], grad_gamma_stern_j[1])
															-grad_gamma_stern_j[2]* vektor));
				  


				   if(!f_){

						Quad5SpaceTimeInterfaceCL<double> res(oberflaechentrafo* dot(grad_gamma_i, grad_gamma_j));
						// WICHTIGE Änderung
						//if(0.5*absdet < 1.0e-10) coup[i][j]= 0.0;
						/*else*/	coup[i][j]= res.quad(absdet);
						
				   
				   }else{
						
					   Quad5SpaceTimeInterfaceCL<double> koeff(triang,f_);
					   Quad5SpaceTimeInterfaceCL<double> res(oberflaechentrafo* koeff* dot(grad_gamma_i, grad_gamma_j));
					   // WICHTIGE Änderung
					   //if(0.5*absdet < 1.0e-10) coup[i][j]= 0.0;
					   /*else*/ 	coup[i][j]= res.quad(absdet);
				   
				   }
			   
			   
			   
		   
		   
		   }
	   
	   
	   } 
	   
   
   
   }else{ // ToDo: Was passiert, falls surfactant==false?
   
   
   
   }
	// spiegel untere Dreieckshälfte auf obere Dreieckshälfte
	for(int i=0; i<3; i++){
		   for(int j=0; j<i; j++){
			   
			   coup[j][i]= coup[i][j];
		   }
	}

   // Update Mapping lokale Nummerierung <-> globale Nummerierung
   for(int k=0; k<3; k++){
		 Numb[k]= triang.getVertex(k)->Unknowns(idx);
   }
   /*
   if(detklein){
	  
	   for(int i=0; i<3; i++){
			for(int j=0; j<3; j++){
				std::cout << "coup[" << Numb[i] << "][" << Numb[j]<<"]= " << coup[i][j] << "\n";
			}
	   }
   
   }
   */
   

}





//*************************************************************
//
//	Klasse: DerivAccumulatorP1CL
//
//*************************************************************

const Point2DCL DerivAccumulatorP1CL::GradRef[3] = {MakePoint2D(-1.0,-1.0), MakePoint2D(1.0, 0.0), MakePoint2D(0.0, 1.0)};

DerivAccumulatorP1CL::DerivAccumulatorP1CL(const MySpaceTimeGridCL & stg,  MatrixCL* Umat, 
        SpaceTimeInterfaceP1CL& fe, scalar_fct_ptr f, vec2D_fct_ptr w, bool surfactant, const double t_):
        stg_(stg), Umat_(Umat),  fe_(fe), f_(f), w_(w),  surfactant_(surfactant),
        U_(0),
        t(t_)
{
    for(int i=0; i<3; i++)
    {
        phi[i][i]=1.;
       
    }
}

void DerivAccumulatorP1CL::begin_accumulation ()
{
   
    if (Umat_)
        U_ = new MatrixBuilderCL( Umat_, fe_.GetNumUnknowns(), fe_.GetNumUnknowns());

}

void DerivAccumulatorP1CL::finalize_accumulation ()
{
    if (U_ != 0){
        U_->Build(); // danach ist Umat_ gesetzt in CRS-Format
        delete U_;
    }
}

void DerivAccumulatorP1CL::visit ( MySpaceTimeGridCL& stg, const MyTriangleCL& triang)
{
    local_setup(stg, triang);
    if (U_ != 0)
        update_global_matrix(stg, triang);
    
}

void DerivAccumulatorP1CL::update_global_matrix(MySpaceTimeGridCL& stg, const MyTriangleCL& triang)
{
	
    for(int i=0; i<3; ++i){      
		
		if(IsOnBoundary(stg, *(triang.getVertex(i)))){ // Knoten i ist auf Rand -> i-te Zeile: 0,0,...0,1,0,..,0
			for(int j=0; j<3; j++){
				if(i==j) {(*U_)(Numb[i], Numb[j])= 1.0;}
				else{	(*U_)(Numb[i], Numb[j])= 0.0;}
			}

		}else{ // Knoten i ist nicht auf Rand
			for(int j=0; j<3; j++){
				(*U_)(Numb[i], Numb[j])+= coup[i][j];
			}
		}
		
		/*
		for(int j=0; j<3; j++){// nur für Dot-Matrix-Test so, will 1 weghaben
		(*U_)(Numb[i], Numb[j])+= coup[i][j]; // nur für Dot-Matrix-Test so, will 1 weghaben
		}
		*/
    }
        
}


void DerivAccumulatorP1CL::local_setup(MySpaceTimeGridCL& stg, const MyTriangleCL& triang)
{
 
   absdet= GetAbsDet(triang);

   //für (1+(V_h)^2)^(-1/2)   
   
   double vh=0.0; //wird in GetN_gamma richtig gesetzt
   Point2DCL vhngamma; // ist nur Dummy, wird hier nicht gebraucht
   Point2DCL ngamma; // wird in GetN_gamma richtig gesetzt
   GetN_gamma(triang, ngamma, vh, vhngamma);
   Quad5SpaceTimeInterfaceCL<> oberflaechentrafo(std::pow(1.0 + vh*vh, -1./2.));
   Quad5SpaceTimeInterfaceCL<Point2DCL> w(triang, w_); // brauche ich für dot(chi_j)
   // Teste dot chi_j mit w_h statt w
   Point2DCL wh;
   GetW_h(triang, wh);
   Quad5SpaceTimeInterfaceCL<Point2DCL> w_h(wh);
   // Teste mit wneu= (w\cdot n_Gamma)n_gamma -> nur Normalenanteile
   Quad5SpaceTimeInterfaceCL<Point2DCL> ngamma_(ngamma);
   Quad5SpaceTimeInterfaceCL<Point2DCL> wneu(dot(w,ngamma_)*ngamma_);
   
   
   //ÄNDERUNG IST AB HIER FALSCH !!!!!!!!!!!!!!!!!!!!!!!!!!
   /*
   double dummy1;
   Point2DCL ngammahelp;
   Point2DCL dummy2;
   GetN_gamma(triang, ngammahelp, dummy1, dummy2); // Ngamma wird gesetzt
   Quad5SpaceTimeInterfaceCL<Point2DCL> ngamma(ngammahelp);
   Quad5SpaceTimeInterfaceCL<Point2DCL> w(triang, w_);
   Quad5SpaceTimeInterfaceCL<> oberflaechentrafo( std::pow(1.0 + dot(w,ngamma)*dot(w,ngamma), -1./2.) );
   */
   //ÄNDERUNG ENDE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   if(surfactant_){

	   SMatrixCL<3,2> M= GetGradTrafoMatrix(triang);

	   for(int i=0; i<3; i++){
		   for(int j=0; j<3; j++){

			   
				   // für Chi_i
				   Quad5SpaceTimeInterfaceCL<double> chi_i(phi[i]);

				   // für grad_\Gamma_stern(Chi_j)
				   
				   Point3DCL grad_gamma_stern_j(M*GradRef[j]);

				   // für dot(chi_j)
				   
				   Quad5SpaceTimeInterfaceCL<Point2DCL> vec (MakePoint2D(grad_gamma_stern_j[0], grad_gamma_stern_j[1]));
				   Quad5SpaceTimeInterfaceCL<double> u_2(grad_gamma_stern_j[2]);
				   // Teste dot chi_j mit w_h statt w
				   //Quad5SpaceTimeInterfaceCL<double> dot_chi_j( u_2 + dot(vec, w_h));
				   // hier dot chi_j mit w statt w_h
				   Quad5SpaceTimeInterfaceCL<double> dot_chi_j( u_2 +  dot(vec, w));
				   // hier dot chi_j mit wneu (=projektion auf normalenrichtung)
				   //Quad5SpaceTimeInterfaceCL<double> dot_chi_j (u_2 + dot(vec, wneu));


				   if(!f_){

						Quad5SpaceTimeInterfaceCL<double> res(oberflaechentrafo* dot_chi_j* chi_i);
						//WICHTIGE ÄNDERUNG
						//if(0.5*absdet < 1.0e-10) coup[i][j]= 0.0;
						/*else*/	coup[i][j]= res.quad(absdet);
				   
				   }else{
						
					   Quad5SpaceTimeInterfaceCL<double> koeff(triang,f_);
					   Quad5SpaceTimeInterfaceCL<double> res(oberflaechentrafo* koeff* dot_chi_j* chi_i);
					   //WICHTIGE ÄNDERUNG
					   //if(0.5*absdet < 1.0e-10) coup[i][j]= 0.0;
					   /*else*/ 	coup[i][j]= res.quad(absdet);
				   
				   }
			   
			   
			   
		   
		   
		   }
	   
	   
	   }
   
   
   
   }else{ // ToDo: Was passiert, falls surfactant==false?
   
   
   
   }

   // Update Mapping lokale Nummerierung <-> globale Nummerierung
   for(int k=0; k<3; k++){
		 Numb[k]= triang.getVertex(k)->Unknowns(idx);
   }

}


}// end of namespace DROPS