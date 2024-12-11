// Definitionen aus mygeometry.h




#include "SarahTest/mygeometry.h"
#include "num/interfacePatch.h"
#include "num/bndData.h"
	// muss ich das alles noch includen?
//#include "geom/multigrid.h"
//#include "misc/problem.h"
//#include "misc/container.h"
#include <map>
#include <list>
#include <utility>
#include <cmath>

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

namespace DROPS{



// Funktion, die in Viereck größten Innenwinkel ermittelt und Index zu Knoten mit größtem Innenwinkel zurückgibt
int MaxAngle (InterfacePatchCL& patch)
{
	Point3DCL P= patch.GetPoint(0);
	Point3DCL Q= patch.GetPoint(1);
	Point3DCL R= patch.GetPoint(2);
	Point3DCL S= patch.GetPoint(3);

	// bestimme Kantenlängen des Vierecks PQR + QRS (Diagonale QR)
	double pq= norm(P-Q);
	double pr= norm(P-R);
	double qr= norm(Q-R); // Diagonale
	double qs= norm(Q-S);
	double rs= norm(R-S);
	
	// Bestimme 4 Innenwinkel
	double innenwinkel[4]=
	{
		std::acos((pq*pq+pr*pr-qr*qr)/(2.*pq*pr)), // phiP (innenwinkel bei P)
		std::acos((pq*pq+qr*qr-pr*pr)/(2.*pq*qr)) + std::acos((qs*qs+qr*qr-rs*rs)/(2.*qs*qr)), // phiQ
		std::acos((pr*pr+qr*qr-pq*pq)/(2.*pr*qr)) + std::acos((rs*rs+qr*qr-qs*qs)/(2.*rs*qr)), // phiR
		std::acos((qs*qs+rs*rs-qr*qr)/(2.*qs*rs)) // phiS
	
	};

	// bestimme maximalen innenwinkel
	double max= 0.;
	int maxindex= 0;
	for(int i=0; i<4; i++){
		if(innenwinkel[i] >= max){ max= innenwinkel[i]; maxindex=i; }
	}
	
	return maxindex;
	
}

// Funktion, mit der ich meine Geometrie mit Hilfe von einem MuliGrid und einer levelset Funktion aufbaue


void MySpaceTimeGridCL::BuildMyGeometry ( const MultiGridCL& mg, const VecDescCL& levelset  ) // leveset wird nachher mein lset.Phi sein
{
	
	MultiGridCL::const_TriangTetraIteratorCL	it=mg.GetTriangTetraBegin();
	std::map <  Point3DCL,  MyVertexCL* , compare >  helper;  //Hilfsmap, über die ich nachher iteriere, um keine knoten mehrfach zu haben

	std::list<double> lastcomponent; // Hilfsliste, in der ich letzte Komponenten der Knoten sammel um t_n-1 und t_n zu finden
	double eps = 1e-10;

	
	//folgende Zeilen sind aus surfactant.cpp main
	//const BndCondT bcls[6]= { DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC };
    //const LsetBndDataCL::bnd_val_fun bfunls[6]= { 0,0,0,0,0,0};
    //LsetBndDataCL lsbnd( 6, bcls, bfunls);
	const BndDataCL<>& lsetbnd = BndDataCL<>(0); //wird bei Init einfach übergeben (keine Leveset-Randdaden)
												 // = NoBndDataCL<> oder = NoBndCondCL<> klappt nicht
	
	

	for(it; it!=mg.GetTriangTetraEnd(); ++it) {

		InterfacePatchCL patch;
		patch.Init(*it, levelset, lsetbnd);   // übergibt levelset-werte und initialisiert patch

		if(patch.Intersects()) { // falls es Intersection mit Tetra *it gibt...

			for (Uint i=0; i<8; i++) { // ...gehe 8 reguläre Kinder von *it durch und berechne Intersection mit Kindern
				
				if( patch.ComputeVerticesOfCut(i,true) ) { // .ComputeVerticesOfCut == true, falls es 2D-patch gibt
					
					// compute_PQRS=true
					
					
					for(Uint j=0; j<patch.GetNumPoints(); j++) { // es werden die Knoten P,Q,R,(S) betrachtet

						// er findet viele knoten
						//std::cout << "gefundener knoten: " << patch.GetPoint(j) << "\n";
						


						
	
						std::map<Point3DCL, MyVertexCL*, compare>::iterator maperator= helper.find(patch.GetPoint(j));

						if(maperator==helper.end()) { // Schlüssel Point3DCL patch.GetPoint(j) ist noch nicht in map

							
							
							vertices_.push_back( MyVertexCL(patch.GetPoint(j)) );// Knoten wird in VertexListe hinzugefügt

							lastcomponent.push_back((&vertices_.back())->GetCoordinates()[2]); // letzte Koordinate des obigen Knoten wird gespeichert
							
							//std::cout << "mapsize vorher: " << helper.size() << "\n";
							
							
							
							helper[patch.GetPoint(j)] = &vertices_.back();
							
								//  Zeiger auf Knoten (zeigt an die Stelle, wo Knoten in VertexList liegt) 
							//std::cout << "mapsize nachher: " << helper.size() << "\n";												//  (braucht man um nachher Dreiecke mit Zeigern aufzubauen,
																			//  die auf Knoten zeigen, die in VertexListe drin sind)
																			//  .back() funktioniert, da vertices_ nicht leer
																			// erzeugt erst Default-wert für Zeiger (NULL?) und weist dann &vertices_.back() zu
						}	
						
						
						
					
					} // end for(uint j=0..patch.GetNumPoints())
					
					if(!patch.IsQuadrilateral()) { // Schnitt des Kindes ist Dreieck

						const MyTriangleCL T1(helper[patch.GetPoint(0)], helper[patch.GetPoint(1)], helper[patch.GetPoint(2)]); // Dreieck T1=(P,Q,R) wird angelegt und zu TriangList hinzugefügt
						triangles_.push_back(T1); 
					}

					else {  // Schnitt des Kindes ist Viereck
							// teile Viereck in 2 Dreiecke, indem Knoten mit größtem Innenwinkel mit 
							// gegenüberliegenden Knoten verbunden wird
							// es gilt: Viereck = PQR + QRS mit Diagonale QR
					
					  int winkelindex= MaxAngle(patch); // index des Punkt mit größtem Innenwinkel
						if(winkelindex==1 || winkelindex==2) { // bei Q oder R liegt größter Innenwinkel
							//std::cout << "erster fall\n";
							const MyTriangleCL T1(helper[patch.GetPoint(0)], helper[patch.GetPoint(1)], helper[patch.GetPoint(2)]); // Dreieck T1=(P,Q,R) wird angelegt und zu TriangList hinzugefügt
							triangles_.push_back(T1);
							const MyTriangleCL T2(helper[patch.GetPoint(1)], helper[patch.GetPoint(2)], helper[patch.GetPoint(3)]);// Dreieck T2=(Q,R,S) wird angelegt und zu TriangList hinzugefügt
							triangles_.push_back(T2);
						}else{ // bei P oder S liegt größter Innenwinkel
							//std::cout << "zweiter fall\n";
							const MyTriangleCL T1(helper[patch.GetPoint(1)], helper[patch.GetPoint(0)], helper[patch.GetPoint(3)]); // Dreieck T1=(Q,P,S) wird angelegt und zu TriangList hinzugefügt
							triangles_.push_back(T1);
							const MyTriangleCL T2(helper[patch.GetPoint(0)], helper[patch.GetPoint(3)], helper[patch.GetPoint(2)]);// Dreieck T2=(P,S,R) wird angelegt und zu TriangList hinzugefügt
							triangles_.push_back(T2);						
						}

					}
					
				
				}// ende von if(patch.IntersectsChild(i))

			
			} // ende von for(Uint i=0; i<8; i++) 


		} // ende von if(patch.Intersects())


	
	}// jetzt sind alle tetras durchgelaufen worden

	


	

	MyVertexIterator iter;
	MyVertexIterator jter;
	
	for(iter=GetVerticesBegin(); iter!=GetVerticesEnd(); ++iter)
	{
		
		for(jter=GetVerticesBegin(); jter!=GetVerticesEnd(); ++jter)
		{
			if(iter!=jter){
				if(iter->GetCoordinates() == jter->GetCoordinates()) {std::cout << "doppelte punkte in vertices_\n";}
			}
		}
	}
	
	// told_ und tnew_ werden gesetzt
	if(lastcomponent.size()!=0){
		lastcomponent.sort();
		told_= lastcomponent.front();
		tnew_= lastcomponent.back();

		
	}


	

}// end of BuildMyGeometry



// baue das örtliche Gitter MySpaceGridCL aus einer Instanz von MySpaceTimeGridCL
void MySpaceGridCL::BuildSpaceGrid(MySpaceTimeGridCL& stg, bool timenew)
{
	double eps= 0.1e-14;
	if(timenew==true){
		t_=stg.Get_t_new();
	}else{
		t_=stg.Get_t_old();
	}
	
	std::map <  Point3DCL,  MyVertexCL* , compare >  helper2;
	for(MySpaceTimeGridCL::MyTriangIterator iter= stg.GetTriangBegin(); iter!=stg.GetTriangEnd(); iter++){ // gehe alle Dreiecke aus SpaceTimeGrid durch
		for(int i=0; i<3; i++){ // (i=0, j=1); (i=1, j=2); (i=2, j=0)
			int j= (i+1)%3;
			// Knoten i und j liegen auf oberem Rand
			if( std::abs((iter->getVertex(i))->GetCoordinates()[2] - t_) < eps && 
						std::abs((iter->getVertex(j))->GetCoordinates()[2] - t_) < eps) {

				std::map<Point3DCL, MyVertexCL*, compare>::iterator maperator0= helper2.find((iter->getVertex(i))->GetCoordinates());
				std::map<Point3DCL, MyVertexCL*, compare>::iterator maperator1= helper2.find((iter->getVertex(j))->GetCoordinates());

				if(maperator0==helper2.end()){ //Knoten i noch nicht in Hilfsmap
					vertices_.push_back(*(iter->getVertex(i)));
					helper2[(iter->getVertex(i))->GetCoordinates()] = &vertices_.back();
				}
				if(maperator1==helper2.end()){ //Knoten j noch nicht in Hilfsmap
					vertices_.push_back(*(iter->getVertex(j)));
					helper2[(iter->getVertex(j))->GetCoordinates()] = &vertices_.back();
				}

				const MyLineCL line(helper2[(iter->getVertex(i))->GetCoordinates()], helper2[(iter->getVertex(j))->GetCoordinates()]);
				lines_.push_back(line); 
			}
		}
		
	}

	
	for(MySpaceGridCL::MyVertexIterator iiter=GetVerticesBegin(); iiter!=GetVerticesEnd(); ++iiter)
	{
		
		for(MySpaceGridCL::MyVertexIterator jjter=GetVerticesBegin(); jjter!=GetVerticesEnd(); ++jjter)
		{
			if(iiter!=jjter){
				if(iiter->GetCoordinates() == jjter->GetCoordinates()) {std::cout << "doppelte punkte in sg- vertices_\n";}
			}
		}
	}



}// end of BuildSpaceGrid


int GetCounter (const Point3DCL& p, std::map<int, Point3DCL>& map)
{
	int size = map.size();
	int res;
	for(int i=0; i<size; i++)
	{
		if(map[i]==p) res=i; // muss ich den vergleich hier auch mit eps machen? p ist ja schon exakt in map drin
	}

  return res;
}

/*
MyTriangleCL Quad2Triang ( int& pp, int& qq, int& rr, int& ss, std::map<int, MyVertexCL*>& pointer )
{
	Point3DCL P = pointer[pp]->GetCoordinates();
	Point3DCL Q = pointer[qq]->GetCoordinates();
	Point3DCL R = pointer[rr]->GetCoordinates();
	Point3DCL S = pointer[ss]->GetCoordinates();

	double length1 = norm(P-Q) + norm(R-S);
	double length2 = norm(P-R) + norm(Q-S);
	double length3 = norm(P-S) + norm(Q-R);

	double max = -1.0;

	if(length1 >= length2) { max = length1; }
	else { max = length2; }
	// max ist jetzt max(length1, length2)

	if(length3 >= max) { max = length3; }
	// max ist jetzt max(length1, length2, length3)

	//if(max==length3) std::cout << "length3\n";
	MyTriangleCL T = MyTriangleCL (pointer[qq], pointer[rr], pointer[ss]); // Annahme: max == length3 -> QR ist Diagonale
	if(max == length1){
		//std::cout << "length1\n";
		 T = MyTriangleCL (pointer[pp], pointer[qq], pointer[ss]); //PQ ist Diagonale
	}else{
			if (max == length2){
				//std::cout << "length2\n";
				 T = MyTriangleCL (pointer[pp], pointer[rr], pointer[ss]); //PR ist Diagonale
			}
	}

	return T; 

}
*/

/*
//**********************************************************
//					 ALTERNATIVE 2
//**********************************************************


void MySpaceTimeGridCL::BuildMyGeometry ( const MultiGridCL& mg, const VecDescCL& levelset  ) // leveset wird nachher mein lset.Phi sein
{
	
	MultiGridCL::const_TriangTetraIteratorCL	it=mg.GetTriangTetraBegin();
	std::map <int, Point3DCL> coordmap;
	std::map <int,  MyVertexCL*>  pointermap;  
	int zaehler=0; // int-zaehler für beide maps
	double eps = 1e-10;
	
	const BndDataCL<>& lsetbnd = BndDataCL<>(0); //wird bei Init einfach übergeben (keine Leveset-Randdaden) (macht man das so mit NoBC?)
												 // = NoBndDataCL<> oder = NoBndCondCL<> klappt nicht
	

	for(it; it!=mg.GetTriangTetraEnd(); ++it) {

		InterfacePatchCL patch;
		patch.Init(*it, levelset, lsetbnd);   // übergibt levelset-werte und initialisiert patch

		if(patch.Intersects()) { // falls es Intersection mit Tetra *it gibt...
			//std::cout << "Tetra wird geschnitten\n";
			for (Uint i=0; i<8; i++) { // ...gehe 8 reguläre Kinder von *it durch und berechne Intersection mit Kindern
				
				if( patch.ComputeVerticesOfCut(i,true) ) { // .ComputeVerticesOfCut() == true, falls es 2D-patch gibt
					
					 // compute_PQRS=true
					//std:: cout << "Viereck? " << patch.IsQuadrilateral() << "\n";
					
					for(Uint j=0; j<patch.GetNumPoints(); j++) { // es werden Knoten P,Q,R,(S) betrachtet
			
						//std::cout << "Gefundener Knoten: " << patch.GetPoint(j) << "\n";
						// eigene find()-Funktion						
						bool found = false;
						std::map<int, Point3DCL>::iterator maperator;
						for (maperator=coordmap.begin(); maperator!=coordmap.end(); ++maperator)
						{
							
							if(norm((maperator->second)-patch.GetPoint(j)) < eps){ found = true;  } // besser, als maperator->second == patch.GetPoint(j)
						}
						

						if(found==false) { //  Point3DCL patch.GetPoint(j) ist noch nicht in coordmap
							
							//std:: cout << "neu gefunden\n";
							//std:: cout << "vertices-size vorher: " << vertices_.size() << "\n";
							vertices_.push_back( MyVertexCL(patch.GetPoint(j)) );// Knoten wird in VertexListe hinzugefügt
							//std:: cout << "vertices-size nachher: " << vertices_.size() << "\n";
							coordmap[zaehler] = patch.GetPoint(j);
							pointermap[zaehler] = &vertices_.back();
							zaehler++;
							
						}	
						
						
						
					} // end for(uint j=0..patch.GetNumPoints())
					
						// hole int
						int p = GetCounter(patch.GetPoint(0), coordmap);
						int q = GetCounter(patch.GetPoint(1), coordmap);
						int r = GetCounter(patch.GetPoint(2), coordmap);
						// baue Dreieck
						const MyTriangleCL T1(pointermap[p], pointermap[q], pointermap[r]); // Dreieck T1=(P,Q,R) wird angelegt und zu TriangList hinzugefügt
						triangles_.push_back(T1); 
					

					if(patch.IsQuadrilateral()) { // Schnitt des Kindes ist Viereck -> Dreieck T2 hinzufügen
						
						std::cout << "in viereck\n";

					  //hole int
					  int s = GetCounter(patch.GetPoint(3), coordmap);

					  //baue Dreieck:

					  const MyTriangleCL T2(pointermap[q],pointermap[r], pointermap[s] );
					  triangles_.push_back(T2);

					}// ende von if(patch.IsQuadrilateral())
					
				
				}// ende von if(patch.IntersectsChild(i))

			
			} // ende von for(Uint i=0; i<8; i++) 


		} // ende von if(patch.Intersects())


	
	}// jetzt sind alle tetras durchgelaufen worden

	// Test, ob es doppelte Knoten in vertices_ gibt: GIBT ES NICHT!
	MyVertexIterator iter1;
	MyVertexIterator jter1;
	
	for(iter1=GetVerticesBegin(); iter1!=GetVerticesEnd(); ++iter1)
	{
		
		for(jter1=GetVerticesBegin(); jter1!=GetVerticesEnd(); ++jter1)
		{
			if(iter1!=jter1){
				if(iter1->GetCoordinates() == jter1->GetCoordinates()) {std::cout << "doppelte punkte in vertices_\n";}
			}
		}
	}
    
	
	// Test, ob es doppelte Dreiecke in triangles_ gibt: GIBT ES NICHT (bei 1-x-y-z)!
	MyTriangIterator iter;
	MyTriangIterator jter;
	
	for(iter=GetTriangBegin(); iter!=GetTriangEnd(); ++iter)
	{
		
		for(jter=GetTriangBegin(); jter!=GetTriangEnd(); ++jter)
		{
			if(iter!=jter){
				if(*iter == *jter) {std::cout << "doppelte dreiecke in triangles_\n";}
			}
		}
	}
	
	
}// end of BuildMyGeometry
*/

}//end of namespace DROPS

