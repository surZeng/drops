// Definitionen aus myfespace.h


#include "SarahTest/myfespace.h"
#include <map>

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


//Konstruiere Nummerierung der Unbekannten an Knoten im Gitter stg für FE-Raum (fortlaufende Nummerierung)
	void SpaceTimeInterfaceP1CL::NumberingVertices(   MySpaceTimeGridCL& stg  )// Aufruf nicht mit const Instanz, da stg verändert wird (es werden Unknowns nummeriert)
	{
		
		const int sysnum= Idx_;
		Uint zaehler= 0;
		for(MySpaceTimeGridCL::MyVertexIterator it=stg.GetVerticesBegin(); it!=stg.GetVerticesEnd(); ++it){
			it->Unknowns.Prepare(sysnum);          //holt Speicher für den Index des FE-Raums 
			it->Unknowns(sysnum)=zaehler;          //setzt den Index an Stelle sysnum auf zaehler
			zaehler++;

			
		}

				
		NumUnknowns_=zaehler; // setze die Anzahl der Unbekannten im FE-Raum
	}


	void SpaceTimeInterfaceP1CL::NumberingVerticesSG(MySpaceGridCL& sg, MySpaceTimeGridCL& stg)
	{
		const int sysnum= Idx_;
		std::map <  Point3DCL,  int , compare >  helper; // speichert Knoten, die auf oberem Rand liegen und Unbekannte dazu
		for(MySpaceTimeGridCL::MyVertexIterator it=stg.GetVerticesBegin(); it!=stg.GetVerticesEnd(); ++it){
			if( std::abs((it->GetCoordinates()[2])- stg.Get_t_new()) < 0.1e-13){ //falls Knoten auf oberem Rand liegt, speicher in helper
				helper[it->GetCoordinates()]= it->Unknowns(sysnum);
			}
		}

		for(MySpaceGridCL::MyVertexIterator oter=sg.GetVerticesBegin(); oter!=sg.GetVerticesEnd(); oter++){
			oter->Unknowns.Prepare(sysnum);
			oter->Unknowns(sysnum)= helper[oter->GetCoordinates()];
			//std::cout << "Knoten hat unbekannte " << oter->Unknowns(sysnum) << "\n";
		}
	}

} // end of namespace DROPS