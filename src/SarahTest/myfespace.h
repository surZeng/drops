// Mein Raum-Zeit-P1 Finite Elemente Raum auf Interface


#ifndef DROPS_MYFESPACE_H
#define DROPS_MYFESPACE_H

#include "SarahTest/mygeometry.h"

namespace DROPS  {

	
// forward declaration
class SpaceTimeInterfaceP1CL;

//************************************************************************************
// SpaceTimeInterfaceP1CL: Klasse f�r meinen Raum-Zeit-FE-Raum auf Interface
//
//                         idx_ ist eindeutiger Index dieses Raums (standardm��ig 0)
//                         NumUnknowns_ die Anzahl der Unbekannte f�r Raum
//************************************************************************************

class SpaceTimeInterfaceP1CL 
{
private: 
	Uint Idx_;  // Index des FE-Raums
	IdxT NumUnknowns_;   // Anzahl der Unbekannte im FE-Raum = Anzahl der Vertices in MySpaceTimeGridCL

public:

	//Default-Konstruktor (standardm��ig werden idx und numUnknowns auf 0 gesetzt)
	SpaceTimeInterfaceP1CL(Uint Idx = 0, IdxT NumUnknowns = 0){ Idx_ = Idx; NumUnknowns_ = NumUnknowns; };

	//Copy-Konstruktor
	SpaceTimeInterfaceP1CL( const SpaceTimeInterfaceP1CL& fe )  { Idx_ = fe.Idx_;   NumUnknowns_ = fe.NumUnknowns_; }

	//member functions

	Uint GetIdx() const { return Idx_; } // gibt Index des FE-Raums zur�ck
	IdxT GetNumUnknowns() const { return NumUnknowns_; } // gibt Anzahl der Unbekannte des FE-Raums zur�ck
	void NumberingVertices ( MySpaceTimeGridCL& stg ); // nummeriert die Knoten aus Gitter stg und setzt NumUnknowns des aufrufenden FE-Raums
	//folgende Funktion macht nur Sinn, wenn f�r stg schon NumberingVertices aufgerufen wurde
	void NumberingVerticesSG (MySpaceGridCL& sg, MySpaceTimeGridCL& stg);



};



} // end of namespace DROPS

#endif


