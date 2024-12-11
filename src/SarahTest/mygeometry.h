//Meine Geometrie: MyVertexCL (Darstellung eines Knoten), 
//				   MyTriangCL (Darstellung eines planaren Dreiecks), 
//				   MySpaceTimeGridCL (Darstellung des Raum-Zeit-Gitters)
//
//				   BuildMyGeometry (baut geometrie aus Multigitter und levelset-Daten auf)
//
//				   MyLineCL (Darstellung eines Liniensegments)
//				   MySpaceGridCL (Darstellung des örtlichen Gitters, hier nur bei t=t_new(=T) )



#ifndef DROPS_MYGEOMETRY_H
#define DROPS_MYGEOMETRY_H

#include "misc/container.h"
#include "num/unknowns.h"
#include "misc/problem.h"
//#include "num/interfacePatch.h"
#include "geom/multigrid.h"
//#include <map>
#include <list>



namespace DROPS
{

//forward declarations 
class MyVertexCL;
class MyTriangleCL;
class MyLineCL;
class MySpaceGridCL;
class MySpaceTimeGridCL;








//*****************************************************************
//	MyVertexCL: Darstellung eines Knoten in 3D (2D-Ort, 1D-Zeit)  *
//*****************************************************************

class MyVertexCL 
{
private: 
	Point3DCL Coordinates_;//3D-Koordinaten des Knoten

public:
	UnknownHandleCL Unknowns;//Unbekannte an Knoten (public, da ich da ran muss)

	//Default-Konstruktor für uninitialisierten Knoten 
	MyVertexCL(){}

	//Konstruktor
	MyVertexCL(const Point3DCL& Coord)
		: Coordinates_(Coord) {}//so wird Coord direkt in Coordinates_ kopiert, ohne dass vorher Coordinates_ angelegt werden muss

	//Copy-Konstruktor
	MyVertexCL(const MyVertexCL& Vertex){Coordinates_=Vertex.Coordinates_; Unknowns=Vertex.Unknowns;}

	//Destruktor (brauche keinen Destruktor)
	


	//member functions

	//Gibt Koordinaten des Knoten
	const Point3DCL& GetCoordinates() const { return Coordinates_; }//warum gebe ich referenz zurück-> damit keine kopie gemacht werden muss
	//Verändert Koordinaten des Knoten
	void ChangeCoordinates(Point3DCL& p) { Coordinates_=p; }// call by reference (s.o.)

};


//********************************************************
// MyTriangleCL: Darstellung eines ebenen Dreiecks       *
//********************************************************

class MyTriangleCL
{

private: 
	SArrayCL<MyVertexCL*, 3> Vertices_;

public:
	//Default-Konstruktor für uninitialisiertes Dreieck (brauche ich nicht)
	//MyTriangleCL(){}

	//Konstruktor
	MyTriangleCL(MyVertexCL* v0, MyVertexCL* v1, MyVertexCL* v2){Vertices_[0]=v0; Vertices_[1]=v1; Vertices_[2]=v2;}

	//Copy-Konstruktor
	MyTriangleCL(const MyTriangleCL& triang){Vertices_=triang.Vertices_;}

	//Destruktor (brauche ich nicht)
	

	//member functions

	//Gibt Knoten i des Dreiecks zurück
	/*const*/ MyVertexCL* getVertex(int i) const { return Vertices_[i]; }
	//Prüft, ob Knoten in Dreieck liegt
	bool HasVertex(const MyVertexCL* vert) const { return (vert==Vertices_[0]||vert==Vertices_[1]||vert==Vertices_[2]); }

	friend bool operator == (MyTriangleCL& T1, MyTriangleCL& T2)
	{
		if(T1.Vertices_[0]->GetCoordinates() == T2.Vertices_[0]->GetCoordinates() 
			&& T1.Vertices_[1]->GetCoordinates() == T2.Vertices_[1]->GetCoordinates() 
			&& T1.Vertices_[2]->GetCoordinates() == T2.Vertices_[2]->GetCoordinates()) return true;

		if(T1.Vertices_[0]->GetCoordinates() == T2.Vertices_[0]->GetCoordinates() 
			&& T1.Vertices_[1]->GetCoordinates() == T2.Vertices_[2]->GetCoordinates() 
			&& T1.Vertices_[2]->GetCoordinates() == T2.Vertices_[1]->GetCoordinates()) return true;

		if(T1.Vertices_[0]->GetCoordinates() == T2.Vertices_[1]->GetCoordinates() 
			&& T1.Vertices_[1]->GetCoordinates() == T2.Vertices_[0]->GetCoordinates() 
			&& T1.Vertices_[2]->GetCoordinates() == T2.Vertices_[2]->GetCoordinates()) return true;

		if(T1.Vertices_[0]->GetCoordinates() == T2.Vertices_[1]->GetCoordinates() 
			&& T1.Vertices_[1]->GetCoordinates() == T2.Vertices_[2]->GetCoordinates() 
			&& T1.Vertices_[2]->GetCoordinates() == T2.Vertices_[0]->GetCoordinates()) return true;

		if(T1.Vertices_[0]->GetCoordinates() == T2.Vertices_[2]->GetCoordinates() 
			&& T1.Vertices_[1]->GetCoordinates() == T2.Vertices_[0]->GetCoordinates() 
			&& T1.Vertices_[2]->GetCoordinates() == T2.Vertices_[1]->GetCoordinates()) return true;

		if(T1.Vertices_[0]->GetCoordinates() == T2.Vertices_[2]->GetCoordinates() 
			&& T1.Vertices_[1]->GetCoordinates() == T2.Vertices_[1]->GetCoordinates() 
			&& T1.Vertices_[2]->GetCoordinates() == T2.Vertices_[0]->GetCoordinates()) return true;

		return false;

	}


};


//*****************************************************************
//	MyLineCL: Darstellung eines Liniensegments					  *
//*****************************************************************

class MyLineCL
{
private: 
	SArrayCL<MyVertexCL*, 2> Vertices_;

public:
	

	//Konstruktor
	MyLineCL(  MyVertexCL* v0,   MyVertexCL* v1){Vertices_[0]=v0; Vertices_[1]=v1;}

	//Copy-Konstruktor
	MyLineCL(const MyLineCL& line){Vertices_=line.Vertices_;}

	
	//member functions

	//Gibt Knoten i der Linie zurück
	const MyVertexCL* getVertex(int i) const { return Vertices_[i]; }



};


//********************************************************************
//	MySpaceGridCL: Darstellung des örtlichen Gitters zum Zeitpunkt t *
//********************************************************************

class MySpaceGridCL
{
public:
	typedef std::list<MyVertexCL> VertexList;
	typedef std::list<MyLineCL> LineList;

	typedef std::list<MyVertexCL>::iterator MyVertexIterator;
	typedef std::list<MyLineCL>::iterator MyLineIterator;

	typedef std::list<MyVertexCL>::const_iterator const_MyVertexIterator;
	typedef std::list<MyLineCL>::const_iterator const_MyLineIterator;
	

private:
	VertexList vertices_;
	LineList lines_;
	double t_; // örtliches Gitter zum Zeitpunkt t

	friend class MySpaceTimeGridCL; // MySpaceTimeGridCL ist friend von MySpaceGridCL und kann deshalb in BuildMyGeometry die private-Member von MySpaceGridCL nutzen
	

public:
	//Default-Konstruktor
	MySpaceGridCL() {}
	//Konstruktor
	MySpaceGridCL( const VertexList& verts, const LineList& lines, const double t)
		: vertices_(verts), lines_(lines), t_(t) {}
	//Copy-Konstruktor
	MySpaceGridCL(const MySpaceGridCL& sg){vertices_=sg.vertices_; lines_=sg.lines_; t_=sg.t_; }
	

	//Iteratoren
	MyVertexIterator GetVerticesBegin()	{ return vertices_.begin(); }
	MyVertexIterator GetVerticesEnd()	{ return vertices_.end();   }
	MyLineIterator GetLineBegin()	{ return lines_.begin();}
	MyLineIterator GetLineEnd()		{ return lines_.end();  }


	//const_Iteratoren (wird benutzt, wenn Instanz der Klasse const ist)
	const_MyVertexIterator GetVerticesBegin()	const { return vertices_.begin(); }
	const_MyVertexIterator GetVerticesEnd()	 const { return vertices_.end();   }
	const_MyLineIterator GetLineBegin()	const { return lines_.begin();}
	const_MyLineIterator GetLineEnd()	const	{ return lines_.end();  }


	//member functions

	//füge Knoten hinzu
	void addVertex( const MyVertexCL& vert) { vertices_.push_back(vert); }
	//füge Linie hinzu
	void addLine(const MyLineCL& line) { lines_.push_back(line); }
	//baut MySpaceGridCL aus MySpaceTimeGridCL, timenew==true -> Gitter bei t=t_new, else -> Gitter bei t=t_old
	void BuildSpaceGrid( MySpaceTimeGridCL& stg, bool timenew);


	//folgende Funktionen machen nur Sinn, wenn BuildSpaceGridCL schon aufgerufen wurde
	
	//Gibt Anzahl der Knoten im SpaceGitter zurück
	int NumVertices() { return vertices_.size(); }
	//Gibt Anzahl der Linien im SpaceGitter zurück
	int NumLines() { return lines_.size(); }
	//Gibt t_ zurück
	double Get_t() { return t_; }

};


//***********************************************************
// MySpaceTimeGridCL: Darstellung des Raum-Zeit-Gitters     *
//***********************************************************


class MySpaceTimeGridCL
{
public:
	typedef std::list<MyVertexCL> VertexList;
	typedef std::list<MyTriangleCL> TriangList;

	typedef std::list<MyVertexCL>::iterator MyVertexIterator;
	typedef std::list<MyTriangleCL>::iterator MyTriangIterator;

	typedef std::list<MyVertexCL>::const_iterator const_MyVertexIterator;
	typedef std::list<MyTriangleCL>::const_iterator const_MyTriangIterator;
	

private:
	VertexList vertices_;
	TriangList triangles_;
	//MySpaceGridCL sg_; // örtliches Gitter zum Zeitpunkt t_n
	double told_; // told_= t_n-1
	double tnew_; // tnew_= t_n
	//int NumDirVerts_; // Anzahl der Dirichlet-Knoten im Gitter

public:
	//Default-Konstruktor
	MySpaceTimeGridCL() {}
	//Konstruktor
	MySpaceTimeGridCL( const VertexList& verts, const TriangList& triang, const MySpaceGridCL sg, const double told, const double tnew)
		: vertices_(verts), triangles_(triang), /*sg_(sg),*/ told_(told), tnew_(tnew) {}
	//Copy-Konstruktor
	MySpaceTimeGridCL(const MySpaceTimeGridCL& stg){vertices_=stg.vertices_; triangles_=stg.triangles_; /*sg_=stg.sg_;*/ told_=stg.told_; tnew_=stg.tnew_; }
	//Destruktor (brauche ich nicht)
	

	//Iteratoren
	MyVertexIterator GetVerticesBegin()	{ return vertices_.begin(); }
	MyVertexIterator GetVerticesEnd()	{ return vertices_.end();   }
	MyTriangIterator GetTriangBegin()	{ return triangles_.begin();}
	MyTriangIterator GetTriangEnd()		{ return triangles_.end();  }


	//const_Iteratoren (wird benutzt, wenn Instanz der Klasse const ist)
	const_MyVertexIterator GetVerticesBegin()	const { return vertices_.begin(); }
	const_MyVertexIterator GetVerticesEnd()	 const { return vertices_.end();   }
	const_MyTriangIterator GetTriangBegin()	const { return triangles_.begin();}
	const_MyTriangIterator GetTriangEnd()	const	{ return triangles_.end();  }


	//member functions

	//füge Knoten hinzu
	void addVertex( const MyVertexCL& vert) { vertices_.push_back(vert); }
	//füge Dreieck hinzu
	void addTriang(const MyTriangleCL& triang) { triangles_.push_back(triang); }
	// Funktion mit der ich meine Geometrie aus Multigitter und levelset-Daten aufbaue, diese Funktion setzt erst sg_, told_ , tnew_ 
	void BuildMyGeometry ( const MultiGridCL& mg, const VecDescCL& levelset  ) ;

	//folgende Funktionen machen nur Sinn, wenn BuildMyGeometry schon aufgerufen wurde
	//Gibt Anzahl der Knoten im Gitter zurück
	int NumVertices() { return vertices_.size(); }
	//Gibt Anzahl der Dreiecke im Gitter zurück
	int NumTriangles() { return triangles_.size(); }
	//Gibt sg_ zurück
	//MySpaceGridCL GetSpaceGrid() {return sg_;}
	//Gibt told_ zurück
	double Get_t_old() { return told_; }
	//Gibt tnew_ zurück
	double Get_t_new() { return tnew_; }
	

};







int GetCounter (const Point3DCL& p, std::map<int, Point3DCL>& map);
//MyTriangleCL Quad2Triang ( int& pp, int& qq, int& rr, int& ss, std::map<int, MyVertexCL*>& pointer );



}// end of namespace DROPS

#endif