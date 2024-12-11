// Klasse für die Ausgabe der Geometrie mit VTK (Paraview)

#ifndef DROPS_VTKOUTPUT_H
#define DROPS_VTKOUTPUT_H

#include "SarahTest/myfespace.h"
#include "misc/problem.h"
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>

namespace DROPS
{

// forward declaration
class MyVTKOutCL;


class MyVTKOutCL
{

private:
	const MySpaceTimeGridCL& stg_; // Referenz auf mein Gitter
	std::string dirname_;			// Verzeichnisname
	std::string filename_;			// Dateiname
	std::ofstream file_;			 // Datei, in die Daten reingeschrieben werden
	
	Uint numVerts_;					// Anzahl der Knoten im Gitter
	Uint numTriangs_;				// Anzahl der Dreiecke im Gitter
	std::map<const MyVertexCL*, Uint> vertexmap_; // Speichert zu einer Adresse eines Knoten einen Int (aufsteigend)
	VectorBaseCL<float> coord_;		// Koordinaten der Knoten im Gitter
	VectorBaseCL<Uint> triang_;		// Konnektivität für die Dreiecke (welche coord gehören zu welchen Dreieck)

	VectorCL solution_; // numerische Daten der skalaren Lösung
	VectorCL residuum_; // numerische Daten des Residuums (zum Debuggen)

public:
	// Konstruktor
	MyVTKOutCL(const MySpaceTimeGridCL& stg, const VectorCL& solution, const VectorCL& residuum, std::string& dirname, std::string& filename)
		: stg_(stg), solution_(solution), residuum_(residuum), dirname_(dirname), filename_(filename) {}

	// Destruktor
	//~MyVTKOutCL();

	// member functions
	void NewFile();
	void PutHeader();
	void PutFooter();
	void Clear();
	void CheckFile( const std::ofstream& os) const;
	void GatherCoord();
	void WriteCoord();
	void GatherTriangles();
	void WriteTriangles();
	void WriteScalarSolution();
	void WriteScalarSolutionResiduum();
	void PutGeometry();
	void PutGeometryAndSolution();
	void PutGeometrySolutionResiduum(); // Residuum zum Debuggen






};


} // end of namespace DROPS

#endif