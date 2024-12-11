// Definitionen aus vtkoutput.h

#include "SarahTest/vtkoutput.h"
//#include "SarahTest/mygeometry.h"

namespace DROPS
{


void MyVTKOutCL::NewFile()
{
	std::string filename(filename_);
	filename+= ".vtu";
    file_.open((dirname_+"/"+filename).c_str());
	if(!file_){
		CreateDirectory( dirname_);
        file_.open((dirname_+filename).c_str());
	}
	CheckFile(file_);
	PutHeader();
	
}

void MyVTKOutCL::PutHeader()
/** Writes the header into the VTK file*/
{
    file_ << "<?xml version=\"1.0\"?>\n"     // this is just the XML declaration, it's unnecessary for the actual VTK file
             "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
             "<UnstructuredGrid>\n";
}

void MyVTKOutCL::PutFooter()
/** Closes the file XML conform*/
{
    file_ <<"\n</Piece>" // ohne PointData
            "\n</UnstructuredGrid>"
            "\n</VTKFile>";
}

void MyVTKOutCL::Clear()
{
    if (numVerts_>0){

        vertexmap_.clear();
        coord_.resize(0);
        triang_.resize(0);
    }
}


void MyVTKOutCL::GatherCoord() // Hier werden Vertices gesammelt. Es wird auch die Anzahl der Dreiecke im Gitter ermittelt
{
	Uint numpoints = 0;
	Uint numtriang = 0;
	for(MySpaceTimeGridCL::const_MyVertexIterator it=stg_.GetVerticesBegin(); it!=stg_.GetVerticesEnd(); it++){
		numpoints++;
	}
	for(MySpaceTimeGridCL::const_MyTriangIterator it=stg_.GetTriangBegin(); it!=stg_.GetTriangEnd(); it++){
		numtriang++;
	}

	numTriangs_=numtriang;
	numVerts_=numpoints;
	//std:: cout << "NumVerts= " << numVerts_ << "\n";
	//std:: cout << "NumTriang= " << numTriangs_ << "\n";
	coord_.resize(3*numVerts_);
	Uint counter = 0;
	for(MySpaceTimeGridCL::const_MyVertexIterator it=stg_.GetVerticesBegin(); it!=stg_.GetVerticesEnd(); it++){
		// speicher aufsteigende Nummer für Vertex
		vertexmap_[&(*it)]=counter;
		for(int i=0; i<3; i++){
			// Schreibe Koordinaten des Vertex in coord_
			coord_[3*counter + i] = (float)(it->GetCoordinates()[i]);	
			//std:: cout << "in datei: "<< coord_[3*counter + i] << "\n";
		}
		++counter;
	}

}


void MyVTKOutCL::WriteCoord()
{
	//if(file_.is_open()) {std::cout << "file offen\n";}
 file_<< "<Piece NumberOfPoints=\""<<numVerts_<<"\" NumberOfCells=\""<<numTriangs_<<"\">"
            "\n\t<Points>"
            "\n\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n\t\t"; 

    
        for (Uint i=0; i<numVerts_; ++i){
            file_<< coord_[3*i+0] << ' ' << coord_[3*i+1] << ' ' << coord_[3*i+2]<< ' ';
		}

    file_<< "\n\t\t</DataArray> \n"
            "\t</Points>\n";

}

void MyVTKOutCL::GatherTriangles()
{
 
    triang_.resize(3*numTriangs_);      //  3 * Anzahl der Dreiecke in Gitter

    // Sammel Konnektivität
    Uint counter=0;
    for (MySpaceTimeGridCL::const_MyTriangIterator it= stg_.GetTriangBegin(); it!=stg_.GetTriangEnd(); ++it){ 
        for (int vert= 0; vert<3; ++vert){
            triang_[counter] = vertexmap_[it->getVertex(vert)]; // warum triang_[counter++] und nicht triang_[counter]...; counter++;
			counter++;
		}
       
    }

}

void MyVTKOutCL::WriteTriangles()
{

file_   << "\t<Cells>\n"
	"\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n\t\t";

        
        // Write out connectivities
        
            for (Uint i=0; i<numTriangs_; ++i)
            {
                file_ << triang_[3*i+0] << ' '<< triang_[3*i+1] << ' '<< triang_[3*i+2] << " ";
            }

file_ << "\n\t\t</DataArray>\n"
          "\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n\t\t";

    // Write out offsets
      for(Uint i=1; i<=numTriangs_; ++i) file_ << i*3<<" "; 

	  //Write Cell Type (vtk-triangle = 5, vtk-triangle-strip = 6)
	  file_ << "\n\t\t</DataArray>"
             "\n\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n\t\t";
    const char* TriangType= "5";
    for(Uint i=1; i<=numTriangs_; ++i)
        file_ << TriangType << " ";

    file_ << "\n\t\t</DataArray>"
             "\n\t</Cells>";
    
   


}

void MyVTKOutCL::WriteScalarSolution(){
 
	file_<< "\n\t<PointData Scalars= \"solution\">"
		"\n\t\t<DataArray type=\"Float32\" Name=\"solution\" NumberOfComponents=\"1\" format=\"ascii\">"
		"\n\t\t";

	for(int i=0; i<solution_.size(); i++){
		file_<< solution_[i] << " ";
	}

	file_<< "\n\t\t</DataArray>"
		"\n\t</PointData>";


}

void MyVTKOutCL::WriteScalarSolutionResiduum(){
 
	file_<< "\n\t<PointData Scalars= \"solution, residuum\">"
		"\n\t\t<DataArray type=\"Float32\" Name=\"solution\" NumberOfComponents=\"1\" format=\"ascii\">"
		"\n\t\t";

	for(int i=0; i<solution_.size(); i++){
		file_<< solution_[i] << " ";
	}

	file_<< "\n\t\t</DataArray>";
	
	file_<< "\n\t\t<DataArray type=\"Float32\" Name=\"residuum\" NumberOfComponents=\"1\" format=\"ascii\">"
		"\n\t\t";
	for(int i=0; i<residuum_.size(); i++){
		file_<< residuum_[i] << " ";
	}
	file_<< "\n\t\t</DataArray>";
	file_<<	"\n\t</PointData>";



	


}

void MyVTKOutCL::CheckFile( const std::ofstream& os) const
/** Checks if a file is open*/
{
    if (!os.is_open())
        throw DROPSErrCL( "MyVTKOutCL: error while opening file!");
}


void MyVTKOutCL::PutGeometry()
/** At first the geometry is put into the VTK file. Therefore this procedure
    opens the file and writes description into the file.
    
*/
{
	//std::cout << "NewFile()\n";
    NewFile();
	//std::cout << "Clear()\n";
    Clear();
	//std::cout << "GatherCoord()\n";
    GatherCoord();
	//std::cout << "GatherTriangles()\n";
    GatherTriangles();
	//file_.close();
	//std::cout << "WriteCoord()\n";
	WriteCoord();
	//std::cout << "WriteTriangles()\n";
	WriteTriangles();
	//std::cout << "PutFooter()\n";
	PutFooter();
}


void MyVTKOutCL::PutGeometryAndSolution()
{

    NewFile();
	
    Clear();
	
    GatherCoord();
	
    GatherTriangles();
	
	WriteCoord();
	
	WriteTriangles();
	
	WriteScalarSolution();
	PutFooter();


}



void MyVTKOutCL::PutGeometrySolutionResiduum()
{

    NewFile();
	
    Clear();
	
    GatherCoord();
	
    GatherTriangles();
	
	WriteCoord();
	
	WriteTriangles();
	
	WriteScalarSolutionResiduum();
	PutFooter();


}


}// end of namespace DROPS