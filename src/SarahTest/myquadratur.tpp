// Definitionen der Template-Sachen und inline-Funktionen aus myquadratur.h

namespace DROPS
{

//*****************************************
// Klasse: Quad5SpaceTimeInterfaceCL
//*****************************************

template<class T>
  inline Quad5SpaceTimeInterfaceCL<T>&
  Quad5SpaceTimeInterfaceCL<T>::assign(const MyTriangleCL& s, instat_fun_ptr f , double t)
{
	SArrayCL<Point3DCL, 3> M;
	M[0]= (s.getVertex(0)->GetCoordinates()); M[1]= (s.getVertex(1)->GetCoordinates()); M[2]= (s.getVertex(2)->GetCoordinates());
	Point3DCL* p= QuadraturData.TransformNodes( M,0 );

    for (Uint i= 0; i<Quad5SpaceTimeInterfaceDataCL::NumNodes; ++i)
	{
        (*this)[i]= f( (*(p+i)), t);
	}

    return *this;
}


template<class T>
  inline Quad5SpaceTimeInterfaceCL<T>&
  Quad5SpaceTimeInterfaceCL<T>::assign(const MyLocalP1CL<T>& func)
{
	for(Uint i= 0; i<Quad5SpaceTimeInterfaceDataCL::NumNodes; ++i){
		(*this)[i]= func( Quad5SpaceTimeInterfaceDataCL::Node[i] );	
	}

    return *this;
}

template<class T>
  inline Quad5SpaceTimeInterfaceCL<T>&
  Quad5SpaceTimeInterfaceCL<T>::assign(const MyLocalP2CL<T>& func)
{
	for(Uint i= 0; i<Quad5SpaceTimeInterfaceDataCL::NumNodes; ++i){
		(*this)[i]= func( Quad5SpaceTimeInterfaceDataCL::Node[i] );	
	}

    return *this;
}

template<class T>
  Quad5SpaceTimeInterfaceCL<T>::Quad5SpaceTimeInterfaceCL(const MyTriangleCL& s,
      instat_fun_ptr f, double t)
  : base_type( value_type(), Quad5SpaceTimeInterfaceDataCL::NumNodes)
{
    this->assign( s, f, t);
}



template<class T>
  Quad5SpaceTimeInterfaceCL<T>::Quad5SpaceTimeInterfaceCL(const MyLocalP1CL<T>& func)
  : base_type( value_type(), Quad5SpaceTimeInterfaceDataCL::NumNodes)
{
    this->assign( func);
}

template<class T>
  Quad5SpaceTimeInterfaceCL<T>::Quad5SpaceTimeInterfaceCL(const MyLocalP2CL<T>& func)
  : base_type( value_type(), Quad5SpaceTimeInterfaceDataCL::NumNodes)
{
    this->assign( func);
}


//*****************************************
// Klasse: Quad9Space1DCL
//*****************************************

template<class T>
  inline Quad9Space1DCL<T>&
  Quad9Space1DCL<T>::assign(const MyLineCL& s, instat_fun_ptr f , double t)
{
	SArrayCL<Point3DCL, 2> M;
	M[0]= (s.getVertex(0)->GetCoordinates()); M[1]= (s.getVertex(1)->GetCoordinates()); 
	Point3DCL* p= QuadraturData.TransformNodes( M,0 );

    for (Uint i= 0; i<Quad9Space1DDataCL::NumNodes; ++i)
	{
        (*this)[i]= f( (*(p+i)), t);
	}

    return *this;
}


template<class T>
  inline Quad9Space1DCL<T>&
  Quad9Space1DCL<T>::assign(value_type values[2])
{
	for(Uint i= 0; i<Quad9Space1DDataCL::NumNodes; ++i){
		// Funktionswerte values werden linear interpoliert
		(*this)[i]= values[0]*Quad9Space1DDataCL::Node[i][0] + values[1]*Quad9Space1DDataCL::Node[i][1];	
	}

    return *this;
}



template<class T>
  Quad9Space1DCL<T>::Quad9Space1DCL(const MyLineCL& s,
      instat_fun_ptr f, double t)
  : base_type( value_type(), Quad9Space1DDataCL::NumNodes)
{
    this->assign( s, f, t);
}



template<class T>
  Quad9Space1DCL<T>::Quad9Space1DCL(value_type values[2])
  : base_type( value_type(), Quad9Space1DDataCL::NumNodes)
{
    this->assign( values);
}




//******************************************
// Klasse: MyLocalP1CL
//******************************************

template<class T>
  inline MyLocalP1CL<T>&
  MyLocalP1CL<T>::assign(const MyTriangleCL& s, instat_fun_ptr f, double t)
{
    for (Uint i= 0; i< 3; ++i)
        (*this)[i]= f( s.getVertex(i)->GetCoordinates(), t);
    return *this;
}

template<class T>
  inline MyLocalP1CL<T>&
  MyLocalP1CL<T>::assign(const MyTriangleCL& s, T a[3])
{
    for (Uint i= 0; i< 3; ++i)
        (*this)[i]= a[i];
    return *this;
}

template<class T>
  MyLocalP1CL<T>::MyLocalP1CL(const MyTriangleCL& s, instat_fun_ptr f , double t)
  : base_type( value_type(), 3)
{
    this->assign( s, f, t);
}

template<class T>
  MyLocalP1CL<T>::MyLocalP1CL(const MyTriangleCL& s, T a[3])
  : base_type( value_type(), 3)
{
    this->assign( s, a);
}


template<class T>
  inline typename MyLocalP1CL<T>::value_type
  MyLocalP1CL<T>::operator() (const Point3DCL& p) const
{
	return p[0]*(*this)[0]+p[1]*(*this)[1]+p[2]*(*this)[2]; 
}



//************************************************
//  Klasse: MyLocalP2CL
//************************************************

template<class T>
  inline MyLocalP2CL<T>&
  MyLocalP2CL<T>::assign(const MyTriangleCL& s, instat_fun_ptr f, double t)
{
    for (Uint i= 0; i< 3; ++i){
        (*this)[i]= f( s.getVertex( i)->GetCoordinates(), t);
	}
	Point3DCL edge0= 0.5*s.getVertex( 0)->GetCoordinates() + 0.5*s.getVertex( 1)->GetCoordinates();
	Point3DCL edge1= 0.5*s.getVertex( 1)->GetCoordinates() + 0.5*s.getVertex( 2)->GetCoordinates();
	Point3DCL edge2= 0.5*s.getVertex( 2)->GetCoordinates() + 0.5*s.getVertex( 0)->GetCoordinates();
	(*this)[3]= f(edge0,t);
	(*this)[4]= f(edge1,t);
	(*this)[5]= f(edge2,t);
    
    return *this;
}

template<class T>
  inline MyLocalP2CL<T>&
  MyLocalP2CL<T>::assign(const MyTriangleCL& s, T a[6])
{
    for (Uint i= 0; i< 6; ++i)
        (*this)[i]= a[i];
    return *this;
}


template<class T>
  inline MyLocalP2CL<T>&
  MyLocalP2CL<T>::assign(const MyLocalP1CL<T>& p1)
{
    for (size_t i= 0; i < 3; ++i){
        (*this)[i]= p1[i];
    }
	(*this)[3]= 0.5* (p1[0]+p1[1]);
	(*this)[4]= 0.5* (p1[1]+p1[2]);
	(*this)[5]= 0.5* (p1[2]+p1[0]);
    return *this;
}

template<class T>
  MyLocalP2CL<T>::MyLocalP2CL(const MyTriangleCL& s, instat_fun_ptr f , double t)
  : base_type( value_type(), 6)
{
    this->assign( s, f, t);
}

template<class T>
  MyLocalP2CL<T>::MyLocalP2CL(const MyTriangleCL& s, T a[6])
  : base_type( value_type(), 6)
{
    this->assign( s, a);
}

template<class T>
  MyLocalP2CL<T>::MyLocalP2CL (const MyLocalP1CL<T>& p1)
    : base_type( value_type(), 6)
{
    this->assign( p1);
}

template<class T>
  inline typename MyLocalP2CL<T>::value_type
  MyLocalP2CL<T>::operator() (const Point3DCL& p) const
{
    return (*this)[0]*H0(p)+ (*this)[1]*H1(p)+(*this)[2]*H2(p)+(*this)[3]*H3(p)+(*this)[4]*H4(p)+(*this)[5]*H5(p);
}

}// end of namespace DROPS


    
  
