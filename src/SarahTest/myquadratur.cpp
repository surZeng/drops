// Definitionen aus myquadratur.h

#include "SarahTest/myquadratur.h"

namespace DROPS
{

//------------------------------ Quad5SpaceTimeInterfaceDataCL------------------------------------------------------------



Point3DCL*
Quad5SpaceTimeInterfaceDataCL::TransformNodes (const SArrayCL<Point3DCL,3>& M, Point3DCL* p)
{
    if (!p) p= new Point3DCL[NumNodes];
    for (Uint i=0; i < NumNodes; ++i)
        //p[i]=M*Node[i]; M (als Matrix) ist spaltenweise gespeichert!
        for (Uint k= 0; k < 3; ++k)
            p[i][k]= M[0][k]*Node[i][0] + M[1][k]*Node[i][1] + M[2][k]*Node[i][2];
    return p;
}

const double Quad5SpaceTimeInterfaceDataCL::Weight[NumNodes]= {
    9./80.,

    (155.+std::sqrt(15.))/2400.,
    (155.+std::sqrt(15.))/2400.,
    (155.+std::sqrt(15.))/2400.,

    (155.-std::sqrt(15.))/2400.,
	(155.-std::sqrt(15.))/2400.,
	(155.-std::sqrt(15.))/2400.
};

// muss man vordeklarieren, weil Node statisch ist und in Quad5SpaceTimeInterfaceDataCL::Quad5SpaceTimeInterfaceDataCL() benutzt wird
Point3DCL Quad5SpaceTimeInterfaceDataCL::Node[NumNodes];

Quad5SpaceTimeInterfaceDataCL::Quad5SpaceTimeInterfaceDataCL()
{
	double a= (9.- 2. * std::sqrt(15.))/ 21.;
	double b= (6. + std::sqrt(15.))/ 21.;
	double c= (9. + 2. * std::sqrt(15.))/ 21.;
	double d= (6. - std::sqrt(15.))/ 21.;

    Node[0]= MakePoint3D(1./3., 1./3., 1./3. );

    Node[1]= MakePoint3D(a, b, b);
    Node[2]= MakePoint3D(b, a, b);
    Node[3]= MakePoint3D(b, b, a);

    Node[4]= MakePoint3D(c, d, d);
	Node[5]= MakePoint3D(d, c, d);
	Node[6]= MakePoint3D(d, d, c);
}




//-----------------------Quad9Space1DDataCL-----------------------------------------------


Point3DCL*
Quad9Space1DDataCL::TransformNodes (const SArrayCL<Point3DCL,2>& M, Point3DCL* p)
{
    if (!p) p= new Point3DCL[NumNodes];
    for (Uint i=0; i < NumNodes; ++i)
        //p[i]=M*Node[i]; M (als Matrix) ist spaltenweise gespeichert!
        for (Uint k= 0; k < 3; ++k)
            p[i][k]= M[0][k]*Node[i][0] + M[1][k]*Node[i][1];
    return p;
}

const double Quad9Space1DDataCL::Weight[NumNodes]= {
    0.5*(322.-13.*std::sqrt(70.))/900.,
	0.5*(322.+13.*std::sqrt(70.))/900.,
    
	0.5*128./225.,

	0.5*(322.+13.*std::sqrt(70.))/900.,
	0.5*(322.-13.*std::sqrt(70.))/900.,
};

// muss man vordeklarieren, weil Node statisch ist und in Quad9Space1DDataCL::Quad9Space1DDataCL() benutzt wird
Point2DCL Quad9Space1DDataCL::Node[NumNodes];

Quad9Space1DDataCL::Quad9Space1DDataCL()
{
	double a= -1./6.* std::sqrt(5.+2.*std::sqrt(10./7.)) + 1./2.;
	double b= -1./6.* std::sqrt(5.-2.*std::sqrt(10./7.)) + 1./2.;
	double c= 1./2.;
	double d= -b + 1.;
	double e= -a + 1.;

    Node[0]= MakePoint2D(1.-a, a );
	Node[1]= MakePoint2D(1.-b, b);
	
    Node[2]= MakePoint2D(1.-c, c);

    Node[3]= MakePoint2D(1.-d, d);
    Node[4]= MakePoint2D(1.-e, e);
	
}














}// end of namespace DROPS