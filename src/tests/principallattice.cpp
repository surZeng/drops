/// \file principallattice.cpp
/// \brief tests the PrincipalLattice-class
/// \author LNM RWTH Aachen: Joerg Grande, Liang Zhang; SC RWTH Aachen:

/*
 * This file is part of DROPS.
 *
 * DROPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DROPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Copyright 2011 LNM/SC RWTH Aachen, Germany
*/

#include "geom/principallattice.h"
#include "geom/subtriangulation.h"
#include "num/quadrature.h"
#include "misc/container.h"
#include "num/discretize.h"
#include "num/lattice-eval.h"
#include "geom/multigrid.h"

#include <iostream>
#include <sstream>
#include <tr1/unordered_map>


void test_tetra_cut ()
{
    std::cout<<"=========================TetraPartition test: \n"
             <<"all 81 level-set sign patterns will be prescribed, then the corresponding tetra partition visualization files will be created."<<std::endl;
    DROPS::GridFunctionCL<> ls( 4);
    ls[0]= -1.; ls[1]= 0.; ls[2]= 0.; ls[3]= 0.;
    DROPS::TetraPartitionCL tet;
    // tet.partition_principal_lattice<DROPS::SortedVertexPolicyCL, DROPS::MergeCutPolicyCL> ( 1, ls);
    // std::cerr << tet;
    int c= 0;
    for (int i= -1; i <= 1; ++i)
      for (int j= -1; j <= 1; ++j)
        for (int k= -1; k <= 1; ++k)
          for (int l= -1; l <= 1; ++l, ++c) {
              if (i == 0 && j == 0 && k == 0 && l == 0) continue;
              ls[0]= i; ls[1]= j; ls[2]= k; ls[3]= l;
              std::cout << "c: " << c << " ls: " << ls[0] << ' ' << ls[1] << ' ' << ls[2] << ' ' << ls[3] << std::endl;
              DROPS::RefTetraPartitionCL cut( static_cast<double*>(&ls[0]));
//              DROPS::SignPatternTraitCL comb_cut( static_cast<double*>(&ls[0]));
              tet.make_partition<DROPS::SortedVertexPolicyCL, DROPS::MergeCutPolicyCL> ( DROPS::PrincipalLatticeCL::instance( 1), ls);
//              if (c == 5) {
//                  std::cerr << comb_cut << std::endl;
//                  std:: cerr << cut << std::endl;
//                  std::cerr << tet << std::endl << std::endl;
//              }
              std::ostringstream name;
              name << "hallo" << c << ".vtu";
              std::ofstream file( name.str().c_str());
              DROPS::write_paraview_vtu( file, tet);
          }
}

void test_cut_surface ()
{
    std::cout<<"=========================Surface patch test: \n"
             <<"all 81 level-set sign patterns will be prescribed in a tetra, then the corresponding interface patch visualization files will be created."<<std::endl;
    DROPS::GridFunctionCL<> ls( 4);
    ls[0]= -1.; ls[1]= 0.; ls[2]= 0.; ls[3]= 0.;
    DROPS::SurfacePatchCL tet;
    // tet.partition_principal_lattice ( 1, ls);
    // std::cerr << tet;
    int c= 0;
    for (int i= -1; i <= 1; ++i)
      for (int j= -1; j <= 1; ++j)
        for (int k= -1; k <= 1; ++k)
          for (int l= -1; l <= 1; ++l, ++c) {
              if (i == 0 && j == 0 && k == 0 && l == 0) continue;
              ls[0]= i; ls[1]= j; ls[2]= k; ls[3]= l;
              std::cout << "c: " << c << " ls: " << ls[0] << ' ' << ls[1] << ' ' << ls[2] << ' ' << ls[3] << std::endl;
              //DROPS::RefTetraPartitionCL cut( static_cast<double*>(&ls[0]));
              //DROPS::SignPatternTraitCL comb_cut( static_cast<double*>(&ls[0]));
              tet.make_patch<DROPS::MergeCutPolicyCL>( DROPS::PrincipalLatticeCL::instance( 1), ls);
              std::ostringstream name;
              name << "hallo_surf" << c << ".vtu";
              std::ofstream file( name.str().c_str());
              DROPS::write_paraview_vtu( file, tet);
          }
}

void test_principal_lattice ()
{
    std::cout<<"=========================PrincipalLatticeCL test: \n"
             <<"4 principal lattice with 1, 2, 3, 4 intervals will be created respectively, and key information will be showed."<<std::endl;
    for (int i= 1; i <= 4; ++i) {
        const DROPS::PrincipalLatticeCL& lat= DROPS::PrincipalLatticeCL::instance( i);
        std::cout << "======================================= \n" 
                  << "Number of intervals: "<<lat.num_intervals() 
                  << "| " << "Number of vertices: "<<lat.vertex_size() 
                  << "| " << "Number of tetra: "   <<lat.tetra_size() << std::endl;
        std::cout << "=======================================Barycentric coordinates of vertices: "<<std::endl; 
        for (DROPS::PrincipalLatticeCL::const_vertex_iterator v= lat.vertex_begin(), end= lat.vertex_end(); v != end; ++v) {
            std::cout << /*lat.num_intervals()*/(*v) << std::endl;
        }
        std::cout << "=======================================Indices of each sub-tetra:" << std::endl;
        for (DROPS::PrincipalLatticeCL::const_tetra_iterator v= lat.tetra_begin(), end= lat.tetra_end(); v != end; ++v) {
            std::cout << (*v)[0] << ' '  << (*v)[1] << ' ' << (*v)[2] << ' ' << (*v)[3] << ' ' << std::endl;
        }
    }
}

inline double sphere (const DROPS::Point3DCL& p)
{
    return p.norm() - 0.5;
}

inline double sphere_instat (const DROPS::Point3DCL& p, double)
{
    return sphere( p);
}

void test_sphere_cut ()
{
    std::cout<<"=========================Sphere cut test: \n"
             <<"A principal lattice with 10 intervals will be cut by a sphere with 0.5 radius;\n"
             <<"Negative tetras and interface patch will be showed in two vtu files respectively."<<std::endl;

    DROPS::TetraBuilderCL tetrabuilder( 0);
    DROPS::MultiGridCL mg( tetrabuilder);

    const DROPS::PrincipalLatticeCL& lat= DROPS::PrincipalLatticeCL::instance( 10);
    DROPS::GridFunctionCL<> ls( lat.vertex_size());
    evaluate_on_vertexes( &sphere_instat, *mg.GetAllTetraBegin(), lat, 0., Addr( ls));
    DROPS::TetraPartitionCL tet;
    tet.make_partition<DROPS::SortedVertexPolicyCL, DROPS::MergeCutPolicyCL>( lat, ls);
    std::ostringstream name;
    name << "sphere.vtu";
    std::ofstream file( name.str().c_str());
    DROPS::write_paraview_vtu( file, tet, DROPS::NegTetraC);
    file.close();

    DROPS::SurfacePatchCL surf;
    surf.make_patch<DROPS::MergeCutPolicyCL>( lat, ls);
    name.str( "");
    name << "sphere_surf.vtu";
    file.open( name.str().c_str());
    DROPS::write_paraview_vtu( file, surf);
}


void test_sphere_integral ()
{
    std::cout<<"=========================Volume integral test: \n"
             <<"A 2x2x2 cubic is cut by a sphere with 0.5 radius;\n"
             <<"Negative part and positive part will be integrated seperately by using QuadDomainCL."<<std::endl;
    //std::cout << "Enter the number of subdivisions of the cube: ";
    DROPS::Uint num_sub = 8;
    //std::cin >> num_sub;
    //std::cout << "Enter the number of subdivisions of the principal lattice: ";
    DROPS::Uint num_sub_lattice = 2;
    //std::cin >> num_sub_lattice;
    DROPS::BrickBuilderCL brick(DROPS::Point3DCL( -1.),
                                2.*DROPS::std_basis<3>(1),
                                2.*DROPS::std_basis<3>(2),
                                2.*DROPS::std_basis<3>(3),
                                num_sub, num_sub, num_sub);
    DROPS::MultiGridCL mg( brick);
    const DROPS::PrincipalLatticeCL& lat= DROPS::PrincipalLatticeCL::instance( num_sub_lattice);
    DROPS::GridFunctionCL<> ls( lat.vertex_size());
    DROPS::TetraPartitionCL tet;
    // DROPS::SurfacePatchCL patch;
    double vol_neg= 0., vol_pos= 0.;
    // double surf= 0.;
    DROPS::QuadDomainCL qdom;

    DROPS_FOR_TRIANG_TETRA( mg, 0, it) {
        evaluate_on_vertexes( sphere_instat, *it, lat, 0., Addr( ls));
        // tet.make_partition<DROPS::UnorderedVertexPolicyCL, DROPS::MergeCutPolicyCL>( lat, ls);
        // tet.make_partition<DROPS::SortedVertexPolicyCL, DROPS::MergeCutPolicyCL>( lat, ls);
        tet.make_partition<DROPS::PartitionedVertexPolicyCL, DROPS::MergeCutPolicyCL>( lat, ls);
        // patch.make_partition<DROPS::SortedVertexPolicyCL, DROPS::MergeCutPolicyCL>( lat, ls);
        DROPS::make_CompositeQuad5Domain( qdom, tet);
        DROPS::GridFunctionCL<> integrand( 1., qdom.vertex_size());
        double tmp_neg, tmp_pos;
        quad( integrand, it->GetVolume()*6., qdom, tmp_neg, tmp_pos);
        vol_neg+= tmp_neg; vol_pos+= tmp_pos;
        // q5_2d.assign( patch);
        // DROPS::GridFunctionCL<> surf_integrand( 1., q5_2d.size());
        // surf+= q5_2d.quad( surf_integrand);

    }
    std::cout << "Volume of the negative part: " << vol_neg << ", volume of the positive part: " << vol_pos << std::endl;
}

void test_extrapolated_sphere_integral ()
{
    std::cout << "Enter the number of subdivisions of the cube: ";
    DROPS::Uint num_sub;
    std::cin >> num_sub;
    std::cout << "Enter the number of extrapolation levels: ";
    DROPS::Uint num_level;
    std::cin >> num_level;
    DROPS::BrickBuilderCL brick(DROPS::Point3DCL( -1.),
                                2.*DROPS::std_basis<3>(1),
                                2.*DROPS::std_basis<3>(2),
                                2.*DROPS::std_basis<3>(3),
                                num_sub, num_sub, num_sub);
    DROPS::MultiGridCL mg( brick);
    double vol_neg= 0., vol_pos= 0.;
    DROPS::QuadDomainCL qdom;
    DROPS::ExtrapolationToZeroCL extra( num_level, DROPS::RombergSubdivisionCL());
    // DROPS::ExtrapolationToZeroCL extra( num_level, DROPS::HarmonicSubdivisionCL());

    DROPS_FOR_TRIANG_TETRA( mg, 0, it) {
        DROPS::LocalP2CL<> ls_loc( *it, &sphere_instat);
        make_ExtrapolatedQuad5Domain( qdom, ls_loc, extra);
        DROPS::GridFunctionCL<> integrand( 1., qdom.vertex_size());
        double tmp_neg, tmp_pos;
        quad( integrand, it->GetVolume()*6., qdom, tmp_neg, tmp_pos);
        vol_neg+= tmp_neg; vol_pos+= tmp_pos;
    }
    std::cout << "Volume of the negative part: " << vol_neg << ", volume of the positive part: " << vol_pos << std::endl;
}

void test_sphere_surface_integral ()
{
    std::cout<<"=========================Interface area integral test: \n"
             <<"A 2x2x2 cubic is cut by a sphere with 0.5 radius;\n"
             <<"Interface area will be computed by using QuadDomain2DCL"<<std::endl;
    //std::cout << "Enter the number of subdivisions of the cube: ";
    DROPS::Uint num_sub = 8;
    //std::cin >> num_sub;
    //std::cout << "Enter the number of subdivisions of the principal lattice: ";
    DROPS::Uint num_sub_lattice =2;
    //std::cin >> num_sub_lattice;
    DROPS::BrickBuilderCL brick(DROPS::Point3DCL( -1.),
                                2.*DROPS::std_basis<3>(1),
                                2.*DROPS::std_basis<3>(2),
                                2.*DROPS::std_basis<3>(3),
                                num_sub, num_sub, num_sub);
    DROPS::MultiGridCL mg( brick);
    const DROPS::PrincipalLatticeCL& lat= DROPS::PrincipalLatticeCL::instance( num_sub_lattice);
    DROPS::GridFunctionCL<> ls( lat.vertex_size());
    DROPS::SurfacePatchCL patch;
    double surf= 0.;
    DROPS::QuadDomain2DCL qdom;

    DROPS_FOR_TRIANG_TETRA( mg, 0, it) {
        evaluate_on_vertexes( &sphere_instat, *it, lat, 0., Addr( ls));
        patch.make_patch<DROPS::MergeCutPolicyCL>( lat, ls);
        DROPS::make_CompositeQuad5Domain2D( qdom, patch, *it);
        DROPS::GridFunctionCL<> integrand( 1., qdom.vertex_size());
        surf+= quad_2D( integrand, qdom);
    }
    std::cout << "Surface: " <<surf << std::endl;
}

void test_extrapolated_sphere_surface_integral ()
{
    std::cout << "Enter the number of subdivisions of the cube: ";
    DROPS::Uint num_sub;
    std::cin >> num_sub;
    std::cout << "Enter the number of extrapolation levels: ";
    DROPS::Uint num_level;
    std::cin >> num_level;
    DROPS::BrickBuilderCL brick(DROPS::Point3DCL( -1.),
                                2.*DROPS::std_basis<3>(1),
                                2.*DROPS::std_basis<3>(2),
                                2.*DROPS::std_basis<3>(3),
                                num_sub, num_sub, num_sub);
    DROPS::MultiGridCL mg( brick);
    double surf= 0.;
    DROPS::QuadDomain2DCL qdom;
    DROPS::ExtrapolationToZeroCL extra( num_level, DROPS::RombergSubdivisionCL());
    // DROPS::ExtrapolationToZeroCL extra( num_level, DROPS::HarmonicSubdivisionCL());

    DROPS_FOR_TRIANG_TETRA( mg, 0, it) {
        DROPS::LocalP2CL<> ls_loc( *it, &sphere_instat);
        make_ExtrapolatedQuad5Domain2D( qdom, ls_loc, *it, extra);
        DROPS::GridFunctionCL<> integrand( 1., qdom.vertex_size());
        surf+= quad_2D( integrand, qdom);
    }
    std::cout << "Surface: " << surf << std::endl;
}

int main()
{
    try {
        test_tetra_cut();
        test_cut_surface();
        test_principal_lattice();
        test_sphere_cut();
        test_sphere_integral();
        //test_extrapolated_sphere_integral();
        test_sphere_surface_integral();
        //test_extrapolated_sphere_surface_integral();
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }
    return 0;
}
