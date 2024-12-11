// Klasse für meine adaptive Verfeinerung mit Krümmung

#ifndef DROPS_MYADAP_H
#define DROPS_MYADAP_H

#include "geom/multigrid.h"
#include "levelset/levelset.h"

namespace DROPS
{

	

	class MyAdapCL
	{
	private:
		MultiGridCL& mg_;
		double width_;                                          ///< width of the refined grid
		int c_level_, f1_level_, f2_level_;                     ///< coarsest level: c_level_, 
																///< finest level at interface without high curvature: f1_level_
																///< finest level at interface with high curvature: f2_level_
																/// c_level < f1_level < f2_level
		bool modified_;
		double konst_;		// konst_ < 1								/// wähle h so, dass h<= konst_ * 1/Kappa_max
		double initialmeshsize_; // initialmeshsize_ = Gitterweite des Ausgangsgitters

		// Wertet Levelsetfunktion und kappafunktion aus
		double GetValue( instat_scalar_fun_ptr dist, const VertexCL& v, double t=0.)  { return dist( v.GetCoord(), t); }
		double GetValue( instat_scalar_fun_ptr dist, const EdgeCL& e, double t=0.)    { return dist( GetBaryCenter( e), t); }
		double GetValue( instat_scalar_fun_ptr dist, const TetraCL& tet, double t=0.) { return dist( GetBaryCenter( tet), t); }

		//mache einen Schritt der adaptiven Verfeinerung
		bool ModifyGridStep( instat_scalar_fun_ptr Dist, instat_scalar_fun_ptr Kappa);

	public:
		MyAdapCL(  MultiGridCL& mg, double width, int c_level, int f1_level, int f2_level, double konst, double initialmeshsize)
			: mg_( mg), width_(width), c_level_(c_level), f1_level_(f1_level), f2_level_(f2_level), modified_(false), konst_(konst), initialmeshsize_(initialmeshsize) {}

		// mache Triangulierung bezüglich levelsetfunktion und krümmung kappa
		void MakeInitialTriang (instat_scalar_fun_ptr Dist, instat_scalar_fun_ptr Kappa);
	
	
	};



}

#endif