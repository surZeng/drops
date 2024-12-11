// Definitionen aus myadap.h

#include "SarahTest/myadap.h"

namespace DROPS
{

	bool MyAdapCL::ModifyGridStep( instat_scalar_fun_ptr Dist, instat_scalar_fun_ptr Kappa)
	{
	
		bool modified= false;
		for (MultiGridCL::TriangTetraIteratorCL it= mg_.GetTriangTetraBegin(),
			 end= mg_.GetTriangTetraEnd(); it!=end; ++it)
		{
			double d= 1e99;
			double k= 0;
			int num_pos= 0;
			for (Uint j=0; j<4; ++j)
			{
				const double dist= GetValue( Dist, *it->GetVertex( j));
				const double kappa= GetValue(Kappa, *it->GetVertex( j));
				if (dist>=0) ++num_pos;
				d= std::min( d, std::abs( dist));
				k= std::max( k, std::abs( kappa));
			}
			for (Uint j=0; j<6; ++j)
			{
				const double dist= GetValue( Dist, *it->GetEdge( j));
				const double kappa= GetValue(Kappa, *it->GetEdge( j));
				if (dist>=0) ++num_pos;
				d= std::min( d, std::abs( dist));
				k= std::max( k, std::abs( kappa));
			}
			d= std::min( d, std::abs( GetValue( Dist, *it)));
			k= std::max( k, std::abs( GetValue( Kappa, *it)));

			double ref= konst_/k; // finde h, so dass h<= ref und dann soll-level aus diesem h berechnen


			const bool vzw= num_pos!=0 && num_pos!=10; // change of sign
			const Uint l= it->GetLevel();
			// In the shell:      level should be f_level1_, if curvature is not high, f_level2, else.
			// Outside the shell: level should be c_level_.
			Uint soll_level;
				//= (d<=width_ || vzw) ? f_level_ : c_level_;
			if(d<=width_ || vzw){
				// h=initial_/2^(soll-level) <= ref <=> soll-level= obere_gaußklammer(log_2(initial_/ref))
				double help= std::ceil( std::log(initialmeshsize_/ref)/std::log(2.0) );
				if(help >=0){ soll_level= help;
				}else{ soll_level= 0;}// konvertierung von double zu int
				if(soll_level < f1_level_) soll_level= f1_level_;
				if(soll_level > f2_level_) soll_level= f2_level_;
				//std::cout << "soll_level= " << soll_level << "\n";
			}else{
				soll_level= c_level_;
			}

			if (l !=  soll_level || (l == soll_level && !it->IsRegular()) )
			{ // tetra will be marked for refinement/removement
				modified= true;
				if (l <= soll_level){
					//if(soll_level==6 ) std::cout << "l kleiner gleich 6\n";
					it->SetRegRefMark();
				}
				else {// l > soll_level
					it->SetRemoveMark();
					//if(soll_level==5 ) std::cout << "l größer als 5\n";
				}
			}
		}

		if (modified) {
			mg_.Refine();
		}
       

		return modified;
	
	}




	void MyAdapCL::MakeInitialTriang (instat_scalar_fun_ptr Dist, instat_scalar_fun_ptr Kappa)
	{
		
		TimerCL time;
		time.Reset();
		time.Start();

		const Uint min_ref_num= f2_level_;
		Uint i;
		bool modified= true;
		for (i=0; i<2*min_ref_num && modified; ++i)
			modified=ModifyGridStep( Dist, Kappa);

		time.Stop();
		const double duration=time.GetTime();

		std::cout << "MakeInitialTriang: " << i
                << " refinements in " << duration << " seconds\n"
                << "last level: " << mg_.GetLastLevel() << '\n';
		mg_.SizeInfo( std::cout);
	
	
	}




}