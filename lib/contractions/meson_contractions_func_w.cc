/*! 
 *  Functions that do meson spin and color contractions.
 *  Color is nested inside of spin.
 *  Authors:
 *  Arjun Gambhir
 *  David Brantley
 */

/* CAVEAT EMPTOR: The tensor elements WILL switch sign under a different convention.
 * For example, the FH routine as it is uses the conventions s_mu_nu = 1./2 [g_mu,g_nu],
 * while Chroma uses the convention s_mu_nu = i/2 [g_mu,g_nu]. The CMDM operator in glue
 * uses the convention s_mu_nu = i./2 [g_mu,g_nu]. The different conventions incur a different
 * sign under transformation g5 Gamma dagger g4
 */

#include "chromabase.h"
#include "util/spin_basis.h"
#include "meson_contractions_func_w.h"

namespace Chroma
{

  void anti_quark_FH_prop(LatticePropagator & FH_quark_prop,
                          LatticePropagator & FH_antiquark_prop,
                          std::string& cur)
  {
    int which_gamma_five = Ns*Ns-1;

    if(cur == "S")
    {
    FH_antiquark_prop = adj(Gamma(which_gamma_five)*FH_quark_prop*Gamma(which_gamma_five));
    }
    if(cur == "P")
    {
    FH_antiquark_prop = adj(Gamma(which_gamma_five)*FH_quark_prop*Gamma(which_gamma_five));
    }
    if(cur == "V1")
    {
    FH_antiquark_prop = -adj(Gamma(which_gamma_five)*FH_quark_prop*Gamma(which_gamma_five));
    }
    if(cur == "V2")
    {
    FH_antiquark_prop = -adj(Gamma(which_gamma_five)*FH_quark_prop*Gamma(which_gamma_five));
    }
    if(cur == "V3")
    {
    FH_antiquark_prop = -adj(Gamma(which_gamma_five)*FH_quark_prop*Gamma(which_gamma_five));
    }
    if(cur == "V4")
    {
    FH_antiquark_prop = -adj(Gamma(which_gamma_five)*FH_quark_prop*Gamma(which_gamma_five));
    }
    if(cur == "A1")
    {
    FH_antiquark_prop = adj(Gamma(which_gamma_five)*FH_quark_prop*Gamma(which_gamma_five));
    }
    if(cur == "A2")
    {
    FH_antiquark_prop = adj(Gamma(which_gamma_five)*FH_quark_prop*Gamma(which_gamma_five));
    }
    if(cur == "A3")
    {
    FH_antiquark_prop = adj(Gamma(which_gamma_five)*FH_quark_prop*Gamma(which_gamma_five));
    }
    if(cur == "A4")
    {
    FH_antiquark_prop = adj(Gamma(which_gamma_five)*FH_quark_prop*Gamma(which_gamma_five));
    }
    if(cur == "T12")
    {
    FH_antiquark_prop = adj(Gamma(which_gamma_five)*FH_quark_prop*Gamma(which_gamma_five));
    }
    if(cur == "T13")
    {
    FH_antiquark_prop = adj(Gamma(which_gamma_five)*FH_quark_prop*Gamma(which_gamma_five));
    }
    if(cur == "T14")
    {
    FH_antiquark_prop = adj(Gamma(which_gamma_five)*FH_quark_prop*Gamma(which_gamma_five));
    }
    if(cur == "T23")
    {
    FH_antiquark_prop = adj(Gamma(which_gamma_five)*FH_quark_prop*Gamma(which_gamma_five));
    }
    if(cur == "T24")
    {
    FH_antiquark_prop = adj(Gamma(which_gamma_five)*FH_quark_prop*Gamma(which_gamma_five));
    }
    if(cur == "T34")
    {
    FH_antiquark_prop = adj(Gamma(which_gamma_five)*FH_quark_prop*Gamma(which_gamma_five));
    }
    if(cur == "CHROMO_MAG")
    {
    FH_antiquark_prop = -adj(Gamma(which_gamma_five)*FH_quark_prop*Gamma(which_gamma_five));
    }
  }


  void FH_I_one_Iz_pm_one_contract(LatticePropagator & quark_1,
                             LatticePropagator & quark_2,
                             LatticeComplex & contracted,
                             bool& is_FH_antiquark,
                             std::string& cur)
  {
    //quark_1 is the forward (possibly FH) prop, quark_2 is the antiquark.

       int which_gamma_five = Ns*Ns-1;

	   if(is_FH_antiquark == false)
       {
        LatticePropagator anti_quark_prop = adj(Gamma(which_gamma_five)*quark_2*Gamma(which_gamma_five));
        contracted = trace(Gamma(which_gamma_five)*quark_1*Gamma(which_gamma_five)*anti_quark_prop);
       }
       else if (is_FH_antiquark == true)
       {
        LatticePropagator anti_quark_prop;

        anti_quark_FH_prop(quark_2,anti_quark_prop,cur);

        contracted = trace(Gamma(which_gamma_five)*quark_1*Gamma(which_gamma_five)*anti_quark_prop);
       }

  }


  void I_one_Iz_pm_one_contract(LatticePropagator & quark_1,
                             LatticePropagator & quark_2,
                             LatticeComplex & contracted)
  {
    //quark_1 is the forward prop, quark_2 is the antiquark.

       int which_gamma_five = Ns*Ns-1;


        LatticePropagator anti_quark_prop = adj(Gamma(which_gamma_five)*quark_2*Gamma(which_gamma_five));
        contracted = trace(Gamma(which_gamma_five)*quark_1*Gamma(which_gamma_five)*anti_quark_prop);

  }


  void write_correlator(bool full_correlator,
			std::string meson_name,
#ifdef BUILD_HDF5
			std::string path,
			HDF5Writer & h5writer,
			HDF5Base::writemode & h5mode,
#endif
			int t_0,
			int Nt,
			multi1d<int> & source_coords,
			LalibeSftMom & FT,
			LatticeComplex & meson)
  {
    if(full_correlator == true)
    {
#ifdef BUILD_HDF5
      std::string correlator_path = path+"/"+meson_name+"/4D_correlator";
      h5writer.push(correlator_path);
      correlator_path = correlator_path+"/x"+std::to_string(source_coords[0])+"_y"+std::to_string(source_coords[1])+"_z"+std::to_string(source_coords[2])+"_t"+std::to_string(source_coords[3]);
      h5writer.write(correlator_path, meson, h5mode);
      h5writer.writeAttribute(correlator_path, "is_shifted", 0, h5mode);
      h5writer.cd("/");
#endif
    }
    else
    {
      //Momentum loop
      multi2d<Complex> FTed_meson = FT.sft(meson);
      //Temp variable for writing below.
      Complex temp_element;
      //Move the h5 pushing here, since all momentum keys will be written in the same general path.
#ifdef BUILD_HDF5
      std::string correlator_path = path+"/"+meson_name+"/x"+std::to_string(source_coords[0])+"_y"+std::to_string(source_coords[1])+"_z"+std::to_string(source_coords[2])+"_t"+std::to_string(source_coords[3]);
      h5writer.push(correlator_path);
#else
      std::string correlator_path = meson_name+"_x"+std::to_string(source_coords[0])+"_y"+std::to_string(source_coords[1])+"_z"+std::to_string(source_coords[2])+"_t"+std::to_string(source_coords[3]);
#endif
      for(int mom = 0; mom < FT.numMom(); mom++)
      {
	//One more temp variable instanited inside loop (once again for writing.)
	multi1d<Complex> meson_correlator;
	meson_correlator.resize(Nt);
	multi1d<int> momenta = FT.numToMom(mom);
#ifndef BUILD_HDF5
	std::string correlator_path_mom = correlator_path+"_px"+std::to_string(momenta[0])+"_py"+std::to_string(momenta[1])+"_pz"+std::to_string(momenta[2]);
	TextFileWriter file_out(correlator_path_mom);
#endif
	for(int t = 0; t < Nt; t++)
	{
	  temp_element = FTed_meson[mom][t];
	  int t_relative = t - t_0;
	  if(t_relative < 0)
	    t_relative += Nt;

#ifndef BUILD_HDF5
	  file_out<<temp_element<<"\n";
#endif
	  meson_correlator[t_relative] = temp_element;
	}
#ifndef BUILD_HDF5
	file_out.close();
#else
	//Change the name of string compred to 4d output so general correlator path is the same.
	std::string correlator_path_mom = correlator_path+"/px"+std::to_string(momenta[0])+"_py"+std::to_string(momenta[1])+"_pz"+std::to_string(momenta[2]);
	h5writer.write(correlator_path_mom, meson_correlator, h5mode);
	h5writer.writeAttribute(correlator_path_mom, "is_shifted", 1, h5mode);
	h5writer.cd("/");
#endif
      }
    }
  }



} // End namespace Chroma
