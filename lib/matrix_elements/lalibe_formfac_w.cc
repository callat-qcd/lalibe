/*
 * Form factor code from chroma, modified to use lalibe's fourier transform routine.
 * Arjun Gambhir
 */

//CHROMA STUFF
#include "chromabase.h"

//LALIBE STUFF
#include "../momentum/lalibe_sftmom.h"
#include "lalibe_formfac_w.h"
#include "bilinear_gamma.h"

namespace Chroma
{

  /*!
   * Structures for hadron parts
   *
   * \ingroup hadron
   *
   * @{
   */

  // Read a momenta struct
  void read(BinaryReader& bin, LalibeFormFac_momenta_t& mom)
  {
    read(bin, mom.magic);
    if (mom.magic != 20301)
    {
      QDPIO::cerr << "read(FormFac_momenta_t): magic number invalid" << std::endl;
      QDP_abort(1);
    }
    read(bin, mom.inser_mom);
    read(bin, mom.local_current);
    read(bin, mom.nonlocal_current);
  }

  //
  void read(BinaryReader& bin, LalibeFormFac_insertion_t& mes)
  {
    read(bin, mes.gamma_value);
    read(bin, mes.momenta);
  }

  //
  void read(BinaryReader& bin, LalibeFormFac_insertions_t& form)
  {
    read(bin, form.output_version);
    read(bin, form.formFac);
  }

  // Write a momenta struct
  void write(BinaryWriter& bin, const LalibeFormFac_momenta_t& mom)
  {
    int magic = 20301;
    write(bin, magic);
    write(bin, mom.inser_mom);
    write(bin, mom.local_current);
    write(bin, mom.nonlocal_current);
  }

  //
  void write(BinaryWriter& bin, const LalibeFormFac_insertion_t& mes)
  {
    write(bin, mes.gamma_value);
    write(bin, mes.momenta);
  }

  //
  void write(BinaryWriter& bin, const LalibeFormFac_insertions_t& form)
  {
    write(bin, form.output_version);
    write(bin, form.formFac);
  }

  /*! @} */  // end of group hadron

  /* Compute the spin and color array sum for the 3pt contractions.
   * \param correlator         final correlation function after spin color sum.( Write )
   * \param gamma_prop         propagator left multiplied by the matrix element insertion. ( Read )
   * \param seq_quark_prop     sequential quark propagator ( Read )
   *
   * This routine computes the correlation function:
   *
   * Sum_ab Sum_kl  Gamma_ac P_cb^kl(z,x) * Sigma_cb^kl(z,x)
   * where P is the forward propagator and Sigma is the sequential propagator.
  */

/*
  LatticeComplex SpinColorArraySum(LatticePropagator& gamma_prop,
                         LatticePropagator& sequential_prop)
   {
    LatticeComplex correlator = zero;
    LatticeColorMatrix tmp1;
    LatticeColorMatrix tmp2;
    for(int spin_ind1 = 0; spin_ind1 < Ns; spin_ind1++)
    {

       for(int spin_ind2 = 0; spin_ind2 < Ns; spin_ind2++)
       {
          for(int color_ind1 = 0; color_ind1 < Nc; color_ind1++)
          {

             for(int color_ind2 = 0; color_ind2 < Nc; color_ind2++)
             {
             tmp1 = peekSpin(gamma_prop,spin_ind1,spin_ind2);
       	     tmp2 = peekSpin(sequential_prop,spin_ind1,spin_ind2);

             correlator += peekColor(tmp1,color_ind1,color_ind2)*peekColor(tmp2,color_ind1,color_ind2);

             }//color_ind2


          }//color_ind1

        }//spin_ind2

    }//spin_ind1
   return correlator;
   } //function
*/

  //! Compute contractions for current insertion 3-point functions.
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions!
   *
   * \param form               structures holding formfactors ( Write )
   * \param u                  gauge fields (used for non-local currents) ( Read )
   * \param quark_propagator   quark propagator ( Read )
   * \param seq_quark_prop     sequential quark propagator ( Read )
   * \param gamma_insertion    extra gamma insertion at source ( Read )
   * \param phases             fourier transform phase factors ( Read )
   * \param t0                 cartesian coordinates of the source ( Read )
   */


  void FormFac(LalibeFormFac_insertions_t& form,
	       const multi1d<LatticeColorMatrix>& u,
	       const LatticePropagator& quark_propagator,
	       const LatticePropagator& seq_quark_prop,
	       int gamma_insertion,
	       const LalibeSftMom& phases,
	       bool full_correlator,
	       multi1d<int> & source_coords,
	       multi1d<std::string> bilinears,
#ifdef BUILD_HDF5
	       std::string path,
	       std::string particle,
	       HDF5Writer & h5writer,
	       HDF5Base::writemode & h5mode,
#endif
	       int t0)
  {
    START_CODE();

    // Length of lattice in j_decay direction and 3pt correlations fcns
    int length = phases.numSubsets();

    int G5 = Ns*Ns-1;

    // Construct the anti-quark propagator from the seq. quark prop.


    LatticePropagator anti_quark_prop = Gamma(G5) * seq_quark_prop * Gamma(G5);
    /*
        NOTE: The following '-' sign is added to account for another '-' we are
        not sure where it comes from - but checking against known results, we
        know we are off by an overall sign, so we add it here
    */
    anti_quark_prop = -anti_quark_prop;

    // Rough timings (arbitrary units):
    //   Variant 1: 120
    //   Variant 2: 140
    // See previous cvs versions (before 1.10) for Variant 2 - only keeping Variant 1

    //form.formFac.resize(Nd*Nd);
    //Now, this only goes up to the number of currents that are specified.
    form.formFac.resize(bilinears.size());

    // Loop over currents
    QDPIO::cout << "LALIBE_FORM_FACTOR: N_currents " << bilinears.size() << std::endl;
    //This preserves the older I/O, but will only index/write the currents that are requested.

    int gamma_value = 0;

    //These are needed in the loop, but never change.
    LatticePropagator gamma_propagator = zero;
    LatticePropagator quark = quark_propagator;

    multi1d<LatticeColorMatrix> gfield = u;

    for(int current_index = 0; current_index < bilinears.size(); current_index++)
    {
      //  For the case where the gamma value indicates we are evaluating either
      //  the std::vector or axial std::vector currents, we will also evaluate
      //  the non-local currents.  The non-local std::vector current is the conserved
      //  current.  The non-local axial std::vector current would be partially
      //  conserved but for the Wilson term.  In these cases we will set
      //  mu = corresponding direction.  In all other cases, we will set mu = -1.

      bool compute_nonlocal = false;
      int mu;

      std::string present_current = bilinears[current_index];

      QDPIO::cout << "LALIBE_FORM_FACTOR: current " << present_current << std::endl;

      if (present_current == "A1" || present_current == "V1")
      {
	mu = 0;
	// AWL - 2019-01-12
	// turning off all nonlocal currents - we don't use them
	//compute_nonlocal = true;
      }
      else if (present_current == "A2" || present_current == "V2")
      {
	mu = 1;
	//compute_nonlocal = true;
      }
      else if (present_current == "A3" || present_current == "V3")
      {
	mu = 2;
	//compute_nonlocal = true;
      }
      else if (present_current == "A4" || present_current == "V4")
      {
	mu = 3;
	//compute_nonlocal = true;
      }
      else
      {
	mu = -1;
	compute_nonlocal = false;
      }

      // The local non-conserved std::vector-current matrix element
      //Use lalibe functions for gamma insertions.
      Bilinear_Gamma(present_current, gamma_propagator, quark, u);

      LatticeComplex local_current = trace(adj(anti_quark_prop) * gamma_propagator * Gamma(gamma_insertion));

       multi2d<DComplex> hsum, hsum_nonlocal;

      //We only do momentum injection if we are not dumping out the full, 4-d correlator.
      if(full_correlator == true)
      {
#ifdef BUILD_HDF5
	std::string correlator_path = path+particle+"/"+present_current+"/x"+std::to_string(source_coords[0])+"_y"+std::to_string(source_coords[1])+"_z"+std::to_string(source_coords[2])+"_t"+std::to_string(source_coords[3])+"/4D_correlator";
	h5writer.push(correlator_path);
	correlator_path = correlator_path+"/local_current";
	h5writer.write(correlator_path, local_current, h5mode);
	h5writer.writeAttribute(correlator_path, "is_shifted", 0, h5mode);
	h5writer.cd("/");
#endif
      }
      else
	hsum = phases.sft(local_current);


      // Construct the non-local current matrix element
      //
      // The form of J_mu = (1/2)*[psibar(x+mu)*U^dag_mu*(1+gamma_mu)*psi(x) -
      //                           psibar(x)*U_mu*(1-gamma_mu)*psi(x+mu)]
      // NOTE: the 1/2  is included down below in the sumMulti stuff
      LatticeComplex non_local_current;


      if(compute_nonlocal){
	Bilinear_Gamma(present_current, gamma_propagator, quark, u);
/*
        LatticePropagator tmp_prop1 = adj(gfield[mu])*(quark_propagator + gamma_propagator);
        LatticePropagator tmp_prop2 = shift(seq_prop,FORWARD,mu);

        non_local_current = 0.5*SpinColorArraySum(tmp_prop1,tmp_prop2);

        tmp_prop1 = shift(quark_propagator,FORWARD,mu);

        tmp_prop1 -= shift(gamma_propagator,FORWARD,mu);

        tmp_prop1 = gfield[mu]*tmp_prop1;
        non_local_current -= 0.5*SpinColorArraySum(tmp_prop1,seq_prop);
*/
	non_local_current = trace(adj(u[mu] * shift(anti_quark_prop, FORWARD, mu)) *
		(quark_propagator + gamma_propagator) * Gamma(gamma_insertion));
	LatticePropagator tmp_prop1 = u[mu] *
	  shift(quark_propagator, FORWARD, mu);
	Bilinear_Gamma(present_current, gamma_propagator, tmp_prop1, u);
	non_local_current -= trace(adj(anti_quark_prop) *
				  (tmp_prop1 - gamma_propagator) * Gamma(gamma_insertion));

        non_local_current = 0.5*non_local_current;

	if(full_correlator == true)
	{
#ifdef BUILD_HDF5
	  std::string correlator_path = path+particle+"/"+present_current+"/x"+std::to_string(source_coords[0])+"_y"+std::to_string(source_coords[1])+"_z"+std::to_string(source_coords[2])+"_t"+std::to_string(source_coords[3])+"/4D_correlator";
	  h5writer.push(correlator_path);
	  correlator_path = correlator_path+"/non_local_current";
	  h5writer.write(correlator_path, non_local_current, h5mode);
	  h5writer.writeAttribute(correlator_path, "is_shifted", 0, h5mode);
	  h5writer.cd("/");
#endif
	}
	else
	  hsum_nonlocal = phases.sft(non_local_current);
      }


      form.formFac[gamma_value].gamma_value = gamma_value;
      if(full_correlator == false)
      {
	form.formFac[gamma_value].momenta.resize(phases.numMom());  // hold momenta output

	// Loop over insertion momenta and print out results
	for(int inser_mom_num=0; inser_mom_num<phases.numMom(); ++inser_mom_num)
	{
	  form.formFac[gamma_value].momenta[inser_mom_num].inser_mom = phases.numToMom(inser_mom_num);

	  multi1d<ComplexF> local_cur3ptfn(length); // always compute
          multi1d<ComplexF> nonlocal_cur3ptfn;

          if(compute_nonlocal)
	    nonlocal_cur3ptfn.resize(length);      // possibly compute

	  for (int t=0; t < length; ++t)
	  {
	    int t_eff = (t - t0 + length) % length;

	    local_cur3ptfn[t_eff] = Complex(hsum[inser_mom_num][t]);
	    if (compute_nonlocal)
        {
	      nonlocal_cur3ptfn[t_eff] = Complex(hsum_nonlocal[inser_mom_num][t]);
	    }

	  } // end for(t)

	  //Write FTed result in hdf5 if it's turned on.
	  multi1d<int> momenta = phases.numToMom(inser_mom_num);
#ifdef BUILD_HDF5
	  std::string correlator_path = path+particle+"/"+present_current+"/x"+std::to_string(source_coords[0])+"_y"+std::to_string(source_coords[1])+"_z"+std::to_string(source_coords[2])+"_t"+std::to_string(source_coords[3])+"/px"+std::to_string(momenta[0])+"_py"+std::to_string(momenta[1])+"_pz"+std::to_string(momenta[2]);
	  h5writer.push(correlator_path);
	  std::string correlator_path_c3pt = correlator_path+"/local_current";;
	  h5writer.write(correlator_path_c3pt, local_cur3ptfn, h5mode);
	  h5writer.writeAttribute(correlator_path_c3pt, "is_shifted", 1, h5mode);
	  h5writer.cd("/");
#endif
	  if(compute_nonlocal)
	  {
#ifdef BUILD_HDF5
	    std::string correlator_path = path+particle+"/"+present_current+"/x"+std::to_string(source_coords[0])+"_y"+std::to_string(source_coords[1])+"_z"+std::to_string(source_coords[2])+"_t"+std::to_string(source_coords[3])+"/px"+std::to_string(momenta[0])+"_py"+std::to_string(momenta[1])+"_pz"+std::to_string(momenta[2]);
	    h5writer.push(correlator_path);
	    std::string correlator_path_c3pt = correlator_path+"/non_local_current";;
	    h5writer.write(correlator_path_c3pt, nonlocal_cur3ptfn, h5mode);
	    h5writer.writeAttribute(correlator_path_c3pt, "is_shifted", 1, h5mode);
	    h5writer.cd("/");
#endif
	  }
	  form.formFac[gamma_value].momenta[inser_mom_num].local_current    = local_cur3ptfn;
	  form.formFac[gamma_value].momenta[inser_mom_num].nonlocal_current = nonlocal_cur3ptfn;

	} // end for(inser_mom_num)
	gamma_value++;
      } // end for(gamma_value)
    }// end if stattement to do mom
    END_CODE();
  }

}  // end namespace Chroma
