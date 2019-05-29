/*
 * Form factor code from chroma, modified to use lalibe's fourier transform routine.
 * Arjun Gambhir
 */

#ifndef __lalibe_formfac_h__
#define __lalibe_formfac_h__

#include "../momentum/lalibe_sftmom.h"

namespace Chroma 
{

 /*
  * Structures for hadron parts
  *
  * \ingroup hadron
  *
  * @{
  */
  struct LalibeFormFac_momenta_t
  {
    int               magic;     // magic number for sanity checks
    multi1d<int>      inser_mom;
    multi1d<ComplexF> local_current;
    multi1d<ComplexF> nonlocal_current;
  };

  struct LalibeFormFac_insertion_t
  {
    int              gamma_value;
    multi1d<LalibeFormFac_momenta_t> momenta;
  };

  struct LalibeFormFac_insertions_t
  {
    int  output_version;   // Unique id for each output version of the structures
    multi1d<LalibeFormFac_insertion_t>  formFac;
  };


  // Readers and writers
  void read(BinaryReader& bin, LalibeFormFac_momenta_t& mom);
  void read(BinaryReader& bin, LalibeFormFac_insertion_t& mes);
  void read(BinaryReader& bin, LalibeFormFac_insertions_t& form);
  void write(BinaryWriter& bin, const LalibeFormFac_momenta_t& mom);
  void write(BinaryWriter& bin, const LalibeFormFac_insertion_t& mes);
  void write(BinaryWriter& bin, const LalibeFormFac_insertions_t& form);


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


  LatticeComplex SpinColorArraySum(LatticePropagator& gamma_prop,
                         LatticePropagator& sequential_prop);

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
	       int t0);

}  // end namespace Chroma

#endif
