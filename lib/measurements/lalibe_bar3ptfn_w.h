/*! 
 * Baryon c3pt calculator, first version based off chroma.
 * Uses lalibe's fourier transform routine.
 * Writes correlators in hdf5.
 * This is likely to change a lot in the coming weeks.
 * Arjun Gambhir
 */

#ifndef __lalibe_bar3ptfn_h__
#define __lalibe_bar3ptfn_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace LalibeBar3ptfnEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct LalibeBar3ptfnParams 
  {
    LalibeBar3ptfnParams();
    LalibeBar3ptfnParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long      frequency;

    struct Param_t
    {
      multi1d<std::string> currents;        //list of currents
      int              j_decay;
      bool output_full_correlator;          //If no momentum is specified, we output the full correlator.
      bool is_mom_max;                      //keeps track of which momentum mode we are using
      int p2_max;                           //max of momentum transfer squared, optional
      multi1d<multi1d<int>> mom_list;       //list of momenta insertions, optional
      multi2d<int> p_list;                  //momentum list the slow fourier transform needs
#ifdef BUILD_HDF5
      std::string file_name;
      std::string obj_path;
#endif
    } param;

    struct SeqProp_t
    {
      std::string      seqprop_id;
      int              gamma_insertion;    /*!< second gamma insertion */
    };

    struct NamedObject_t
    {
      std::string           gauge_id;
      std::string           prop_id;
      multi1d<SeqProp_t>    seqprops;
#ifndef BUILD_HDF5
      std::string           bar3ptfn_file;
#endif
    } named_obj;
  };


  //! Inline measurement of 3pt functions
  /*! \ingroup inlinehadron */
  class LalibeBar3ptfn : public AbsInlineMeasurement 
  {
  public:
    ~LalibeBar3ptfn() {}
    LalibeBar3ptfn(const LalibeBar3ptfnParams& p) : params(p) {}
    LalibeBar3ptfn(const LalibeBar3ptfn& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    LalibeBar3ptfnParams params;
  };

};

#endif
