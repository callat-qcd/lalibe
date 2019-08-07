/*
Authors
Arjun Gambhir

First pass at porting the linear combo executable of latscat over to lalibe.
After initial prototype working, lots of optimizations are likely to follow.
Memory management when boosting needs to be looked into.
*/

#ifndef __lalibe_NN_LC_prop_w_h__
#define __lalibe_NN_LC_prop_w_h__

//Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace LalibeNucleonNucleonLinearComboPropagatorEnv
  {
    extern const std::string name;
    bool registerAll();

    struct NNLCPropParams 
    {
      NNLCPropParams();
      NNLCPropParams(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long      frequency;

      // only same flavor bilinear operators are used
      struct NNLCProp_t{
        multi1d<std::string> currents;        //list of currents
        //int  t0;                            //t0 of input prop
        //int j_decay;                        //orthogonal direction of FT
        //Apparently I can read these things from a src prop.
        bool is_mom_max;                      //keeps track of which momentum mode we are using
        int p2_max;                           //max of momentum transfer squared, optional
        multi1d<multi1d<int>> mom_list;       //list of momenta insertions, optional
        multi2d<int> p_list;                  //momentum list the slow fourier transform needs
        ChromaProp_t prop_param;              //params for next lin solve
        std::string contractions_filename;    //filename of hdf5 file that contains contraction terms
        int contractions_n_sq;                //FIXME What comment do I put for this?
      } nnlcparam ;

      struct NamedObject_t
      {
        std::string                     gauge_id;
        std::string                     src_prop_id;
        multi1d<std::string>            fh_prop_id;

      } named_obj;
    };
   
    //Stuff below is ported from latscat.
    //LatticePars is read and set up by lalibe's main function.
    /*struct LatticePars{
      multi1d<int> nrow;
      int tDir;
      int tLength;
      bool gfix;
    };*/

    struct FFTPars{
      unsigned int chunksize_block;
    };

    struct InverterPars{
      GroupXML_t invParamUpLo;
      GroupXML_t invParamUpHi;
      GroupXML_t invParamStrangeLo;
      GroupXML_t invParamStrangeHi;
      std::string file;
      unsigned int highPrecFrequency;
      bool do_high;
    };

    //Changed Ullong to ullong. 
    struct SourcePars{
      unsigned int nsources;
      ullong startseed;
      Real sigmasq;
      Real b;
      multi1d<int> boost;
      multi2d<int> displacements;
      std::string type;
    };

    struct ContractionPars{
      std::string filename;
      int nsqmax;
    };

    struct Outputpars{
      int stripesize;
      std::string path;
    };

    struct FermPars{
      std::string file;
      multi1d<int> boundary;
      GroupXML_t strange_pars;
      GroupXML_t up_pars;
    };
   
    //Linkage chroma hack was here in the original, not needed now.
    
    const std::string dirlist[4]={"x","y","z","t"};
    const std::string contterms[12]={"PP_SING0_loc","PP_SING0","PP_TRIPP","PP_TRIP0","PP_TRIPM","PN_TRIPP_loc","PN_TRIP0_loc","PN_TRIPM_loc","PN_SING0","PN_TRIPP","PN_TRIP0","PN_TRIPM"};

    const Real gausspar=1.4;
    const Real oneoversqrt6=0.4082482904638631;
    const Real oneoversqrt2=0.7071067811865475;

    //Read and write propagator functions from chroma were here, this is already built in now.

    //! Compute a spacetime dependent sequential source
    /*! \ingroup inlinehadron */
    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const NNLCPropParams& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the measurement
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    private:
      NNLCPropParams params;
    };

  }

}

#endif
