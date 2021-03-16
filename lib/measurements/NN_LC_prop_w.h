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
        std::string contractions_filename;    //filename of hdf5 file that contains contraction terms
        int contractions_n_sq;                //FIXME What comment do I put for this?
        unsigned int fft_chunksize;           //originally the only parameter in FFTPar struct
        bool fft_tune;                        //tune the fft?
        //multi1d<int> boosts;                //boosts; this will only give one boost...not what we want
        multi1d< multi1d<int> > boosts;       //boosts
        std::string output_filename;          //output file
        int output_stripesize;                //output stripesize; default recommended
        bool dirac_basis;                     //specifies props in dirac basis, this is false by default
        multi1d<GroupXML_t> sink_xml;         //used for Thorsten's sink construction, TODO: replace this with Chroma's SINK_SMEAR
        //FIXME: Not sure if the way I read sink_xml will work for a multi1d. 
      } nnlcparam ;

      struct NamedObject_t
      {
        //Gauge field likely needed for latscat from this section.
        std::string                     gauge_id;
        std::string                     prop0_id;
        std::string                     prop1_id;
        //For now, this will only support 2 propagators, but in the future this could be extended.

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

    /*struct FFTPars{
      unsigned int chunksize_block;
    };*/
    //We'll read this in NN_LC_Prop instead having a single member struct.

    /*struct InverterPars{
      GroupXML_t invParamUpLo;
      GroupXML_t invParamUpHi;
      GroupXML_t invParamStrangeLo;
      GroupXML_t invParamStrangeHi;
      std::string file;
      unsigned int highPrecFrequency;
      bool do_high;
    };*/
    //I don't think we'll need invert params, because that's not done in this measurement.

    //Changed Ullong to ullong. 
    /*struct SourcePars{
      unsigned int nsources;
      ullong startseed;
      Real sigmasq;
      Real b;
      multi1d<int> boost;
      multi2d<int> displacements;
      std::string type;
    };*/
    //Only boost is needed for reading, the other parms will be instantiated as needed.

    /*struct ContractionPars{
      std::string filename;
      int nsqmax;
    };*/
    //Don't need this. This is read in NNLCProp. 

    /*struct Outputpars{
      int stripesize;
      std::string path;
    };*/
    //Read these in the main struct as well.

    //This struct was apparently never used???
    /*struct FermPars{
      std::string file;
      multi1d<int> boundary;
      GroupXML_t strange_pars;
      GroupXML_t up_pars;
    };*/
   
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