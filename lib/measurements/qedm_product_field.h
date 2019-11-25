#ifndef __inline_qedm_product_field_h__
#define __inline_qedm_product_field_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace LalibeQEDMProductFieldEnv
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineProductFieldParams 
  {
    InlineProductFieldParams();
    InlineProductFieldParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {
      std::string quark_type;
      std::string ext_field_filename;
      Real coupling;
      //bool debug;
      //std::string debug_outfile;
    } param;

    struct NamedObject_t
    {
      //std::string external_id;
      std::string gauge_in_id;
      std::string gauge_out_id;
    } named_obj;

    //std::string xml_file;  // Alternate XML file pattern
  };


  //! Inline measurement of Wilson loops
  /*! \ingroup inlinehadron */
  class InlineProductField : public AbsInlineMeasurement 
  {
  public:
    ~InlineProductField() {}
    InlineProductField(const InlineProductFieldParams& p) : params(p) {}
    InlineProductField(const InlineProductField& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no, XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no, XMLWriter& xml_out); 
    //! Read in external field
    multi1d<LatticeReal> importExternalField(const std::string& ext_field_filename);

  private:
    InlineProductFieldParams params;
  };

};

#endif
