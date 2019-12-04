//Arjun Singh Gambhir
//Real authors: Ben Hoerz, John Bulava, Colin Morningstar, ??
//Code directly ported from chroma_laph, more acknowledgements coming later.
#ifndef LAPH_GAUGE_CONFIGURATION_HANDLER_H
#define LAPH_GAUGE_CONFIGURATION_HANDLER_H

#include "meas/inline/io/named_objmap.h"
#include "gauge_configuration_info.h"
#include "chromabase.h"


namespace Chroma {
  namespace LaphEnv {


// *******************************************************************
// *                                                                 *
// *  Class "GaugeConfigurationHandler" manages information about    *
// *  and access to the gauge configuration in the Laph environment. *
// *  It contains a "GaugeConfigurationInfo" that holds the gauge    *
// *  header, as well as a reference to the actual cfg.              *
// *                                                                 *
// *    Basic usage:                                                 *
// *                                                                 *
// *      (a) declare a GaugeConfigurationHandler                    *
// *      (b) set the info  -- via constructor or setInfo(..)        *
// *      (c) set the data  -- via setData(..)                       *
// *      (d) use getData(..) for access to gauge configuration      *
// *                                                                 *
// *     -- getInfo(..) provides access to configuration info        *
// *     -- can be cleared and info/data reset                       *
// *                                                                 *
// *    Example:                                                     *
// *       GaugeConfigurationInfo Gin;                               *
// *       GaugeConfigurationHandler uHandler;                       *
// *       uHandler.setInfo(Gin);                                    *
// *       uHandler.setData();                                       *
// *       const multi1d<LatticeColorMatrix>& U=uHandler.getData();  *
// *    --> access to config through U                               *
// *                                                                 *
// *******************************************************************


class GaugeConfigurationHandler
{

      // pointer to the info about the gauge config (internally
      // managed by this handler)
      
    const GaugeConfigurationInfo* gauge_info;
    std::string objmap_gauge_id;

      // pointer to the gauge field (external: in NamedObjMap)
 
    const multi1d<LatticeColorMatrix>* cfg;
    
    
      // prevent copying
    GaugeConfigurationHandler(const GaugeConfigurationHandler& u);
    GaugeConfigurationHandler& operator=(const GaugeConfigurationHandler& u);


  public:

    GaugeConfigurationHandler();
    
    GaugeConfigurationHandler(const GaugeConfigurationInfo& gaugeinfo,
                              const std::string& gauge_id = "default_gauge_field" );

    void setInfo(const GaugeConfigurationInfo& gaugeinfo,
                 const std::string& gauge_id = "default_gauge_field" );

    ~GaugeConfigurationHandler();
    
    void clear();
    
    void setData();
    
    
      // access to the gauge configuration and its info

    const multi1d<LatticeColorMatrix>& getData();

    const GaugeConfigurationInfo& getGaugeConfigurationInfo() const;

    bool isInfoSet() const { return (gauge_info!=0);}

    bool isDataSet() const { return (cfg!=0);}

    std::string getGaugeId() const {return objmap_gauge_id;}

    std::string getFullRecordXML() const;


  private:

    void set_info(const GaugeConfigurationInfo& gaugeinfo,
                  const std::string& gauge_id);  
    void check_info_set(const std::string& name) const;
 
};

// *************************************************************************
  }
}
#endif
