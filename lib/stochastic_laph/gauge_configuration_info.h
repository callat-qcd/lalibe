//Arjun Singh Gambhir
//Real authors: Ben Hoerz, John Bulava, Colin Morningstar, ??
//Code directly ported from chroma_laph, more acknowledgements coming later.
#ifndef LAPH_GAUGE_CONFIGURATION_INFO_H
#define LAPH_GAUGE_CONFIGURATION_INFO_H

#include "chromabase.h"
#include "xml_help.h"

namespace Chroma {
  namespace LaphEnv {

// ********************************************************************
// *                                                                  *
// *  Class "GaugeConfigurationInfo" holds information about a        *
// *  gauge configuration.  It also checks that its config info       *
// *  is the same when compared against another object of class       *
// *  GaugeConfigurationInfo. This will also return the trajectory    *
// *  number as well as some other info about the gauge               *
// *  configuration.                                                  *
// *                                                                  *
// *  The constructor takes either an XmlReader or a string.  The     *
// *  XmlReader version is called as follows:                         *
// *                                                                  *
// *    XmlReader xml_in(...);                                        *
// *    GaugeConfigurationInfo U(xml_in);                             *
// *                                                                  *
// *  There are three ways of constructing a GaugeConfigurationInfo   *
// *  depending on whether a "<gauge_id>" or a "<getFromFile>" is     *
// *  given as a tag in the XmlReader:                                *
// *                                                                  *
// *    (1) If "<gauge_id>" is given as a tag:                        *
// *                                                                  *
// *      --> This constructor expects XML of the form                *
// *               <GaugeConfigurationInfo>                           *
// *                  <gauge_id>....</gauge_id>                       *
// *               </GaugeConfigurationInfo>                          *
// *          then information about the configuration is obtained    *
// *          from TheNamedObjMap.                                    *
// *                                                                  *
// *    (2) If "<getFromFile>" is given as a tag:                     *
// *                                                                  *
// *      --> This constructor expects XML of the form                *
// *               <GaugeConfigurationInfo>                           *
// *                  <getFromFile> gaugefilename </getFromFile>      *
// *               </GaugeConfigurationInfo>                          *
// *          then information about the configuration is obtained    *
// *          from xml header in the gauge configuration file named   *
// *          in the tag.                                             *
// *                                                                  *
// *    (3) If neither "<gauge_id>" nor "<getFromFile>" is given as   *
// *        a tag:                                                    *
// *                                                                  *
// *      --> This version of the constructor expects all of the      *
// *          information to be given in the XML itself.              *
// *            <GaugeConfigurationInfo>                              *
// *               <HMCTrajectoryNumber>...<HMCTrajectoryNumber>      *
// *               <TimeDir> .... <TimeDir>                           *
// *               <ConfigType> .... <ConfigType>                     *
// *               <FileName> .... <FileName>                         *
// *               <TimeExtent> .... <TimeExtent>                     *
// *               <NumberOfDir> .... <NumberOfDir>                   *
// *               <LatticeExtents> .... <LatticeExtents>             *
// *            </GaugeConfigurationInfo>                             *
// *                                                                  *
// *                                                                  *
// *    string headerInfo(...);                                       *
// *    GaugeConfigurationInfo U(headerInfo);                         *
// *                                                                  *
// *                                                                  *
// *  Example usage:                                                  *
// *                                                                  *
// *    GaugeConfigurationInfo u2(xml_in);                            *
// *    u.checkEqual(u2); --> checks that u2 and u have same          *
// *                           content; throws exception if not       *
// *                                                                  *
// *    string out = u.output();   // xml output                      *
// *    string out = u.output(2);  // indented xml output             *
// *    int j = u.getTrajNum();    // returns  number of RHMC         *
// *                               // trajectory for this config      *
// *    int nt = u.getTimeExtent(); // time extent of lattice         *
// *    int tdir = u.getTimeDir();  // index of time                  *
// *                                                                  *
// ********************************************************************


class GaugeConfigurationInfo
{

  std::string file_name;
  std::string config_type;
  int traj_num;
  int time_dir;
  int time_extent;
  multi1d<int> extents;
  int number_dir;

 public:

  GaugeConfigurationInfo(XmlReader& xml_rdr);
         
  GaugeConfigurationInfo(XmlReader& xml_rdr, std::string& gaugeFileXML);

  GaugeConfigurationInfo(const GaugeConfigurationInfo& rhs);
                 
  GaugeConfigurationInfo& operator=(const GaugeConfigurationInfo& rhs);

  ~GaugeConfigurationInfo(){}
  
  void checkEqual(const GaugeConfigurationInfo& rhs) const;

  bool operator==(const GaugeConfigurationInfo& rhs) const;




  std::string output(int indent = 0) const;

  void output(XmlWriter& xmlout) const;

  int getTrajNum() const { return traj_num; }

  std::string getFileName() const { return file_name; }

  std::string getGaugeConfigHeader() const { return output(0); }

  void getGaugeConfigHeader(XmlWriter& xmlout) const { output(xmlout); }

  int getTimeExtent() const { return time_extent; }

  int getTimeDir() const { return time_dir; }
  
  int getNumberOfDirections() const { return number_dir; }

  std::string getFileXML() const;

  const multi1d<int> getExtents() const { return extents; }


 private:

  void set_info(XmlReader& xmlg, std::string& gauge_xml,
                bool get_gauge_xml=false);
  void set_from_named_obj_map(const std::string& gauge_id,
                              std::string& gauge_xml);
  void set_from_file_name(const std::string& configfile,
                          std::string& gauge_xml);
  void set_from_file_xml(XmlReader& xmlg);
  void set_from_xml(XmlReader& xmlg);
  void checkDimensions();
  void readGaugeXML(const std::string& config_file_name,
                    std::string& gauge_xml) const;

};


// **********************************************************
  }
}
#endif
