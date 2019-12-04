//Arjun Singh Gambhir
//Real authors: Ben Hoerz, John Bulava, Colin Morningstar, ??
//Code directly ported from chroma_laph, more acknowledgements coming later.
#ifndef FIELD_SMEARING_H
#define FIELD_SMEARING_H

#include "qdp.h"
#include "chromabase.h"
#include "xml_help.h"

namespace Chroma {
  namespace LaphEnv {


// *******************************************************************
// *                                                                 *
// *   Objects of class "GluonSmearingInfo" and "QuarkSmearingInfo"  *
// *   store identifying info about both the link and quark field    *
// *   smearing used.  Stout link smearing and quark Laph smearing   *
// *   are enforced.  The XML input must have the format             *
// *                                                                 *
// *       <GluonStoutSmearingInfo>                                  *
// *           <LinkIterations> 4 </LinkIterations>                  *
// *           <LinkStapleWeight>  0.25 </LinkStapleWeight>          *
// *       </GluonStoutSmearingInfo>                                 *
// *                                                                 *
// *       <QuarkLaphSmearingInfo>                                   *
// *           <LaphSigmaCutoff> 0.76 </LaphSigmaCutoff>             *
// *           <NumberLaphEigvecs> 32 </NumberLaphEigvecs>           *
// *       </QuarkLaphSmearingInfo>                                  *
// *                                                                 *
// *   Remember that the Laplacian Heaviside smearing (Laph) is      *
// *   defined by                                                    *
// *              theta( sigma_s^2 + Delta )                         *
// *                                                                 *
// *   where Delta = covariant Laplacian and sigma_s is specified    *
// *   in the <LaphSigmaCutoff> tag.                                 *
// *                                                                 *
// *   Example usage:                                                *
// *                                                                 *
// *     XmlReader xml_in(...);                                      *
// *     GluonSmearingInfo Gsmear(xml_in);                           *
// *     QuarkSmearingInfo Qsmear(xml_in);                           *
// *                                                                 *
// *     GluonSmearingInfo smear2(....);                             *
// *     Gsmear.checkEqual(smear2); // throws string exception       *
// *                                //   if smear2 != smear          *
// *     if (smear==smear2) ...     // returns boolean               *
// *                                                                 *
// *     int ival = Gsmear.getNumberOfLinkIterations();              *
// *     double dval = Gsmear.getLinkStapleWeight();                 *
// *     int jval = Qsmear.getNumberOfLaplacianEigenvectors();       *
// *     double dval = Qsmear.getLaphSigmaCutoff();                  *
// *     string sval = Gsmear.getSmearingType();                     *
// *     string ssval = Qsmear.getSmearingType();                    *
// *                                                                 *
// *     string out = smear.output();    // xml output               *
// *     string out = smear.output(2);   // indented xml output      *
// *                                                                 *
// *******************************************************************


class GluonSmearingInfo
{

  int linkIterations;
  double linkStapleWeight;

 public:  

  GluonSmearingInfo(XmlReader& xml_in);

  GluonSmearingInfo(const GluonSmearingInfo& in);

  GluonSmearingInfo& operator=(const GluonSmearingInfo& in);

  ~GluonSmearingInfo(){}

  void checkEqual(const GluonSmearingInfo& in) const;

  bool operator==(const GluonSmearingInfo& in) const;


    // output functions

  int getNumberOfLinkIterations() const { return linkIterations; }

  double getLinkStapleWeight() const { return linkStapleWeight; }

  std::string getSmearingType() const { return "STOUT"; }

  std::string output(int indent = 0) const;

  void output(XmlWriter& xmlout) const;

  std::string getHeader() const { return output(0); }

  void getHeader(XmlWriter& xmlout) const { return output(xmlout); }

 private:

  void extract_info_from_reader(XmlReader& xml_in);

};



class QuarkSmearingInfo
{

  int laphNumEigvecs;
  double laphSigma;

 public:  

  QuarkSmearingInfo(XmlReader& xml_in);

  QuarkSmearingInfo(const QuarkSmearingInfo& in);

  QuarkSmearingInfo& operator=(const QuarkSmearingInfo& in);

  void increaseUpdate(const QuarkSmearingInfo& in);

  ~QuarkSmearingInfo(){}

  void checkEqual(const QuarkSmearingInfo& in) const;

  bool operator==(const QuarkSmearingInfo& in) const;

    // same as checkEqual, but in.laphNumEigvecs >= laphNumEigvecs 
    // is needed only

  void checkOK(const QuarkSmearingInfo& in) const;


    // output functions

  int getNumberOfLaplacianEigenvectors() const {return laphNumEigvecs;}
  
  double getLaphSigmaCutoff() const { return laphSigma; }

  std::string getSmearingType() const { return "LAPH"; }

  std::string output(int indent = 0) const;

  void output(XmlWriter& xmlout) const;

  std::string getHeader() const { return output(0); }

  void getHeader(XmlWriter& xmlout) const { return output(xmlout); }

 private:

  void extract_info_from_reader(XmlReader& xml_in);

};


// **************************************************
  }
}
#endif
