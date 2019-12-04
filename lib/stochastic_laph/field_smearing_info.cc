//Arjun Singh Gambhir
//Real authors: Ben Hoerz, John Bulava, Colin Morningstar, ??
//Code directly ported from chroma_laph, more acknowledgements coming later.
#include "field_smearing_info.h"
#include "xml_help.h"
using namespace std;

namespace Chroma {
  namespace LaphEnv {


// *************************************************************

   // XmlReader constructor

GluonSmearingInfo::GluonSmearingInfo(XmlReader& xml_in)
{
 xml_tag_assert(xml_in,"GluonStoutSmearingInfo","GluonSmearingInfo");
 XmlReader xmlr(xml_in, "./descendant-or-self::GluonStoutSmearingInfo");
 extract_info_from_reader(xmlr);
}


void GluonSmearingInfo::extract_info_from_reader(XmlReader& xml_in)
{
 xmlread(xml_in,"LinkIterations", linkIterations, "GluonSmearingInfo");
 xmlread(xml_in,"LinkStapleWeight", linkStapleWeight, "GluonSmearingInfo");
 if ((linkIterations<0)||(linkStapleWeight<0.0)){ 
    xml_cerr(xml_in,"invalid smearing scheme parameters in GluonSmearingInfo");
    xmlreadfail(xml_in,"GluonSmearingInfo");}
}


 // *************************************************************

    // copy constructor

GluonSmearingInfo::GluonSmearingInfo(const GluonSmearingInfo& in) 
            : linkIterations(in.linkIterations),
              linkStapleWeight(in.linkStapleWeight) {}

GluonSmearingInfo& GluonSmearingInfo::operator=(const GluonSmearingInfo& in)
{
 linkIterations=in.linkIterations;
 linkStapleWeight=in.linkStapleWeight;
 return *this;
}

void GluonSmearingInfo::checkEqual(const GluonSmearingInfo& in) const
{
 if  ((linkIterations!=in.linkIterations)
    ||(abs(linkStapleWeight-in.linkStapleWeight)>1e-12)){
    std::cerr << "GluonSmearingInfo checkEqual failed"<<std::endl;
    std::cerr << "LHS:"<<std::endl<<output()<<std::endl<<"RHS:"<<std::endl<<in.output()<<std::endl;
    throw string("GluonSmearingInfo checkEqual failed...");}
}

bool GluonSmearingInfo::operator==(const GluonSmearingInfo& in) const
{
 return ((linkIterations==in.linkIterations)
       &&(abs(linkStapleWeight-in.linkStapleWeight)<1e-12));
}


string GluonSmearingInfo::output(int indent) const
{
 string pad(3*indent,' ');
 ostringstream oss;
 oss << pad << "<GluonStoutSmearingInfo>"<<endl;
 oss << pad << "   <LinkIterations>" << linkIterations 
     << "</LinkIterations>"<<endl;
 oss << pad << "   <LinkStapleWeight>" << linkStapleWeight 
     << "</LinkStapleWeight>"<<endl;
 oss << pad << "</GluonStoutSmearingInfo>"<<endl;
 return oss.str();
}

void GluonSmearingInfo::output(XmlWriter& xmlout) const
{
 push(xmlout,"GluonStoutSmearingInfo");
 write(xmlout,"LinkIterations",linkIterations);
 write(xmlout,"LinkStapleWeight",linkStapleWeight);
 pop(xmlout);
}


// *************************************************************



   // XmlReader constructor

QuarkSmearingInfo::QuarkSmearingInfo(XmlReader& xml_in)
{
 xml_tag_assert(xml_in,"QuarkLaphSmearingInfo","QuarkSmearingInfo");
 XmlReader xmlr(xml_in, "./descendant-or-self::QuarkLaphSmearingInfo");
 extract_info_from_reader(xmlr);
}


void QuarkSmearingInfo::extract_info_from_reader(XmlReader& xml_in)
{
 xmlread(xml_in,"LaphSigmaCutoff", laphSigma, "QuarkSmearingInfo");
 xmlread(xml_in,"NumberLaphEigvecs", laphNumEigvecs, "QuarkSmearingInfo");
 if ((laphNumEigvecs<1)||(laphSigma<=0)){
    xml_cerr(xml_in,"invalid smearing scheme parameters in QuarkSmearingInfo");
    throw(string("error"));}
}


  // ************************************************************

    // copy constructor

QuarkSmearingInfo::QuarkSmearingInfo(const QuarkSmearingInfo& in) 
            : laphNumEigvecs(in.laphNumEigvecs),
              laphSigma(in.laphSigma) {}

QuarkSmearingInfo& QuarkSmearingInfo::operator=(const QuarkSmearingInfo& in)
{
 laphNumEigvecs=in.laphNumEigvecs;
 laphSigma=in.laphSigma;
 return *this;
}
  
void QuarkSmearingInfo::increaseUpdate(const QuarkSmearingInfo& in)
{
 if (in.laphNumEigvecs>laphNumEigvecs){
    laphNumEigvecs=in.laphNumEigvecs;
    laphSigma=in.laphSigma;}
}


void QuarkSmearingInfo::checkEqual(const QuarkSmearingInfo& in) const
{
 if  ((abs(laphSigma-in.laphSigma)>1e-12)
    ||(laphNumEigvecs!=in.laphNumEigvecs)){
    std::cerr << "QuarkSmearingInfo checkEqual failed"<<endl;
    std::cerr << "LHS:"<<endl<<output()<<endl<<"RHS:"<<endl<<in.output()<<endl;
    throw string("QuarkSmearingInfo checkEqual failed...");}
}

bool QuarkSmearingInfo::operator==(const QuarkSmearingInfo& in) const
{
 return ((abs(laphSigma-in.laphSigma)<1e-12)
       &&(laphNumEigvecs==in.laphNumEigvecs));
}


void QuarkSmearingInfo::checkOK(const QuarkSmearingInfo& in) const
{
 if  (laphNumEigvecs>in.laphNumEigvecs){
    std::cerr << "QuarkSmearingInfo checkOK failed"<<endl;
    std::cerr << "LHS:"<<endl<<output()<<endl<<"RHS:"<<endl<<in.output()<<endl;
    throw string("QuarkSmearingInfo checkOK failed...");}
}

string QuarkSmearingInfo::output(int indent) const
{
 string pad(3*indent,' ');
 ostringstream oss;
 oss << pad << "<QuarkLaphSmearingInfo>"<<endl;
 oss << pad << "   <LaphSigmaCutoff>" << laphSigma 
     << "</LaphSigmaCutoff>"<<endl;
 oss << pad << "   <NumberLaphEigvecs>" << laphNumEigvecs 
     << "</NumberLaphEigvecs>"<<endl;
 oss << pad << "</QuarkLaphSmearingInfo>"<<endl;
 return oss.str();
}

void QuarkSmearingInfo::output(XmlWriter& xmlout) const
{
 push(xmlout,"QuarkLaphSmearingInfo");
 write(xmlout,"LaphSigmaCutoff",laphSigma);
 write(xmlout,"NumberLaphEigvecs",laphNumEigvecs);
 pop(xmlout);
}


// *************************************************************
  }
}
