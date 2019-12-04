//Arjun Singh Gambhir
//Real authors: Ben Hoerz, John Bulava, Colin Morningstar, ??
//Code directly ported from chroma_laph, more acknowledgements coming later.
#include "laph_eigen_info.h"
#include "xml_help.h"
using namespace std;

namespace Chroma {
  namespace LaphEnv {


// *************************************************************

   // XmlReader constructor

LaphEigenSolverInfo::LaphEigenSolverInfo(XmlReader& xml_in)
{
 xml_tag_assert(xml_in,"LaphEigenSolverInfo");
 XmlReader xmlr(xml_in, "./descendant-or-self::LaphEigenSolverInfo");
 extract_info_from_reader(xmlr);
}


void LaphEigenSolverInfo::extract_info_from_reader(XmlReader& xmlr)
{
 xmlread(xmlr,"ResidualTolerance", tolerance, "LaphEigenSolverInfo");
 xmlread(xmlr,"MaxIterations", maxIterations, "LaphEigenSolverInfo");
 xmlread(xmlr,"KrylovDimension", dimKrylov, "LaphEigenSolverInfo");

 chebyshevOrder=1;
 if (xml_tag_count(xmlr,"ChebyshevOrder")==1)
    xmlread(xmlr,"ChebyshevOrder", chebyshevOrder, "LaphEigenSolverInfo");
 if (chebyshevOrder<1) chebyshevOrder=1;

 maxEigenvalue=15.0;
 if (xml_tag_count(xmlr,"MaxEigenvalue")==1){
    xmlread(xmlr,"MaxEigenvalue",maxEigenvalue,"LaphEigenSolverInfo");
    if (maxEigenvalue<12.0){
       xml_cerr(xmlr,"invalid MaxEigenvalue in LaphEigenSolverInfo: should exceed 12");
    xmlreadfail(xmlr,"LaphEigenSolverInfo");}}
 else{
    if (chebyshevOrder>1){
       xml_cerr(xmlr,"when Chebyshev acceleration used, must");
       xml_cerr(xmlr," supply a maximum eigenvalue");
       QDP_abort(1);}}

 cutoffEigenvalue=0.5;
 if (xml_tag_count(xmlr,"CutoffEigenvalue")==1){
    xmlread(xmlr,"CutoffEigenvalue",cutoffEigenvalue,"LaphEigenSolverInfo");
    if ((cutoffEigenvalue<0.03)||(cutoffEigenvalue>(maxEigenvalue-0.9))){
       xml_cerr(xmlr,"invalid CutoffEigenvalue in LaphEigenSolverInfo");
    xmlreadfail(xmlr,"LaphEigenSolverInfo");}}
 else{
    if (chebyshevOrder>1){
       xml_cerr(xmlr,"when Chebyshev acceleration used, must");
       xml_cerr(xmlr," supply a cutoff eigenvalue");
       QDP_abort(1);}}

 startVector="equal_components";
 if (xml_tag_count(xmlr,"StartingVectorType")==1){
    string svread;
    xmlread(xmlr,"StartingVectorType",svread,"LaphEigenSolverInfo");
    svread=tidyString(svread);
    if (svread=="random") startVector="random";
    else if (svread=="equal_components") startVector="equal_components";
    else{
       xml_cerr(xmlr,"Invalid <StartingVectorType> tag in LaphEigenSolverInfo");
       throw(string("error"));}}

 if (xml_tag_count(xmlr,"OutputVerbosity")==1){
    xmlread(xmlr,"OutputVerbosity",outputVerbosity,"LaphEigenSolverInfo");
    if (outputVerbosity<0) outputVerbosity=0;
    if (outputVerbosity>2) outputVerbosity=2;}
 else
    outputVerbosity=0;

 if ((maxIterations<4)||(dimKrylov<6)||(tolerance<=0)){
    xml_cerr(xmlr,"invalid input parameters in LaphEigenSolverInfo");
    xmlreadfail(xmlr,"LaphEigenSolverInfo");}

}


 // *************************************************************

    // copy constructor

LaphEigenSolverInfo::LaphEigenSolverInfo(const LaphEigenSolverInfo& in) 
            :  maxIterations(in.maxIterations),
               dimKrylov(in.dimKrylov),
               tolerance(in.tolerance),
               chebyshevOrder(in.chebyshevOrder),
               maxEigenvalue(in.maxEigenvalue),
               cutoffEigenvalue(in.cutoffEigenvalue),
               startVector(in.startVector),
               outputVerbosity(in.outputVerbosity) {}

LaphEigenSolverInfo& LaphEigenSolverInfo::operator=(
               const LaphEigenSolverInfo& in)
{
 maxIterations=in.maxIterations;
 startVector=in.startVector;
 dimKrylov=in.dimKrylov;
 tolerance=in.tolerance;
 chebyshevOrder=in.chebyshevOrder;
 maxEigenvalue=in.maxEigenvalue;
 cutoffEigenvalue=in.cutoffEigenvalue;
 outputVerbosity=in.outputVerbosity;
 return *this;
}


string LaphEigenSolverInfo::output(int indent) const
{
 string pad(3*indent,' ');
 ostringstream oss;
 oss << pad << "<LaphEigenSolverInfo>"<<endl;
 oss << pad << "   <ResidualTolerance>" << tolerance 
            << "</ResidualTolerance>" << endl;
 oss << pad << "   <MaxIterations>" << maxIterations 
            << "</MaxIterations>" << endl;
 oss << pad << "   <KrylovDimension>" << dimKrylov 
            << "</KrylovDimension>" <<endl;
 if (chebyshevOrder>1){
   oss << pad << "   <ChebyshevOrder>" << chebyshevOrder 
              << "</ChebyshevOrder>" <<endl;
   oss << pad << "   <MaxEigenvalue>" << maxEigenvalue 
              << "</MaxEigenvalue>"<<endl;
   oss << pad << "   <CutoffEigenvalue>" << cutoffEigenvalue 
              << "</CutoffEigenvalue>"<<endl;}
 oss << pad << "   <StartingVectorType>" << startVector 
            << "</StartingVectorType>" << endl;
 oss << pad << "   <OutputVerbosity>" << outputVerbosity
            << "</OutputVerbosity>" << endl;
 oss << pad << "</LaphEigenSolverInfo>"<<endl;
 return oss.str();
}

void LaphEigenSolverInfo::output(XmlWriter& xmlout) const
{
 push(xmlout,"LaphEigenSolverInfo");
 write(xmlout,"ResidualTolerance", tolerance);
 write(xmlout,"MaxIterations", maxIterations);
 write(xmlout,"KrylovDimension", dimKrylov);
 if (chebyshevOrder>1){
    write(xmlout,"ChebyshevOrder", chebyshevOrder);
    write(xmlout,"MaxEigenvalue",maxEigenvalue);
    write(xmlout,"CutoffEigenvalue",cutoffEigenvalue);}
 write(xmlout,"StartingVectorType", startVector);
 write(xmlout,"OutputVerbosity", outputVerbosity);
 pop(xmlout);
}


// *************************************************************
  }
}
