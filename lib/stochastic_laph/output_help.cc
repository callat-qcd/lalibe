//Arjun Singh Gambhir
//Real authors: Ben Hoerz, John Bulava, Colin Morningstar, ??
//Code directly ported from chroma_laph, more acknowledgements coming later.
#include "output_help.h"
using namespace std;

namespace Chroma {
  namespace LaphEnv {


// *********************************************************

void output(const Complex& data, TextFileWriter& tout)
{
 tout <<"("<<data.elem().elem().elem().real()<<", "
           <<data.elem().elem().elem().imag()<<")";
}

void output(const multi1d<Complex>& data, TextFileWriter& tout)
{
 int n=data.size();
 for (int k=0;k<n;k++){
    tout << "elem["<<k<<"] = ";
    output(data[k],tout);
    tout << "\n";}
}


void output(const multi2d<Complex>& data, TextFileWriter& tout)
{
 int n1=data.size1();
 int n2=data.size2();
 for (int k1=0;k1<n1;k1++)
 for (int k2=0;k2<n2;k2++){
    tout << "elem["<<k2<<"]["<<k1<<"] = ";
    output(data(k2,k1),tout);
    tout << "\n";}
}


void output(const multi3d<Complex>& data, TextFileWriter& tout)
{
 int n1=data.size1();
 int n2=data.size2();
 int n3=data.size3();
 for (int k1=0;k1<n1;k1++)
 for (int k2=0;k2<n2;k2++)
 for (int k3=0;k3<n3;k3++){
    tout << "elem["<<k3<<"]["<<k2<<"]["<<k1<<"] = ";
    output(data(k3,k2,k1),tout);
    tout << "\n";}
}


// **********************************************************
  }
}
