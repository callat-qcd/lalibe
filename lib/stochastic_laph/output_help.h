//Arjun Singh Gambhir
//Real authors: Ben Hoerz, John Bulava, Colin Morningstar, ??
//Code directly ported from chroma_laph, more acknowledgements coming later.
#ifndef OUTPUT_HELP_H
#define OUTPUT_HELP_H

#include "qdp.h"
#include "chromabase.h"

namespace Chroma {
  namespace LaphEnv {

// *********************************************************

void output(const Complex& data, TextFileWriter& tout);
void output(const multi1d<Complex>& data, TextFileWriter& tout);
void output(const multi2d<Complex>& data, TextFileWriter& tout);
void output(const multi3d<Complex>& data, TextFileWriter& tout);

// *********************************************************
  }
}
#endif
