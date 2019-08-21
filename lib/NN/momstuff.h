//Arjun Singh Gambhir
//Code ported from latscat.
#ifndef _MOMSTUFF
#define _MOMSTUFF
#include "chromabase.h"
//We include chroma types this way.

//Better use chroma namespace for this stuff.

namespace Chroma {

	struct momentum{
	private:
	  int comp[3];
	public:
	  momentum(){}
	  momentum(const int& x, const int& y, const int& z){
		  comp[0]=x;
		  comp[1]=y;
		  comp[2]=z;
	  }
	  
	  momentum(const multi1d<int>& mom){
		  if(mom.size()==3){
			  comp[0]=mom[0];
			  comp[1]=mom[1];
			  comp[2]=mom[2];
		  }
		  else{
			  comp[0]=0;
			  comp[1]=0;
			  comp[2]=0;
		  }
	  }
	  
	  int operator[](const int index) const;
	  int& operator[](const int index);
	  int magsq() const;
	};

	//io-utils
	void read(XMLReader&,const std::string&,momentum&);
	std::ostream& operator<<(std::ostream &os,const momentum &obj);
	StandardOutputStream& operator<<(StandardOutputStream &os,const momentum &obj);

	//projection routines:
	LatticeComplex get_phases(momentum mom, int time_dir, Real fact=1.);
	LatticeComplex get_sin_phases(momentum mom, int time_dir, Real fact=1.);
	LatticeComplex get_cos_phases(momentum mom, int time_dir, Real fact=1.);

	//other fields, relevant for convoluting into propagators:
	LatticeComplex getPhasesA1(int displacement, int time_dir, Real fact);
	LatticeComplex getPhasesT1d1(int component, int displacement, int time_dir, Real fact);

}

#endif
