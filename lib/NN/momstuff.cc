//Arjun Singh Gambhir
//Code ported from latscat.
//#include "progs.h"
//progs.h was included here, but we are trying to avoid that.

//Let's include the header.
#include "momstuff.h"

//this file contains utility routines for momentum sinks and sources:
//using namespace Chroma;
//Lalibe way.

namespace Chroma {

	//operators:
	int& momentum::operator[](const int index){
	  return comp[index%3];
	}

	int momentum::operator[](const int index) const{
	  return comp[index%3];
	}

	//members:
	int momentum::magsq()const{
		return comp[0]*comp[0]+comp[1]*comp[1]+comp[2]*comp[2];
	}

	//io-routines:
	void read(XMLReader& xml, const std::string& path, momentum &mom){
	  XMLReader top(xml,path);
	  multi1d<int> tmp;

	  try {
	    read(top,"components",tmp);
	  }
	  catch( const std::string& e ) {
	    QDPIO::cerr << "Caught Exception : " << e << std::endl;
	    QDP_abort(1);
	  }
	  if(tmp.size()!=3){
	    QDPIO::cerr << "Error reading components for momentum vector: found " << tmp.size() << " components, only 3 are supported." << std::endl;
	    QDP_abort(1);
	  }

	  for(unsigned int i=0; i<3; i++) mom[i]=tmp[i];
	}

	std::ostream& operator<<(std::ostream &os,const momentum &obj){
	  os << "(" << obj[0] << "," << obj[1] << "," << obj[2] << ")";
	  return os;
	}

	StandardOutputStream& operator<<(StandardOutputStream &os,const momentum &obj){
	  os << "(" << obj[0] << "," << obj[1] << "," << obj[2] << ")";
	  return os;
	}

	//projection routines:
	LatticeComplex get_phases(momentum mom, int time_dir, Real fact){
	  START_CODE();

	  Real twopi = Chroma::twopi;
	  LatticeReal p_dot_x(float(0.0));

	  int j = 0;
	  for(int nu = 0; nu < Nd; ++nu){
	    //skip time coordinate:
	    if (nu == time_dir) continue;

	    p_dot_x += LatticeReal(fact * Layout::latticeCoordinate(nu)) * twopi * static_cast<Real>(mom[j]) / Layout::lattSize()[nu];
	    j++;
	  }

	  END_CODE();

	  // The complex phases  exp(i p.x )
	  return LatticeComplex(cmplx(cos(p_dot_x), sin(p_dot_x)));
	}

	//get real and imaginary parts separately
	LatticeComplex get_sin_phases(momentum mom, int time_dir, Real fact){
	  START_CODE();

	  Real twopi = Chroma::twopi;
	  LatticeReal p_dot_x(float(0.0));

	  int j = 0;
	  for(int nu = 0; nu < Nd; ++nu){
	    //skip time coordinate:
	    if (nu == time_dir) continue;

	    p_dot_x += LatticeReal(fact * Layout::latticeCoordinate(nu)) * twopi * static_cast<Real>(mom[j]) / Layout::lattSize()[nu];
	    j++;
	  }

	  END_CODE();

	  // The real phases  sin( p.x )
	  return LatticeComplex(cmplx(sin(p_dot_x),0.));
	}

	//get real and imaginary parts separately
	LatticeComplex get_cos_phases(momentum mom, int time_dir, Real fact){
	  START_CODE();

	  Real twopi = Chroma::twopi;
	  LatticeReal p_dot_x(float(0.0));

	  int j = 0;
	  for(int nu = 0; nu < Nd; ++nu){
	    //skip time coordinate:
	    if (nu == time_dir) continue;

	    p_dot_x += LatticeReal(fact * Layout::latticeCoordinate(nu)) * twopi * static_cast<Real>(mom[j]) / Layout::lattSize()[nu];
	    j++;
	  }

	  END_CODE();

	  // The real phases cos( p.x )
	  return LatticeComplex(cmplx(cos(p_dot_x),0.));
	}


	//create non-local A1 symmetric sink: 1/(sqrt(6)*sqrt(V)) * sum_{i=1}^3 cos( p * r * e_i ), where r is the displacement and e_i the i-th unit vector
	LatticeComplex getPhasesA1(int displacement, int time_dir, Real fact){
	  START_CODE();

	  LatticeComplex result=zero;
	  momentum shift;
	  Real norm=0.816496580927726;
	  //fourier factor:
	  int vol3=1;
	  for(int nu = 0; nu < Nd; ++nu){
	    //skip time coordinate:
	    if (nu == time_dir) continue;
	    vol3*=Layout::lattSize()[nu];

	    for(int mu=0; mu<Nd; mu++) shift[mu]=0;
	    shift[nu]=displacement;
	    result+=get_cos_phases(shift,time_dir,fact);
	  }
	  result*=norm/sqrt(static_cast<Real>(vol3));

	  END_CODE();

	  return result;
	}

	LatticeComplex getPhasesT1d1(int component, int displacement, int time_dir, Real fact){
	  START_CODE();

	  LatticeComplex result=zero;
	  Complex imath=cmplx(0.,static_cast<Real>(1.));
	  momentum shift;
	  Real norm;
	  for(int mu=0; mu<Nd; mu++) shift[mu]=0;

	  int vol3=1;
	  for(int nu = 0; nu < Nd; ++nu){
	    if (nu == time_dir) continue;
	    vol3*=Layout::lattSize()[nu];
	  }

	  switch(component){
	  case 1:
	    //iD_x-D_y at source, so iD_x+D_y at sink, note that exp(ipx)-exp(-ipx)= 2i sin(px)
	    norm=0.5/sqrt(static_cast<Real>(vol3));
	    shift[0]=displacement;
	    result-=2.*get_sin_phases(shift,time_dir,fact);
	    shift[0]=0;
	    shift[1]=displacement;
	    result+=2.*imath*get_sin_phases(shift,time_dir,fact);
	    result*=norm;
	    break;
	  case -1:
	    //-iD_x-D_y at source, so -iD_x+D_y at sink, note that exp(ipx)-exp(-ipx)= 2i sin(px)
	    norm=0.5/sqrt(static_cast<Real>(vol3));
	    shift[0]=displacement;
	    result+=2.*get_sin_phases(shift,time_dir,fact);
	    shift[0]=0;
	    shift[1]=displacement;
	    result+=2.*imath*get_sin_phases(shift,time_dir,fact);
	    result*=norm;
	    break;
	  case 0:
	    //-iD_z at source, so -iD_z at sink
	    norm=0.7071067811865475/sqrt(static_cast<Real>(vol3));
	    shift[2]=displacement;
	    result+=2.*get_sin_phases(shift,time_dir,fact);
	    result*=norm;
	    break;
	  }
	  END_CODE();

	  return result;
	}
}
