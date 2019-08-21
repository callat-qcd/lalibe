//Arjun Singh Gambhir
//Code ported from latscat.
#ifndef _TRANSFORM
#define _TRANSFORM

//Include chroma stuff this way
#include "chromabase.h"
//Also, we need to put this stuff inside the chroma namespace.

namespace Chroma {

  //shifting coordinates:
  multi1d<int> shiftcoord(const multi1d<int>& pos, const unsigned int& dir, const int& distance);
  multi1d<int> shiftcoord(const multi1d<int>& pos, const multi1d<int>& shift, const int& sign=+1);

  //mapping function for shifting output:
  struct Timeshift : public MapFunc
  {
    int t0, tDir, tLength;

    Timeshift(){}
    Timeshift(int t00, int tDirr, int tLengthh) : t0(t00), tDir(tDirr), tLength(tLengthh) {}

    void set(int t00, int tDirr, int tLengthh){
      t0=t00;
      tDir=tDirr;
      tLength=tLengthh;
    }

    multi1d<int> operator()(const multi1d<int>& x, int sign)const{
      multi1d<int> result(x);
      if(sign<0) result[tDir]=(result[tDir]+tLength-t0)%tLength;
      else result[tDir]=(result[tDir]+t0)%tLength;
      return result;
    }
  };

  struct Timeshiftmap{
    int t0, t1, t2, tDir, tLength;
    unsigned int nsites_t;
    Timeshift tshift1, tshift2;
    Map map1, map2;
    
    Timeshiftmap(const int& t00, const int& tDirr, const int& tLengthh) : t0(t00), tDir(tDirr), tLength(tLengthh){
      nsites_t=Layout::subgridLattSize()[tDir];
      
      //create a map to perform a shift in order to align with node boundaries
      t1=t0%nsites_t;
      tshift1.set(t1,tDir,tLength);
      map1.make(tshift1);
      
      t2=(t0-t1);
      tshift2.set(t2,tDir,tLength);
      map2.make(tshift2);
    }
    
    template<class T>
    OLattice<T> operator()(const OLattice<T>& field, const bool& antiper=false){
      OLattice<T> fieldcpy;
      OLattice<T> result;
      
      //multiply with sign for anti-periodic quantities
      if(antiper){
        fieldcpy=where( Layout::latticeCoordinate(tDir)>=t0, field, (-1)*field);
      }
      else fieldcpy=field;
      
      if(t1>0)result=map1(fieldcpy);
      else result=fieldcpy;
      if(t2>0) return map2(result);
      else return result;
    }
  };

  //arbitrary shift
  template<class T>
  OLattice<T> vecshift(const OLattice<T>& field, const multi1d<int>& shiftvec, const int& sign=+1){
    if(shiftvec.size()!=Nd){
      QDPIO::cout << "Error, please specify a valid shit vector (all 4 dimensions must be specified)!" << std::endl;
      return field;
    }
    
    OLattice<T> result(field);
    for(unsigned int dd=0; dd<Nd; dd++){
      unsigned int numshifts=abs(shiftvec[dd]);
      
      if(numshifts>0){
        //this can be confusing: it actually defines an "inverse" shift on f, such that ftilde(x)=f(x-shift), or ftilde(shift)=f(0):
        int direction=( shiftvec[dd]*sign>0 ? FORWARD : BACKWARD); 
        
        for(unsigned int s=0; s<numshifts; s++){
          result=shift(result, direction, dd);
        }
      }
    }
    return result;
  }
}
#endif
