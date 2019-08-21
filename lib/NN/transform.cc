//Arjun Singh Gambhir
//Code ported from latscat.
//#include "progs.h"
//No progs included this way...

//this file is for generating source vectors:
//using namespace Chroma;
//Let's do this the lalibe way, also let's actually include the declarations, just like last time...
#include "transform.h"

namespace Chroma {

  //specialized shift in specified direction:
  multi1d<int> shiftcoord(const multi1d<int>& pos, const unsigned int& dir, const int& distance){
    multi1d<int> result(pos);
    
    if(abs(static_cast<int>(dir))>=Nd) return result;
    result[dir]=(result[dir]+Layout::lattSize()[dir]+distance)%Layout::lattSize()[dir];
    
    return result;
  }

  //general shift
  multi1d<int> shiftcoord(const multi1d<int>& pos, const multi1d<int>& shift, const int& sign){
    multi1d<int> result(pos);
    
    if(shift.size()>=Nd) return result;
    //allows for 3D and 4D shifts
    for(unsigned int d=0; d<shift.size(); d++){
      result[d]=(result[d]+Layout::lattSize()[d]+sign*shift[d])%Layout::lattSize()[d];
    }
    
    return result;
  }
}
