//Arjun Singh Gambhir
//Code ported from latscat.
//#include "progs.h"
//No progs, but we include the header here for consistency.
#include "checkpointstuff.h"

//create directory
int checkpoint::create_directory(const std::string& dirname){
	if(!isopen) open(HDF5Base::ate);
	
	QDPIO::cout << "create directory " << dirname << std::flush;
	file.mkdir(dirname);
	QDPIO::cout << " done!" << std::endl;
	return EXIT_SUCCESS;
}

//template specializations for checking parameters:
template<>
bool checkpoint::check_parameter<std::string>(const std::string& paramname, const std::string& val){
  if(!isopen) open(HDF5Base::ate);

  QDPIO::cout << "checking parameter " << paramname << std::flush;
  std::string readval;
  file.read(mainpath+paramname,readval);
  if(readval.compare(val.c_str())!=0){
    QDPIO::cout << " failed!" << std::endl;
    return false;
  }
  else{
    QDPIO::cout << " done!" << std::endl;
    return true;
  }
}

//template specializations for setting parameters
template<>
int checkpoint::set_parameter<Complex>(const std::string& paramname, const multi1d<Complex>& arrayval){
  if(!isopen) open(HDF5Base::ate);

  QDPIO::cout << "writing complex array parameter " << paramname << std::flush;
  file.write(mainpath+paramname,arrayval,HDF5Base::trunc);
  QDPIO::cout << " done!" << std::endl;
  return EXIT_SUCCESS;
}

template<>
int checkpoint::set_parameter<int>(const std::string& paramname, const int& val){
	if(!isopen) open(HDF5Base::ate);

	QDPIO::cout<< "writing parameter " << paramname << " with value= " << val << std::flush;
	file.write(mainpath+paramname,val,HDF5Base::trunc);
	QDPIO::cout << " done!" << std::endl;
	return EXIT_SUCCESS;
}

template<>
int checkpoint::set_parameter<unsigned int>(const std::string& paramname, const unsigned int& val){
	if(!isopen) open(HDF5Base::ate);

	QDPIO::cout<< "writing parameter " << paramname << " with value= " << val << std::flush;
	file.write(mainpath+paramname,val,HDF5Base::trunc);
	QDPIO::cout << " done!" << std::endl;
	return EXIT_SUCCESS;
}

template<>
int checkpoint::set_parameter<LatticePropagatorD3>(const std::string& paramname, const LatticePropagatorD3& val){
    if(!isopen) open(HDF5Base::ate);
    
    QDPIO::cout << "writing lattice propagator " << paramname << std::flush;
    file.write(mainpath+paramname,val,HDF5Base::trunc);
    QDPIO::cout << " done!" << std::endl;
    return EXIT_SUCCESS;
}

template<>
int checkpoint::set_parameter<LatticePropagatorF3>(const std::string& paramname, const LatticePropagatorF3& val){
    if(!isopen) open(HDF5Base::ate);
    
    QDPIO::cout << "writing lattice propagator " << paramname << std::flush;
    file.write(mainpath+paramname,val,HDF5Base::trunc);
    QDPIO::cout << " done!" << std::endl;
    return EXIT_SUCCESS;
}


//template specializations for getting parameters
template<>
int checkpoint::get_parameter<Complex>(const std::string& paramname, multi1d<Complex>& arrayval){
  if(!isopen) open(HDF5Base::ate);
  
  QDPIO::cout << "reading complex array parameter " << paramname << std::flush;
  file.read(mainpath+paramname,arrayval);
  QDPIO::cout << " done!" << std::endl;
  return EXIT_SUCCESS;
}
