//Arjun Singh Gambhir
//Code ported from latscat.
#ifndef _CHECKPOINTSTUFF
#define _CHECKPOINTSTUFF

#include "chromabase.h"
//This should bring a lot of stuff with it related to qdp objects, hdf5, and std types.

class checkpoint{
private:
	std::string filename, mainpath;
	HDF5Writer file;
	bool isopen;
	int stripesize;

public:
	//constructors:
	checkpoint(const std::string& filenamee) : isopen(false), stripesize(-1), mainpath(""){
		filename=filenamee;
#ifdef PROFILE
		file.set_profiling(true);
#endif
	}
  
	checkpoint(const std::string& filenamee, const int& stripesizee) : isopen(false), stripesize(stripesizee), mainpath(""){
		filename=filenamee;
#ifdef PROFILE
		file.set_profiling(true);
#endif
	}
  
	//access members:
	bool check_exist(){
		return file.check_exists(filename);
	}

	void set_stripesize(const int& stripesizee){stripesize=stripesizee;}

	void open(const HDF5Base::writemode& mode=HDF5Base::ate){
		if(!isopen){
			file.open(filename,mode);
			file.set_stripesize(stripesize);
		}
		isopen=true;
	}

	void close(){
		if(isopen) file.close();
		isopen=false;
	}

	void remove(){
		close();
		std::remove(filename.c_str());
	}

	//data manipulation members
	void set_mainpath(const std::string& pathname){
		mainpath=pathname;
	}

	//state handling
	void set_consistency(const bool& state){
		if(!isopen) open(HDF5Base::ate);
		file.write(mainpath+".consistent",static_cast<int>(state),HDF5Base::trunc);
	}
	
	bool check_consistency(){
		if(!isopen) open(HDF5Base::ate);
		int state;
		QDPIO::cout << "checking consistency";
		file.read(mainpath+".consistent",state);
		if(state==false){
			QDPIO::cout << " failed!" << std::endl;
		}
		else QDPIO::cout << " done!" << std::endl;
		return static_cast<bool>(state);
	}

	//checking parameters:
	template<class T>
	bool check_parameter(const std::string& paramname, const T& val){
		if(!isopen) open(HDF5Base::ate);

		QDPIO::cout << "checking parameter " << paramname << std::flush;
		T readval;
		file.read(mainpath+paramname,readval);

		if(readval!=val){
			QDPIO::cout << " failed!" << std::endl;
			return false;
		}
		else{
			QDPIO::cout << " done!" << std::endl;
			return true;
		}
	}

	template<class T>
	bool check_parameter(const std::string& paramname, const OScalar<T>& val){
		if(!isopen) open(HDF5Base::ate);

		QDPIO::cout<< "checking parameter " << paramname << std::flush;
		OScalar<T> readval;
		file.read(mainpath+paramname,readval);
		if( (toFloat(fabs(readval-val))>std::numeric_limits<REAL>::epsilon()) ){
			QDPIO::cout << " failed!" << std::endl;
			return false;
		}
		else{
			QDPIO::cout << " done!" << std::endl;
			return true;
		}
	}
	
	template<class T>
	bool check_parameter(const std::string& paramname, const multi2d<T>& val){
		if(!isopen) open(HDF5Base::ate);

		QDPIO::cout << "checking parameter " << paramname << std::flush;
		multi2d<T> readval;
		file.read(mainpath+paramname,readval);

		if( (readval.nrows()!=val.nrows()) || (readval.ncols()!=val.ncols()) ){
			QDPIO::cout << " failed!" << std::endl;
			return false;
		}

		for(unsigned int n=0; n<readval.nrows(); n++){
			for(unsigned int m=0; m<readval.ncols(); m++){
				if(readval(n,m)!=val(n,m)){
					QDPIO::cout << " failed at component (" << n << "," << m << ")!" << std::endl;
					return false;
				}
			}
		}
		QDPIO::cout << " done!" << std::endl;
		return true;
	}

	//create directory:
	int create_directory(const std::string& dirname);

	//setting parameters:
	template<class T>
	int set_parameter(const std::string& paramname, const T& val){
		if(!isopen) open(HDF5Base::ate);
    
		QDPIO::cout<< "writing parameter " << paramname << std::flush;
		file.write(mainpath+paramname,val,HDF5Base::trunc);
		QDPIO::cout << " done!" << std::endl;
		return EXIT_SUCCESS;
	}

	template<class T>
	int set_parameter(const std::string& paramname, const multi1d<T>& arrayval){
		if(!isopen) open(HDF5Base::ate);

		QDPIO::cout<< "writing parameter " << paramname << std::flush;
		file.write(mainpath+paramname,arrayval,HDF5Base::trunc);
		QDPIO::cout << " done!" << std::endl;
		return EXIT_SUCCESS;
	}
	
	template<class T>
	int set_parameter(const std::string& paramname, const multi2d<T>& array2val){
		if(!isopen) open(HDF5Base::ate);

		QDPIO::cout<< "writing parameter " << paramname << std::flush;
		file.write(mainpath+paramname,array2val,HDF5Base::trunc);
		QDPIO::cout << " done!" << std::endl;
		return EXIT_SUCCESS;
	}
	
	template<class T>
	int set_parameter(const std::string& paramname, const PSeed<T>& seed){
		if(!isopen) open(HDF5Base::ate);
		multi1d<T> readseed(4);
		for(unsigned int i=0; i<4; i++) readseed[i]=seed.elem(i);
		file.write(mainpath+paramname,readseed,HDF5Base::trunc);
		return EXIT_SUCCESS;
	}

	//getting parameters:
	template<class T>
	int get_parameter(const std::string& paramname, T& val){
		if(!isopen) open(HDF5Base::ate);

		QDPIO::cout<< "reading parameter " << paramname << std::flush;
		file.read(mainpath+paramname,val);
		QDPIO::cout << " done!" << std::endl;
		return EXIT_SUCCESS;
	}

	template<class T>
	int get_parameter(const std::string& paramname, multi1d<T>& arrayval){
		if(!isopen) open(HDF5Base::ate);
		QDPIO::cout<< "reading parameter " << paramname << std::flush;
		file.read(mainpath+paramname,arrayval);
		QDPIO::cout << " done!" << std::endl;
		return EXIT_SUCCESS;
	}
	
	template<class T>
	int get_parameter(const std::string& paramname, multi2d<T>& array2val){
		if(!isopen) open(HDF5Base::ate);
		QDPIO::cout<< "reading parameter " << paramname << std::flush;
		file.read(mainpath+paramname,array2val);
		QDPIO::cout << " done!" << std::endl;
		return EXIT_SUCCESS;
	}
	
	template<class T>
	int get_parameter(const std::string& paramname, PSeed<T>& seed){
		if(!isopen) open(HDF5Base::ate);
		multi1d<T> readseed(4);
		file.read(mainpath+paramname,readseed);
		for(unsigned int i=0; i<4; i++) seed.elem(i)=readseed[i];
		return EXIT_SUCCESS;
	}
};

//specializations
//checking
template<>
bool checkpoint::check_parameter<std::string>(const std::string& paramname, const std::string& val);
//setting
template<>
int checkpoint::set_parameter<Complex>(const std::string& paramname, const multi1d<Complex>& arrayval);
template<>
int checkpoint::set_parameter<int>(const std::string& paramname, const int& val);
template<>
int checkpoint::set_parameter<unsigned int>(const std::string& paramname, const unsigned int& val);
template<>
int checkpoint::set_parameter<LatticePropagatorD3>(const std::string& paramname, const LatticePropagatorD3& val);
template<>
int checkpoint::set_parameter<LatticePropagatorF3>(const std::string& paramname, const LatticePropagatorF3& val);
//getting
template<>
int checkpoint::get_parameter<Complex>(const std::string& paramname, multi1d<Complex>& arrayval);
#endif
