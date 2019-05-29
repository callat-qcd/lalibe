// -*- C++ -*-
/*! This file reads an object written from the totally snazzy hdf5 writer and puts in named object storage.
 *  If you find bugs, blame Arjun Singh Gambhir
 */

//Protect everything with a preprocessing directive.
#ifdef BUILD_HDF5

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/named_objmap.h"

#include "util/ferm/map_obj/map_obj_factory_w.h"
#include "util/ferm/map_obj/map_obj_aggregate_w.h"

#include "handle.h"
#include "actions/ferm/invert/containers.h"

#include "util/ferm/transf.h"

//LALIBE stuff
#include "hdf5_read_obj.h"

namespace Chroma 
{ 
  namespace LalibeHDF5ReadNamedObjEnv 
  { 
    namespace HDF5ReadObjectEnv
    { 
      class HDF5ReadObject
      {
      public:
	virtual ~HDF5ReadObject() {}

	virtual void operator()() = 0;
      };


      namespace
      {
	typedef SingletonHolder< 
	  ObjectFactory<HDF5ReadObject, 
			std::string,
			TYPELIST_1(const Params&),
			HDF5ReadObject* (*)(const Params&), 
			StringFactoryError> >
	TheHDF5ReadObjectFactory;
    
      }

      // Anonymous namespaces are so useful, get outta here static...
      namespace
      {
	class HDF5ReadLatProp : public HDF5ReadObject
	{
	private:
	  Params params;

	public:
	  HDF5ReadLatProp(const Params& p) : params(p) {}

	  void operator()() {
	    LatticePropagator obj;
	    //XMLReader file_xml, record_xml;

	    HDF5Reader reader;
	    reader.open(params.file.file_name);
	    reader.cd(params.file.path);
	    reader.read(params.file.obj_name,obj);
	    std::string file;
	    std::string record;
	    reader.readAttribute(params.file.obj_name, "file_xml", file);
	    reader.readAttribute(params.file.obj_name, "record_xml", record);
	    reader.cd("/");
	    reader.close();
	
	    /*std::string cut_name = params.file.file_name.substr(0, params.file.file_name.find(".h5", 0));
	    //Assume no conflicts with object names.
	    cut_name = cut_name+"_"+params.file.obj_name;
	    std::string file_xml_name = cut_name+"_file.xml";
	    std::string record_xml_name = cut_name+"_record.xml";
	    //Making new text files is a little janky, perhaps there's a way to do this internally?
	    TextFileWriter file_out(file_xml_name);
	    file_out<<file;
	    file_out.close();
	    TextFileWriter record_out(record_xml_name);
	    record_out<<record;
	    record_out.close();
	    XMLReader file_xml(file_xml_name);
	    XMLReader record_xml(record_xml_name);*/
	    //Above we were writing temporary files, which was annoying.
	    
	    std::istringstream  file_xml_stream(file);
	    std::istringstream  record_xml_stream(record);
	    XMLReader  file_xml(file_xml_stream);
	    XMLReader  record_xml(record_xml_stream);

	    TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.object_id);
	    TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.object_id) = obj;
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);
	  }
	};

	HDF5ReadObject* hdf5ReadLatProp(const Params& p)
	{
	  return new HDF5ReadLatProp(p);
	}


	class HDF5ReadLatPropF : public HDF5ReadObject
	{
	private:
	  Params params;

	public:
	  HDF5ReadLatPropF(const Params& p) : params(p) {}

	  void operator()() {
	    LatticePropagatorF obj;
	    //XMLReader file_xml, record_xml;
	   
	    HDF5Reader reader;
	    reader.open(params.file.file_name);
	    reader.cd(params.file.path);
	    reader.read(params.file.obj_name,obj);
	    std::string file;
	    std::string record;
	    reader.readAttribute(params.file.obj_name, "file_xml", file);
	    reader.readAttribute(params.file.obj_name, "record_xml", record);
	    reader.cd("/");
	    reader.close();
	
	    /*std::string cut_name = params.file.file_name.substr(0, params.file.file_name.find(".h5", 0));
	    //Assume no conflicts with object names.
	    cut_name = cut_name+"_"+params.file.obj_name;
	    std::string file_xml_name = cut_name+"_file.xml";
	    std::string record_xml_name = cut_name+"_record.xml";
	    //Making new text files is a little janky, perhaps there's a way to do this internally?
	    TextFileWriter file_out(file_xml_name);
	    file_out<<file;
	    file_out.close();
	    TextFileWriter record_out(record_xml_name);
	    record_out<<record;
	    record_out.close();
	    XMLReader file_xml(file_xml_name);
	    XMLReader record_xml(record_xml_name);*/
	    //Above we were writing temporary files, which was annoying.
	    
	    std::istringstream  file_xml_stream(file);
	    std::istringstream  record_xml_stream(record);
	    XMLReader  file_xml(file_xml_stream);
	    XMLReader  record_xml(record_xml_stream);

	    TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.object_id);
	    TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.object_id) = obj;
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);

	  }
	};

	HDF5ReadObject* hdf5ReadLatPropF(const Params& p)
	{
	  return new HDF5ReadLatPropF(p);
	}


	class HDF5ReadLatPropD : public HDF5ReadObject
	{
	private:
	  Params params;

	public:
	  HDF5ReadLatPropD(const Params& p) : params(p) {}

	  void operator()() {
	    LatticePropagatorD obj;
	    //XMLReader file_xml, record_xml;
	   
	    HDF5Reader reader;
	    reader.open(params.file.file_name);
	    reader.cd(params.file.path);
	    reader.read(params.file.obj_name,obj);
	    std::string file;
	    std::string record;
	    reader.readAttribute(params.file.obj_name, "file_xml", file);
	    reader.readAttribute(params.file.obj_name, "record_xml", record);
	    reader.cd("/");
	    reader.close();
	
	    /*std::string cut_name = params.file.file_name.substr(0, params.file.file_name.find(".h5", 0));
	    //Assume no conflicts with object names.
	    cut_name = cut_name+"_"+params.file.obj_name;
	    std::string file_xml_name = cut_name+"_file.xml";
	    std::string record_xml_name = cut_name+"_record.xml";
	    //Making new text files is a little janky, perhaps there's a way to do this internally?
	    TextFileWriter file_out(file_xml_name);
	    file_out<<file;
	    file_out.close();
	    TextFileWriter record_out(record_xml_name);
	    record_out<<record;
	    record_out.close();
	    XMLReader file_xml(file_xml_name);
	    XMLReader record_xml(record_xml_name);*/
	    //Above we were writing temporary files, which was annoying.
	    
	    std::istringstream  file_xml_stream(file);
	    std::istringstream  record_xml_stream(record);
	    XMLReader  file_xml(file_xml_stream);
	    XMLReader  record_xml(record_xml_stream);

	    TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.object_id);
	    TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.object_id) = obj;
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);

	  }
	};

	HDF5ReadObject* hdf5ReadLatPropD(const Params& p)
	{
	  return new HDF5ReadLatPropD(p);
	}
	
	class HDF5ReadLatUpperProp : public HDF5ReadObject
	{
	private:
	  Params params;

	public:
	  HDF5ReadLatUpperProp(const Params& p) : params(p) {}

	  void operator()() {
	    LatticePropagator obj;

	    HDF5Reader reader;
	    reader.open(params.file.file_name);
	    //reader.cd(params.file.path);
	    //We need to cd to the full path specified, since multiple fermions are stored inside.
	    std::string propagator_path=params.file.path+"/"+params.file.obj_name;
	    reader.cd(propagator_path);
	    std::string fermion_path;
	    //Now, we do the reading in this loop.
	    for(int color_source = 0; color_source < Nc ; ++color_source) 
	    {
	      for(int spin_source = Ns/2; spin_source < Ns; ++spin_source) 
	      { 
		int original_spin = spin_source - Ns/2;
		LatticeFermion psi;
	
		QDPIO::cout<<"Reading color: "<<color_source<<" spin: "<<original_spin<<" and copying into spin "<<spin_source<<std::endl;
		fermion_path = "color_"+std::to_string(color_source)+"_spin_"+std::to_string(original_spin);
		reader.read(fermion_path,psi);
		FermToProp(psi, obj, color_source, original_spin);
		FermToProp(LatticeFermion(-psi), obj, color_source, spin_source);
	      }
	    }
	    std::string file;
	    std::string record;
	    reader.readAttribute(fermion_path, "file_xml", file);
	    reader.readAttribute(fermion_path, "record_xml", record);
	    reader.cd("/");
	    reader.close();
	
	    std::istringstream  file_xml_stream(file);
	    std::istringstream  record_xml_stream(record);
	    XMLReader  file_xml(file_xml_stream);
	    XMLReader  record_xml(record_xml_stream);

	    TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.object_id);
	    TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.object_id) = obj;
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);
	  }
	};

	HDF5ReadObject* hdf5ReadLatUpperProp(const Params& p)
	{
	  return new HDF5ReadLatUpperProp(p);
	}
	
	class HDF5ReadLatUpperPropF : public HDF5ReadObject
	{
	private:
	  Params params;

	public:
	  HDF5ReadLatUpperPropF(const Params& p) : params(p) {}

	  void operator()() {
	    LatticePropagatorF obj;

	    HDF5Reader reader;
	    reader.open(params.file.file_name);
	    //reader.cd(params.file.path);
	    //We need to cd to the full path specified, since multiple fermions are stored inside.
	    std::string propagator_path=params.file.path+"/"+params.file.obj_name;
	    reader.cd(propagator_path);
	    std::string fermion_path;
	    //Now, we do the reading in this loop.
	    for(int color_source = 0; color_source < Nc ; ++color_source) 
	    {
	      for(int spin_source = Ns/2; spin_source < Ns; ++spin_source) 
	      { 
		int original_spin = spin_source - Ns/2;
		LatticeFermionF psi;
	
		QDPIO::cout<<"Reading color: "<<color_source<<" spin: "<<original_spin<<" and copying into spin "<<spin_source<<std::endl;
		fermion_path = "color_"+std::to_string(color_source)+"_spin_"+std::to_string(original_spin);
		reader.read(fermion_path,psi);
		FermToProp(psi, obj, color_source, original_spin);
		FermToProp(LatticeFermion(-psi), obj, color_source, spin_source);
	      }
	    }
	    std::string file;
	    std::string record;
	    reader.readAttribute(fermion_path, "file_xml", file);
	    reader.readAttribute(fermion_path, "record_xml", record);
	    reader.cd("/");
	    reader.close();
	
	    std::istringstream  file_xml_stream(file);
	    std::istringstream  record_xml_stream(record);
	    XMLReader  file_xml(file_xml_stream);
	    XMLReader  record_xml(record_xml_stream);

	    TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.object_id);
	    TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.object_id) = obj;
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);
	  }
	};

	HDF5ReadObject* hdf5ReadLatUpperPropF(const Params& p)
	{
	  return new HDF5ReadLatUpperPropF(p);
	}
	
	class HDF5ReadLatUpperPropD : public HDF5ReadObject
	{
	private:
	  Params params;

	public:
	  HDF5ReadLatUpperPropD(const Params& p) : params(p) {}

	  void operator()() {
	    LatticePropagatorD obj;

	    HDF5Reader reader;
	    reader.open(params.file.file_name);
	    //reader.cd(params.file.path);
	    //We need to cd to the full path specified, since multiple fermions are stored inside.
	    std::string propagator_path=params.file.path+"/"+params.file.obj_name;
	    reader.cd(propagator_path);
	    std::string fermion_path;
	    //Now, we do the reading in this loop.
	    for(int color_source = 0; color_source < Nc ; ++color_source) 
	    {
	      for(int spin_source = Ns/2; spin_source < Ns; ++spin_source) 
	      { 
		int original_spin = spin_source - Ns/2;
		LatticeFermionD psi;
	
		QDPIO::cout<<"Reading color: "<<color_source<<" spin: "<<original_spin<<" and copying into spin "<<spin_source<<std::endl;
		fermion_path = "color_"+std::to_string(color_source)+"_spin_"+std::to_string(original_spin);
		reader.read(fermion_path,psi);
		FermToProp(psi, obj, color_source, original_spin);
		FermToProp(LatticeFermion(-psi), obj, color_source, spin_source);
	      }
	    }
	    std::string file;
	    std::string record;
	    reader.readAttribute(fermion_path, "file_xml", file);
	    reader.readAttribute(fermion_path, "record_xml", record);
	    reader.cd("/");
	    reader.close();
	
	    std::istringstream  file_xml_stream(file);
	    std::istringstream  record_xml_stream(record);
	    XMLReader  file_xml(file_xml_stream);
	    XMLReader  record_xml(record_xml_stream);

	    TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.object_id);
	    TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.object_id) = obj;
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);
	  }
	};

	HDF5ReadObject* hdf5ReadLatUpperPropD(const Params& p)
	{
	  return new HDF5ReadLatUpperPropD(p);
	}
	
	class HDF5ReadLatLowerProp : public HDF5ReadObject
	{
	private:
	  Params params;

	public:
	  HDF5ReadLatLowerProp(const Params& p) : params(p) {}

	  void operator()() {
	    LatticePropagator obj;

	    HDF5Reader reader;
	    reader.open(params.file.file_name);
	    //reader.cd(params.file.path);
	    //We need to cd to the full path specified, since multiple fermions are stored inside.
	    std::string propagator_path=params.file.path+"/"+params.file.obj_name;
	    reader.cd(propagator_path);
	    std::string fermion_path;
	    //Now, we do the reading in this loop.
	    for(int color_source = 0; color_source < Nc ; ++color_source) 
	    {
	      for(int spin_source = 0; spin_source < Ns/2; ++spin_source) 
	      { 
		int original_spin = spin_source + Ns/2;
		LatticeFermion psi;
	
		QDPIO::cout<<"Reading color: "<<color_source<<" spin: "<<original_spin<<" and copying into spin "<<spin_source<<std::endl;
		fermion_path = "color_"+std::to_string(color_source)+"_spin_"+std::to_string(original_spin);
		reader.read(fermion_path,psi);
		FermToProp(psi, obj, color_source, original_spin);
		FermToProp(psi, obj, color_source, spin_source);
	      }
	    }
	    std::string file;
	    std::string record;
	    reader.readAttribute(fermion_path, "file_xml", file);
	    reader.readAttribute(fermion_path, "record_xml", record);
	    reader.cd("/");
	    reader.close();
	
	    std::istringstream  file_xml_stream(file);
	    std::istringstream  record_xml_stream(record);
	    XMLReader  file_xml(file_xml_stream);
	    XMLReader  record_xml(record_xml_stream);

	    TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.object_id);
	    TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.object_id) = obj;
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);
	  }
	};

	HDF5ReadObject* hdf5ReadLatLowerProp(const Params& p)
	{
	  return new HDF5ReadLatLowerProp(p);
	}
	
	class HDF5ReadLatLowerPropF : public HDF5ReadObject
	{
	private:
	  Params params;

	public:
	  HDF5ReadLatLowerPropF(const Params& p) : params(p) {}

	  void operator()() {
	    LatticePropagatorF obj;

	    HDF5Reader reader;
	    reader.open(params.file.file_name);
	    //reader.cd(params.file.path);
	    //We need to cd to the full path specified, since multiple fermions are stored inside.
	    std::string propagator_path=params.file.path+"/"+params.file.obj_name;
	    reader.cd(propagator_path);
	    std::string fermion_path;
	    //Now, we do the reading in this loop.
	    for(int color_source = 0; color_source < Nc ; ++color_source) 
	    {
	      for(int spin_source = 0; spin_source < Ns/2; ++spin_source) 
	      { 
		int original_spin = spin_source + Ns/2;
		LatticeFermionF psi;
	
		QDPIO::cout<<"Reading color: "<<color_source<<" spin: "<<original_spin<<" and copying into spin "<<spin_source<<std::endl;
		fermion_path = "color_"+std::to_string(color_source)+"_spin_"+std::to_string(original_spin);
		reader.read(fermion_path,psi);
		FermToProp(psi, obj, color_source, original_spin);
		FermToProp(psi, obj, color_source, spin_source);
	      }
	    }
	    std::string file;
	    std::string record;
	    reader.readAttribute(fermion_path, "file_xml", file);
	    reader.readAttribute(fermion_path, "record_xml", record);
	    reader.cd("/");
	    reader.close();
	
	    std::istringstream  file_xml_stream(file);
	    std::istringstream  record_xml_stream(record);
	    XMLReader  file_xml(file_xml_stream);
	    XMLReader  record_xml(record_xml_stream);

	    TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.object_id);
	    TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.object_id) = obj;
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);
	  }
	};

	HDF5ReadObject* hdf5ReadLatLowerPropF(const Params& p)
	{
	  return new HDF5ReadLatLowerPropF(p);
	}
	
	class HDF5ReadLatLowerPropD : public HDF5ReadObject
	{
	private:
	  Params params;

	public:
	  HDF5ReadLatLowerPropD(const Params& p) : params(p) {}

	  void operator()() {
	    LatticePropagatorD obj;

	    HDF5Reader reader;
	    reader.open(params.file.file_name);
	    //reader.cd(params.file.path);
	    //We need to cd to the full path specified, since multiple fermions are stored inside.
	    std::string propagator_path=params.file.path+"/"+params.file.obj_name;
	    reader.cd(propagator_path);
	    std::string fermion_path;
	    //Now, we do the reading in this loop.
	    for(int color_source = 0; color_source < Nc ; ++color_source) 
	    {
	      for(int spin_source = 0; spin_source < Ns/2; ++spin_source) 
	      { 
		int original_spin = spin_source + Ns/2;
		LatticeFermionD psi;
	
		QDPIO::cout<<"Reading color: "<<color_source<<" spin: "<<original_spin<<" and copying into spin "<<spin_source<<std::endl;
		fermion_path = "color_"+std::to_string(color_source)+"_spin_"+std::to_string(original_spin);
		reader.read(fermion_path,psi);
		FermToProp(psi, obj, color_source, original_spin);
		FermToProp(psi, obj, color_source, spin_source);
	      }
	    }
	    std::string file;
	    std::string record;
	    reader.readAttribute(fermion_path, "file_xml", file);
	    reader.readAttribute(fermion_path, "record_xml", record);
	    reader.cd("/");
	    reader.close();
	
	    std::istringstream  file_xml_stream(file);
	    std::istringstream  record_xml_stream(record);
	    XMLReader  file_xml(file_xml_stream);
	    XMLReader  record_xml(record_xml_stream);

	    TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.object_id);
	    TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.object_id) = obj;
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);
	  }
	};

	HDF5ReadObject* hdf5ReadLatLowerPropD(const Params& p)
	{
	  return new HDF5ReadLatLowerPropD(p);
	}



	class HDF5ReadStagLatProp : public HDF5ReadObject
	{
	private:
	  Params params;

	public:
	  HDF5ReadStagLatProp(const Params& p) : params(p) {}

	  void operator()() {
	    LatticeStaggeredPropagator obj;
	    //XMLReader file_xml, record_xml;
	   
	    HDF5Reader reader;
	    reader.open(params.file.file_name);
	    reader.cd(params.file.path);
	    reader.read(params.file.obj_name,obj);
	    std::string file;
	    std::string record;
	    reader.readAttribute(params.file.obj_name, "file_xml", file);
	    reader.readAttribute(params.file.obj_name, "record_xml", record);
	    reader.cd("/");
	    reader.close();
	
	    /*std::string cut_name = params.file.file_name.substr(0, params.file.file_name.find(".h5", 0));
	    //Assume no conflicts with object names.
	    cut_name = cut_name+"_"+params.file.obj_name;
	    std::string file_xml_name = cut_name+"_file.xml";
	    std::string record_xml_name = cut_name+"_record.xml";
	    //Making new text files is a little janky, perhaps there's a way to do this internally?
	    TextFileWriter file_out(file_xml_name);
	    file_out<<file;
	    file_out.close();
	    TextFileWriter record_out(record_xml_name);
	    record_out<<record;
	    record_out.close();
	    XMLReader file_xml(file_xml_name);
	    XMLReader record_xml(record_xml_name);*/
	    //Above we were writing temporary files, which was annoying.
	    
	    std::istringstream  file_xml_stream(file);
	    std::istringstream  record_xml_stream(record);
	    XMLReader  file_xml(file_xml_stream);
	    XMLReader  record_xml(record_xml_stream);

	    TheNamedObjMap::Instance().create<LatticeStaggeredPropagator>(params.named_obj.object_id);
	    TheNamedObjMap::Instance().getData<LatticeStaggeredPropagator>(params.named_obj.object_id) = obj;
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);

	  }
	};
	

	HDF5ReadObject* hdf5ReadStagLatProp(const Params& p)
	{
	  return new HDF5ReadStagLatProp(p);
	}


	class HDF5ReadStagLatPropF : public HDF5ReadObject
	{
	private:
	  Params params;

	public:
	  HDF5ReadStagLatPropF(const Params& p) : params(p) {}

	  void operator()() {
	    LatticeStaggeredPropagatorF obj;
	    //XMLReader file_xml, record_xml;
	   
	    HDF5Reader reader;
	    reader.open(params.file.file_name);
	    reader.cd(params.file.path);
	    reader.read(params.file.obj_name,obj);
	    std::string file;
	    std::string record;
	    reader.readAttribute(params.file.obj_name, "file_xml", file);
	    reader.readAttribute(params.file.obj_name, "record_xml", record);
	    reader.cd("/");
	    reader.close();
	
	    /*std::string cut_name = params.file.file_name.substr(0, params.file.file_name.find(".h5", 0));
	    //Assume no conflicts with object names.
	    cut_name = cut_name+"_"+params.file.obj_name;
	    std::string file_xml_name = cut_name+"_file.xml";
	    std::string record_xml_name = cut_name+"_record.xml";
	    //Making new text files is a little janky, perhaps there's a way to do this internally?
	    TextFileWriter file_out(file_xml_name);
	    file_out<<file;
	    file_out.close();
	    TextFileWriter record_out(record_xml_name);
	    record_out<<record;
	    record_out.close();
	    XMLReader file_xml(file_xml_name);
	    XMLReader record_xml(record_xml_name);*/
	    //Above we were writing temporary files, which was annoying.
	    
	    std::istringstream  file_xml_stream(file);
	    std::istringstream  record_xml_stream(record);
	    XMLReader  file_xml(file_xml_stream);
	    XMLReader  record_xml(record_xml_stream);

	    TheNamedObjMap::Instance().create<LatticeStaggeredPropagator>(params.named_obj.object_id);
	    TheNamedObjMap::Instance().getData<LatticeStaggeredPropagator>(params.named_obj.object_id) = obj;
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);

	  }
	};

	HDF5ReadObject* hdf5ReadStagLatPropF(const Params& p)
	{
	  return new HDF5ReadStagLatPropF(p);
	}

	class HDF5ReadStagLatPropD : public HDF5ReadObject
	{
	private:
	  Params params;

	public:
	  HDF5ReadStagLatPropD(const Params& p) : params(p) {}

	  void operator()() {
	    LatticeStaggeredPropagatorD obj;
	    //XMLReader file_xml, record_xml;
	   
	    HDF5Reader reader;
	    reader.open(params.file.file_name);
	    reader.cd(params.file.path);
	    reader.read(params.file.obj_name,obj);
	    std::string file;
	    std::string record;
	    reader.readAttribute(params.file.obj_name, "file_xml", file);
	    reader.readAttribute(params.file.obj_name, "record_xml", record);
	    reader.cd("/");
	    reader.close();
	
	    /*std::string cut_name = params.file.file_name.substr(0, params.file.file_name.find(".h5", 0));
	    //Assume no conflicts with object names.
	    cut_name = cut_name+"_"+params.file.obj_name;
	    std::string file_xml_name = cut_name+"_file.xml";
	    std::string record_xml_name = cut_name+"_record.xml";
	    //Making new text files is a little janky, perhaps there's a way to do this internally?
	    TextFileWriter file_out(file_xml_name);
	    file_out<<file;
	    file_out.close();
	    TextFileWriter record_out(record_xml_name);
	    record_out<<record;
	    record_out.close();
	    XMLReader file_xml(file_xml_name);
	    XMLReader record_xml(record_xml_name);*/
	    //Above we were writing temporary files, which was annoying.
	    
	    std::istringstream  file_xml_stream(file);
	    std::istringstream  record_xml_stream(record);
	    XMLReader  file_xml(file_xml_stream);
	    XMLReader  record_xml(record_xml_stream);

	    TheNamedObjMap::Instance().create<LatticeStaggeredPropagator>(params.named_obj.object_id);
	    TheNamedObjMap::Instance().getData<LatticeStaggeredPropagator>(params.named_obj.object_id) = obj;
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);

	  }
	};

	HDF5ReadObject* hdf5ReadStagLatPropD(const Params& p)
	{
	  return new HDF5ReadStagLatPropD(p);
	}



	class HDF5ReadLatFerm : public HDF5ReadObject
	{
	private:
	  Params params;

	public:
	  HDF5ReadLatFerm(const Params& p) : params(p) {}

	  //! Read a propagator
	  void operator()() {
	    LatticeFermion obj;
	    //XMLReader file_xml, record_xml;
	   
	    HDF5Reader reader;
	    reader.open(params.file.file_name);
	    reader.cd(params.file.path);
	    reader.read(params.file.obj_name,obj);
	    std::string file;
	    std::string record;
	    reader.readAttribute(params.file.obj_name, "file_xml", file);
	    reader.readAttribute(params.file.obj_name, "record_xml", record);
	    reader.cd("/");
	    reader.close();
	
	    /*std::string cut_name = params.file.file_name.substr(0, params.file.file_name.find(".h5", 0));
	    //Assume no conflicts with object names.
	    cut_name = cut_name+"_"+params.file.obj_name;
	    std::string file_xml_name = cut_name+"_file.xml";
	    std::string record_xml_name = cut_name+"_record.xml";
	    //Making new text files is a little janky, perhaps there's a way to do this internally?
	    TextFileWriter file_out(file_xml_name);
	    file_out<<file;
	    file_out.close();
	    TextFileWriter record_out(record_xml_name);
	    record_out<<record;
	    record_out.close();
	    XMLReader file_xml(file_xml_name);
	    XMLReader record_xml(record_xml_name);*/
	    //Above we were writing temporary files, which was annoying.
	    
	    std::istringstream  file_xml_stream(file);
	    std::istringstream  record_xml_stream(record);
	    XMLReader  file_xml(file_xml_stream);
	    XMLReader  record_xml(record_xml_stream);

	    TheNamedObjMap::Instance().create<LatticeFermion>(params.named_obj.object_id);
	    TheNamedObjMap::Instance().getData<LatticeFermion>(params.named_obj.object_id) = obj;
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);

	  }
	};

	HDF5ReadObject* hdf5ReadLatFerm(const Params& p)
	{
	  return new HDF5ReadLatFerm(p);
	}

//There are problems with casting a LatticeFermion between precisions, not sure what's going on here so I only do the "floating-precision" here.

	class HDF5ReadArrayLatColMat : public HDF5ReadObject
	{
	private:
	  Params params;

	public:
	  HDF5ReadArrayLatColMat(const Params& p) : params(p) {}

	  void operator()() {
	    multi1d<LatticeColorMatrix> obj;
	    //XMLReader file_xml, record_xml;
	   
	    HDF5Reader reader;
	    reader.open(params.file.file_name);
	    reader.cd(params.file.path);
	    //This needs to be fixed.
	    reader.read(params.file.obj_name,obj);
	    std::string file;
	    std::string record;
	    reader.readAttribute(params.file.obj_name, "file_xml", file);
	    reader.readAttribute(params.file.obj_name, "record_xml", record);
	    reader.cd("/");
	    reader.close();
	
	    /*std::string cut_name = params.file.file_name.substr(0, params.file.file_name.find(".h5", 0));
	    //Assume no conflicts with object names.
	    cut_name = cut_name+"_"+params.file.obj_name;
	    std::string file_xml_name = cut_name+"_file.xml";
	    std::string record_xml_name = cut_name+"_record.xml";
	    //Making new text files is a little janky, perhaps there's a way to do this internally?
	    TextFileWriter file_out(file_xml_name);
	    file_out<<file;
	    file_out.close();
	    TextFileWriter record_out(record_xml_name);
	    record_out<<record;
	    record_out.close();
	    XMLReader file_xml(file_xml_name);
	    XMLReader record_xml(record_xml_name);*/
	    //Above we were writing temporary files, which was annoying.
	    
	    std::istringstream  file_xml_stream(file);
	    std::istringstream  record_xml_stream(record);
	    XMLReader  file_xml(file_xml_stream);
	    XMLReader  record_xml(record_xml_stream);

	    TheNamedObjMap::Instance().create<multi1d<LatticeColorMatrix>>(params.named_obj.object_id);
	    TheNamedObjMap::Instance().getData<multi1d<LatticeColorMatrix>>(params.named_obj.object_id) = obj;
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);

	  }
	};

	HDF5ReadObject* hdf5ReadArrayLatColMat(const Params& p)
	{
	  return new HDF5ReadArrayLatColMat(p);
	}


/*	class HDF5ReadArrayLatColMatF : public HDF5ReadObject
	{
	private:
	  Params params;

	public:
	  HDF5ReadArrayLatColMatF(const Params& p) : params(p) {}

	  void operator()() {
	    multi1d<LatticeColorMatrixF> obj;
	    //XMLReader file_xml, record_xml;
	   
	    HDF5Reader reader;
	    reader.open(params.file.file_name);
	    reader.cd(params.file.path);
	    //This needs to be fixed.
	    //reader.read(params.file.obj_name,obj);
	    std::string file;
	    std::string record;
	    reader.readAttribute(params.file.obj_name, "file_xml", file);
	    reader.readAttribute(params.file.obj_name, "record_xml", record);
	    reader.cd("/");
	    reader.close();
	
	    std::string cut_name = params.file.file_name.substr(0, params.file.file_name.find(".h5", 0));
	    std::string file_xml_name = cut_name+"_file.xml";
	    std::string record_xml_name = cut_name+"_record.xml";
	    //Making new text files is a little janky, perhaps there's a way to do this internally?
	    TextFileWriter file_out(file_xml_name);
	    file_out<<file;
	    file_out.close();
	    TextFileWriter record_out(record_xml_name);
	    record_out<<record;
	    record_out.close();
	    XMLReader file_xml(file_xml_name);
	    XMLReader record_xml(record_xml_name);

	    TheNamedObjMap::Instance().create<multi1d<LatticeColorMatrix>>(params.named_obj.object_id);
	    //This needs to be fixed.
	    //TheNamedObjMap::Instance().getData<multi1d<LatticeColorMatrix>>(params.named_obj.object_id) = obj;
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);*/

	    /*TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(params.named_obj.object_id);
	    (TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.object_id)).resize(Nd);
	    for(int i=0; i < Nd; i++) {
	      (TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.object_id))[i] = obj[i];
	    }
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);*/
	  /*}
	};

	HDF5ReadObject* hdf5ReadArrayLatColMatF(const Params& p)
	{
	  return new HDF5ReadArrayLatColMatF(p);
	}


	class HDF5ReadArrayLatColMatD : public HDF5ReadObject
	{
	private:
	  Params params;

	public:
	  HDF5ReadArrayLatColMatD(const Params& p) : params(p) {}

	  void operator()() {
	    multi1d<LatticeColorMatrixD> obj;
	    //XMLReader file_xml, record_xml;
	   
	    HDF5Reader reader;
	    reader.open(params.file.file_name);
	    reader.cd(params.file.path);
	    //This needs to be fixed.
	    //reader.read(params.file.obj_name,obj);
	    std::string file;
	    std::string record;
	    reader.readAttribute(params.file.obj_name, "file_xml", file);
	    reader.readAttribute(params.file.obj_name, "record_xml", record);
	    reader.cd("/");
	    reader.close();
	
	    std::string cut_name = params.file.file_name.substr(0, params.file.file_name.find(".h5", 0));
	    std::string file_xml_name = cut_name+"_file.xml";
	    std::string record_xml_name = cut_name+"_record.xml";
	    //Making new text files is a little janky, perhaps there's a way to do this internally?
	    TextFileWriter file_out(file_xml_name);
	    file_out<<file;
	    file_out.close();
	    TextFileWriter record_out(record_xml_name);
	    record_out<<record;
	    record_out.close();
	    XMLReader file_xml(file_xml_name);
	    XMLReader record_xml(record_xml_name);

	    TheNamedObjMap::Instance().create<multi1d<LatticeColorMatrix>>(params.named_obj.object_id);
	    //This needs to be fixed.
	    //TheNamedObjMap::Instance().getData<multi1d<LatticeColorMatrix>>(params.named_obj.object_id) = obj;
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);*/

	    /*TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(params.named_obj.object_id);
	    (TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.object_id)).resize(Nd);
	    for(int i=0; i < Nd;i++) {
	      (TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.object_id))[i] = obj[i];
	    }

	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);*/
	  /*}
	};


	HDF5ReadObject* hdf5ReadArrayLatColMatD(const Params& p)
	{
	  return new HDF5ReadArrayLatColMatD(p);
	}*/


	bool registered = false;

      }  // ends the namespace


      bool registerAll() 
      {
	bool success = true; 
	if (! registered)
	{
	  success &= TheHDF5ReadObjectFactory::Instance().registerObject(std::string("LatticePropagator"),   
									hdf5ReadLatProp);
	  success &= TheHDF5ReadObjectFactory::Instance().registerObject(std::string("LatticePropagatorF"), 
									hdf5ReadLatPropF);
	  success &= TheHDF5ReadObjectFactory::Instance().registerObject(std::string("LatticePropagatorD"), 
									hdf5ReadLatPropD);

	  success &= TheHDF5ReadObjectFactory::Instance().registerObject(std::string("LatticeUpperPropagator"),   
									hdf5ReadLatUpperProp);
	  success &= TheHDF5ReadObjectFactory::Instance().registerObject(std::string("LatticeUpperPropagatorF"),   
									hdf5ReadLatUpperPropF);
	  success &= TheHDF5ReadObjectFactory::Instance().registerObject(std::string("LatticeUpperPropagatorD"),   
									hdf5ReadLatUpperPropD);
	  success &= TheHDF5ReadObjectFactory::Instance().registerObject(std::string("LatticeLowerPropagator"),   
									hdf5ReadLatLowerProp);
	  success &= TheHDF5ReadObjectFactory::Instance().registerObject(std::string("LatticeLowerPropagatorF"),   
									hdf5ReadLatLowerPropF);
	  success &= TheHDF5ReadObjectFactory::Instance().registerObject(std::string("LatticeLowerPropagatorD"),   
									hdf5ReadLatLowerPropD);

	  success &= TheHDF5ReadObjectFactory::Instance().registerObject(std::string("LatticeStaggeredPropagator"),   
									hdf5ReadStagLatProp);
	  success &= TheHDF5ReadObjectFactory::Instance().registerObject(std::string("LatticeStaggeredPropagatorF"), 
									hdf5ReadStagLatPropF);
	  success &= TheHDF5ReadObjectFactory::Instance().registerObject(std::string("LatticeStaggeredPropagatorD"), 
									hdf5ReadStagLatPropD);

	  success &= TheHDF5ReadObjectFactory::Instance().registerObject(std::string("LatticeFermion"), 
									hdf5ReadLatFerm);

//      success &= TheHDF5ReadObjectFactory::Instance().registerObject(std::string("LatticeFermionF"), 
//								   hdf5ReadLatFermF);
//      success &= TheHDF5ReadObjectFactory::Instance().registerObject(std::string("LatticeFermionD"), 
//								   hdf5ReadLatFermD);

	  success &= TheHDF5ReadObjectFactory::Instance().registerObject(std::string("Multi1dLatticeColorMatrix"), 
									hdf5ReadArrayLatColMat);

//	  success &= TheHDF5ReadObjectFactory::Instance().registerObject(std::string("Multi1dLatticeColorMatrixF"), 
//									hdf5ReadArrayLatColMatF);

//	  success &= TheHDF5ReadObjectFactory::Instance().registerObject(std::string("Multi1dLatticeColorMatrixD"), 
//									hdf5ReadArrayLatColMatD);

	  registered = true;
	}
	return success;
      }
    } // namespace HDF5ReadObjectEnv


    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }

      bool registered = false;

      const std::string name = "HDF5_READ_NAMED_OBJECT";
    }
    
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	success &= HDF5ReadObjectEnv::registerAll();
	success &= MapObjectWilson4DEnv::registerAll();
	registered = true;
      }
      return success;
    }


    //------------------------------------------------------------------------
    struct NamedObject_t
    {
      std::string   object_id;
      std::string   object_type;
    };


    void read(XMLReader& xml, const std::string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "object_id", input.object_id);
      read(inputtop, "object_type", input.object_type);
    }

    void read(XMLReader& xml, const std::string& path, Params::File_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "file_name", input.file_name);
      read(inputtop, "path", input.path);
      read(inputtop, "obj_name", input.obj_name);
    }


    Params::Params() { frequency = 0; }

    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      try 
      {
	XMLReader paramtop(xml_in, path);

	if (paramtop.count("Frequency") == 1)
	  read(paramtop, "Frequency", frequency);
	else
	  frequency = 1;

	named_obj_xml = readXMLGroup(paramtop, "NamedObject", "object_type");

	read(paramtop, "NamedObject", named_obj);

	read(paramtop, "File", file);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << __func__ << ": caught Exception reading XML: " << e << std::endl;
	QDP_abort(1);
      }
    }


    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      START_CODE();

      push(xml_out, "hdf5_read_named_obj");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": object reader" << std::endl;
      StopWatch swatch;

      try
      {	
	QDPIO::cout << "Attempt to read object name = " << params.named_obj.object_id << std::endl;
	write(xml_out, "object_id", params.named_obj.object_id);

	// Create the hdf5 reader (generic)
	Handle<HDF5ReadObjectEnv::HDF5ReadObject> hdf5ReadObject(
	  HDF5ReadObjectEnv::TheHDF5ReadObjectFactory::Instance().createObject(params.named_obj.object_type,
									     params));

	// Actually do the reading
	swatch.start();

	(*hdf5ReadObject)();

	swatch.stop();

	QDPIO::cout << "Object successfully read: time= " 
		    << swatch.getTimeInSeconds() 
		    << " secs" << std::endl;
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << name << ": caught dynamic cast error" 
		    << std::endl;
	QDP_abort(1);
      }
      catch (const std::string& e) 
      {
	QDPIO::cerr << name << ": std::map call failed: " << e 
		    << std::endl;
	QDP_abort(1);
      }
    
      QDPIO::cout << name << ": ran successfully" << std::endl;

      pop(xml_out);  // popping hdf5_read_named_obj xml

      END_CODE();
    } 

  }

}


#endif
