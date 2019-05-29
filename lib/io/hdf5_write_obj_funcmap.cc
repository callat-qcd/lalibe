// -*- C++ -*-
/*! Arjun Singh Gambhir
 *  The meat of the hdf5 writing is being done here.
 */

//Protect everything with a preprocessing directive.
#ifdef BUILD_HDF5

#include "named_obj.h"
#include "meas/inline/io/named_objmap.h"

#include "handle.h"
#include "qdp_map_obj_memory.h"

#include "actions/ferm/invert/containers.h"
#include "meas/inline/make_xml_file.h"

#include "util/ferm/transf.h"

//LALIBE stuff...
#include "hdf5_write_obj_funcmap.h"

namespace Chroma
{
 
  //! IO function std::map environment
  /*! \ingroup inlineio */
  namespace HDF5WriteObjCallMapEnv
  { 
    // Anonymous namespace
    namespace
    {
      //------------------------------------------------------------------------
      //! Write a propagator
      void HDF5WriteLatProp(const std::string& buffer_id,
			   const std::string& outputfile,
			   const std::string& obj_name,
			   const std::string& path, HDF5Base::writemode wmode)
      {
	LatticePropagator obj;
	//XMLReader file_xml, record_xml;
	XMLBufferWriter file_xml, record_xml;

	//obj = TheNamedObjMap::Instance().getData<LatticePropagatorD3>(buffer_id);
	obj = TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);

	HDF5Writer h5out(outputfile);
	h5out.push(path);
	std::string propagator_path = path+"/"+obj_name;
	h5out.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB
	h5out.write(propagator_path, obj, wmode);
	std::string record = record_xml.str();
	std::string file = file_xml.str();
	//This string needs to be broadcasted, right now only the head node has it.
	QDPInternal::broadcast_str(file);
	QDPInternal::broadcast_str(record);
	std::string file_formatted = file;
	std::string record_formatted = record;
	//This is cutting xml record/file_xml and formatting them to be only slightly more human readable.
	for(int i = 0; i < record_formatted.length(); i++)
	{
	  if(record_formatted[i] == '\n'|| record_formatted[i] == ' ' || record_formatted[i] == '\"') 
	  {
	    record_formatted.erase(i, 1); 
	    i--;
	  }
	}
	for(int i = 0; i < file_formatted.length(); i++)
	{
	  if(file_formatted[i] == '\n'|| file_formatted[i] == ' ' || file_formatted[i] == '\"') 
	  {
	    file_formatted.erase(i, 1); 
	    i--;
	  }
	}
	//QDPIO::cout<<"This is the formatted xml_record that is going to be written: "<<std::endl<<record<<std::endl;
	h5out.writeAttribute(propagator_path, "file_xml", file, wmode);
	h5out.writeAttribute(propagator_path, "file_xml_formatted", file_formatted, wmode);
	h5out.writeAttribute(propagator_path, "record_xml", record, wmode);
	h5out.writeAttribute(propagator_path, "record_xml_formatted", record_formatted, wmode);

	h5out.cd("/");
	h5out.close();

      }

      //! Write a single prec propagator
      void HDF5WriteLatPropF(const std::string& buffer_id,
			   const std::string& outputfile,
			   const std::string& obj_name,
			   const std::string& path, HDF5Base::writemode wmode)
      {
    //LatticeDiracPropagatorF3 obj;
	LatticePropagatorF obj;
	XMLBufferWriter file_xml, record_xml;

	obj = TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);
	
	HDF5Writer h5out(outputfile);
	h5out.push(path);
	std::string propagator_path = path+"/"+obj_name;
	h5out.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB
	h5out.write(propagator_path, obj, wmode);
	std::string record = record_xml.str();
	std::string file = file_xml.str();
	//This string needs to be broadcasted, right now only the head node has it.
	QDPInternal::broadcast_str(file);
	QDPInternal::broadcast_str(record);
	std::string file_formatted = file;
	std::string record_formatted = record;
	//This is cutting xml record/file_xml and formatting them to be only slightly more human readable.
	for(int i = 0; i < record_formatted.length(); i++)
	{
	  if(record_formatted[i] == '\n'|| record_formatted[i] == ' ' || record_formatted[i] == '\"') 
	  {
	    record_formatted.erase(i, 1); 
	    i--;
	  }
	}
	for(int i = 0; i < file_formatted.length(); i++)
	{
	  if(file_formatted[i] == '\n'|| file_formatted[i] == ' ' || file_formatted[i] == '\"') 
	  {
	    file_formatted.erase(i, 1); 
	    i--;
	  }
	}
	//QDPIO::cout<<"This is the formatted xml_record that is going to be written: "<<std::endl<<record<<std::endl;
	h5out.writeAttribute(propagator_path, "file_xml", file, wmode);
	h5out.writeAttribute(propagator_path, "file_xml_formatted", file_formatted, wmode);
	h5out.writeAttribute(propagator_path, "record_xml", record, wmode);
	h5out.writeAttribute(propagator_path, "record_xml_formatted", record_formatted, wmode);

	h5out.cd("/");
	h5out.close();

      }


      //! Write a double prec propagator
      void HDF5WriteLatPropD(const std::string& buffer_id,
			   const std::string& outputfile,
			   const std::string& obj_name,
			   const std::string& path, HDF5Base::writemode wmode)
      {
	LatticePropagatorD obj;
	XMLBufferWriter file_xml, record_xml;

	obj = TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);

	HDF5Writer h5out(outputfile);
	h5out.push(path);
	std::string propagator_path = path+"/"+obj_name;
	h5out.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB
	h5out.write(propagator_path, obj, wmode);
	std::string record = record_xml.str();
	std::string file = file_xml.str();
	//This string needs to be broadcasted, right now only the head node has it.
	QDPInternal::broadcast_str(file);
	QDPInternal::broadcast_str(record);
	std::string file_formatted = file;
	std::string record_formatted = record;
	//This is cutting xml record/file_xml and formatting them to be only slightly more human readable.
	for(int i = 0; i < record_formatted.length(); i++)
	{
	  if(record_formatted[i] == '\n'|| record_formatted[i] == ' ' || record_formatted[i] == '\"') 
	  {
	    record_formatted.erase(i, 1); 
	    i--;
	  }
	}
	for(int i = 0; i < file_formatted.length(); i++)
	{
	  if(file_formatted[i] == '\n'|| file_formatted[i] == ' ' || file_formatted[i] == '\"') 
	  {
	    file_formatted.erase(i, 1); 
	    i--;
	  }
	}
	//QDPIO::cout<<"This is the formatted xml_record that is going to be written: "<<std::endl<<record<<std::endl;
	h5out.writeAttribute(propagator_path, "file_xml", file, wmode);
	h5out.writeAttribute(propagator_path, "file_xml_formatted", file_formatted, wmode);
	h5out.writeAttribute(propagator_path, "record_xml", record, wmode);
	h5out.writeAttribute(propagator_path, "record_xml_formatted", record_formatted, wmode);

	h5out.cd("/");
	h5out.close();

      }
      //! Write the upper two components of a non-relativistic half-spin propagator
      void HDF5WriteUpperLatProp(const std::string& buffer_id,
			   const std::string& outputfile,
			   const std::string& obj_name,
			   const std::string& path, HDF5Base::writemode wmode)
      {
	LatticePropagator obj;
	XMLBufferWriter file_xml, record_xml;

	obj = TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);

	//Check that this is a non-relativistic propagator before any writing.
	Double epsilon = 1e-8; //Maybe this should be more strict. 
	for(int color_source = 0; color_source < Nc ; ++color_source) 
	{
	  for(int spin_source = Ns/2; spin_source < Ns; ++spin_source) 
	  { 
	    int original_spin = spin_source - Ns/2;
	    LatticeFermion psi;
	    LatticeFermion chi;

	    PropToFerm(obj, psi, color_source, original_spin);
	    PropToFerm(obj, chi, color_source, spin_source);

	    QDPIO::cout<<"Checking color: "<<color_source<<" spin: "<<original_spin<<" against spin: "<<spin_source<<std::endl;
	    Double difference_norm = norm2(psi + chi);
	    QDPIO::cout<<"Difference norm is: "<<difference_norm<<std::endl;
	    if (toBool(difference_norm > epsilon))
	    {
	      QDPIO::cout<<"ERROR, based on checking the spin components, you are not really writing an upper non-relativistic propagator; exiting..."<<std::endl;
	      QDP_abort(1);
	    }
	  }
	}

	HDF5Writer h5out(outputfile);
	//h5out.push(path);
	std::string propagator_path = path+"/"+obj_name;
	h5out.push(propagator_path); //We want to push here since each spin component will be written separately.
	h5out.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB
	//Loop over the the upper spin components and write them.
	for(int color_source = 0; color_source < Nc ; ++color_source) 
	{
	  for(int spin_source = Ns/2; spin_source < Ns; ++spin_source) 
	  { 
	    int original_spin = spin_source - Ns/2;
	    LatticeFermion psi;

	    PropToFerm(obj, psi, color_source, original_spin);

	    std::string fermion_path = propagator_path+"/color_"+std::to_string(color_source)+"_spin_"+std::to_string(original_spin);
	    QDPIO::cout<<"Writing color: "<<color_source<<" spin: "<<original_spin<<std::endl;
	    //Now all the usual stuff happens, only difference is we are writing a fermion.
	    h5out.write(fermion_path, psi, wmode);
	    std::string record = record_xml.str();
	    std::string file = file_xml.str();
	    //This string needs to be broadcasted, right now only the head node has it.
	    QDPInternal::broadcast_str(file);
	    QDPInternal::broadcast_str(record);
	    std::string file_formatted = file;
	    std::string record_formatted = record;
	    //This is cutting xml record/file_xml and formatting them to be only slightly more human readable.
	    for(int i = 0; i < record_formatted.length(); i++)
	    {
	      if(record_formatted[i] == '\n'|| record_formatted[i] == ' ' || record_formatted[i] == '\"') 
	      {
		record_formatted.erase(i, 1); 
		i--;
	      }
	    }
	    for(int i = 0; i < file_formatted.length(); i++)
	    {
	      if(file_formatted[i] == '\n'|| file_formatted[i] == ' ' || file_formatted[i] == '\"') 
	      {
		file_formatted.erase(i, 1); 
		i--;
	      }
	    }
	    h5out.writeAttribute(fermion_path, "file_xml", file, wmode);
	    h5out.writeAttribute(fermion_path, "file_xml_formatted", file_formatted, wmode);
	    h5out.writeAttribute(fermion_path, "record_xml", record, wmode);
	    h5out.writeAttribute(fermion_path, "record_xml_formatted", record_formatted, wmode);
	  }
	}
	//The h5 file needs to be closed outside of this loop.
	h5out.cd("/");
	h5out.close();
      }
      
      //! Write the upper two components of a non-relativistic half-spin propagator
      void HDF5WriteUpperLatPropF(const std::string& buffer_id,
			   const std::string& outputfile,
			   const std::string& obj_name,
			   const std::string& path, HDF5Base::writemode wmode)
      {
	LatticePropagatorF obj;
	XMLBufferWriter file_xml, record_xml;

	obj = TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);

	//Check that this is a non-relativistic propagator before any writing.
	Double epsilon = 1e-8; //Maybe this should be more strict. 
	for(int color_source = 0; color_source < Nc ; ++color_source) 
	{
	  for(int spin_source = Ns/2; spin_source < Ns; ++spin_source) 
	  { 
	    int original_spin = spin_source - Ns/2;
	    LatticeFermionF psi;
	    LatticeFermionF chi;

	    PropToFerm(obj, psi, color_source, original_spin);
	    PropToFerm(obj, chi, color_source, spin_source);

	    QDPIO::cout<<"Checking color: "<<color_source<<" spin: "<<original_spin<<" against spin: "<<spin_source<<std::endl;
	    Double difference_norm = norm2(psi + chi);
	    QDPIO::cout<<"Difference norm is: "<<difference_norm<<std::endl;
	    if (toBool(difference_norm > epsilon))
	    {
	      QDPIO::cout<<"ERROR, based on checking the spin components, you are not really writing an upper non-relativistic propagator; exiting..."<<std::endl;
	      QDP_abort(1);
	    }
	  }
	}

	HDF5Writer h5out(outputfile);
	//h5out.push(path);
	std::string propagator_path = path+"/"+obj_name;
	h5out.push(propagator_path); //We want to push here since each spin component will be written separately.
	h5out.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB
	//Loop over the the upper spin components and write them.
	for(int color_source = 0; color_source < Nc ; ++color_source) 
	{
	  for(int spin_source = Ns/2; spin_source < Ns; ++spin_source) 
	  { 
	    int original_spin = spin_source - Ns/2;
	    LatticeFermionF psi;

	    PropToFerm(obj, psi, color_source, original_spin);

	    std::string fermion_path = propagator_path+"/color_"+std::to_string(color_source)+"_spin_"+std::to_string(original_spin);
	    QDPIO::cout<<"Writing color: "<<color_source<<" spin: "<<original_spin<<std::endl;
	    //Now all the usual stuff happens, only difference is we are writing a fermion.
	    h5out.write(fermion_path, psi, wmode);
	    std::string record = record_xml.str();
	    std::string file = file_xml.str();
	    //This string needs to be broadcasted, right now only the head node has it.
	    QDPInternal::broadcast_str(file);
	    QDPInternal::broadcast_str(record);
	    std::string file_formatted = file;
	    std::string record_formatted = record;
	    //This is cutting xml record/file_xml and formatting them to be only slightly more human readable.
	    for(int i = 0; i < record_formatted.length(); i++)
	    {
	      if(record_formatted[i] == '\n'|| record_formatted[i] == ' ' || record_formatted[i] == '\"') 
	      {
		record_formatted.erase(i, 1); 
		i--;
	      }
	    }
	    for(int i = 0; i < file_formatted.length(); i++)
	    {
	      if(file_formatted[i] == '\n'|| file_formatted[i] == ' ' || file_formatted[i] == '\"') 
	      {
		file_formatted.erase(i, 1); 
		i--;
	      }
	    }
	    h5out.writeAttribute(fermion_path, "file_xml", file, wmode);
	    h5out.writeAttribute(fermion_path, "file_xml_formatted", file_formatted, wmode);
	    h5out.writeAttribute(fermion_path, "record_xml", record, wmode);
	    h5out.writeAttribute(fermion_path, "record_xml_formatted", record_formatted, wmode);
	  }
	}
	//The h5 file needs to be closed outside of this loop.
	h5out.cd("/");
	h5out.close();
      }
      
      //! Write the upper two components of a non-relativistic half-spin propagator
      void HDF5WriteUpperLatPropD(const std::string& buffer_id,
			   const std::string& outputfile,
			   const std::string& obj_name,
			   const std::string& path, HDF5Base::writemode wmode)
      {
	LatticePropagatorD obj;
	XMLBufferWriter file_xml, record_xml;

	obj = TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);

	//Check that this is a non-relativistic propagator before any writing.
	Double epsilon = 1e-8; //Maybe this should be more strict. 
	for(int color_source = 0; color_source < Nc ; ++color_source) 
	{
	  for(int spin_source = Ns/2; spin_source < Ns; ++spin_source) 
	  { 
	    int original_spin = spin_source - Ns/2;
	    LatticeFermionD psi;
	    LatticeFermionD chi;

	    PropToFerm(obj, psi, color_source, original_spin);
	    PropToFerm(obj, chi, color_source, spin_source);

	    QDPIO::cout<<"Checking color: "<<color_source<<" spin: "<<original_spin<<" against spin: "<<spin_source<<std::endl;
	    Double difference_norm = norm2(psi + chi);
	    QDPIO::cout<<"Difference norm is: "<<difference_norm<<std::endl;
	    if (toBool(difference_norm > epsilon))
	    {
	      QDPIO::cout<<"ERROR, based on checking the spin components, you are not really writing an upper non-relativistic propagator; exiting..."<<std::endl;
	      QDP_abort(1);
	    }
	  }
	}

	HDF5Writer h5out(outputfile);
	//h5out.push(path);
	std::string propagator_path = path+"/"+obj_name;
	h5out.push(propagator_path); //We want to push here since each spin component will be written separately.
	h5out.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB
	//Loop over the the upper spin components and write them.
	for(int color_source = 0; color_source < Nc ; ++color_source) 
	{
	  for(int spin_source = Ns/2; spin_source < Ns; ++spin_source) 
	  { 
	    int original_spin = spin_source - Ns/2;
	    LatticeFermionD psi;

	    PropToFerm(obj, psi, color_source, original_spin);

	    std::string fermion_path = propagator_path+"/color_"+std::to_string(color_source)+"_spin_"+std::to_string(original_spin);
	    QDPIO::cout<<"Writing color: "<<color_source<<" spin: "<<original_spin<<std::endl;
	    //Now all the usual stuff happens, only difference is we are writing a fermion.
	    h5out.write(fermion_path, psi, wmode);
	    std::string record = record_xml.str();
	    std::string file = file_xml.str();
	    //This string needs to be broadcasted, right now only the head node has it.
	    QDPInternal::broadcast_str(file);
	    QDPInternal::broadcast_str(record);
	    std::string file_formatted = file;
	    std::string record_formatted = record;
	    //This is cutting xml record/file_xml and formatting them to be only slightly more human readable.
	    for(int i = 0; i < record_formatted.length(); i++)
	    {
	      if(record_formatted[i] == '\n'|| record_formatted[i] == ' ' || record_formatted[i] == '\"') 
	      {
		record_formatted.erase(i, 1); 
		i--;
	      }
	    }
	    for(int i = 0; i < file_formatted.length(); i++)
	    {
	      if(file_formatted[i] == '\n'|| file_formatted[i] == ' ' || file_formatted[i] == '\"') 
	      {
		file_formatted.erase(i, 1); 
		i--;
	      }
	    }
	    h5out.writeAttribute(fermion_path, "file_xml", file, wmode);
	    h5out.writeAttribute(fermion_path, "file_xml_formatted", file_formatted, wmode);
	    h5out.writeAttribute(fermion_path, "record_xml", record, wmode);
	    h5out.writeAttribute(fermion_path, "record_xml_formatted", record_formatted, wmode);
	  }
	}
	//The h5 file needs to be closed outside of this loop.
	h5out.cd("/");
	h5out.close();
      }
      
      //! Write the lower two components of a non-relativistic half-spin propagator
      void HDF5WriteLowerLatProp(const std::string& buffer_id,
			   const std::string& outputfile,
			   const std::string& obj_name,
			   const std::string& path, HDF5Base::writemode wmode)
      {
	LatticePropagator obj;
	XMLBufferWriter file_xml, record_xml;

	obj = TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);

	//Check that this is a non-relativistic propagator before any writing.
	Double epsilon = 1e-8; //Maybe this should be more strict. 
	for(int color_source = 0; color_source < Nc ; ++color_source) 
	{
	  for(int spin_source = 0; spin_source < Ns/2; ++spin_source) 
	  { 
	    int original_spin = spin_source + Ns/2;
	    LatticeFermion psi;
	    LatticeFermion chi;

	    PropToFerm(obj, psi, color_source, original_spin);
	    PropToFerm(obj, chi, color_source, spin_source);

	    QDPIO::cout<<"Checking color: "<<color_source<<" spin: "<<original_spin<<" against spin: "<<spin_source<<std::endl;
	    Double difference_norm = norm2(psi - chi);
	    QDPIO::cout<<"Difference norm is: "<<difference_norm<<std::endl;
	    if (toBool(difference_norm > epsilon))
	    {
	      QDPIO::cout<<"ERROR, based on checking the spin components, you are not really writing a lower non-relativistic propagator; exiting..."<<std::endl;
	      QDP_abort(1);
	    }
	  }
	}

	HDF5Writer h5out(outputfile);
	//h5out.push(path);
	std::string propagator_path = path+"/"+obj_name;
	h5out.push(propagator_path); //We want to push here since each spin component will be written separately.
	h5out.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB
	//Loop over the the upper spin components and write them.
	for(int color_source = 0; color_source < Nc ; ++color_source) 
	{
	  for(int spin_source = 0; spin_source < Ns/2; ++spin_source) 
	  { 
	    int original_spin = spin_source + Ns/2;
	    LatticeFermion psi;

	    PropToFerm(obj, psi, color_source, original_spin);

	    std::string fermion_path = propagator_path+"/color_"+std::to_string(color_source)+"_spin_"+std::to_string(original_spin);
	    QDPIO::cout<<"Writing color: "<<color_source<<" spin: "<<original_spin<<std::endl;
	    //Now all the usual stuff happens, only difference is we are writing a fermion.
	    h5out.write(fermion_path, psi, wmode);
	    std::string record = record_xml.str();
	    std::string file = file_xml.str();
	    //This string needs to be broadcasted, right now only the head node has it.
	    QDPInternal::broadcast_str(file);
	    QDPInternal::broadcast_str(record);
	    std::string file_formatted = file;
	    std::string record_formatted = record;
	    //This is cutting xml record/file_xml and formatting them to be only slightly more human readable.
	    for(int i = 0; i < record_formatted.length(); i++)
	    {
	      if(record_formatted[i] == '\n'|| record_formatted[i] == ' ' || record_formatted[i] == '\"') 
	      {
		record_formatted.erase(i, 1); 
		i--;
	      }
	    }
	    for(int i = 0; i < file_formatted.length(); i++)
	    {
	      if(file_formatted[i] == '\n'|| file_formatted[i] == ' ' || file_formatted[i] == '\"') 
	      {
		file_formatted.erase(i, 1); 
		i--;
	      }
	    }
	    h5out.writeAttribute(fermion_path, "file_xml", file, wmode);
	    h5out.writeAttribute(fermion_path, "file_xml_formatted", file_formatted, wmode);
	    h5out.writeAttribute(fermion_path, "record_xml", record, wmode);
	    h5out.writeAttribute(fermion_path, "record_xml_formatted", record_formatted, wmode);
	  }
	}
	//The h5 file needs to be closed outside of this loop.
	h5out.cd("/");
	h5out.close();
      }
      
      //! Write the lower two components of a non-relativistic half-spin propagator
      void HDF5WriteLowerLatPropF(const std::string& buffer_id,
			   const std::string& outputfile,
			   const std::string& obj_name,
			   const std::string& path, HDF5Base::writemode wmode)
      {
	LatticePropagatorF obj;
	XMLBufferWriter file_xml, record_xml;

	obj = TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);

	//Check that this is a non-relativistic propagator before any writing.
	Double epsilon = 1e-8; //Maybe this should be more strict. 
	for(int color_source = 0; color_source < Nc ; ++color_source) 
	{
	  for(int spin_source = 0; spin_source < Ns/2; ++spin_source) 
	  { 
	    int original_spin = spin_source + Ns/2;
	    LatticeFermionF psi;
	    LatticeFermionF chi;

	    PropToFerm(obj, psi, color_source, original_spin);
	    PropToFerm(obj, chi, color_source, spin_source);

	    QDPIO::cout<<"Checking color: "<<color_source<<" spin: "<<original_spin<<" against spin: "<<spin_source<<std::endl;
	    Double difference_norm = norm2(psi - chi);
	    QDPIO::cout<<"Difference norm is: "<<difference_norm<<std::endl;
	    if (toBool(difference_norm > epsilon))
	    {
	      QDPIO::cout<<"ERROR, based on checking the spin components, you are not really writing a lower non-relativistic propagator; exiting..."<<std::endl;
	      QDP_abort(1);
	    }
	  }
	}

	HDF5Writer h5out(outputfile);
	//h5out.push(path);
	std::string propagator_path = path+"/"+obj_name;
	h5out.push(propagator_path); //We want to push here since each spin component will be written separately.
	h5out.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB
	//Loop over the the upper spin components and write them.
	for(int color_source = 0; color_source < Nc ; ++color_source) 
	{
	  for(int spin_source = 0; spin_source < Ns/2; ++spin_source) 
	  { 
	    int original_spin = spin_source + Ns/2;
	    LatticeFermionF psi;

	    PropToFerm(obj, psi, color_source, original_spin);

	    std::string fermion_path = propagator_path+"/color_"+std::to_string(color_source)+"_spin_"+std::to_string(original_spin);
	    QDPIO::cout<<"Writing color: "<<color_source<<" spin: "<<original_spin<<std::endl;
	    //Now all the usual stuff happens, only difference is we are writing a fermion.
	    h5out.write(fermion_path, psi, wmode);
	    std::string record = record_xml.str();
	    std::string file = file_xml.str();
	    //This string needs to be broadcasted, right now only the head node has it.
	    QDPInternal::broadcast_str(file);
	    QDPInternal::broadcast_str(record);
	    std::string file_formatted = file;
	    std::string record_formatted = record;
	    //This is cutting xml record/file_xml and formatting them to be only slightly more human readable.
	    for(int i = 0; i < record_formatted.length(); i++)
	    {
	      if(record_formatted[i] == '\n'|| record_formatted[i] == ' ' || record_formatted[i] == '\"') 
	      {
		record_formatted.erase(i, 1); 
		i--;
	      }
	    }
	    for(int i = 0; i < file_formatted.length(); i++)
	    {
	      if(file_formatted[i] == '\n'|| file_formatted[i] == ' ' || file_formatted[i] == '\"') 
	      {
		file_formatted.erase(i, 1); 
		i--;
	      }
	    }
	    h5out.writeAttribute(fermion_path, "file_xml", file, wmode);
	    h5out.writeAttribute(fermion_path, "file_xml_formatted", file_formatted, wmode);
	    h5out.writeAttribute(fermion_path, "record_xml", record, wmode);
	    h5out.writeAttribute(fermion_path, "record_xml_formatted", record_formatted, wmode);
	  }
	}
	//The h5 file needs to be closed outside of this loop.
	h5out.cd("/");
	h5out.close();
      }
      
      //! Write the lower two components of a non-relativistic half-spin propagator
      void HDF5WriteLowerLatPropD(const std::string& buffer_id,
			   const std::string& outputfile,
			   const std::string& obj_name,
			   const std::string& path, HDF5Base::writemode wmode)
      {
	LatticePropagatorD obj;
	XMLBufferWriter file_xml, record_xml;

	obj = TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);

	//Check that this is a non-relativistic propagator before any writing.
	Double epsilon = 1e-8; //Maybe this should be more strict. 
	for(int color_source = 0; color_source < Nc ; ++color_source) 
	{
	  for(int spin_source = 0; spin_source < Ns/2; ++spin_source) 
	  { 
	    int original_spin = spin_source + Ns/2;
	    LatticeFermionD psi;
	    LatticeFermionD chi;

	    PropToFerm(obj, psi, color_source, original_spin);
	    PropToFerm(obj, chi, color_source, spin_source);

	    QDPIO::cout<<"Checking color: "<<color_source<<" spin: "<<original_spin<<" against spin: "<<spin_source<<std::endl;
	    Double difference_norm = norm2(psi - chi);
	    QDPIO::cout<<"Difference norm is: "<<difference_norm<<std::endl;
	    if (toBool(difference_norm > epsilon))
	    {
	      QDPIO::cout<<"ERROR, based on checking the spin components, you are not really writing a lower non-relativistic propagator; exiting..."<<std::endl;
	      QDP_abort(1);
	    }
	  }
	}

	HDF5Writer h5out(outputfile);
	//h5out.push(path);
	std::string propagator_path = path+"/"+obj_name;
	h5out.push(propagator_path); //We want to push here since each spin component will be written separately.
	h5out.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB
	//Loop over the the upper spin components and write them.
	for(int color_source = 0; color_source < Nc ; ++color_source) 
	{
	  for(int spin_source = 0; spin_source < Ns/2; ++spin_source) 
	  { 
	    int original_spin = spin_source + Ns/2;
	    LatticeFermionD psi;

	    PropToFerm(obj, psi, color_source, original_spin);

	    std::string fermion_path = propagator_path+"/color_"+std::to_string(color_source)+"_spin_"+std::to_string(original_spin);
	    QDPIO::cout<<"Writing color: "<<color_source<<" spin: "<<original_spin<<std::endl;
	    //Now all the usual stuff happens, only difference is we are writing a fermion.
	    h5out.write(fermion_path, psi, wmode);
	    std::string record = record_xml.str();
	    std::string file = file_xml.str();
	    //This string needs to be broadcasted, right now only the head node has it.
	    QDPInternal::broadcast_str(file);
	    QDPInternal::broadcast_str(record);
	    std::string file_formatted = file;
	    std::string record_formatted = record;
	    //This is cutting xml record/file_xml and formatting them to be only slightly more human readable.
	    for(int i = 0; i < record_formatted.length(); i++)
	    {
	      if(record_formatted[i] == '\n'|| record_formatted[i] == ' ' || record_formatted[i] == '\"') 
	      {
		record_formatted.erase(i, 1); 
		i--;
	      }
	    }
	    for(int i = 0; i < file_formatted.length(); i++)
	    {
	      if(file_formatted[i] == '\n'|| file_formatted[i] == ' ' || file_formatted[i] == '\"') 
	      {
		file_formatted.erase(i, 1); 
		i--;
	      }
	    }
	    h5out.writeAttribute(fermion_path, "file_xml", file, wmode);
	    h5out.writeAttribute(fermion_path, "file_xml_formatted", file_formatted, wmode);
	    h5out.writeAttribute(fermion_path, "record_xml", record, wmode);
	    h5out.writeAttribute(fermion_path, "record_xml_formatted", record_formatted, wmode);
	  }
	}
	//The h5 file needs to be closed outside of this loop.
	h5out.cd("/");
	h5out.close();
      }


      //! Write a fermion
      void HDF5WriteLatFerm(const std::string& buffer_id,
			   const std::string& outputfile,
			   const std::string& obj_name,
			   const std::string& path, HDF5Base::writemode wmode)
      {
	LatticeFermion obj;
	XMLBufferWriter file_xml, record_xml;

	obj = TheNamedObjMap::Instance().getData<LatticeFermion>(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);

	HDF5Writer h5out(outputfile);
	h5out.push(path);
	std::string propagator_path = path+"/"+obj_name;
	h5out.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB
	h5out.write(propagator_path, obj, wmode);
	std::string record = record_xml.str();
	std::string file = file_xml.str();
	//This string needs to be broadcasted, right now only the head node has it.
	QDPInternal::broadcast_str(file);
	QDPInternal::broadcast_str(record);
	std::string file_formatted = file;
	std::string record_formatted = record;
	//This is cutting xml record/file_xml and formatting them to be only slightly more human readable.
	for(int i = 0; i < record_formatted.length(); i++)
	{
	  if(record_formatted[i] == '\n'|| record_formatted[i] == ' ' || record_formatted[i] == '\"') 
	  {
	    record_formatted.erase(i, 1); 
	    i--;
	  }
	}
	for(int i = 0; i < file_formatted.length(); i++)
	{
	  if(file_formatted[i] == '\n'|| file_formatted[i] == ' ' || file_formatted[i] == '\"') 
	  {
	    file_formatted.erase(i, 1); 
	    i--;
	  }
	}
	//QDPIO::cout<<"This is the formatted xml_record that is going to be written: "<<std::endl<<record<<std::endl;
	h5out.writeAttribute(propagator_path, "file_xml", file, wmode);
	h5out.writeAttribute(propagator_path, "file_xml_formatted", file_formatted, wmode);
	h5out.writeAttribute(propagator_path, "record_xml", record, wmode);
	h5out.writeAttribute(propagator_path, "record_xml_formatted", record_formatted, wmode);

	h5out.cd("/");
	h5out.close();

      }


      //! Write a single prec fermion
      void HDF5WriteLatFermF(const std::string& buffer_id,
			   const std::string& outputfile,
			   const std::string& obj_name,
			   const std::string& path, HDF5Base::writemode wmode)
      {
	LatticeFermionF obj;
	XMLBufferWriter file_xml, record_xml;

	obj = TheNamedObjMap::Instance().getData<LatticeFermion>(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);

	HDF5Writer h5out(outputfile);
	h5out.push(path);
	std::string propagator_path = path+"/"+obj_name;
	h5out.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB
	h5out.write(propagator_path, obj, wmode);
	std::string record = record_xml.str();
	std::string file = file_xml.str();
	//This string needs to be broadcasted, right now only the head node has it.
	QDPInternal::broadcast_str(file);
	QDPInternal::broadcast_str(record);
	std::string file_formatted = file;
	std::string record_formatted = record;
	//This is cutting xml record/file_xml and formatting them to be only slightly more human readable.
	for(int i = 0; i < record_formatted.length(); i++)
	{
	  if(record_formatted[i] == '\n'|| record_formatted[i] == ' ' || record_formatted[i] == '\"') 
	  {
	    record_formatted.erase(i, 1); 
	    i--;
	  }
	}
	for(int i = 0; i < file_formatted.length(); i++)
	{
	  if(file_formatted[i] == '\n'|| file_formatted[i] == ' ' || file_formatted[i] == '\"') 
	  {
	    file_formatted.erase(i, 1); 
	    i--;
	  }
	}
	//QDPIO::cout<<"This is the formatted xml_record that is going to be written: "<<std::endl<<record<<std::endl;
	h5out.writeAttribute(propagator_path, "file_xml", file, wmode);
	h5out.writeAttribute(propagator_path, "file_xml_formatted", file_formatted, wmode);
	h5out.writeAttribute(propagator_path, "record_xml", record, wmode);
	h5out.writeAttribute(propagator_path, "record_xml_formatted", record_formatted, wmode);

	h5out.cd("/");
	h5out.close();

      }


      //! Write a double prec fermion
      void HDF5WriteLatFermD(const std::string& buffer_id,
			   const std::string& outputfile,
			   const std::string& obj_name,
			   const std::string& path, HDF5Base::writemode wmode)
      {
	LatticeFermionD obj;
	XMLBufferWriter file_xml, record_xml;

	obj = TheNamedObjMap::Instance().getData<LatticeFermion>(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);

	HDF5Writer h5out(outputfile);
	h5out.push(path);
	std::string propagator_path = path+"/"+obj_name;
	h5out.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB
	h5out.write(propagator_path, obj, wmode);
	std::string record = record_xml.str();
	std::string file = file_xml.str();
	//This string needs to be broadcasted, right now only the head node has it.
	QDPInternal::broadcast_str(file);
	QDPInternal::broadcast_str(record);
	std::string file_formatted = file;
	std::string record_formatted = record;
	//This is cutting xml record/file_xml and formatting them to be only slightly more human readable.
	for(int i = 0; i < record_formatted.length(); i++)
	{
	  if(record_formatted[i] == '\n'|| record_formatted[i] == ' ' || record_formatted[i] == '\"') 
	  {
	    record_formatted.erase(i, 1); 
	    i--;
	  }
	}
	for(int i = 0; i < file_formatted.length(); i++)
	{
	  if(file_formatted[i] == '\n'|| file_formatted[i] == ' ' || file_formatted[i] == '\"') 
	  {
	    file_formatted.erase(i, 1); 
	    i--;
	  }
	}
	//QDPIO::cout<<"This is the formatted xml_record that is going to be written: "<<std::endl<<record<<std::endl;
	h5out.writeAttribute(propagator_path, "file_xml", file, wmode);
	h5out.writeAttribute(propagator_path, "file_xml_formatted", file_formatted, wmode);
	h5out.writeAttribute(propagator_path, "record_xml", record, wmode);
	h5out.writeAttribute(propagator_path, "record_xml_formatted", record_formatted, wmode);

	h5out.cd("/");
	h5out.close();

      }

      //! Write a staggered propagator
      void HDF5WriteLatStagProp(const std::string& buffer_id,
			   const std::string& outputfile,
			   const std::string& obj_name,
			   const std::string& path, HDF5Base::writemode wmode)
      {
	LatticeStaggeredPropagator obj;
	XMLBufferWriter file_xml, record_xml;

	obj = TheNamedObjMap::Instance().getData<LatticeStaggeredPropagator>(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);

	HDF5Writer h5out(outputfile);
	h5out.push(path);
	std::string propagator_path = path+"/"+obj_name;
	h5out.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB
	h5out.write(propagator_path, obj, wmode);
	std::string record = record_xml.str();
	std::string file = file_xml.str();
	//This string needs to be broadcasted, right now only the head node has it.
	QDPInternal::broadcast_str(file);
	QDPInternal::broadcast_str(record);
	std::string file_formatted = file;
	std::string record_formatted = record;
	//This is cutting xml record/file_xml and formatting them to be only slightly more human readable.
	for(int i = 0; i < record_formatted.length(); i++)
	{
	  if(record_formatted[i] == '\n'|| record_formatted[i] == ' ' || record_formatted[i] == '\"') 
	  {
	    record_formatted.erase(i, 1); 
	    i--;
	  }
	}
	for(int i = 0; i < file_formatted.length(); i++)
	{
	  if(file_formatted[i] == '\n'|| file_formatted[i] == ' ' || file_formatted[i] == '\"') 
	  {
	    file_formatted.erase(i, 1); 
	    i--;
	  }
	}
	//QDPIO::cout<<"This is the formatted xml_record that is going to be written: "<<std::endl<<record<<std::endl;
	h5out.writeAttribute(propagator_path, "file_xml", file, wmode);
	h5out.writeAttribute(propagator_path, "file_xml_formatted", file_formatted, wmode);
	h5out.writeAttribute(propagator_path, "record_xml", record, wmode);
	h5out.writeAttribute(propagator_path, "record_xml_formatted", record_formatted, wmode);

	h5out.cd("/");
	h5out.close();

      }


      //! Write a single prec staggered propagator
      void HDF5WriteLatStagPropF(const std::string& buffer_id,
			   const std::string& outputfile,
			   const std::string& obj_name,
			   const std::string& path, HDF5Base::writemode wmode)
      {
	LatticeStaggeredPropagatorF obj;
	XMLBufferWriter file_xml, record_xml;

	obj = TheNamedObjMap::Instance().getData<LatticeStaggeredPropagator>(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);

	HDF5Writer h5out(outputfile);
	h5out.push(path);
	std::string propagator_path = path+"/"+obj_name;
	h5out.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB
	h5out.write(propagator_path, obj, wmode);
	std::string record = record_xml.str();
	std::string file = file_xml.str();
	//This string needs to be broadcasted, right now only the head node has it.
	QDPInternal::broadcast_str(file);
	QDPInternal::broadcast_str(record);
	std::string file_formatted = file;
	std::string record_formatted = record;
	//This is cutting xml record/file_xml and formatting them to be only slightly more human readable.
	for(int i = 0; i < record_formatted.length(); i++)
	{
	  if(record_formatted[i] == '\n'|| record_formatted[i] == ' ' || record_formatted[i] == '\"') 
	  {
	    record_formatted.erase(i, 1); 
	    i--;
	  }
	}
	for(int i = 0; i < file_formatted.length(); i++)
	{
	  if(file_formatted[i] == '\n'|| file_formatted[i] == ' ' || file_formatted[i] == '\"') 
	  {
	    file_formatted.erase(i, 1); 
	    i--;
	  }
	}
	//QDPIO::cout<<"This is the formatted xml_record that is going to be written: "<<std::endl<<record<<std::endl;
	h5out.writeAttribute(propagator_path, "file_xml", file, wmode);
	h5out.writeAttribute(propagator_path, "file_xml_formatted", file_formatted, wmode);
	h5out.writeAttribute(propagator_path, "record_xml", record, wmode);
	h5out.writeAttribute(propagator_path, "record_xml_formatted", record_formatted, wmode);

	h5out.cd("/");
	h5out.close();
   
      }

      //! Write a double prec staggered propagator
      void HDF5WriteLatStagPropD(const std::string& buffer_id,
			   const std::string& outputfile,
			   const std::string& obj_name,
			   const std::string& path, HDF5Base::writemode wmode)
      {
	LatticeStaggeredPropagatorD obj;
	XMLBufferWriter file_xml, record_xml;

	obj = TheNamedObjMap::Instance().getData<LatticeStaggeredPropagator>(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);

	HDF5Writer h5out(outputfile);
	h5out.push(path);
	std::string propagator_path = path+"/"+obj_name;
	h5out.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB
	h5out.write(propagator_path, obj, wmode);
	std::string record = record_xml.str();
	std::string file = file_xml.str();
	//This string needs to be broadcasted, right now only the head node has it.
	QDPInternal::broadcast_str(file);
	QDPInternal::broadcast_str(record);
	std::string file_formatted = file;
	std::string record_formatted = record;
	//This is cutting xml record/file_xml and formatting them to be only slightly more human readable.
	for(int i = 0; i < record_formatted.length(); i++)
	{
	  if(record_formatted[i] == '\n'|| record_formatted[i] == ' ' || record_formatted[i] == '\"') 
	  {
	    record_formatted.erase(i, 1); 
	    i--;
	  }
	}
	for(int i = 0; i < file_formatted.length(); i++)
	{
	  if(file_formatted[i] == '\n'|| file_formatted[i] == ' ' || file_formatted[i] == '\"') 
	  {
	    file_formatted.erase(i, 1); 
	    i--;
	  }
	}
	//QDPIO::cout<<"This is the formatted xml_record that is going to be written: "<<std::endl<<record<<std::endl;
	h5out.writeAttribute(propagator_path, "file_xml", file, wmode);
	h5out.writeAttribute(propagator_path, "file_xml_formatted", file_formatted, wmode);
	h5out.writeAttribute(propagator_path, "record_xml", record, wmode);
	h5out.writeAttribute(propagator_path, "record_xml_formatted", record_formatted, wmode);

	h5out.cd("/");
	h5out.close();
   
      }


      //! Write a gauge field
      void HDF5WriteArrayLatColMat(const std::string& buffer_id,
			   const std::string& outputfile,
			   const std::string& obj_name,
			   const std::string& path, HDF5Base::writemode wmode)
      {
	multi1d<LatticeColorMatrix> obj;
	XMLBufferWriter file_xml, record_xml;

	obj = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);

	HDF5Writer h5out(outputfile);
	h5out.push(path);
	std::string propagator_path = path+"/"+obj_name;
	h5out.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB
	h5out.write(propagator_path, obj, wmode);
	std::string record = record_xml.str();
	std::string file = file_xml.str();
	//This string needs to be broadcasted, right now only the head node has it.
	QDPInternal::broadcast_str(file);
	QDPInternal::broadcast_str(record);
	std::string file_formatted = file;
	std::string record_formatted = record;
	//This is cutting xml record/file_xml and formatting them to be only slightly more human readable.
	for(int i = 0; i < record_formatted.length(); i++)
	{
	  if(record_formatted[i] == '\n'|| record_formatted[i] == ' ' || record_formatted[i] == '\"') 
	  {
	    record_formatted.erase(i, 1); 
	    i--;
	  }
	}
	for(int i = 0; i < file_formatted.length(); i++)
	{
	  if(file_formatted[i] == '\n'|| file_formatted[i] == ' ' || file_formatted[i] == '\"') 
	  {
	    file_formatted.erase(i, 1); 
	    i--;
	  }
	}
	//QDPIO::cout<<"This is the formatted xml_record that is going to be written: "<<std::endl<<record<<std::endl;
	h5out.writeAttribute(propagator_path, "file_xml", file, wmode);
	h5out.writeAttribute(propagator_path, "file_xml_formatted", file_formatted, wmode);
	h5out.writeAttribute(propagator_path, "record_xml", record, wmode);
	h5out.writeAttribute(propagator_path, "record_xml_formatted", record_formatted, wmode);

	h5out.cd("/");
	h5out.close();
    
      }


      //! Write single prec a gauge field
      void HDF5WriteArrayLatColMatF(const std::string& buffer_id,
			   const std::string& outputfile,
			   const std::string& obj_name,
			   const std::string& path, HDF5Base::writemode wmode)
      {
	XMLBufferWriter file_xml, record_xml;

	multi1d<LatticeColorMatrix>& temp 
	  = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(buffer_id);
	multi1d<LatticeColorMatrixF> obj(temp.size());
	for(int mu=0; mu < obj.size(); ++mu)
	  obj[mu] = temp[mu];

	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);

	HDF5Writer h5out(outputfile);
	h5out.push(path);
	std::string propagator_path = path+"/"+obj_name;
	h5out.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB
	h5out.write(propagator_path, obj, wmode);
	std::string record = record_xml.str();
	std::string file = file_xml.str();
	//This string needs to be broadcasted, right now only the head node has it.
	QDPInternal::broadcast_str(file);
	QDPInternal::broadcast_str(record);
	std::string file_formatted = file;
	std::string record_formatted = record;
	//This is cutting xml record/file_xml and formatting them to be only slightly more human readable.
	for(int i = 0; i < record_formatted.length(); i++)
	{
	  if(record_formatted[i] == '\n'|| record_formatted[i] == ' ' || record_formatted[i] == '\"') 
	  {
	    record_formatted.erase(i, 1); 
	    i--;
	  }
	}
	for(int i = 0; i < file_formatted.length(); i++)
	{
	  if(file_formatted[i] == '\n'|| file_formatted[i] == ' ' || file_formatted[i] == '\"') 
	  {
	    file_formatted.erase(i, 1); 
	    i--;
	  }
	}
	//QDPIO::cout<<"This is the formatted xml_record that is going to be written: "<<std::endl<<record<<std::endl;
	h5out.writeAttribute(propagator_path, "file_xml", file, wmode);
	h5out.writeAttribute(propagator_path, "file_xml_formatted", file_formatted, wmode);
	h5out.writeAttribute(propagator_path, "record_xml", record, wmode);
	h5out.writeAttribute(propagator_path, "record_xml_formatted", record_formatted, wmode);

	h5out.cd("/");
	h5out.close();
    
      }

      //! Write a double prec gauge field
      void HDF5WriteArrayLatColMatD(const std::string& buffer_id,
			   const std::string& outputfile,
			   const std::string& obj_name,
			   const std::string& path, HDF5Base::writemode wmode)
      {
	XMLBufferWriter file_xml, record_xml;

	multi1d<LatticeColorMatrix>& temp 
	  = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(buffer_id);
	multi1d<LatticeColorMatrixD> obj(temp.size());
	for(int mu=0; mu < obj.size(); ++mu)
	  obj[mu] = temp[mu];

	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);

	HDF5Writer h5out(outputfile);
	h5out.push(path);
	std::string propagator_path = path+"/"+obj_name;
	h5out.set_stripesize(1048576); // MAGIC NUMBER ALERT!  1048576 = 1024^2 = 1MB
	h5out.write(propagator_path, obj, wmode);
	std::string record = record_xml.str();
	std::string file = file_xml.str();
	//This string needs to be broadcasted, right now only the head node has it.
	QDPInternal::broadcast_str(file);
	QDPInternal::broadcast_str(record);
	std::string file_formatted = file;
	std::string record_formatted = record;
	//This is cutting xml record/file_xml and formatting them to be only slightly more human readable.
	for(int i = 0; i < record_formatted.length(); i++)
	{
	  if(record_formatted[i] == '\n'|| record_formatted[i] == ' ' || record_formatted[i] == '\"') 
	  {
	    record_formatted.erase(i, 1); 
	    i--;
	  }
	}
	for(int i = 0; i < file_formatted.length(); i++)
	{
	  if(file_formatted[i] == '\n'|| file_formatted[i] == ' ' || file_formatted[i] == '\"') 
	  {
	    file_formatted.erase(i, 1); 
	    i--;
	  }
	}
	//QDPIO::cout<<"This is the formatted xml_record that is going to be written: "<<std::endl<<record<<std::endl;
	h5out.writeAttribute(propagator_path, "file_xml", file, wmode);
	h5out.writeAttribute(propagator_path, "file_xml_formatted", file_formatted, wmode);
	h5out.writeAttribute(propagator_path, "record_xml", record, wmode);
	h5out.writeAttribute(propagator_path, "record_xml_formatted", record_formatted, wmode);

	h5out.cd("/");
	h5out.close();
    
      }
      //! Local registration flag
      bool registered = false;

    }  // end namespace WriteObjCallMap


    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheHDF5WriteObjFuncMap::Instance().registerFunction(std::string("LatticePropagator"), 
								      HDF5WriteLatProp);
	success &= TheHDF5WriteObjFuncMap::Instance().registerFunction(std::string("LatticePropagatorF"), 
								      HDF5WriteLatPropF);
	success &= TheHDF5WriteObjFuncMap::Instance().registerFunction(std::string("LatticePropagatorD"), 
								      HDF5WriteLatPropD);
	success &= TheHDF5WriteObjFuncMap::Instance().registerFunction(std::string("LatticeUpperPropagator"), 
								      HDF5WriteUpperLatProp);
	success &= TheHDF5WriteObjFuncMap::Instance().registerFunction(std::string("LatticeUpperPropagatorF"), 
								      HDF5WriteUpperLatPropF);
	success &= TheHDF5WriteObjFuncMap::Instance().registerFunction(std::string("LatticeUpperPropagatorD"), 
								      HDF5WriteUpperLatPropD);
	success &= TheHDF5WriteObjFuncMap::Instance().registerFunction(std::string("LatticeLowerPropagator"), 
								      HDF5WriteLowerLatProp);
	success &= TheHDF5WriteObjFuncMap::Instance().registerFunction(std::string("LatticeLowerPropagatorF"), 
								      HDF5WriteLowerLatPropF);
	success &= TheHDF5WriteObjFuncMap::Instance().registerFunction(std::string("LatticeLowerPropagatorD"), 
								      HDF5WriteLowerLatPropD);
	success &= TheHDF5WriteObjFuncMap::Instance().registerFunction(std::string("LatticeFermion"), 
								      HDF5WriteLatFerm);
        success &= TheHDF5WriteObjFuncMap::Instance().registerFunction(std::string("LatticeFermionF"), 
								    HDF5WriteLatFermF);
        success &= TheHDF5WriteObjFuncMap::Instance().registerFunction(std::string("LatticeFermionD"), 
								 HDF5WriteLatFermD);
	success &= TheHDF5WriteObjFuncMap::Instance().registerFunction(std::string("LatticeStaggeredPropagator"), 
								      HDF5WriteLatStagProp);
	success &= TheHDF5WriteObjFuncMap::Instance().registerFunction(std::string("LatticeStaggeredPropagatorF"), 
								      HDF5WriteLatStagPropF);
	success &= TheHDF5WriteObjFuncMap::Instance().registerFunction(std::string("LatticeStaggeredPropagatorD"), 
								      HDF5WriteLatStagPropD);
	success &= TheHDF5WriteObjFuncMap::Instance().registerFunction(std::string("Multi1dLatticeColorMatrix"), 
								      HDF5WriteArrayLatColMat);
	success &= TheHDF5WriteObjFuncMap::Instance().registerFunction(std::string("Multi1dLatticeColorMatrixF"), 
								      HDF5WriteArrayLatColMatF);
	success &= TheHDF5WriteObjFuncMap::Instance().registerFunction(std::string("Multi1dLatticeColorMatrixD"), 
								      HDF5WriteArrayLatColMatD);

	registered = true;
      }
      return success;
    }
  }

}

#endif
