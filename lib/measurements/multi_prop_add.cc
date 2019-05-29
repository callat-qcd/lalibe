/*
Authors
Arjun Gambhir
Andre Walker-Loud

Add multiple propagators
INPUT
    List of Propagators
    optional List of weights
    - if no weights are provide, do simple average
OUTPUT
    added propagator
*/

// Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
//#include "util/info/unique_id.h"

// Lalibe Stuff
#include "multi_prop_add.h"

namespace Chroma
{
    namespace LalibeMultiPropagatorAddEnv
    {
        namespace
        {
            AbsInlineMeasurement* createMeasurement(XMLReader& xml_in,
                const std::string& path)
            {
                return new InlineMeas(PropWeights(xml_in, path));
            }
            //! Local registration flag
            bool registered = false;
        }
        const std::string name = "MULTI_PROP_ADD";

        //! Register all the factories
        bool registerAll()
        {
            bool success = true;
            if (! registered)
            {
                success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
                registered = true;
            }
            return success;
        }

        void read(XMLReader& xml, const std::string& path, PropWeights::weights_t& par)
        {
            XMLReader paramtop(xml, path);
            if (paramtop.count("weights") != 0)
            {
                read(paramtop, "weights", par.weights); // list of weights
                par.have_weights = true;
  	            QDPIO::cout<<"Using user specified weights "<<std::endl;
            }
            else
            {
                par.have_weights = false;
                QDPIO::cout<<"No weights specified: doing simple average with 1/N_props "<<std::endl;
            }
	    if (paramtop.count("delete_props") != 0)
	    {
		read(paramtop, "delete_props", par.delete_props); // list of weights
	        QDPIO::cout<<"User specified to turn delete_props to "<<par.delete_props<<std::endl;
	    }
	    else
	    {
		par.delete_props = false;
		QDPIO::cout<<"By default, no props will be deleted after sum "<<std::endl;
	    }
        }

        void write(XMLWriter& xml, const std::string& path, PropWeights::weights_t& par)
        {
            push(xml, path);
            write(xml, "weights" ,par.weights); //list of weights
	        write(xml, "delete_props" ,par.delete_props);
            pop(xml);
        }

        //! NamedObject input
        void read(XMLReader& xml, const std::string& path, PropWeights::NamedObject_t& input)
        {
            XMLReader inputtop(xml, path);
            read(inputtop, "prop_ids", input.prop_ids);
            read(inputtop, "result_prop", input.result_prop);
        }

        //! NamedObject output
        void write(XMLWriter& xml, const std::string& path, const PropWeights::NamedObject_t& input)
        {
            push(xml, path);
            write(xml, "prop_ids", input.prop_ids);
            write(xml, "result_prop", input.result_prop);
            pop(xml);
        }

        // Param stuff
        PropWeights::PropWeights()
        {
            frequency = 0;
        }

        PropWeights::PropWeights(XMLReader& xml_in, const std::string& path)
        {
            try
            {
                XMLReader paramtop(xml_in, path);
                if (paramtop.count("Frequency") == 1)
                    read(paramtop, "Frequency", frequency);
                else
                    frequency = 1;

                // Read weights
                read(paramtop, "PropWeights", weights_lst);

                // Read in the NamedObject info
                read(paramtop, "NamedObject", named_obj);
            }
            catch(const std::string& e)
            {
                QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
                QDP_abort(1);
            }
        }

        void PropWeights::writeXML(XMLWriter& xml_out, const std::string& path)
        {
            push(xml_out, path);
            write(xml_out, "PropWeights", weights_lst);
            write(xml_out, "NamedObject", named_obj);
            pop(xml_out);
        }

        // Function call
        void  InlineMeas::operator()(unsigned long update_no, XMLWriter& xml_out)
        {
            START_CODE();

            StopWatch snoop;
            snoop.reset();
            snoop.start();
            QDPIO::cout << "MULTI_PROP_ADD: start" << std::endl;
            QDPIO::cout << "  ASSUMING User EITHER has all props from the same source OR has gauged fixed the lattice..." << std::endl;

            int N_props = params.named_obj.prop_ids.size();
            QDPIO::cout << "MULTI_PROP_ADD: N_props " << N_props << std::endl;
            if (params.weights_lst.have_weights)
            {
                int N_weights = params.weights_lst.weights.size();
                if (N_weights != N_props)
                {
                    QDPIO::cout << "MULTI_PROP_ADD ERROR: N_weights must = N_props" << std::endl;
                    QDP_abort(1);
                }
            }
            else
            {
                params.weights_lst.weights.resize(N_props);
                for(int ni = 0; ni < N_props; ni++)
                    params.weights_lst.weights[ni] = 1./N_props;
            }

            // Instantiate xml stuff and LatticePropagator; Grab xml record from first prop being added
            XMLReader file_xml, record_xml;
            LatticePropagator prop_to_add;
            try
            {
	            // Try the cast to see if this is a valid source
	            prop_to_add = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_ids[0]);

	            TheNamedObjMap::Instance().get(params.named_obj.prop_ids[0]).getFileXML(file_xml);
                TheNamedObjMap::Instance().get(params.named_obj.prop_ids[0]).getRecordXML(record_xml);
            }
            catch (std::bad_cast)
            {
	            QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
	            QDP_abort(1);
            }
            catch (const std::string& e)
            {
	           QDPIO::cerr << name << ": error extracting source_header: " << e << std::endl;
	           QDP_abort(1);
            }
            // Start adding props to result
            LatticePropagator multi_prop_result = prop_to_add * params.weights_lst.weights[0];
            for(int ni = 1; ni < N_props; ni++)
            {
                try
                {
                    prop_to_add = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_ids[ni]);
                    // ADD to result prop
                    multi_prop_result += prop_to_add * params.weights_lst.weights[ni];
                }
                catch (std::bad_cast)
                {
    	            QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
    	            QDP_abort(1);
                }
                catch (const std::string& e)
                {
    	           QDPIO::cerr << name << ": error extracting source_header: " << e << std::endl;
    	           QDP_abort(1);
                }
		if (params.weights_lst.delete_props)
		{
		    // Deleting the object.
		    TheNamedObjMap::Instance().erase(params.named_obj.prop_ids[ni]);
		}
            }

            /*
            *  Write the a source out to a named buffer
            */
            try
            {
                QDPIO::cout << "Attempt to store result prop" << std::endl;

                // Store the seqsource
                TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.result_prop);
                TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.result_prop) = multi_prop_result ;
                TheNamedObjMap::Instance().get(params.named_obj.result_prop).setFileXML(file_xml);
                TheNamedObjMap::Instance().get(params.named_obj.result_prop).setRecordXML(record_xml);

                QDPIO::cout << "Propagator sum successfully stored"  << std::endl;
            }
            catch (std::bad_cast)
            {
                QDPIO::cerr << name << ": dynamic cast error" << std::endl;
                QDP_abort(1);
            }
            catch (const std::string& e)
            {
                QDPIO::cerr << name << ": error storing seqsource: " << e << std::endl;
                QDP_abort(1);
            }

	    //This was causing problems, since nothing is written here.
            //pop(xml_out);   // multi_prop_add

	    //Redundant
            //QDPIO::cout << "MULTI_PROP_ADD ran successfully" << std::endl;

            snoop.stop();
            QDPIO::cout << LalibeMultiPropagatorAddEnv::name << ": total time = " << snoop.getTimeInSeconds() << " secs" << std::endl;
            QDPIO::cout << LalibeMultiPropagatorAddEnv::name<< ": ran successfully" << std::endl;
            END_CODE();
        }
    }// LalibeMultiPropagatorAddEnv
};
