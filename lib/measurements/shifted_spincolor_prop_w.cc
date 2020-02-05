/*
Authors
Arjun Gambhir

INPUT
    Source
OUTPUT
    Propagator solved by shifting spin/color indices of first solution.
*/

// Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "fermact.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "util/info/unique_id.h"
//We need this to insert Fermions to a Prop, this is imperitive for dilution.
#include "util/ferm/transf.h"

// Lalibe Stuff
#include "shifted_spincolor_prop_w.h"

namespace Chroma
{
    namespace LalibeShiftedSpinColorPropagatorEnv
    {
        namespace
        {
            AbsInlineMeasurement* createMeasurement(XMLReader& xml_in,
                const std::string& path)
            {
                return new InlineMeas(ShiftedSpinColorParams(xml_in, path));
            }
            //! Local registration flag
            bool registered = false;
        }
        const std::string name = "SHIFTED_SPINCOLOR_PROPAGATOR";

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

        void read(XMLReader& xml, const std::string& path, ShiftedSpinColorParams::ShiftedSpinColorProp_t& par)
        {
            XMLReader paramtop(xml, path);
	    read(paramtop, "PropagatorParam" ,par.prop_param ); //params for next lin solve
        }

        void write(XMLWriter& xml, const std::string& path, ShiftedSpinColorParams::ShiftedSpinColorProp_t& par)
        {
            push(xml, path);
	    write(xml, "PropagatorParam" ,par.prop_param); //params for next lin solve
            pop(xml);

        }

        //! NamedObject input
        void read(XMLReader& xml, const std::string& path, ShiftedSpinColorParams::NamedObject_t& input)
        {
            XMLReader inputtop(xml, path);
            read(inputtop, "gauge_id"     , input.gauge_id);
            read(inputtop, "source_id"  , input.src_id);
            read(inputtop, "prop_id"   , input.prop_id);
        }

        //! NamedObject output
        void write(XMLWriter& xml, const std::string& path, const ShiftedSpinColorParams::NamedObject_t& input)
        {
            push(xml, path);
            write(xml, "gauge_id"     , input.gauge_id    );
            write(xml, "source_id"  , input.src_id     );
            write(xml, "prop_id"   , input.prop_id);
            pop(xml);
        }

        // Param stuff
        ShiftedSpinColorParams::ShiftedSpinColorParams()
        {
            frequency = 0;
        }

        ShiftedSpinColorParams::ShiftedSpinColorParams(XMLReader& xml_in, const std::string& path)
        {
            try
            {
                XMLReader paramtop(xml_in, path);
                if (paramtop.count("Frequency") == 1)
                    read(paramtop, "Frequency", frequency);
                else
                    frequency = 1;

                // Parameters for source construction
                read(paramtop, "ShiftedSpinColorParams", shiftedscparam);

                // Read in the NamedObject info
                read(paramtop, "NamedObject", named_obj);
            }
            catch(const std::string& e)
            {
                QDPIO::cerr << __func__ << ": Caught Exception reading XML: "
                                                            << e << std::endl;
                QDP_abort(1);
            }
        }

        void ShiftedSpinColorParams::writeXML(XMLWriter& xml_out, const std::string& path)
        {
            push(xml_out, path);
            write(xml_out, "ShiftedSpinColorParams", shiftedscparam);
            write(xml_out, "NamedObject", named_obj);
            pop(xml_out);
        }

        void  InlineMeas::operator()(unsigned long update_no, XMLWriter& xml_out)
        {
            START_CODE();

            StopWatch snoop;
            snoop.reset();
            snoop.start();
            QDPIO::cout << "SHIFTED_SPINCOLOR_PROPAGATOR: start" << std::endl;

            // Test and grab a reference to the gauge field
            XMLBufferWriter gauge_xml;
            try
            {
                TheNamedObjMap::Instance().getData
                    <multi1d <LatticeColorMatrix> >(params.named_obj.gauge_id);
                TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
            }
            catch( std::bad_cast )
            {
                QDPIO::cerr << LalibeShiftedSpinColorPropagatorEnv::name
                    << ": caught dynamic cast error" << std::endl;
                QDP_abort(1);
            }
            catch (const std::string& e)
            {
                QDPIO::cerr << LalibeShiftedSpinColorPropagatorEnv::name
                    << ": map call failed: " << e << std::endl;
                QDP_abort(1);
            }
            const multi1d<LatticeColorMatrix>& u =
                TheNamedObjMap::Instance().getData
                    <multi1d <LatticeColorMatrix> >(params.named_obj.gauge_id);

	    // Read src quark propagator
            XMLReader src_file_xml, src_record_xml;
            LatticePropagator quark_source;
	    int t0;
	    int j_decay;
	    // Record the type of header
	    bool make_sourceP = false;
	    bool seqsourceP = false;
	    //Need origin for fourier transform!
	    multi1d<int> origin;
            QDPIO::cout << "Attempt to read forward propagator" << std::endl;
            try
            {
                quark_source = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.src_id);
                TheNamedObjMap::Instance().get(params.named_obj.src_id).getFileXML(src_file_xml);
                TheNamedObjMap::Instance().get(params.named_obj.src_id).getRecordXML(src_record_xml);
		// Try to invert this record XML into a source struct
		// First identify what kind of source might be here
		if (src_record_xml.count("/MakeSource") != 0)
		{
		  make_sourceP = true;
		  MakeSourceProp_t  orig_header;
		  read(src_record_xml, "/MakeSource", orig_header);

		  j_decay = orig_header.source_header.j_decay;
		  t0      = orig_header.source_header.t_source;
		}
		else if (src_record_xml.count("/SequentialSource") != 0)
		{
		  seqsourceP = true;
		  SequentialSource_t   orig_header;
		  read(src_record_xml, "/SequentialSource", orig_header);

		  j_decay = orig_header.seqsource_header.j_decay;
		  t0      = orig_header.seqsource_header.t_sink;   // funky, but probably not really needed
		}
		else
		{
		  throw std::string("No appropriate header found");
		}
	    }
	    catch( std::bad_cast )
	    {
		QDPIO::cerr << LalibeShiftedSpinColorPropagatorEnv::name
		    << ": caught dynamic cast error" << std::endl;
		QDP_abort(1);
	    }
	    catch (const std::string& e)
	    {
		QDPIO::cerr << LalibeShiftedSpinColorPropagatorEnv::name
		    << ": map call failed: " << e << std::endl;
		QDP_abort(1);
	    }
	    
	    //
	    // Initialize fermion action
	    //
	    std::istringstream  xml_s(params.shiftedscparam.prop_param.fermact.xml);
	    XMLReader  fermacttop(xml_s);
	    QDPIO::cout << "FermAct = " << params.shiftedscparam.prop_param.fermact.id << std::endl;


            StopWatch swatch;
            swatch.reset();
            QDPIO::cout << "Try the various factories" << std::endl;

            // Typedefs to save typing
            typedef LatticeFermion               T;
            typedef multi1d<LatticeColorMatrix>  P;
            typedef multi1d<LatticeColorMatrix>  Q;

            // Generic Wilson-Type stuff
            Handle< FermionAction<T,P,Q> >
	    S_f(TheFermionActionFactory::Instance().createObject(params.shiftedscparam.prop_param.fermact.id,
							         fermacttop,
							         params.shiftedscparam.prop_param.fermact.path));

           Handle< FermState<T,P,Q> > state(S_f->createState(u));

           Handle< SystemSolver<LatticeFermion> > PP = S_f->qprop(state,
							          params.shiftedscparam.prop_param.invParam);
      
           QDPIO::cout << "Suitable factory found: compute the quark prop" << std::endl;

           int ncg_had = 0;
           LatticeFermion init_guess = zero;
           LatticePropagator quark_soln = zero;

           for(int color_source(0);color_source<Nc;color_source++){
             for(int spin_source=0; spin_source < Ns; spin_source++){
               StopWatch swatch;
               swatch.reset();
               LatticeFermion src_ferm = zero;
               LatticeFermion soln_ferm = init_guess;
               PropToFerm(quark_source, src_ferm, color_source, spin_source);
               //Do the solve.
               swatch.start();

               SystemSolverResults_t res = (*PP)(soln_ferm, src_ferm);
               ncg_had = res.n_count;

               swatch.stop();
               QDPIO::cout << "Propagator computed: time= " 
             	      << swatch.getTimeInSeconds() 
             	      << " secs" << std::endl;
	       if (color_source == 0 && spin_source == 0)
		 init_guess = soln_ferm; 
               FermToProp(soln_ferm, quark_soln, color_source, spin_source);
             }
           }
	      
           // Save the propagator info
           try
           {
             QDPIO::cout << "Start writing propagator info" << std::endl;

             XMLBufferWriter file_xml;
             push(file_xml, "propagator");
             write(file_xml, "id", uniqueId());  // NOTE: new ID form
             pop(file_xml);

             XMLBufferWriter record_xml;
             if (make_sourceP)
             {
               MakeSourceProp_t  orig_header;
               read(src_record_xml, "/MakeSource", orig_header);

               Propagator_t  new_header;   // note, abandoning state_info
               new_header.prop_header   = params.shiftedscparam.prop_param;
               new_header.source_header = orig_header.source_header;
               new_header.gauge_header  = orig_header.gauge_header;
               write(record_xml, "Propagator", new_header);  
             } 
             else if (seqsourceP)
             {
               SequentialSource_t  orig_header;
               read(src_record_xml, "/SequentialSource", orig_header);

               SequentialProp_t  new_header;   // note, abandoning state_info
               new_header.seqprop_header   = params.shiftedscparam.prop_param;
               new_header.sink_header      = orig_header.sink_header;
               new_header.seqsource_header = orig_header.seqsource_header;
               new_header.forward_props    = orig_header.forward_props;
               new_header.gauge_header     = orig_header.gauge_header;
               write(record_xml, "SequentialProp", new_header);  
             }

	     LatticePropagator& quark_propagator = 
	       TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_id);
             // Write the propagator xml info
             TheNamedObjMap::Instance().get(params.named_obj.prop_id).setFileXML(file_xml);
             TheNamedObjMap::Instance().get(params.named_obj.prop_id).setRecordXML(record_xml);

             QDPIO::cout << "Propagator successfully updated" << std::endl;
           }
           catch (std::bad_cast)
           {
             QDPIO::cerr << LalibeShiftedSpinColorPropagatorEnv::name << ": caught dynamic cast error" 
             	    << std::endl;
             QDP_abort(1);
           }
           catch (const std::string& e) 
           {
             QDPIO::cerr << LalibeShiftedSpinColorPropagatorEnv::name << ": error extracting prop_header: " << e << std::endl;
             QDP_abort(1);
           }

           snoop.stop();
           QDPIO::cout << LalibeShiftedSpinColorPropagatorEnv::name << ": total time = " << snoop.getTimeInSeconds() << " secs" << std::endl;
           QDPIO::cout << LalibeShiftedSpinColorPropagatorEnv::name<< ": ran successfully" << std::endl;
           END_CODE();

	}
    }// LalibeShiftedSpinColorPropagatorEnv
  };
