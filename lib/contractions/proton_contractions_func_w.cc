/*!
 *  Functions that do baryon spin and color contractions.
 *  Color is nested inside of spin.
 *  Authors:
 *  Arjun Gambhir
 */

#include "proton_contractions_func_w.h"
#include "chromabase.h"
#include "util/ferm/diractodr.h"


namespace Chroma
{
namespace LegacyProton {

    void rotate_to_Dirac_Basis(LatticePropagator & quark_to_be_rotated)
    {
        //I am lazy so I copy Robert to rotate this stuff...
        SpinMatrix U = DiracToDRMat();
        quark_to_be_rotated = adj(U)*quark_to_be_rotated*U ;
    }

    void get_spin_wavefunctions(
        multi2d<int> & src_spins,
        multi2d<int> & snk_spins,
        multi1d<Real> & src_weights,
        multi1d<Real> & snk_weights,
        std::string baryon_name,
        std::string spin,
        int parity
    )
    {
        //Names should all be self-explanatory, parity=0 is positive parity, while parity=1 is negative.
        if (baryon_name == "proton")
        {
            src_spins.resize(Ns/2, Nc);
            snk_spins.resize(Ns, Nc);
            src_weights.resize(Ns/2);
            snk_weights.resize(Ns);
            src_weights[0] = 1/sqrt(2); src_weights[1] = -1/sqrt(2);
            snk_weights[0] = 1/sqrt(2); snk_weights[1] = -1/sqrt(2);
            snk_weights[2] = 1/sqrt(2); snk_weights[3] = -1/sqrt(2);
            if(parity == 0)
            {
                if(spin == "up")
                {
                    //spin up to spin up piece:
                    QDPIO::cout<<"Computing spin-up to spin-up piece..."<<std::endl;
                    //source spins:
                    src_spins[0][0] = 0; src_spins[0][1] = 0; src_spins[0][2] = 1;
                    src_spins[1][0] = 0; src_spins[1][1] = 1; src_spins[1][2] = 0;
	                //sink spins:
                    snk_spins[0][0] = 0; snk_spins[0][1] = 0; snk_spins[0][2] = 1;
                    snk_spins[1][0] = 0; snk_spins[1][1] = 1; snk_spins[1][2] = 0;
                    snk_spins[2][0] = 0; snk_spins[2][1] = 0; snk_spins[2][2] = 1;
                    snk_spins[3][0] = 1; snk_spins[3][1] = 0; snk_spins[3][2] = 0;
                }
                else if(spin == "dn")
                {
                    QDPIO::cout<<"Computing spin-down to spin-down piece..."<<std::endl;
                    //spin down to spin down piece:
                    //sink spins:
                    src_spins[0][0] = 1; src_spins[0][1] = 0; src_spins[0][2] = 1;
                    src_spins[1][0] = 1; src_spins[1][1] = 1; src_spins[1][2] = 0;
                    //source spins:
                    snk_spins[0][0] = 1; snk_spins[0][1] = 0; snk_spins[0][2] = 1;
                    snk_spins[1][0] = 1; snk_spins[1][1] = 1; snk_spins[1][2] = 0;
                    snk_spins[2][0] = 0; snk_spins[2][1] = 1; snk_spins[2][2] = 1;
                    snk_spins[3][0] = 1; snk_spins[3][1] = 1; snk_spins[3][2] = 0;
                }
            }
            else if(parity == 1)
            {
                if(spin == "up")
                {
                    //spin up to spin up piece:
                    QDPIO::cout<<"Computing negative parity spin-up to spin-up piece..."<<std::endl;
                    //source spins:
                    src_spins[0][0] = 2; src_spins[0][1] = 2; src_spins[0][2] = 3;
                    src_spins[1][0] = 2; src_spins[1][1] = 3; src_spins[1][2] = 2;
                    //sink spins:
                    snk_spins[0][0] = 2, snk_spins[0][1] = 2, snk_spins[0][2] = 3;
                    snk_spins[1][0] = 2, snk_spins[1][1] = 3, snk_spins[1][2] = 2;
                    snk_spins[2][0] = 2, snk_spins[2][1] = 2, snk_spins[2][2] = 3;
                    snk_spins[3][0] = 3, snk_spins[3][1] = 2, snk_spins[3][2] = 2;
                }
                else if(spin == "dn")
                {
                    QDPIO::cout<<"Computing negative parity spin-down to spin-down piece..."<<std::endl;
                    //spin down to spin down piece:
                    //sink spins:
                    src_spins[0][0] = 3; src_spins[0][1] = 2; src_spins[0][2] = 3;
                    src_spins[1][0] = 3; src_spins[1][1] = 3; src_spins[1][2] = 2;
                    //source spins:
                    snk_spins[0][0] = 3; snk_spins[0][1] = 2; snk_spins[0][2] = 3;
                    snk_spins[1][0] = 3; snk_spins[1][1] = 3; snk_spins[1][2] = 2;
                    snk_spins[2][0] = 2; snk_spins[2][1] = 3; snk_spins[2][2] = 3;
                    snk_spins[3][0] = 3; snk_spins[3][1] = 3; snk_spins[3][2] = 2;
                }
            }
        }
    }

    void color_contraction(
        LatticeColorMatrix & quark_1,
        LatticeColorMatrix & quark_2,
        LatticeColorMatrix & quark_3,
        LatticeComplex & spin_color_contracted_thing
    )
    {
        //This does the double epsilon color contractions, there are 6 combinations per epislon, so 36 total.
        //Antisymmetry makes half of these terms negative, positives ones below.
        spin_color_contracted_thing  = peekColor(quark_1,0,0)*peekColor(quark_2,1,1)*peekColor(quark_3,2,2);
        spin_color_contracted_thing += peekColor(quark_1,0,0)*peekColor(quark_2,2,2)*peekColor(quark_3,1,1);
        spin_color_contracted_thing += peekColor(quark_1,0,1)*peekColor(quark_2,1,2)*peekColor(quark_3,2,0);
        spin_color_contracted_thing += peekColor(quark_1,0,1)*peekColor(quark_2,2,0)*peekColor(quark_3,1,2);
        spin_color_contracted_thing += peekColor(quark_1,0,2)*peekColor(quark_2,1,0)*peekColor(quark_3,2,1);
        spin_color_contracted_thing += peekColor(quark_1,0,2)*peekColor(quark_2,2,1)*peekColor(quark_3,1,0);
        spin_color_contracted_thing += peekColor(quark_1,1,0)*peekColor(quark_2,0,2)*peekColor(quark_3,2,1);
        spin_color_contracted_thing += peekColor(quark_1,1,0)*peekColor(quark_2,2,1)*peekColor(quark_3,0,2);
        spin_color_contracted_thing += peekColor(quark_1,1,1)*peekColor(quark_2,0,0)*peekColor(quark_3,2,2);
        spin_color_contracted_thing += peekColor(quark_1,1,1)*peekColor(quark_2,2,2)*peekColor(quark_3,0,0);
        spin_color_contracted_thing += peekColor(quark_1,1,2)*peekColor(quark_2,0,1)*peekColor(quark_3,2,0);
        spin_color_contracted_thing += peekColor(quark_1,1,2)*peekColor(quark_2,2,0)*peekColor(quark_3,0,1);
        spin_color_contracted_thing += peekColor(quark_1,2,0)*peekColor(quark_2,0,1)*peekColor(quark_3,1,2);
        spin_color_contracted_thing += peekColor(quark_1,2,0)*peekColor(quark_2,1,2)*peekColor(quark_3,0,1);
        spin_color_contracted_thing += peekColor(quark_1,2,1)*peekColor(quark_2,0,2)*peekColor(quark_3,1,0);
        spin_color_contracted_thing += peekColor(quark_1,2,1)*peekColor(quark_2,1,0)*peekColor(quark_3,0,2);
        spin_color_contracted_thing += peekColor(quark_1,2,2)*peekColor(quark_2,0,0)*peekColor(quark_3,1,1);
        spin_color_contracted_thing += peekColor(quark_1,2,2)*peekColor(quark_2,1,1)*peekColor(quark_3,0,0);
        //Negative terms below.
        spin_color_contracted_thing -= peekColor(quark_1,0,0)*peekColor(quark_2,1,2)*peekColor(quark_3,2,1);
        spin_color_contracted_thing -= peekColor(quark_1,0,0)*peekColor(quark_2,2,1)*peekColor(quark_3,1,2);
        spin_color_contracted_thing -= peekColor(quark_1,0,1)*peekColor(quark_2,1,0)*peekColor(quark_3,2,2);
        spin_color_contracted_thing -= peekColor(quark_1,0,1)*peekColor(quark_2,2,2)*peekColor(quark_3,1,0);
        spin_color_contracted_thing -= peekColor(quark_1,0,2)*peekColor(quark_2,1,1)*peekColor(quark_3,2,0);
        spin_color_contracted_thing -= peekColor(quark_1,0,2)*peekColor(quark_2,2,0)*peekColor(quark_3,1,1);
        spin_color_contracted_thing -= peekColor(quark_1,1,0)*peekColor(quark_2,0,1)*peekColor(quark_3,2,2);
        spin_color_contracted_thing -= peekColor(quark_1,1,0)*peekColor(quark_2,2,2)*peekColor(quark_3,0,1);
        spin_color_contracted_thing -= peekColor(quark_1,1,1)*peekColor(quark_2,0,2)*peekColor(quark_3,2,0);
        spin_color_contracted_thing -= peekColor(quark_1,1,1)*peekColor(quark_2,2,0)*peekColor(quark_3,0,2);
        spin_color_contracted_thing -= peekColor(quark_1,1,2)*peekColor(quark_2,0,0)*peekColor(quark_3,2,1);
        spin_color_contracted_thing -= peekColor(quark_1,2,0)*peekColor(quark_2,0,2)*peekColor(quark_3,1,1);
        spin_color_contracted_thing -= peekColor(quark_1,2,0)*peekColor(quark_2,1,1)*peekColor(quark_3,0,2);
        spin_color_contracted_thing -= peekColor(quark_1,1,2)*peekColor(quark_2,2,1)*peekColor(quark_3,0,0);
        spin_color_contracted_thing -= peekColor(quark_1,2,1)*peekColor(quark_2,0,0)*peekColor(quark_3,1,2);
        spin_color_contracted_thing -= peekColor(quark_1,2,1)*peekColor(quark_2,1,2)*peekColor(quark_3,0,0);
        spin_color_contracted_thing -= peekColor(quark_1,2,2)*peekColor(quark_2,0,1)*peekColor(quark_3,1,0);
        spin_color_contracted_thing -= peekColor(quark_1,2,2)*peekColor(quark_2,1,0)*peekColor(quark_3,0,1);
    }

    void spin_contraction(
        LatticePropagator & quark_1,
        LatticePropagator & quark_2,
        LatticePropagator & quark_3,
        multi2d<int> & src_spins,
        multi2d<int> & snk_spins,
        multi1d<Real> & src_weights,
        multi1d<Real> & snk_weights,
        LatticeComplex & baryon_contracted_thing
    )
    {
        //This should be passed in as zero, but I'll do it here too just to be safe.
        baryon_contracted_thing = zero;
        for(int src_index = 0; src_index < src_weights.size(); src_index++)
        {
            for(int snk_index = 0; snk_index < snk_weights.size(); snk_index++)
            {
                LatticeComplex temp_contraction = zero;
                LatticeColorMatrix q1 = peekSpin(quark_1, snk_spins[snk_index][0], src_spins[src_index][0]);
                LatticeColorMatrix q2 = peekSpin(quark_2, snk_spins[snk_index][1], src_spins[src_index][1]);
                LatticeColorMatrix q3 = peekSpin(quark_3, snk_spins[snk_index][2], src_spins[src_index][2]);
                color_contraction(q1, q2, q3, temp_contraction);
                baryon_contracted_thing += src_weights[src_index] * snk_weights[snk_index] * temp_contraction;
            }
        }
    }

    void write_correlator(
        bool full_correlator,
        bool antiperiodic,
        std::string baryon_name,
        std::string spin,
#ifdef BUILD_HDF5
        std::string path,
        HDF5Writer & h5writer,
        HDF5Base::writemode & h5mode,
#endif
        int t_0,
        int Nt,
        multi1d<int> & source_coords,
        LalibeSftMom & FT,
        LatticeComplex & baryon
    )
    {
        if(full_correlator == true)
        {
#ifdef BUILD_HDF5
            std::string correlator_path = path+"/"+baryon_name+"/spin_"+spin+"/4D_correlator";
            h5writer.push(correlator_path);
            correlator_path = correlator_path   +"/x"+std::to_string(source_coords[0])
                                                +"_y"+std::to_string(source_coords[1])
                                                +"_z"+std::to_string(source_coords[2])
                                                +"_t"+std::to_string(source_coords[3]);
            h5writer.write(correlator_path, baryon, h5mode);
            h5writer.writeAttribute(correlator_path, "is_shifted", 0, h5mode);
            h5writer.cd("/");
#endif
        }
        else
        {
            //Momentum loop
            multi2d<Complex> FTed_baryon = FT.sft(baryon);
            //Temp variable for writing below.
            Complex temp_element;
            //Move the h5 pushing here, since all momentum keys will be written in the same general path.
#ifdef BUILD_HDF5
            std::string correlator_path = path+"/"+baryon_name+"/spin_"+spin
                                            +"/x"+std::to_string(source_coords[0])
                                            +"_y"+std::to_string(source_coords[1])
                                            +"_z"+std::to_string(source_coords[2])
                                            +"_t"+std::to_string(source_coords[3]);
            h5writer.push(correlator_path);
#else
            std::string correlator_path = baryon_name+"_spin-"+spin
                                            +"_x"+std::to_string(source_coords[0])
                                            +"_y"+std::to_string(source_coords[1])
                                            +"_z"+std::to_string(source_coords[2])
                                            +"_t"+std::to_string(source_coords[3]);
#endif
            for(int mom = 0; mom < FT.numMom(); mom++)
            {
                //One more temp variable instanited inside loop (once again for writing.)
                multi1d<Complex> baryon_correlator;
                baryon_correlator.resize(Nt);
                multi1d<int> momenta = FT.numToMom(mom);
#ifndef BUILD_HDF5
                std::string correlator_path_mom = correlator_path
                                                    +"_px"+std::to_string(momenta[0])
                                                    +"_py"+std::to_string(momenta[1])
                                                    +"_pz"+std::to_string(momenta[2]);
                TextFileWriter file_out(correlator_path_mom);
#endif
                for(int t = 0; t < Nt; t++)
                {
                    temp_element = FTed_baryon[mom][t];
                    int t_relative = t - t_0;
                    if(t_relative < 0)
                        t_relative += Nt;
                    if((t_relative >= (Nt - t_0)) && antiperiodic == true)
                        temp_element = -temp_element;
#ifndef BUILD_HDF5
                    file_out<<temp_element<<"\n";
#endif
                    baryon_correlator[t_relative] = temp_element;
                }
#ifndef BUILD_HDF5
                file_out.close();
#else
                //Change the name of string compred to 4d output so general correlator path is the same.
                std::string correlator_path_mom = correlator_path
                                                    +"/px"+std::to_string(momenta[0])
                                                    +"_py"+std::to_string(momenta[1])
                                                    +"_pz"+std::to_string(momenta[2]);
                h5writer.write(correlator_path_mom, baryon_correlator, h5mode);
                h5writer.writeAttribute(correlator_path_mom, "is_shifted", 1, h5mode);
                h5writer.cd("/");
#endif
            }
        }
    }
} // End namespace LegacyProton
} // End namespace Chroma
