/*! 
 *  Functions that do baryon spin and color contractions.
 *  Color is nested inside of spin.
 *  Authors:
 *  Ben Hoerz
 *  Andre Walker-Loud
 */

#include "baryon_contractions_func_w.h"
#include "chromabase.h"
#include "util/spin_basis.h"


namespace Chroma 
{ 

    // all the spin elementals required for various baryons in Dirac-Pauli basis,
    // each iPair corresponds to a pair of (sinkSpin, sourceSpin) indices
  namespace {

    typedef unsigned int uint;
    typedef std::pair<uint, uint> iPair;
    typedef std::vector<std::pair<std::tuple<iPair, iPair, iPair>, double>> SpinElemListType;

    //std::map<std::pair<std::string, std::string>, SpinElemListType> elemMap;


    std::map<std::string, std::tuple<char, char, char>> flavMap {
        // Octet
        { "proton", std::make_tuple('d', 'u', 'u')},
        { "neutron", std::make_tuple('d', 'd', 'u')},
        { "lambda_z", std::make_tuple('d', 's', 'u')},
        { "lambda_to_sigma", std::make_tuple('d', 's', 'u')},
        { "sigma_to_lambda", std::make_tuple('d', 's', 'u')},
        { "sigma_p", std::make_tuple('s', 'u', 'u')},
        { "sigma_z", std::make_tuple('d', 's', 'u')},
        { "sigma_m", std::make_tuple('d', 'd', 's')},
        { "xi_z", std::make_tuple('s', 's', 'u')},
        { "xi_m", std::make_tuple('d', 's', 's')},
        // Decuplet
        { "delta_pp", std::make_tuple('u', 'u', 'u')},
        { "delta_p", std::make_tuple('d', 'u', 'u')},
        { "delta_z", std::make_tuple('d', 'd', 'u')},
        { "delta_m", std::make_tuple('d', 'd', 'd')},
        { "sigma_star_p", std::make_tuple('s', 'u', 'u')},
        { "sigma_star_z", std::make_tuple('d', 's', 'u')},
        { "sigma_star_m", std::make_tuple('d', 'd', 's')},
        { "xi_star_z", std::make_tuple('s', 's', 'u')},
        { "xi_star_m", std::make_tuple('d', 's', 's')},
        { "omega_m", std::make_tuple('s', 's', 's')},
        // Octet negative parity
        { "proton_np", std::make_tuple('d', 'u', 'u')},
        { "neutron_np", std::make_tuple('d', 'd', 'u')},
        { "lambda_z_np", std::make_tuple('d', 's', 'u')},
        { "lambda_to_sigma_np", std::make_tuple('d', 's', 'u')},
        { "sigma_to_lambda_np", std::make_tuple('d', 's', 'u')},
        { "sigma_p_np", std::make_tuple('s', 'u', 'u')},
        { "sigma_z_np", std::make_tuple('d', 's', 'u')},
        { "sigma_m_np", std::make_tuple('d', 'd', 's')},
        { "xi_z_np", std::make_tuple('s', 's', 'u')},
        { "xi_m_np", std::make_tuple('d', 's', 's')},
        // Decuplet negative parity
        { "delta_pp_np", std::make_tuple('u', 'u', 'u')},
        { "delta_p_np", std::make_tuple('d', 'u', 'u')},
        { "delta_z_np", std::make_tuple('d', 'd', 'u')},
        { "delta_m_np", std::make_tuple('d', 'd', 'd')},
        { "sigma_star_p_np", std::make_tuple('s', 'u', 'u')},
        { "sigma_star_z_np", std::make_tuple('d', 's', 'u')},
        { "sigma_star_m_np", std::make_tuple('d', 'd', 's')},
        { "xi_star_z_np", std::make_tuple('s', 's', 'u')},
        { "xi_star_m_np", std::make_tuple('d', 's', 's')},
        { "omega_m_np", std::make_tuple('s', 's', 's')},
    };

    std::map<std::string, std::map<std::string, SpinElemListType>> elemMap = {
	{ "proton", {
	  { "up", {
	    { std::make_tuple(iPair(1,1), iPair(0,0), iPair(0,0)), 1.0 },
	    { std::make_tuple(iPair(0,1), iPair(0,0), iPair(1,0)), -1.0 },
	    { std::make_tuple(iPair(1,0), iPair(0,0), iPair(0,1)), -1.0 },
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(1,1)), 0.5 },
	    { std::make_tuple(iPair(0,0), iPair(0,1), iPair(1,0)), 0.5 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(1,1), iPair(0,0), iPair(1,1)), 0.5 },
	    { std::make_tuple(iPair(1,1), iPair(0,1), iPair(1,0)), 0.5 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(1,1)), -1.0 },
	    { std::make_tuple(iPair(1,0), iPair(0,1), iPair(1,1)), -1.0 },
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(1,1)), 1.0 },
	  } },
	} },
	{ "neutron", {
	  { "up", {
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(0,0)), 0.5 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(0,0)), 0.5 },
	    { std::make_tuple(iPair(0,0), iPair(0,1), iPair(1,0)), -1.0 },
	    { std::make_tuple(iPair(0,0), iPair(1,0), iPair(0,1)), -1.0 },
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(1,1)), 1.0 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(1,1), iPair(1,1), iPair(0,0)), 1.0 },
	    { std::make_tuple(iPair(0,1), iPair(1,1), iPair(1,0)), -1.0 },
	    { std::make_tuple(iPair(1,0), iPair(1,1), iPair(0,1)), -1.0 },
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(1,1)), 0.5 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(1,1)), 0.5 },
	  } },
	} },
	{ "omega_m", {
	  { "upup", {
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(0,0)), 6.0 },
	  } },
	  { "up", {
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(1,1)), 6.0 },
	    { std::make_tuple(iPair(0,0), iPair(0,1), iPair(1,0)), 12.0 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(1,1)), 6.0 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(1,1)), 12.0 },
	  } },
	  { "dndn", {
	    { std::make_tuple(iPair(1,1), iPair(1,1), iPair(1,1)), 6.0 },
	  } },
	} },
	{ "delta_pp", {
	  { "upup", {
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(0,0)), 6.0 },
	  } },
	  { "up", {
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(1,1)), 6.0 },
	    { std::make_tuple(iPair(0,0), iPair(0,1), iPair(1,0)), 12.0 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(1,1)), 6.0 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(1,1)), 12.0 },
	  } },
	  { "dndn", {
	    { std::make_tuple(iPair(1,1), iPair(1,1), iPair(1,1)), 6.0 },
	  } },
	} },
	{ "delta_p", {
	  { "upup", {
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(0,0)), 6.0 },
	  } },
	  { "up", {
	    { std::make_tuple(iPair(1,1), iPair(0,0), iPair(0,0)), 2.0 },
	    { std::make_tuple(iPair(0,1), iPair(0,0), iPair(1,0)), 4.0 },
	    { std::make_tuple(iPair(1,0), iPair(0,0), iPair(0,1)), 4.0 },
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(1,1)), 4.0 },
	    { std::make_tuple(iPair(0,0), iPair(0,1), iPair(1,0)), 4.0 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(1,1), iPair(0,0), iPair(1,1)), 4.0 },
	    { std::make_tuple(iPair(1,1), iPair(0,1), iPair(1,0)), 4.0 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(1,1)), 4.0 },
	    { std::make_tuple(iPair(1,0), iPair(0,1), iPair(1,1)), 4.0 },
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(1,1)), 2.0 },
	  } },
	  { "dndn", {
	    { std::make_tuple(iPair(1,1), iPair(1,1), iPair(1,1)), 6.0 },
	  } },
	} },
	{ "delta_z", {
	  { "upup", {
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(0,0)), 6.0 },
	  } },
	  { "up", {
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(0,0)), 4.0 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(0,0)), 4.0 },
	    { std::make_tuple(iPair(0,0), iPair(0,1), iPair(1,0)), 4.0 },
	    { std::make_tuple(iPair(0,0), iPair(1,0), iPair(0,1)), 4.0 },
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(1,1)), 2.0 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(1,1), iPair(1,1), iPair(0,0)), 2.0 },
	    { std::make_tuple(iPair(0,1), iPair(1,1), iPair(1,0)), 4.0 },
	    { std::make_tuple(iPair(1,0), iPair(1,1), iPair(0,1)), 4.0 },
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(1,1)), 4.0 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(1,1)), 4.0 },
	  } },
	  { "dndn", {
	    { std::make_tuple(iPair(1,1), iPair(1,1), iPair(1,1)), 6.0 },
	  } },
	} },
	{ "delta_m", {
	  { "upup", {
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(0,0)), 6.0 },
	  } },
	  { "up", {
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(1,1)), 6.0 },
	    { std::make_tuple(iPair(0,0), iPair(0,1), iPair(1,0)), 12.0 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(1,1)), 6.0 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(1,1)), 12.0 },
	  } },
	  { "dndn", {
	    { std::make_tuple(iPair(1,1), iPair(1,1), iPair(1,1)), 6.0 },
	  } },
	} },
	{ "lambda_z", {
	  { "up", {
	    { std::make_tuple(iPair(1,1), iPair(0,0), iPair(0,0)), 0.5 },
	    { std::make_tuple(iPair(0,1), iPair(0,0), iPair(1,0)), -0.5 },
	    { std::make_tuple(iPair(1,0), iPair(0,0), iPair(0,1)), -0.5 },
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(1,1)), 0.5 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(1,1), iPair(1,1), iPair(0,0)), 0.5 },
	    { std::make_tuple(iPair(0,1), iPair(1,1), iPair(1,0)), -0.5 },
	    { std::make_tuple(iPair(1,0), iPair(1,1), iPair(0,1)), -0.5 },
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(1,1)), 0.5 },
	  } },
	} },
	{ "sigma_p", {
	  { "up", {
	    { std::make_tuple(iPair(1,1), iPair(0,0), iPair(0,0)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,1), iPair(0,0), iPair(1,0)), -1.3333333333333333 },
	    { std::make_tuple(iPair(1,0), iPair(0,0), iPair(0,1)), -1.3333333333333333 },
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(1,1)), 0.6666666666666666 },
	    { std::make_tuple(iPair(0,0), iPair(0,1), iPair(1,0)), 0.6666666666666666 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(1,1)), 1.3333333333333333 },
	    { std::make_tuple(iPair(1,0), iPair(0,1), iPair(1,1)), -1.3333333333333333 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(1,1)), -1.3333333333333333 },
	    { std::make_tuple(iPair(1,1), iPair(0,0), iPair(1,1)), 0.6666666666666666 },
	    { std::make_tuple(iPair(1,1), iPair(0,1), iPair(1,0)), 0.6666666666666666 },
	  } },
	} },
	{ "sigma_z", {
	  { "up", {
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(0,0)), 1.3333333333333333 },
	    { std::make_tuple(iPair(1,0), iPair(0,1), iPair(0,0)), -0.6666666666666666 },
	    { std::make_tuple(iPair(0,0), iPair(0,1), iPair(1,0)), -0.6666666666666666 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(0,0)), -0.6666666666666666 },
	    { std::make_tuple(iPair(1,1), iPair(0,0), iPair(0,0)), 0.3333333333333333 },
	    { std::make_tuple(iPair(0,1), iPair(0,0), iPair(1,0)), 0.3333333333333333 },
	    { std::make_tuple(iPair(0,0), iPair(1,0), iPair(0,1)), -0.6666666666666666 },
	    { std::make_tuple(iPair(1,0), iPair(0,0), iPair(0,1)), 0.3333333333333333 },
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(1,1)), 0.3333333333333333 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(1,1), iPair(0,0), iPair(1,1)), 1.3333333333333333 },
	    { std::make_tuple(iPair(1,1), iPair(1,0), iPair(0,1)), -0.6666666666666666 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(1,1)), -0.6666666666666666 },
	    { std::make_tuple(iPair(1,1), iPair(0,1), iPair(1,0)), -0.6666666666666666 },
	    { std::make_tuple(iPair(1,1), iPair(1,1), iPair(0,0)), 0.3333333333333333 },
	    { std::make_tuple(iPair(0,1), iPair(1,1), iPair(1,0)), 0.3333333333333333 },
	    { std::make_tuple(iPair(1,0), iPair(0,1), iPair(1,1)), -0.6666666666666666 },
	    { std::make_tuple(iPair(1,0), iPair(1,1), iPair(0,1)), 0.3333333333333333 },
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(1,1)), 0.3333333333333333 },
	  } },
	} },
	{ "sigma_m", {
	  { "up", {
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(1,1)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,0), iPair(1,0), iPair(0,1)), -1.3333333333333333 },
	    { std::make_tuple(iPair(0,0), iPair(0,1), iPair(1,0)), -1.3333333333333333 },
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(0,0)), 0.6666666666666666 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(0,0)), 0.6666666666666666 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(1,1), iPair(1,1), iPair(0,0)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,1), iPair(1,1), iPair(1,0)), -1.3333333333333333 },
	    { std::make_tuple(iPair(1,0), iPair(1,1), iPair(0,1)), -1.3333333333333333 },
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(1,1)), 0.6666666666666666 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(1,1)), 0.6666666666666666 },
	  } },
	} },
	{ "xi_z", {
	  { "up", {
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(1,1)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,0), iPair(1,0), iPair(0,1)), -1.3333333333333333 },
	    { std::make_tuple(iPair(0,0), iPair(0,1), iPair(1,0)), -1.3333333333333333 },
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(0,0)), 0.6666666666666666 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(0,0)), 0.6666666666666666 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(1,1), iPair(1,1), iPair(0,0)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,1), iPair(1,1), iPair(1,0)), -1.3333333333333333 },
	    { std::make_tuple(iPair(1,0), iPair(1,1), iPair(0,1)), -1.3333333333333333 },
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(1,1)), 0.6666666666666666 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(1,1)), 0.6666666666666666 },
	  } },
	} },
	{ "xi_m", {
	  { "up", {
	    { std::make_tuple(iPair(1,1), iPair(0,0), iPair(0,0)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,1), iPair(0,0), iPair(1,0)), -1.3333333333333333 },
	    { std::make_tuple(iPair(1,0), iPair(0,0), iPair(0,1)), -1.3333333333333333 },
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(1,1)), 0.6666666666666666 },
	    { std::make_tuple(iPair(0,0), iPair(0,1), iPair(1,0)), 0.6666666666666666 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(1,1)), 1.3333333333333333 },
	    { std::make_tuple(iPair(1,0), iPair(0,1), iPair(1,1)), -1.3333333333333333 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(1,1)), -1.3333333333333333 },
	    { std::make_tuple(iPair(1,1), iPair(0,0), iPair(1,1)), 0.6666666666666666 },
	    { std::make_tuple(iPair(1,1), iPair(0,1), iPair(1,0)), 0.6666666666666666 },
	  } },
	} },
	{ "sigma_star_p", {
	  { "upup", {
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(0,0)), 2.0 },
	  } },
	  { "up", {
	    { std::make_tuple(iPair(1,1), iPair(0,0), iPair(0,0)), 0.6666666666666666 },
	    { std::make_tuple(iPair(0,1), iPair(0,0), iPair(1,0)), 1.3333333333333333 },
	    { std::make_tuple(iPair(1,0), iPair(0,0), iPair(0,1)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(1,1)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,0), iPair(0,1), iPair(1,0)), 1.3333333333333333 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(1,1), iPair(0,0), iPair(1,1)), 1.3333333333333333 },
	    { std::make_tuple(iPair(1,1), iPair(0,1), iPair(1,0)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(1,1)), 1.3333333333333333 },
	    { std::make_tuple(iPair(1,0), iPair(0,1), iPair(1,1)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(1,1)), 0.6666666666666666 },
	  } },
	  { "dndn", {
	    { std::make_tuple(iPair(1,1), iPair(1,1), iPair(1,1)), 2.0 },
	  } },
	} },
	{ "sigma_star_z", {
	  { "upup", {
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(0,0)), 2.0 },
	  } },
	  { "up", {
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(0,0)), 0.6666666666666666 },
	    { std::make_tuple(iPair(1,0), iPair(0,1), iPair(0,0)), 0.6666666666666666 },
	    { std::make_tuple(iPair(0,0), iPair(0,1), iPair(1,0)), 0.6666666666666666 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(0,0)), 0.6666666666666666 },
	    { std::make_tuple(iPair(1,1), iPair(0,0), iPair(0,0)), 0.6666666666666666 },
	    { std::make_tuple(iPair(0,1), iPair(0,0), iPair(1,0)), 0.6666666666666666 },
	    { std::make_tuple(iPair(0,0), iPair(1,0), iPair(0,1)), 0.6666666666666666 },
	    { std::make_tuple(iPair(1,0), iPair(0,0), iPair(0,1)), 0.6666666666666666 },
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(1,1)), 0.6666666666666666 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(1,1), iPair(1,1), iPair(0,0)), 0.6666666666666666 },
	    { std::make_tuple(iPair(1,1), iPair(0,1), iPair(1,0)), 0.6666666666666666 },
	    { std::make_tuple(iPair(0,1), iPair(1,1), iPair(1,0)), 0.6666666666666666 },
	    { std::make_tuple(iPair(1,1), iPair(1,0), iPair(0,1)), 0.6666666666666666 },
	    { std::make_tuple(iPair(1,1), iPair(0,0), iPair(1,1)), 0.6666666666666666 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(1,1)), 0.6666666666666666 },
	    { std::make_tuple(iPair(1,0), iPair(1,1), iPair(0,1)), 0.6666666666666666 },
	    { std::make_tuple(iPair(1,0), iPair(0,1), iPair(1,1)), 0.6666666666666666 },
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(1,1)), 0.6666666666666666 },
	  } },
	  { "dndn", {
	    { std::make_tuple(iPair(1,1), iPair(1,1), iPair(1,1)), 2.0 },
	  } },
	} },
	{ "sigma_star_m", {
	  { "upup", {
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(0,0)), 2.0 },
	  } },
	  { "up", {
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(1,1)), 0.6666666666666666 },
	    { std::make_tuple(iPair(0,0), iPair(1,0), iPair(0,1)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,0), iPair(0,1), iPair(1,0)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(0,0)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(0,0)), 1.3333333333333333 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(1,1)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(1,1)), 1.3333333333333333 },
	    { std::make_tuple(iPair(1,0), iPair(1,1), iPair(0,1)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,1), iPair(1,1), iPair(1,0)), 1.3333333333333333 },
	    { std::make_tuple(iPair(1,1), iPair(1,1), iPair(0,0)), 0.6666666666666666 },
	  } },
	  { "dndn", {
	    { std::make_tuple(iPair(1,1), iPair(1,1), iPair(1,1)), 2.0 },
	  } },
	} },
	{ "xi_star_z", {
	  { "upup", {
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(0,0)), 2.0 },
	  } },
	  { "up", {
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(1,1)), 0.6666666666666666 },
	    { std::make_tuple(iPair(0,0), iPair(1,0), iPair(0,1)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,0), iPair(0,1), iPair(1,0)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(0,0)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(0,0)), 1.3333333333333333 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(1,1)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(1,1)), 1.3333333333333333 },
	    { std::make_tuple(iPair(1,0), iPair(1,1), iPair(0,1)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,1), iPair(1,1), iPair(1,0)), 1.3333333333333333 },
	    { std::make_tuple(iPair(1,1), iPair(1,1), iPair(0,0)), 0.6666666666666666 },
	  } },
	  { "dndn", {
	    { std::make_tuple(iPair(1,1), iPair(1,1), iPair(1,1)), 2.0 },
	  } },
	} },
	{ "xi_star_m", {
	  { "upup", {
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(0,0)), 2.0 },
	  } },
	  { "up", {
	    { std::make_tuple(iPair(1,1), iPair(0,0), iPair(0,0)), 0.6666666666666666 },
	    { std::make_tuple(iPair(0,1), iPair(0,0), iPair(1,0)), 1.3333333333333333 },
	    { std::make_tuple(iPair(1,0), iPair(0,0), iPair(0,1)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(1,1)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,0), iPair(0,1), iPair(1,0)), 1.3333333333333333 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(1,1), iPair(0,0), iPair(1,1)), 1.3333333333333333 },
	    { std::make_tuple(iPair(1,1), iPair(0,1), iPair(1,0)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(1,1)), 1.3333333333333333 },
	    { std::make_tuple(iPair(1,0), iPair(0,1), iPair(1,1)), 1.3333333333333333 },
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(1,1)), 0.6666666666666666 },
	  } },
	  { "dndn", {
	    { std::make_tuple(iPair(1,1), iPair(1,1), iPair(1,1)), 2.0 },
	  } },
	} },
	{ "proton_np", {
	  { "up", {
	    { std::make_tuple(iPair(3,3), iPair(2,2), iPair(2,2)), 1.0 },
	    { std::make_tuple(iPair(2,3), iPair(2,2), iPair(3,2)), -1.0 },
	    { std::make_tuple(iPair(3,2), iPair(2,2), iPair(2,3)), -1.0 },
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(3,3)), 0.5 },
	    { std::make_tuple(iPair(2,2), iPair(2,3), iPair(3,2)), 0.5 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(3,3), iPair(2,2), iPair(3,3)), 0.5 },
	    { std::make_tuple(iPair(3,3), iPair(2,3), iPair(3,2)), 0.5 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(3,3)), -1.0 },
	    { std::make_tuple(iPair(3,2), iPair(2,3), iPair(3,3)), -1.0 },
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(3,3)), 1.0 },
	  } },
	} },
	{ "neutron_np", {
	  { "up", {
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(2,2)), 0.5 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(2,2)), 0.5 },
	    { std::make_tuple(iPair(2,2), iPair(2,3), iPair(3,2)), -1.0 },
	    { std::make_tuple(iPair(2,2), iPair(3,2), iPair(2,3)), -1.0 },
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(3,3)), 1.0 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(3,3), iPair(3,3), iPair(2,2)), 1.0 },
	    { std::make_tuple(iPair(2,3), iPair(3,3), iPair(3,2)), -1.0 },
	    { std::make_tuple(iPair(3,2), iPair(3,3), iPair(2,3)), -1.0 },
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(3,3)), 0.5 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(3,3)), 0.5 },
	  } },
	} },
	{ "omega_m_np", {
	  { "upup", {
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(2,2)), 6.0 },
	  } },
	  { "up", {
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(3,3)), 6.0 },
	    { std::make_tuple(iPair(2,2), iPair(2,3), iPair(3,2)), 12.0 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(3,3)), 6.0 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(3,3)), 12.0 },
	  } },
	  { "dndn", {
	    { std::make_tuple(iPair(3,3), iPair(3,3), iPair(3,3)), 6.0 },
	  } },
	} },
	{ "delta_pp_np", {
	  { "upup", {
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(2,2)), 6.0 },
	  } },
	  { "up", {
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(3,3)), 6.0 },
	    { std::make_tuple(iPair(2,2), iPair(2,3), iPair(3,2)), 12.0 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(3,3)), 6.0 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(3,3)), 12.0 },
	  } },
	  { "dndn", {
	    { std::make_tuple(iPair(3,3), iPair(3,3), iPair(3,3)), 6.0 },
	  } },
	} },
	{ "delta_p_np", {
	  { "upup", {
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(2,2)), 6.0 },
	  } },
	  { "up", {
	    { std::make_tuple(iPair(3,3), iPair(2,2), iPair(2,2)), 2.0 },
	    { std::make_tuple(iPair(2,3), iPair(2,2), iPair(3,2)), 4.0 },
	    { std::make_tuple(iPair(3,2), iPair(2,2), iPair(2,3)), 4.0 },
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(3,3)), 4.0 },
	    { std::make_tuple(iPair(2,2), iPair(2,3), iPair(3,2)), 4.0 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(3,3), iPair(2,2), iPair(3,3)), 4.0 },
	    { std::make_tuple(iPair(3,3), iPair(2,3), iPair(3,2)), 4.0 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(3,3)), 4.0 },
	    { std::make_tuple(iPair(3,2), iPair(2,3), iPair(3,3)), 4.0 },
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(3,3)), 2.0 },
	  } },
	  { "dndn", {
	    { std::make_tuple(iPair(3,3), iPair(3,3), iPair(3,3)), 6.0 },
	  } },
	} },
	{ "delta_z_np", {
	  { "upup", {
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(2,2)), 6.0 },
	  } },
	  { "up", {
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(2,2)), 4.0 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(2,2)), 4.0 },
	    { std::make_tuple(iPair(2,2), iPair(2,3), iPair(3,2)), 4.0 },
	    { std::make_tuple(iPair(2,2), iPair(3,2), iPair(2,3)), 4.0 },
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(3,3)), 2.0 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(3,3), iPair(3,3), iPair(2,2)), 2.0 },
	    { std::make_tuple(iPair(2,3), iPair(3,3), iPair(3,2)), 4.0 },
	    { std::make_tuple(iPair(3,2), iPair(3,3), iPair(2,3)), 4.0 },
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(3,3)), 4.0 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(3,3)), 4.0 },
	  } },
	  { "dndn", {
	    { std::make_tuple(iPair(3,3), iPair(3,3), iPair(3,3)), 6.0 },
	  } },
	} },
	{ "delta_m_np", {
	  { "upup", {
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(2,2)), 6.0 },
	  } },
	  { "up", {
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(3,3)), 6.0 },
	    { std::make_tuple(iPair(2,2), iPair(2,3), iPair(3,2)), 12.0 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(3,3)), 6.0 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(3,3)), 12.0 },
	  } },
	  { "dndn", {
	    { std::make_tuple(iPair(3,3), iPair(3,3), iPair(3,3)), 6.0 },
	  } },
	} },
	{ "lambda_z_np", {
	  { "up", {
	    { std::make_tuple(iPair(3,3), iPair(2,2), iPair(2,2)), 0.5 },
	    { std::make_tuple(iPair(2,3), iPair(2,2), iPair(3,2)), -0.5 },
	    { std::make_tuple(iPair(3,2), iPair(2,2), iPair(2,3)), -0.5 },
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(3,3)), 0.5 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(3,3), iPair(3,3), iPair(2,2)), 0.5 },
	    { std::make_tuple(iPair(2,3), iPair(3,3), iPair(3,2)), -0.5 },
	    { std::make_tuple(iPair(3,2), iPair(3,3), iPair(2,3)), -0.5 },
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(3,3)), 0.5 },
	  } },
	} },
	{ "sigma_p_np", {
	  { "up", {
	    { std::make_tuple(iPair(3,3), iPair(2,2), iPair(2,2)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,3), iPair(2,2), iPair(3,2)), -1.3333333333333333 },
	    { std::make_tuple(iPair(3,2), iPair(2,2), iPair(2,3)), -1.3333333333333333 },
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(3,3)), 0.6666666666666666 },
	    { std::make_tuple(iPair(2,2), iPair(2,3), iPair(3,2)), 0.6666666666666666 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(3,3)), 1.3333333333333333 },
	    { std::make_tuple(iPair(3,2), iPair(2,3), iPair(3,3)), -1.3333333333333333 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(3,3)), -1.3333333333333333 },
	    { std::make_tuple(iPair(3,3), iPair(2,2), iPair(3,3)), 0.6666666666666666 },
	    { std::make_tuple(iPair(3,3), iPair(2,3), iPair(3,2)), 0.6666666666666666 },
	  } },
	} },
	{ "sigma_z_np", {
	  { "up", {
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(2,2)), 1.3333333333333333 },
	    { std::make_tuple(iPair(3,2), iPair(2,3), iPair(2,2)), -0.6666666666666666 },
	    { std::make_tuple(iPair(2,2), iPair(2,3), iPair(3,2)), -0.6666666666666666 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(2,2)), -0.6666666666666666 },
	    { std::make_tuple(iPair(3,3), iPair(2,2), iPair(2,2)), 0.3333333333333333 },
	    { std::make_tuple(iPair(2,3), iPair(2,2), iPair(3,2)), 0.3333333333333333 },
	    { std::make_tuple(iPair(2,2), iPair(3,2), iPair(2,3)), -0.6666666666666666 },
	    { std::make_tuple(iPair(3,2), iPair(2,2), iPair(2,3)), 0.3333333333333333 },
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(3,3)), 0.3333333333333333 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(3,3), iPair(2,2), iPair(3,3)), 1.3333333333333333 },
	    { std::make_tuple(iPair(3,3), iPair(3,2), iPair(2,3)), -0.6666666666666666 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(3,3)), -0.6666666666666666 },
	    { std::make_tuple(iPair(3,3), iPair(2,3), iPair(3,2)), -0.6666666666666666 },
	    { std::make_tuple(iPair(3,3), iPair(3,3), iPair(2,2)), 0.3333333333333333 },
	    { std::make_tuple(iPair(2,3), iPair(3,3), iPair(3,2)), 0.3333333333333333 },
	    { std::make_tuple(iPair(3,2), iPair(2,3), iPair(3,3)), -0.6666666666666666 },
	    { std::make_tuple(iPair(3,2), iPair(3,3), iPair(2,3)), 0.3333333333333333 },
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(3,3)), 0.3333333333333333 },
	  } },
	} },
	{ "sigma_m_np", {
	  { "up", {
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(3,3)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,2), iPair(3,2), iPair(2,3)), -1.3333333333333333 },
	    { std::make_tuple(iPair(2,2), iPair(2,3), iPair(3,2)), -1.3333333333333333 },
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(2,2)), 0.6666666666666666 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(2,2)), 0.6666666666666666 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(3,3), iPair(3,3), iPair(2,2)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,3), iPair(3,3), iPair(3,2)), -1.3333333333333333 },
	    { std::make_tuple(iPair(3,2), iPair(3,3), iPair(2,3)), -1.3333333333333333 },
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(3,3)), 0.6666666666666666 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(3,3)), 0.6666666666666666 },
	  } },
	} },
	{ "xi_z_np", {
	  { "up", {
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(3,3)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,2), iPair(3,2), iPair(2,3)), -1.3333333333333333 },
	    { std::make_tuple(iPair(2,2), iPair(2,3), iPair(3,2)), -1.3333333333333333 },
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(2,2)), 0.6666666666666666 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(2,2)), 0.6666666666666666 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(3,3), iPair(3,3), iPair(2,2)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,3), iPair(3,3), iPair(3,2)), -1.3333333333333333 },
	    { std::make_tuple(iPair(3,2), iPair(3,3), iPair(2,3)), -1.3333333333333333 },
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(3,3)), 0.6666666666666666 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(3,3)), 0.6666666666666666 },
	  } },
	} },
	{ "xi_m_np", {
	  { "up", {
	    { std::make_tuple(iPair(3,3), iPair(2,2), iPair(2,2)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,3), iPair(2,2), iPair(3,2)), -1.3333333333333333 },
	    { std::make_tuple(iPair(3,2), iPair(2,2), iPair(2,3)), -1.3333333333333333 },
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(3,3)), 0.6666666666666666 },
	    { std::make_tuple(iPair(2,2), iPair(2,3), iPair(3,2)), 0.6666666666666666 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(3,3)), 1.3333333333333333 },
	    { std::make_tuple(iPair(3,2), iPair(2,3), iPair(3,3)), -1.3333333333333333 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(3,3)), -1.3333333333333333 },
	    { std::make_tuple(iPair(3,3), iPair(2,2), iPair(3,3)), 0.6666666666666666 },
	    { std::make_tuple(iPair(3,3), iPair(2,3), iPair(3,2)), 0.6666666666666666 },
	  } },
	} },
	{ "sigma_star_p_np", {
	  { "upup", {
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(2,2)), 2.0 },
	  } },
	  { "up", {
	    { std::make_tuple(iPair(3,3), iPair(2,2), iPair(2,2)), 0.6666666666666666 },
	    { std::make_tuple(iPair(2,3), iPair(2,2), iPair(3,2)), 1.3333333333333333 },
	    { std::make_tuple(iPair(3,2), iPair(2,2), iPair(2,3)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(3,3)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,2), iPair(2,3), iPair(3,2)), 1.3333333333333333 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(3,3), iPair(2,2), iPair(3,3)), 1.3333333333333333 },
	    { std::make_tuple(iPair(3,3), iPair(2,3), iPair(3,2)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(3,3)), 1.3333333333333333 },
	    { std::make_tuple(iPair(3,2), iPair(2,3), iPair(3,3)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(3,3)), 0.6666666666666666 },
	  } },
	  { "dndn", {
	    { std::make_tuple(iPair(3,3), iPair(3,3), iPair(3,3)), 2.0 },
	  } },
	} },
	{ "sigma_star_z_np", {
	  { "upup", {
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(2,2)), 2.0 },
	  } },
	  { "up", {
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(2,2)), 0.6666666666666666 },
	    { std::make_tuple(iPair(3,2), iPair(2,3), iPair(2,2)), 0.6666666666666666 },
	    { std::make_tuple(iPair(2,2), iPair(2,3), iPair(3,2)), 0.6666666666666666 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(2,2)), 0.6666666666666666 },
	    { std::make_tuple(iPair(3,3), iPair(2,2), iPair(2,2)), 0.6666666666666666 },
	    { std::make_tuple(iPair(2,3), iPair(2,2), iPair(3,2)), 0.6666666666666666 },
	    { std::make_tuple(iPair(2,2), iPair(3,2), iPair(2,3)), 0.6666666666666666 },
	    { std::make_tuple(iPair(3,2), iPair(2,2), iPair(2,3)), 0.6666666666666666 },
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(3,3)), 0.6666666666666666 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(3,3), iPair(3,3), iPair(2,2)), 0.6666666666666666 },
	    { std::make_tuple(iPair(3,3), iPair(2,3), iPair(3,2)), 0.6666666666666666 },
	    { std::make_tuple(iPair(2,3), iPair(3,3), iPair(3,2)), 0.6666666666666666 },
	    { std::make_tuple(iPair(3,3), iPair(3,2), iPair(2,3)), 0.6666666666666666 },
	    { std::make_tuple(iPair(3,3), iPair(2,2), iPair(3,3)), 0.6666666666666666 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(3,3)), 0.6666666666666666 },
	    { std::make_tuple(iPair(3,2), iPair(3,3), iPair(2,3)), 0.6666666666666666 },
	    { std::make_tuple(iPair(3,2), iPair(2,3), iPair(3,3)), 0.6666666666666666 },
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(3,3)), 0.6666666666666666 },
	  } },
	  { "dndn", {
	    { std::make_tuple(iPair(3,3), iPair(3,3), iPair(3,3)), 2.0 },
	  } },
	} },
	{ "sigma_star_m_np", {
	  { "upup", {
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(2,2)), 2.0 },
	  } },
	  { "up", {
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(3,3)), 0.6666666666666666 },
	    { std::make_tuple(iPair(2,2), iPair(3,2), iPair(2,3)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,2), iPair(2,3), iPair(3,2)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(2,2)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(2,2)), 1.3333333333333333 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(3,3)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(3,3)), 1.3333333333333333 },
	    { std::make_tuple(iPair(3,2), iPair(3,3), iPair(2,3)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,3), iPair(3,3), iPair(3,2)), 1.3333333333333333 },
	    { std::make_tuple(iPair(3,3), iPair(3,3), iPair(2,2)), 0.6666666666666666 },
	  } },
	  { "dndn", {
	    { std::make_tuple(iPair(3,3), iPair(3,3), iPair(3,3)), 2.0 },
	  } },
	} },
	{ "xi_star_z_np", {
	  { "upup", {
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(2,2)), 2.0 },
	  } },
	  { "up", {
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(3,3)), 0.6666666666666666 },
	    { std::make_tuple(iPair(2,2), iPair(3,2), iPair(2,3)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,2), iPair(2,3), iPair(3,2)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(2,2)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(2,2)), 1.3333333333333333 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(3,3)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(3,3)), 1.3333333333333333 },
	    { std::make_tuple(iPair(3,2), iPair(3,3), iPair(2,3)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,3), iPair(3,3), iPair(3,2)), 1.3333333333333333 },
	    { std::make_tuple(iPair(3,3), iPair(3,3), iPair(2,2)), 0.6666666666666666 },
	  } },
	  { "dndn", {
	    { std::make_tuple(iPair(3,3), iPair(3,3), iPair(3,3)), 2.0 },
	  } },
	} },
	{ "xi_star_m_np", {
	  { "upup", {
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(2,2)), 2.0 },
	  } },
	  { "up", {
	    { std::make_tuple(iPair(3,3), iPair(2,2), iPair(2,2)), 0.6666666666666666 },
	    { std::make_tuple(iPair(2,3), iPair(2,2), iPair(3,2)), 1.3333333333333333 },
	    { std::make_tuple(iPair(3,2), iPair(2,2), iPair(2,3)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(3,3)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,2), iPair(2,3), iPair(3,2)), 1.3333333333333333 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(3,3), iPair(2,2), iPair(3,3)), 1.3333333333333333 },
	    { std::make_tuple(iPair(3,3), iPair(2,3), iPair(3,2)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(3,3)), 1.3333333333333333 },
	    { std::make_tuple(iPair(3,2), iPair(2,3), iPair(3,3)), 1.3333333333333333 },
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(3,3)), 0.6666666666666666 },
	  } },
	  { "dndn", {
	    { std::make_tuple(iPair(3,3), iPair(3,3), iPair(3,3)), 2.0 },
	  } },
	} },
	{ "sigma_to_lambda", {
	  { "up", {
	    { std::make_tuple(iPair(1,0), iPair(0,1), iPair(0,0)), 0.816496580927726 },
	    { std::make_tuple(iPair(0,0), iPair(0,1), iPair(1,0)), -0.816496580927726 },
	    { std::make_tuple(iPair(1,1), iPair(0,0), iPair(0,0)), -0.408248290463863 },
	    { std::make_tuple(iPair(0,1), iPair(0,0), iPair(1,0)), 0.408248290463863 },
	    { std::make_tuple(iPair(1,0), iPair(0,0), iPair(0,1)), -0.408248290463863 },
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(1,1)), 0.408248290463863 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(1,1), iPair(1,0), iPair(0,1)), -0.816496580927726 },
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(1,1)), 0.816496580927726 },
	    { std::make_tuple(iPair(1,1), iPair(1,1), iPair(0,0)), 0.408248290463863 },
	    { std::make_tuple(iPair(0,1), iPair(1,1), iPair(1,0)), -0.408248290463863 },
	    { std::make_tuple(iPair(1,0), iPair(1,1), iPair(0,1)), 0.408248290463863 },
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(1,1)), -0.408248290463863 },
	  } },
	} },
	{ "lambda_to_sigma", {
	  { "up", {
	    { std::make_tuple(iPair(0,1), iPair(1,0), iPair(0,0)), 0.816496580927726 },
	    { std::make_tuple(iPair(1,1), iPair(0,0), iPair(0,0)), -0.408248290463863 },
	    { std::make_tuple(iPair(0,1), iPair(0,0), iPair(1,0)), -0.408248290463863 },
	    { std::make_tuple(iPair(0,0), iPair(1,0), iPair(0,1)), -0.816496580927726 },
	    { std::make_tuple(iPair(1,0), iPair(0,0), iPair(0,1)), 0.408248290463863 },
	    { std::make_tuple(iPair(0,0), iPair(0,0), iPair(1,1)), 0.408248290463863 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(1,1), iPair(0,1), iPair(1,0)), -0.816496580927726 },
	    { std::make_tuple(iPair(1,1), iPair(1,1), iPair(0,0)), 0.408248290463863 },
	    { std::make_tuple(iPair(0,1), iPair(1,1), iPair(1,0)), 0.408248290463863 },
	    { std::make_tuple(iPair(1,0), iPair(0,1), iPair(1,1)), 0.816496580927726 },
	    { std::make_tuple(iPair(1,0), iPair(1,1), iPair(0,1)), -0.408248290463863 },
	    { std::make_tuple(iPair(0,0), iPair(1,1), iPair(1,1)), -0.408248290463863 },
	  } },
	} },
	{ "sigma_to_lambda_np", {
	  { "up", {
	    { std::make_tuple(iPair(3,2), iPair(2,3), iPair(2,2)), 0.816496580927726 },
	    { std::make_tuple(iPair(2,2), iPair(2,3), iPair(3,2)), -0.816496580927726 },
	    { std::make_tuple(iPair(3,3), iPair(2,2), iPair(2,2)), -0.408248290463863 },
	    { std::make_tuple(iPair(2,3), iPair(2,2), iPair(3,2)), 0.408248290463863 },
	    { std::make_tuple(iPair(3,2), iPair(2,2), iPair(2,3)), -0.408248290463863 },
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(3,3)), 0.408248290463863 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(3,3), iPair(3,2), iPair(2,3)), -0.816496580927726 },
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(3,3)), 0.816496580927726 },
	    { std::make_tuple(iPair(3,3), iPair(3,3), iPair(2,2)), 0.408248290463863 },
	    { std::make_tuple(iPair(2,3), iPair(3,3), iPair(3,2)), -0.408248290463863 },
	    { std::make_tuple(iPair(3,2), iPair(3,3), iPair(2,3)), 0.408248290463863 },
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(3,3)), -0.408248290463863 },
	  } },
	} },
	{ "lambda_to_sigma_np", {
	  { "up", {
	    { std::make_tuple(iPair(2,3), iPair(3,2), iPair(2,2)), 0.816496580927726 },
	    { std::make_tuple(iPair(3,3), iPair(2,2), iPair(2,2)), -0.408248290463863 },
	    { std::make_tuple(iPair(2,3), iPair(2,2), iPair(3,2)), -0.408248290463863 },
	    { std::make_tuple(iPair(2,2), iPair(3,2), iPair(2,3)), -0.816496580927726 },
	    { std::make_tuple(iPair(3,2), iPair(2,2), iPair(2,3)), 0.408248290463863 },
	    { std::make_tuple(iPair(2,2), iPair(2,2), iPair(3,3)), 0.408248290463863 },
	  } },
	  { "dn", {
	    { std::make_tuple(iPair(3,3), iPair(2,3), iPair(3,2)), -0.816496580927726 },
	    { std::make_tuple(iPair(3,3), iPair(3,3), iPair(2,2)), 0.408248290463863 },
	    { std::make_tuple(iPair(2,3), iPair(3,3), iPair(3,2)), 0.408248290463863 },
	    { std::make_tuple(iPair(3,2), iPair(2,3), iPair(3,3)), 0.816496580927726 },
	    { std::make_tuple(iPair(3,2), iPair(3,3), iPair(2,3)), -0.408248290463863 },
	    { std::make_tuple(iPair(2,2), iPair(3,3), iPair(3,3)), -0.408248290463863 },
	  } },
	} },
    };

  }

  std::tuple<char,char,char> get_flavor_code(const std::string& baryon_name) {
    auto fIt = flavMap.find(baryon_name);
    if (fIt != flavMap.end()) {
      return fIt->second;
    }
    else {
      QDPIO::cerr << "Did not find flavor_code for "<<baryon_name << std::endl;
      QDP_abort(1);
    }
  }


  std::vector<std::string> get_spin_components(const std::string& baryon_name) {
    auto mIt = elemMap.find(baryon_name);
    if (mIt == elemMap.end()) {
      QDPIO::cerr << "Did not find spin elementals for " << baryon_name << std::endl;
      QDP_abort(1);
    }

    std::vector<std::string> retVec;
    for (auto sIt : mIt->second)
      retVec.push_back(sIt.first);

    return retVec;
  }


  void color_contraction(const LatticeColorMatrix & quark_1,
			 const LatticeColorMatrix & quark_2,
			 const LatticeColorMatrix & quark_3,
			 LatticeComplex & spin_color_contracted_thing)
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


  void do_contraction(const LatticePropagator & quark_1,
		      const LatticePropagator & quark_2,
		      const LatticePropagator & quark_3,
		      const std::string& baryon_name,
		      const std::string& spin,
		      LatticeComplex & baryon_contracted_thing)
  {
      // get spin elemental list
    auto bIt = elemMap.find(baryon_name);
    if (bIt == elemMap.end()) {
      QDPIO::cerr << "ERROR: Did not find flavor "<< baryon_name<< std::endl;
      QDP_abort(1);
    }

    auto mIt = bIt->second.find(spin);
    if (mIt == bIt->second.end()) {
      QDPIO::cerr << "ERROR: Did not find spin elementals for " << baryon_name << " " <<spin<< std::endl;
      QDP_abort(1);
    }

    //This should be passed in as zero, but I'll do it here too just to be safe.
    baryon_contracted_thing = zero;

    for (auto spinEl : mIt->second) {
      std::tuple<iPair, iPair, iPair> spinComb;
      double coeff;

      std::tie(spinComb, coeff) = spinEl;

      LatticeComplex temp_contraction = zero;
      LatticeColorMatrix q1 = peekSpin(quark_1, std::get<0>(spinComb).first, std::get<0>(spinComb).second);
      LatticeColorMatrix q2 = peekSpin(quark_2, std::get<1>(spinComb).first, std::get<1>(spinComb).second);
      LatticeColorMatrix q3 = peekSpin(quark_3, std::get<2>(spinComb).first, std::get<2>(spinComb).second);
      color_contraction(q1, q2, q3, temp_contraction);
      baryon_contracted_thing += coeff * temp_contraction;
    }
  }

  void write_correlator(bool full_correlator,
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
			LatticeComplex & baryon)
  {
    if(full_correlator == true)
    {
#ifdef BUILD_HDF5
      std::string correlator_path = path+"/"+baryon_name+"/spin_"+spin+"/4D_correlator";
      h5writer.push(correlator_path);
      correlator_path = correlator_path+"/x"+std::to_string(source_coords[0])+"_y"+std::to_string(source_coords[1])+"_z"+std::to_string(source_coords[2])+"_t"+std::to_string(source_coords[3]);
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
      std::string correlator_path = path+"/"+baryon_name+"/spin_"+spin+"/x"+std::to_string(source_coords[0])+"_y"+std::to_string(source_coords[1])+"_z"+std::to_string(source_coords[2])+"_t"+std::to_string(source_coords[3]);
      h5writer.push(correlator_path);
#else
      std::string correlator_path = baryon_name+"_spin-"+spin+"_x"+std::to_string(source_coords[0])+"_y"+std::to_string(source_coords[1])+"_z"+std::to_string(source_coords[2])+"_t"+std::to_string(source_coords[3]);
#endif
      for(int mom = 0; mom < FT.numMom(); mom++)
      {
	//One more temp variable instanited inside loop (once again for writing.)
	multi1d<Complex> baryon_correlator;
	baryon_correlator.resize(Nt);
	multi1d<int> momenta = FT.numToMom(mom);
#ifndef BUILD_HDF5
	std::string correlator_path_mom = correlator_path+"_px"+std::to_string(momenta[0])+"_py"+std::to_string(momenta[1])+"_pz"+std::to_string(momenta[2]);
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
	std::string correlator_path_mom = correlator_path+"/px"+std::to_string(momenta[0])+"_py"+std::to_string(momenta[1])+"_pz"+std::to_string(momenta[2]);
	h5writer.write(correlator_path_mom, baryon_correlator, h5mode);
	h5writer.writeAttribute(correlator_path_mom, "is_shifted", 1, h5mode);
	h5writer.cd("/");
#endif
      }
    }
  }

} // End namespace Chroma
