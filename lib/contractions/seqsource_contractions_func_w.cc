/*
 * Helper routines for constructing sequential sources.
 * David Brantley
 * IMPORTANT NOTE: The quarks must be in the DIRAC BASIS for these routines to work.
   11-Feb-2019 AWL rolling dense spin loops into for loops
   NOTE: in contract1,2,3,4, all the spin sums can be cut in half since we do
         parity projection
   11-Feb-2019 AWL adding support for all time-slice seqsrc

   23-Feb-2019 AWL adding Ben's suggested TimeSliceFunc

 */

//CHROMA STUFF
#include "chromabase.h"

//LALIBE STUFF
#include "seqsource_contractions_func_w.h"
#include "util/ferm/diractodr.h"

namespace Chroma
{
    class TimeSliceSelector : public SetFunc
    {
    public:
        TimeSliceSelector(int dir, int _tslice):
            dir_decay(dir), tslice(_tslice) {}

        // return 1 if the coord is the required tslice, 0 otherwise
        int operator() (const multi1d<int>& coordinate) const
            {
                return coordinate[dir_decay] == tslice;
            }
        int numSubsets() const {return 2;}

        int dir_decay;
        int tslice;

        private:
            TimeSliceSelector() {}  // hide default constructor
    };
    /*
    class TimeSliceFunc : public SetFunc
    {
    public:
        TimeSliceFunc(int dir): dir_decay(dir) {}

        int operator() (const multi1d<int>& coordinate) const {return coordinate[dir_decay];}
        int numSubsets() const {return Layout::lattSize()[dir_decay];}

        int dir_decay;

        private:
            TimeSliceFunc() {}  // hide default constructor
    };
    */

    void rotate_to_Dirac_Basis(LatticePropagator & quark_to_be_rotated)
    {
        //copied from Chroma
        SpinMatrix U = DiracToDRMat();
        quark_to_be_rotated = adj(U)*quark_to_be_rotated*U;
    }

    void rotate_from_Dirac_Basis(LatticePropagator & quark_to_be_rotated)
    {
        SpinMatrix U = DiracToDRMat();
        quark_to_be_rotated = U*quark_to_be_rotated*adj(U);
    }

    SpinMatrix protonDiquarkSpin(int parity)
    {
        SpinMatrix g_one = 0.0;
        Complex one = 1.0;
        Complex srt = 1/sqrt(2);
        //Pick out upper or lower spin components based on whether we have a proton or neg_par proton.
        int spin_one = 0;
        int spin_two = 1;
        if(parity == 1)
        {
	        spin_one = 2;
	        spin_two = 3;
        }

        pokeSpin(g_one,one,spin_one,spin_two);
        pokeSpin(g_one,-one,spin_two,spin_one);
        g_one *= srt;

        return g_one;
    }

    SpinMatrix protonSpinUp(int parity)
    {
        SpinMatrix g_one = 0.0;
        Complex one = 1.0;
        //Pick out upper or lower spin components based on whether we have a proton or neg_par proton.
        int spin_one = 0;
        int spin_two = 1;
        if(parity == 1)
        {
	        spin_one = 2;
	        spin_two = 3;
        }

        pokeSpin(g_one,one,spin_one,spin_one);
        return g_one;
    }

    SpinMatrix protonSpinDn(int parity)
    {
        SpinMatrix g_one = 0.0;
        Complex one = 1.0;
        //Pick out upper or lower spin components based on whether we have a proton or neg_par proton.
        int spin_one = 0;
        int spin_two = 1;
        if(parity == 1)
        {
	        spin_one = 2;
	        spin_two = 3;
        }

        pokeSpin(g_one,one,spin_two,spin_two);
        return g_one;
    }

    /* Helper functions for doing color matrix contractions. */

    /*  Color contract a pair of lattice color matrix objects
        eps_ijk * eps_lmn Q1[j,l ; spin] * Q2[i,m ; spin]
    */
    LatticeColorMatrix dblEpsContract1
    (
        LatticeColorMatrix& cmat1,
        LatticeColorMatrix& cmat2
    )
    {
        LatticeColorMatrix result;

        LatticeComplex tmp1;
        // 0,0
        tmp1  =  peekColor(cmat1,2,1)*peekColor(cmat2,1,2);
        tmp1 += -peekColor(cmat1,2,2)*peekColor(cmat2,1,1);
        tmp1 += -peekColor(cmat1,1,1)*peekColor(cmat2,2,2);
        tmp1 +=  peekColor(cmat1,1,2)*peekColor(cmat2,2,1);
        pokeColor(result,tmp1,0,0);
        // 0,1
        tmp1  = -peekColor(cmat1,2,0)*peekColor(cmat2,1,2);
        tmp1 +=  peekColor(cmat1,2,2)*peekColor(cmat2,1,0);
        tmp1 +=  peekColor(cmat1,1,0)*peekColor(cmat2,2,2);
        tmp1 += -peekColor(cmat1,1,2)*peekColor(cmat2,2,0);
        pokeColor(result,tmp1,0,1);
        // 0,2
        tmp1  =  peekColor(cmat1,2,0)*peekColor(cmat2,1,1);
        tmp1 += -peekColor(cmat1,2,1)*peekColor(cmat2,1,0);
        tmp1 += -peekColor(cmat1,1,0)*peekColor(cmat2,2,1);
        tmp1 +=  peekColor(cmat1,1,1)*peekColor(cmat2,2,0);
        pokeColor(result,tmp1,0,2);
        // 1,0
        tmp1  = -peekColor(cmat1,2,1)*peekColor(cmat2,0,2);
        tmp1 +=  peekColor(cmat1,2,2)*peekColor(cmat2,0,1);
        tmp1 +=  peekColor(cmat1,0,1)*peekColor(cmat2,2,2);
        tmp1 += -peekColor(cmat1,0,2)*peekColor(cmat2,2,1);
        pokeColor(result,tmp1,1,0);
        // 1,1
        tmp1  =  peekColor(cmat1,2,0)*peekColor(cmat2,0,2);
        tmp1 += -peekColor(cmat1,2,2)*peekColor(cmat2,0,0);
        tmp1 += -peekColor(cmat1,0,0)*peekColor(cmat2,2,2);
        tmp1 +=  peekColor(cmat1,0,2)*peekColor(cmat2,2,0);
        pokeColor(result,tmp1,1,1);
        // 1,2
        tmp1  = -peekColor(cmat1,2,0)*peekColor(cmat2,0,1);
        tmp1 +=  peekColor(cmat1,2,1)*peekColor(cmat2,0,0);
        tmp1 +=  peekColor(cmat1,0,0)*peekColor(cmat2,2,1);
        tmp1 += -peekColor(cmat1,0,1)*peekColor(cmat2,2,0);
        pokeColor(result,tmp1,1,2);
        // 2,0
        tmp1  =  peekColor(cmat1,1,1)*peekColor(cmat2,0,2);
        tmp1 += -peekColor(cmat1,1,2)*peekColor(cmat2,0,1);
        tmp1 += -peekColor(cmat1,0,1)*peekColor(cmat2,1,2);
        tmp1 +=  peekColor(cmat1,0,2)*peekColor(cmat2,1,1);
        pokeColor(result,tmp1,2,0);
        // 2,1
        tmp1  = -peekColor(cmat1,1,0)*peekColor(cmat2,0,2);
        tmp1 +=  peekColor(cmat1,1,2)*peekColor(cmat2,0,0);
        tmp1 +=  peekColor(cmat1,0,0)*peekColor(cmat2,1,2);
        tmp1 += -peekColor(cmat1,0,2)*peekColor(cmat2,1,0);
        pokeColor(result,tmp1,2,1);
        // 2,2
        tmp1  =  peekColor(cmat1,1,0)*peekColor(cmat2,0,1);
        tmp1 += -peekColor(cmat1,1,1)*peekColor(cmat2,0,0);
        tmp1 += -peekColor(cmat1,0,0)*peekColor(cmat2,1,1);
        tmp1 +=  peekColor(cmat1,0,1)*peekColor(cmat2,1,0);
        pokeColor(result,tmp1,2,2);

        return result;
    }

    // Color contract a pair of lattice color matrix objects eps_ijk * eps_lmn quark_1^jm * quark_2^il.
    LatticeColorMatrix dblEpsContract2
    (
        LatticeColorMatrix& cmat1,
        LatticeColorMatrix& cmat2
    )
    {
        LatticeColorMatrix result;

        LatticeComplex tmp1 = peekColor(cmat1,2,2)*peekColor(cmat2,1,1);
        tmp1 += -peekColor(cmat1,2,1)*peekColor(cmat2,1,2);
        tmp1 += -peekColor(cmat1,1,2)*peekColor(cmat2,2,1);
        tmp1 += peekColor(cmat1,1,1)*peekColor(cmat2,2,2);

        pokeColor(result,tmp1,0,0);

        tmp1 = -peekColor(cmat1,2,2)*peekColor(cmat2,1,0);
        tmp1 += peekColor(cmat1,2,0)*peekColor(cmat2,1,2);
        tmp1 += peekColor(cmat1,1,2)*peekColor(cmat2,2,0);
        tmp1 += -peekColor(cmat1,1,0)*peekColor(cmat2,2,2);

        pokeColor(result,tmp1,0,1);

        tmp1 = peekColor(cmat1,2,1)*peekColor(cmat2,1,0);
        tmp1 += -peekColor(cmat1,2,0)*peekColor(cmat2,1,1);
        tmp1 += -peekColor(cmat1,1,1)*peekColor(cmat2,2,0);
        tmp1 += peekColor(cmat1,1,0)*peekColor(cmat2,2,1);

        pokeColor(result,tmp1,0,2);

        tmp1 = -peekColor(cmat1,2,2)*peekColor(cmat2,0,1);
        tmp1 += peekColor(cmat1,2,1)*peekColor(cmat2,0,2);
        tmp1 += peekColor(cmat1,0,2)*peekColor(cmat2,2,1);
        tmp1 += -peekColor(cmat1,0,1)*peekColor(cmat2,2,2);

        pokeColor(result,tmp1,1,0);

        tmp1 = peekColor(cmat1,2,2)*peekColor(cmat2,0,0);
        tmp1 += -peekColor(cmat1,2,0)*peekColor(cmat2,0,2);
        tmp1 += -peekColor(cmat1,0,2)*peekColor(cmat2,2,0);
        tmp1 += peekColor(cmat1,0,0)*peekColor(cmat2,2,2);

        pokeColor(result,tmp1,1,1);

        tmp1 = -peekColor(cmat1,2,1)*peekColor(cmat2,0,0);
        tmp1 += peekColor(cmat1,2,0)*peekColor(cmat2,0,1);
        tmp1 += peekColor(cmat1,0,1)*peekColor(cmat2,2,0);
        tmp1 += -peekColor(cmat1,0,0)*peekColor(cmat2,2,1);

        pokeColor(result,tmp1,1,2);

        tmp1 = peekColor(cmat1,1,2)*peekColor(cmat2,0,1);
        tmp1 += -peekColor(cmat1,1,1)*peekColor(cmat2,0,2);
        tmp1 += -peekColor(cmat1,0,2)*peekColor(cmat2,1,1);
        tmp1 += peekColor(cmat1,0,1)*peekColor(cmat2,1,2);

        pokeColor(result,tmp1,2,0);

        tmp1 = -peekColor(cmat1,1,2)*peekColor(cmat2,0,0);
        tmp1 += peekColor(cmat1,1,0)*peekColor(cmat2,0,2);
        tmp1 += peekColor(cmat1,0,2)*peekColor(cmat2,1,0);
        tmp1 += -peekColor(cmat1,0,0)*peekColor(cmat2,1,2);

        pokeColor(result,tmp1,2,1);

        tmp1 = peekColor(cmat1,1,1)*peekColor(cmat2,0,0);
        tmp1 += -peekColor(cmat1,1,0)*peekColor(cmat2,0,1);
        tmp1 += -peekColor(cmat1,0,1)*peekColor(cmat2,1,0);
        tmp1 += peekColor(cmat1,0,0)*peekColor(cmat2,1,1);

        pokeColor(result,tmp1,2,2);

        return result;
    }

    // Color contract a pair of lattice color matrix objects eps_ijk * eps_lmn quark_1^jl * quark_2^kn.
    LatticeColorMatrix dblEpsContract3
    (
        LatticeColorMatrix& cmat1,
        LatticeColorMatrix& cmat2
    )
    {
        LatticeColorMatrix result;

        LatticeComplex tmp1 = -peekColor(cmat1,1,1)*peekColor(cmat2,2,2);
        tmp1 += peekColor(cmat1,1,2)*peekColor(cmat2,2,1);
        tmp1 += peekColor(cmat1,2,1)*peekColor(cmat2,1,2);
        tmp1 += -peekColor(cmat1,2,2)*peekColor(cmat2,1,1);

        pokeColor(result,tmp1,0,0);

        tmp1 = peekColor(cmat1,1,0)*peekColor(cmat2,2,2);
        tmp1 += -peekColor(cmat1,1,2)*peekColor(cmat2,2,0);
        tmp1 += -peekColor(cmat1,2,0)*peekColor(cmat2,1,2);
        tmp1 += peekColor(cmat1,2,2)*peekColor(cmat2,1,0);

        pokeColor(result,tmp1,0,1);

        tmp1 = -peekColor(cmat1,1,0)*peekColor(cmat2,2,1);
        tmp1 += peekColor(cmat1,1,1)*peekColor(cmat2,2,0);
        tmp1 += peekColor(cmat1,2,0)*peekColor(cmat2,1,1);
        tmp1 += -peekColor(cmat1,2,1)*peekColor(cmat2,1,0);

        pokeColor(result,tmp1,0,2);

        tmp1 = peekColor(cmat1,0,1)*peekColor(cmat2,2,2);
        tmp1 += -peekColor(cmat1,0,2)*peekColor(cmat2,2,1);
        tmp1 += -peekColor(cmat1,2,1)*peekColor(cmat2,0,2);
        tmp1 += peekColor(cmat1,2,2)*peekColor(cmat2,0,1);

        pokeColor(result,tmp1,1,0);

        tmp1 = -peekColor(cmat1,0,0)*peekColor(cmat2,2,2);
        tmp1 += peekColor(cmat1,0,2)*peekColor(cmat2,2,0);
        tmp1 += peekColor(cmat1,2,0)*peekColor(cmat2,0,2);
        tmp1 += -peekColor(cmat1,2,2)*peekColor(cmat2,0,0);

        pokeColor(result,tmp1,1,1);

        tmp1 = peekColor(cmat1,0,0)*peekColor(cmat2,2,1);
        tmp1 += -peekColor(cmat1,0,1)*peekColor(cmat2,2,0);
        tmp1 += -peekColor(cmat1,2,0)*peekColor(cmat2,0,1);
        tmp1 += peekColor(cmat1,2,1)*peekColor(cmat2,0,0);

        pokeColor(result,tmp1,1,2);

        tmp1 = -peekColor(cmat1,0,1)*peekColor(cmat2,1,2);
        tmp1 += peekColor(cmat1,0,2)*peekColor(cmat2,1,1);
        tmp1 += peekColor(cmat1,1,1)*peekColor(cmat2,0,2);
        tmp1 += -peekColor(cmat1,1,2)*peekColor(cmat2,0,1);

        pokeColor(result,tmp1,2,0);

        tmp1 = peekColor(cmat1,0,0)*peekColor(cmat2,1,2);
        tmp1 += -peekColor(cmat1,0,2)*peekColor(cmat2,1,0);
        tmp1 += -peekColor(cmat1,1,0)*peekColor(cmat2,0,2);
        tmp1 += peekColor(cmat1,1,2)*peekColor(cmat2,0,0);

        pokeColor(result,tmp1,2,1);

        tmp1 = -peekColor(cmat1,0,0)*peekColor(cmat2,1,1);
        tmp1 += peekColor(cmat1,0,1)*peekColor(cmat2,1,0);
        tmp1 += peekColor(cmat1,1,0)*peekColor(cmat2,0,1);
        tmp1 += -peekColor(cmat1,1,1)*peekColor(cmat2,0,0);

        pokeColor(result,tmp1,2,2);

        return result;
    }

    /*  Contract two quark props as:
        i,j = color; a,b = spin
        DQ[k,n ; a,b] = e_ijk * Q1[j,l ; a,c] * Q2[i,m ; c,b] * e_lmn
        a,b = color; i,j = spin
        DQ[c,f ; i,j] = e_abc * Q1[b,d ; i,k] * Q2[a,e ; k,j] * e_def
    */
    LatticePropagator contract1
    (
        LatticePropagator& quark_1,
        LatticePropagator& quark_2
    )
    {
        LatticePropagator contracted_prop;

        LatticeColorMatrix tmp1;
        LatticeColorMatrix tmp2;
        LatticeColorMatrix tmp3;

        // s_i = spin_i
        // s_j = spin_k
        // s_k = summed spin
        for(int s_i=0; s_i < Ns; ++s_i)
        {
            for(int s_j=0; s_j < Ns; ++s_j)
            {
                tmp3 = 0;
                for(int s_k = 0; s_k < Ns; ++s_k)
                {
                    tmp1  = peekSpin(quark_1,s_i,s_k);
                    tmp2  = peekSpin(quark_2,s_k,s_j);
                    tmp3 += dblEpsContract1(tmp1,tmp2);
                }
                pokeSpin(contracted_prop,tmp3,s_i,s_j);
            }
        }

        return contracted_prop;
    }// Done contract1

    /*  Contract two quark props as:
        i,j = color; a,b = spin
        DQ[k,n ; a,b] = e_ijk * Q1[j,m ; a,b] * Q2[i,l ; c,c] * e_lmn
        a,b = color; i,j = spin
        DQ[c,f ; i,j] = e_abc * Q1[b,e ; i,j] * Q2[a,d ; k,k] * e_def
    */
    LatticePropagator contract2
    (
        LatticePropagator& quark_1,
        LatticePropagator& quark_2
    )
    {
        LatticePropagator contracted_prop;

        LatticeColorMatrix tmp1;
        LatticeColorMatrix tmp2;
        LatticeColorMatrix tmp3;

        tmp2 = traceSpin(quark_2);
        for(int s_i=0; s_i < Ns; ++s_i)
        {
            for(int s_j=0; s_j < Ns; ++s_j)
            {
                tmp1 = peekSpin(quark_1,s_i,s_j);
                tmp3 = dblEpsContract2(tmp1,tmp2);
                pokeSpin(contracted_prop,tmp3,s_i,s_j);
            }
        }

        return contracted_prop;
    } // end contract2

    /*  Contract two quark props as:
        i,j = color; a,b = spin
        DQ[k,n ; a,b] = e_ijk * Q1[j,m ; a,c] * Q2[i,l ; c,b] * e_lmn
        a,b = color; i,j = spin
        DQ[c,f ; i,j] = e_abc * Q1[b,e ; i,k] * Q2[a,d ; k,j] * e_def
    */
    LatticePropagator contract3
    (
        LatticePropagator& quark_1,
        LatticePropagator& quark_2
    )
    {
        LatticePropagator contracted_prop;

        LatticeColorMatrix tmp1;
        LatticeColorMatrix tmp2;
        LatticeColorMatrix tmp3;

        for(int s_i=0; s_i < Ns; ++s_i)
        {
            for(int s_j=0; s_j < Ns; ++s_j)
            {
                tmp3 = 0;
                for(int s_k=0; s_k < Ns; ++s_k)
                {
                    tmp1  = peekSpin(quark_1,s_i,s_k);
                    tmp2  = peekSpin(quark_2,s_k,s_j);
                    tmp3 += dblEpsContract2(tmp1,tmp2);
                }
                pokeSpin(contracted_prop,tmp3,s_i,s_j);
            }
        }
        return contracted_prop;
    }// end contract3

    /*  Contract two quark props as:\
        i,j = color; a,b = spin
        DQ[i,m ; a,b] = e_ijk * Q1[j,l ; a,c] * Q2[k,n ; c,b] * e_lmn
        a,b = color; i,j = spin
        DQ[a,e ; i,j] = e_abc * Q1[b,d ; i,k] * Q2[c,f ; k,j] * e_def
    */
    LatticePropagator contract4
    (
        LatticePropagator& quark_1,
        LatticePropagator& quark_2
    )
    {
        LatticePropagator contracted_prop;

        LatticeColorMatrix tmp1;
        LatticeColorMatrix tmp2;
        LatticeColorMatrix tmp3;

        for(int s_i=0; s_i < Ns; ++s_i)
        {
            for(int s_j=0; s_j < Ns; ++s_j)
            {
                tmp3 = 0;
                for(int s_k=0; s_k < Ns; ++s_k)
                {
                    tmp1  = peekSpin(quark_1,s_i,s_k);
                    tmp2  = peekSpin(quark_2,s_k,s_j);
                    tmp3 += dblEpsContract3(tmp1,tmp2);
                }
                pokeSpin(contracted_prop,tmp3,s_i,s_j);
            }
        }
        return contracted_prop;
    }//end contract4

/* Routines necessary for taking the fully contracted sequential propagator and returning the proper time slice with the proper phase.*/

    LatticeComplex singlePhase
    (
        multi1d<int>& mom,
        multi1d<int>& origin_off,
        int j_decay
    )
    {
        // Taken from sft routine.
	    const Real twopi = 6.283185307179586476925286;

        // Coordinates for sink momenta
        multi1d<LatticeInteger> my_coord(Nd);
        for (int mu=0; mu < Nd; ++mu)
            my_coord[mu] = Layout::latticeCoordinate(mu);

        LatticeReal p_dot_x ;
        p_dot_x = 0. ;

        int j = 0;
        for(int mu = 0; mu < Nd; ++mu)
        {
            if (mu == j_decay) continue ;
            p_dot_x += LatticeReal(my_coord[mu] - origin_off[mu])*twopi*Real(mom[j]) / Layout::lattSize()[mu];
            ++j ;
        } // end for(mu)

        LatticeComplex phase;
        phase = cmplx(cos(p_dot_x),sin(p_dot_x));

        return phase;
    } // end singlePhase

    // Function to ensure proper time ordering.
    void baryonTimeOrder
    (
        LatticeComplex& timeorder,
        int t_source,
        multi1d<int>& BC,
        int j_decay,
        int parity
    )
    {
        Complex one  = 1.0;
        Complex mone = -1.0;
        if (BC[j_decay] < 0)
        {
            QDPIO::cout << "LALIBE_SEQSOURCE: multiplying by (-) sign for terms going around the boundary"<<std::endl;
            if(parity == 0)
            {
                timeorder = where(QDP::Layout::latticeCoordinate(j_decay) >= t_source, one, mone);
            }
            else
            {
                timeorder = where(QDP::Layout::latticeCoordinate(j_decay) <= t_source, one, mone);
            }
        }
        else timeorder = one;

    } // end baryonTimeOrder
    // Project out time slice.

    void projectTimeSlice
    (
        LatticePropagator& source_prop,
        int t_sink,
        int j_decay
    )
    {
        // Make the time slice set.
        Set timeslice;
        //timeslice.make(TimeSliceFunc(j_decay));
        timeslice.make(TimeSliceSelector(j_decay, t_sink));

        // Select only the t_sink components, all others are zero.
        //source_prop[timeslice[t_sink]] = source_prop;
        source_prop[timeslice[0]] = zero;
    }

    LatticePropagator projectBaryonSeqSource
    (
        LatticePropagator& seq_source,
        multi1d<int>& mom,
        multi1d<int>& origin_off,
        int t_sink,
        int j_decay,
        multi1d<int>& BC,
        int parity,
        bool t_all
    )
    {
        int Nt = QDP::Layout::lattSize()[j_decay];
        // We presume the seq_source has already been g5-herm conjugated.
        int t_source = origin_off[j_decay];

        // phase.
        LatticeComplex phase = singlePhase(mom,origin_off,j_decay);

        /*  Since we will add multiple sinks together to form a coherent sink,
            we want to put in the correct (-) sign associated with going through
            the boundary when making these sequential sources, so that they can
            be added together and used to compute 3pt functions without having
            to add the appropriate (-) signs when constructing the correlation
            function
        */
        //Complex one  = 1.0;
        LatticeComplex timeorder;
        baryonTimeOrder(timeorder, t_source, BC, j_decay, parity);

        //int Nt = QDP::Layout::lattSize()[j_decay];
        /*
        QDPIO::cout << "DEBUG timeorder: parity "<<parity<<std::endl;
        QDPIO::cout << "DEBUG timeorder: ";
        for (int t=0; t<Nt; ++t){
            multi1d<int> peek_time(4);
            peek_time[0] = peek_time[1] = peek_time[2] = 0;
            peek_time[3] = t;
            QDPIO::cout<<
                QDP::toDouble(real(peekSite(timeorder,peek_time)))
                <<" ";
        }
        QDPIO::cout<<std::endl;
        */

        LatticePropagator seq_source_projected = timeorder * phase * seq_source;
        if( ! t_all )
        {
            QDPIO::cout << "Selecting only t_sink = "<<t_sink<<std::endl;
            projectTimeSlice(seq_source_projected, t_sink, j_decay);
        }
        return seq_source_projected;
    }

}  // end namespace Chroma
