//Arjun Singh Gambhir
//Code ported from latscat, with many changes included as comments.
//Trying to get rid of this progs include everywhere.
//#include "progs.h"
#include "sinks.h"
//gausSmear function is needed, so we include this stuff from Chroma.
#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/gaus_quark_smearing.h"
#include "meas/smear/gaus_smear.h"

namespace Chroma {

        /*********************************************************
         *********************************************************  GAUSSIAN SINKS
         *********************************************************/


        // To apply a gauss sink in a position-invariant manner, we can take advantage of the fact that
        // smeared(x) = sum_r Gaussian(r) unsmeared(x+r)
        // This must be right, because if s is dirac_delta, we'd get smeared(x) = unsmeared(x), exactly as we expect.
        // Therefore, it can executed as a convolution, which can be accelerated with FFTs.
        // In fact, we can store the Gaussian smearing in Fourier space, doing the FFT only at construction time.
        gauss_sink::gauss_sink(const Real sigmasq, Fourier &FFT, const int& time_direction){ 
            START_CODE();
            
            sigmasquared=sigmasq; 
            fft=&FFT; 

            snk=zero;
            
            // Gaussians are generally written as 
            //    gaussian(r) = e^(- (r-µ)^2 / 2 σ^2) / sqrt(2π σ^2)^(spatial dimensions).
            // Because we're going to momentum space, the central position µ doesn't matter.
            // Then, we just need the displacement r^2 and the width σ^2 to construct everything:
            LatticeReal displacement=zero;
            Real norm = 1./pow(sqrt(Chroma::twopi * sigmasquared),Nd-1);
            
            // To build the gaussian, we need to populate the displacement vector with r^2.
            // We can build up r^2 direction-by-direction,
            LatticeReal r;
            for(int nu = 0; nu < Nd; ++nu){
                // skipping the time direction,
                if ( nu == time_direction ) continue;
                
                // by finding the coordinate distance from the origin,
                r = LatticeReal((LatticeInteger(Layout::latticeCoordinate(nu)) + Layout::lattSize()[nu])%Layout::lattSize()[nu]);
                // accounting for periodicity,
                r=where(r<LatticeReal(Layout::lattSize()[nu]/2),r,LatticeReal(Layout::lattSize()[nu])-r);
                // and squaring.
                displacement += pow(r, 2);
            
                // For convolution normalization, see derivation and NB below.
                norm = norm * sqrt(Real(Layout::lattSize()[nu]));
            }
            
            //Finally, the exponent is -r^2 / 2 σ^2
            LatticeReal exponent;
            exponent = (- 0.5 / sigmasquared ) * displacement ;
            
            // Now we can form the Gaussian
            snk = cmplx( norm * exp(exponent), 0);
            // and store it in momentum space, as its FFT.
            snk = (*fft)(snk,-1); // The sign is explained below.

            END_CODE();
        };

        // Now that we have the FFT of the Gaussian stored, to apply the smearing we just need to convolve it with the propagator.
        // To minimize communication, what's ideal is if the post-FFT things that need to be multiplied wind up on the same node.
        // Then, the best thing to do is to transform unsmeared forward, smearing backward, and then their product backward.
        // The proof follows.

        // Let FFT(f, sign) = 1/sqrt(V) sum_n e^{i sign m n} f(n), the unitary convention for fourier transforms, which is implemented in our Fourier objects.
        // Then, FFT is (up to the sign of the exponent) its own inverse and we can define
        //      s(x) = FFT( s(p), a)        s(p) = FFT( s(x), -a)
        //      φ(x) = FFT( φ(p), b)        φ(p) = FFT( φ(x), -b)
        // where a and b are signs.
        // Now,
        // smeared(x)   = sum_r s(r) φ(x+r)
        //              = 1/V sum_{pqr} e^{i p r a} s(p) e^{i q (x+r) b} φ(q)
        //              = sum_{pq}  s(p) φ(q) e^{i q x b}   1/V sum_r e^{i r (p a + q b)}
        //              = sum_{pq}  s(p) φ(q) e^{i q x b} dirac_delta(p a + q b)
        // Because a and b are signs we can rewrite       dirac_delta(p a + q b) as dirac_delta(p + q a b).
        // Executing the sum on p,
        //              = sum_q s(- a b q) φ(q) e^{i q x b}
        // Since we want the transformed s and φ to be easy to multiply (without communication), we want their arguments to be the same (ie. to be on the same lattice site).
        // So, we must want 1 == - a b.  Let a=1 and b=-1.
        //              = sum_q s(q) φ(q) e^{- i q x}
        //              = sqrt(V) FFT( s(q) φ(q), -1 )
        //              = sqrt(V) FFT( FFT( s(x), -1) FFT( φ(x), +1), -1 ).
        // So, we can execute the convolution trick using our Fourier objects, as long as we remember to multiply by sqrt(V).

        // For further proof that it is a convolution note that
        // sum_r s(r) φ(x+r)    =   sum_{y-x} s(y-x) φ(y)
        // by simply letting y=x+r.  But, since x is fixed and the sum is over all space, sum_{y-x} = sum_y.
        // sum_r s(r) φ(x+r)    =   sum_y s(y-x) φ(y), transparently a convolution.
        LatticePropagator gauss_sink::operator()(const LatticePropagator& prop) const{
            START_CODE();

            QDPIO::cout << "Applying gaussian sink..." << std::flush;

            LatticePropagator temp;
            LatticePropagator result;
            // Transform the propagator forward,
            temp=(*fft)(prop,   1);
            // Multiply by sink, which is stored FFTed backwards.
            temp*=snk;
            // And transform backward to position space.
            result=(*fft)(temp, -1);
            // Noto bene:  what happened to multiplying by sqrt(V) for normalization?
            //             We stuck it into the normalization of the gaussian, to save arithmetic operations.
            //             See the constructor.
            
            QDPIO::cout << "done!" << std::endl;
            
            END_CODE();

            return result; 
        }

        /*********************************************************
         *********************************************************  GAUGE INVARIANT GAUSSIAN SINK SMEARING
         *********************************************************/

        // Without gauge fixing we can still perform Gaussian smearing, but it requires knowledge of the gauge field.
        // The idea is to form the quantity d = σ^2 / N, where σ is the width and N the number of iterations.
        // With the 3D gauge-invariant laplacian ∆, we have
        //      smeared(x) = ( δ_{x,y} - d ∆_{x,y})^N unsmeared(y).
        //
        // Consider a single application of ( delta_{x,y} - d ∆_{x,y}) in the unit-gauge limit.
        // Then,
        //      smeared(x) = ( δ_{x,y} + d sum_µ - 2 δ_{x,y} + δ_{x+µ,y} + δ_{x-µ,y}) unsmeared(y)
        //                 = (1-6d)unsmeared(x) + d sum_{x: nearest neighbors of y} unsmeared(y)
        //                 = (1-6d)( δ_{x,y} + 1/(1/d - 6) sum_{x: nearest neighbors of y}) unsmeared(y)
        //                 = a unsmeared(x) + a alpha sum_{x: nearest neighbors of y} unsmeared(y)
        // where a = 1-6d and alpha = 1/(1/d - 6)

        // Sergey recommends holding alpha approximately fixed while trying to find good values of σ and N.
        // alpha is in some sense the amount of nearest neighbors we add---as we go to the continuum limit, 
        // the fluctuations of the quark fields on neighboring sites will be smoother,
        // so probably smaller alpha is desired nearer to the continuum.
             
        // When you construct, you need to pass σ, N, the time direction (so that you don't smear in that direction), and the gauge field.
        gauge_gauss_sink::gauge_gauss_sink(const Real sgma, const int n, const int time_direction, const multi1d<LatticeColorMatrix>& U){
            gauge_field  = &U;
            t_direction  = time_direction;
            sigma        = sgma;
            N            = n;
        }

        gauge_gauss_sink::gauge_gauss_sink(const double sgma, const int n, const int time_direction, const multi1d<LatticeColorMatrix>& U){
            gauge_field  = &U;
            t_direction  = time_direction;
            sigma        = Real(sgma);
            N            = n;
        }



        // You can construct this sink without passing the gauge field.  Just remember to pass it through set_gauge() later!
        gauge_gauss_sink::gauge_gauss_sink(const Real sgma, const int n, const int time_direction){
            t_direction  = time_direction;
            sigma        = sgma;
            N            = n;
        }

        gauge_gauss_sink::gauge_gauss_sink(const double sgma, const int n, const int time_direction){
            t_direction  = time_direction;
            sigma        = Real(sgma);
            N            = n;
        }


        // You can reset the gauge field, if need be.
        void gauge_gauss_sink::set_gauge(const multi1d<LatticeColorMatrix> &U){
            gauge_field  = &U;
        }

        // And finally, you can apply it!
        LatticePropagator gauge_gauss_sink::operator()(const LatticePropagator& prop) const{
            
            QDPIO::cout << "Applying gauge-invariant gaussian sink..." << std::flush;
            
            // Chroma actually has a built-in!  Handy!
            // You can find it in $CHROMA/lib/meas/smear/gaus_quark_smearing.cc
            //                    $CHROMA/lib/meas/smear/gaus_smear.cc
            LatticePropagator result = prop;
            gausSmear( *gauge_field, result, sigma, N, t_direction );
            
            QDPIO::cout << "done!" << std::endl;
            
            return result;
        }


        /*********************************************************
         *********************************************************  SINK PARSING
         *********************************************************/

        void parseSink(GroupXML_t& input, sink*& snk, Complex& weight, Fourier &FFT, const int& time_direction){
            QDPIO::cout << "Parsing a sink:" << std::endl;
            try{
                std::istringstream xml_sink(input.xml);
                XMLReader sink_top(xml_sink);
                std::string sink_type;
                read(sink_top,"sink_type", sink_type);
                QDPIO::cout << "    type:   " << sink_type << std::endl;
                Complex w;
                read(sink_top,"weight", w);
                QDPIO::cout << "    weight: " << w << std::endl;
                weight = w;
                
                if( "POINT" == sink_type ){
                    snk = new point_sink();
                }
                else if("GAUSS" == sink_type ){
                    Real sigmasq;
                    read(sink_top,"sigma_squared",sigmasq);
                    snk = new gauss_sink(sigmasq, FFT, time_direction);
                }
                else if("GAUGE_INV_GAUSS" == sink_type){
                    Real sigma;
                    int iterations;
                    read(sink_top,"sigma",sigma);
                    read(sink_top,"N",iterations);
                    snk = new gauge_gauss_sink(sigma,iterations,time_direction);
                }
                else{
                    QDPIO::cerr << "Unrecognized sink type " << sink_type <<"." << std::endl;
                    QDPIO::cerr << "Valid types: POINT GAUSS GAUGE_INV_GAUSS." << std::endl;
                    QDP_abort(1);
                }
            }
            catch(std::bad_cast){
                QDPIO::cerr << "parseSink: caught cast error" << std::endl;
                QDP_abort(1);
            }
            catch(std::bad_alloc)
            {
                // This might happen on any node, so report it
                std::cerr << "parseSink: caught bad memory allocation" << std::endl;
                QDP_abort(1);
            }
            catch(const std::string& e)
            {
                QDPIO::cerr << "parseSink: Caught Exception: " << e << std::endl;
                QDP_abort(1);
            }
            catch(std::exception& e)
            {
                QDPIO::cerr << "parseSink: Caught standard library exception: " << e.what() << std::endl;
                QDP_abort(1);
            }
            catch(...)
            {
                // This might happen on any node, so report it
                std::cerr << "parseSink: caught generic exception during gaugeInit" << std::endl;
                //QDP_abort(1);
                throw;
            }
            
        }

}
