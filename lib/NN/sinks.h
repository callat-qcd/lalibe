//Arjun Singh Gambhir
//Code ported from latscat, with many changes included as comments.
#ifndef _SINK_CLASSES
#define _SINK_CLASSES
#include "chromabase.h"
//We include chroma types this way.
#include "io/xml_group_reader.h"
//This is the GroupXML info needed from chroma, done the proper way.
#ifndef CUFFT
#include "fourier_cpu.h"
#else
#include "fourier_cuda.h"
#endif
//FFT needed.

namespace Chroma {

        /*********************************************************
         *********************************************************  SINKS
         *********************************************************/


        // Here we build an infrastructure for describing sinks.
        // Sinks are used for smearing propagators on the side opposite the source.

        class sink{
        public:
            // Sinks are immutable functions that are LatticePropagator --> LatticePropagator:
            virtual LatticePropagator operator()(const LatticePropagator&) const = 0;
            // They also can report their type.
            virtual std::string type() const = 0;
            // Other than that, they can essentially be whatever you want.
            // Of course, some sinks are more useful than others.
        };

        // We'll need a way to read XML and produce a sink of interest.
        void parseSink(GroupXML_t& input, sink*& snk, Complex& weight, Fourier &FFT, const int& time_direction);

        /*********************************************************
         *********************************************************  POINT SINK
         *********************************************************/

        // It's handy to have the ``identity'' sink, often called the point sink.
        // The point sink does nothing when applied --- it leaves the propagator totally alone.
        // ``Point'' stands in contrast to ``smear''.
        class point_sink : public sink {
        public:
            point_sink(){};
           ~point_sink(){};
            std::string type() const { return "POINT"; };
            LatticePropagator operator()(const LatticePropagator& prop) const{
                QDPIO::cout << "Applying point sink...done!" << std::endl;
                return prop; // See?  Do nothing!
                
            }
        };

        /*********************************************************
         *********************************************************  SMEARED SINKS
         *********************************************************/

        // A Gaussian-smeared sink combines the values on the nearby sites according to a (user-specifiable) Gaussian weight.
        class gauss_sink : public sink {
        public:
            gauss_sink(const Real sigmasq, Fourier &FFT, const int& time_direction);
           ~gauss_sink(){ return; };
            std::string type() const { return "GAUSS"; };
            LatticePropagator operator()(const LatticePropagator& prop) const;
        private:
            Real sigmasquared;      // Width of the gaussian.
            Fourier * fft;          // Needed for convolution.
            LatticeComplex snk;     // Implementation detail: store in momentum space for easy convolution.
        };

        // Chroma also has built-in gauge-invariant gauss smearing.
        // Without gauge fixing we can still perform Gaussian smearing, but it requires knowledge of the gauge field.
        // The idea is to form the quantity d = σ^2 / N, where σ is the width and N the number of iterations.
        // With the 3D gauge-invariant laplacian ∆, we have
        //      smeared(x) = ( δ_{x,y} - d ∆_{x,y})^N unsmeared(y).
        // Here we wrap the Chroma built-in in a sink class.
        // It carries more data than previously implemented sinks, in that you need to pass a gauge field at some point.
        // You can reset the gauge field, though.
        class gauge_gauss_sink : public sink {
        public:
            std::string type() const { return "GAUGE_INV_GAUSS"; };
            LatticePropagator operator()(const LatticePropagator& prop) const;
            gauge_gauss_sink(const double sgma, const int n, const int time_direction, const multi1d<LatticeColorMatrix>& U);
            gauge_gauss_sink(const Real sgma, const int n, const int time_direction, const multi1d<LatticeColorMatrix>& U);
            // You can construct without providing a gauge field...
            gauge_gauss_sink(const double sgma, const int n, const int time_direction);
            gauge_gauss_sink(const Real sgma, const int n, const int time_direction);
           ~gauge_gauss_sink(){ return; };
            // ... if you remember to pass a gauge field later.
            void set_gauge(const multi1d<LatticeColorMatrix> &U);

        private:
            Real sigma;                                      // Width of the gaussian.
            int N;                                           // Number of hits.
            const multi1d<LatticeColorMatrix> * gauge_field; // Gauge field for making the 3D laplacian.
            int  t_direction;                                // Direction to not smear.
        };

}
#endif
