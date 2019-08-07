//Arjun Singh Gambhir
//Code ported from latscat.
#ifndef _FOURIER
#define _FOURIER

//necessary stuff
#define MINIMUM(a,b) (a<b ? a : b)

#ifndef CUFFT
#include <fftw3.h>
typedef fftw_complex fft_complex;
typedef fftw_plan fft_plan;
#define FFT_DESTROY(plan) fftw_destroy_plan(plan)
#define FFT_ALLOC(data,size) data=reinterpret_cast<fft_complex*>(fftw_malloc(size))
#define FFT_FREE(data) fftw_free(data)
#define FFT_MEMSET(data,value,count) memset(data,value,count)
#define FFT_FORWARD FFTW_FORWARD
#define FFT_BACKWARD FFTW_BACKWARD
//#define PROFILE
#else
#include "cudautils.h"
#include <cufft.h>
typedef cufftDoubleComplex fft_complex;
typedef cufftHandle fft_plan;
#define FFT_DESTROY(plan) cufftDestroy(plan)
#define FFT_ALLOC(data,size) CU_ALLOC_DEV(data, size)
#define FFT_ALLOC_HOST(data,size,flag) CU_ALLOC_HOST(data,size,flag)
#define FFT_FREE(data) CU_FREE_DEV(data)
#define FFT_MEMSET(data,value,count) CU_MEMSET(data,value,count)

//useful macros
#define FFT_FORWARD CUFFT_FORWARD
#define FFT_BACKWARD CUFFT_INVERSE
#define CUFFT_PLAN_MANY(plan,rank,dims,inembed,istride,idist,onembed,ostride,odist,type,batch){								\
	cufftResult status=cufftPlanMany(plan, rank, dims, inembed, istride, idist, onembed, ostride, odist, type, batch);		\
	if(status!=CUFFT_SUCCESS){																								\
		std::string errmsg="Cannot create cufft plan. Reason:";																\
		switch(status){																										\
			case CUFFT_ALLOC_FAILED: errmsg+="the plan could not be allocated."; break;										\
			case CUFFT_INVALID_VALUE: errmsg+="invalid parameter(s) passed to API."; break;									\
			case CUFFT_INTERNAL_ERROR: errmsg+="internal driver error."; break;												\
			case CUFFT_SETUP_FAILED: errmsg+="failed to initialize cuFFT library."; break;									\
			case CUFFT_INVALID_SIZE: errmsg+="parameter(s) do not have supported size."; break;								\
			default: errmsg+="unknown error.";																				\
		}																													\
		QDP_error_exit(errmsg.c_str());																						\
	}																														\
}																															\

#define CUFFT_EXEC(plan,idata,odata,sign){																					\
	cufftResult status=cufftExecZ2Z(plan, idata, odata, sign);																\
	if(status!=CUFFT_SUCCESS){																								\
		std::string errmsg="Cannot execute FFT. Reason: ";																	\
		switch(status){																										\
			case CUFFT_INVALID_PLAN: errmsg+="invalid plan."; break;														\
			case CUFFT_INVALID_VALUE: errmsg+="invalid parameter(s) passed to API."; break;									\
			case CUFFT_INTERNAL_ERROR: errmsg+="internal driver error."; break;												\
			case CUFFT_SETUP_FAILED: errmsg+="failed to initialize cuFFT library."; break;									\
			case CUFFT_EXEC_FAILED: errmsg+="failed to execute transform."; break;											\
			default: errmsg+="unknown error.";																				\
		}																													\
		QDP_error_exit(errmsg.c_str());																						\
	}																														\
}																															\

#endif

class Fourier{
private:
	//info about this node
	int my_node;
	multi1d<int> my_node_coords;
	int nsites; //number of sites in local volume
	multi1d<int> nsites_dir; //number of sites in local volume for each direction
	int nvsites; //number of volume sites

	//info about global topology
	int num_nodes;
	multi1d<int> nnodes_dir; //number of nodes in each direction
	multi1d<int> ldim; //global number of sites in each direction, lattice extends
	int vol; //lattice volume
	int tdir;
	int float_size;

	//node table used for routing
	multi1d< multi1d<int> > node_id_table;
	multi1d<QMP_comm_t> comms;
	
	//mapping layvout->lexi
	multi1d<int> reordermap;

	//plan data
	fft_plan plans[Nd][2];
	fft_complex* allvec;
	bool tuned, is3dtuned, printflops;
	unsigned int tunedvecsize;

	void init_node_table(){
		if(Nd!=4){
			QDP_error_exit("fourier::fourier: error, the dimension has to be 4!");
		}
		node_id_table.resize(Nd);

		//for all directions, check for the nodes which belong to the same coordinates perpendicular to direction
		for(unsigned int dd=0; dd<Nd; dd++){
			node_id_table[dd].resize(nnodes_dir[dd]);
			for(unsigned int l=0; l<nnodes_dir[dd]; l++){
				multi1d<int> coord=my_node_coords;
				coord[dd]=l;
				node_id_table[dd][l]=Layout::getNodeNumberFrom(coord);
			}
		}
		QMP_barrier();

		//create communicators
		comms.resize(Nd);
		for(unsigned int dd=0; dd<Nd; dd++){

			//determine new rank:
			unsigned int mynewrank=0;
			for(unsigned int l=0; l<nnodes_dir[dd]; l++){
				if(node_id_table[dd][l]==my_node) mynewrank=l;
			}

			//determine new color: all processes perpendicular to the dir direction get a unique color assignment
			unsigned int dira=(dd+1)%Nd;
			unsigned int dirb=(dd+2)%Nd;
			unsigned int dirc=(dd+3)%Nd;
			unsigned int mycolor=my_node_coords[dira]+nnodes_dir[dira]*(my_node_coords[dirb]+nnodes_dir[dirb]*my_node_coords[dirc]);

			//split new communicator from old one:
			QMP_comm_split(QMP_comm_get_default(),mycolor,mynewrank,&comms[dd]);
		}
	}
	
	
	//init reordering map:
	void init_reorder_map(){
                
		//reorder
		reordermap.resize(nsites);
                
		unsigned int run=0;
		for(int site=0; site<vol; ++site){
			multi1d<int> coord = crtesn(site, Layout::lattSize());
			int node = Layout::nodeNumber(coord);

			if(node==my_node){
				reordermap[run]=Layout::linearSiteIndex(coord);
				run++;
			}
		}
	}


	void printNodeTable(){
		QDPIO::cout << "Printing Node Table:" << std::endl;
		for(int n=0; n<num_nodes; n++){
			if(n==my_node){
				std::cout << "Node " << my_node << ":" << std::endl;
				std::cout << "(" << my_node << ")\tcoords (" << my_node_coords[0] << "," << my_node_coords[1] << "," << my_node_coords[2] << "," << my_node_coords[3] << ")" << std::endl;
				for(unsigned int dd=0; dd<Nd; dd++){
					std::cout << "(" << my_node << ")\tdir " << dd << ":";
					for(unsigned int l=0; l<nnodes_dir[dd]; l++){
						std::cout << " " << node_id_table[dd][l];
					}
					std::cout << std::endl;
				}
				std::cout<< std::endl;
			}
			QMP_barrier();
		}
	}

	//mimic MPI_Alltoall (global data transposition) in direction of dir:
	void alltoall_dir(char* dest, char* src, int ncount, int dir){
		QMP_status_t status=QMP_comm_alltoall(comms[dir],dest,src,ncount);
		if(status!=QMP_SUCCESS){
			QDP_error_exit("fourier: error in all to all routine!");
		}
	}

public:
 Fourier(const int& time_dir) : tuned(false), printflops(false){
		tdir=time_dir;
		float_size=sizeof(REAL);

		//set up info about this node:
		my_node=Layout::nodeNumber();
		my_node_coords=Layout::nodeCoord();
		nsites=Layout::sitesOnNode();
		nsites_dir=Layout::subgridLattSize();
		nvsites=1;
		for(unsigned int dd=0; dd<Nd; dd++){
			if(dd!=tdir) nvsites*=nsites_dir[dd];
		}

		//global infos:
		vol=Layout::vol();
		num_nodes=Layout::numNodes();
		nnodes_dir=Layout::logicalSize();
		ldim=Layout::lattSize();

		//OpenMP:
#ifndef CUFFT
		fftw_init_threads();
#endif

		//init node table:
		init_node_table();
		
		//init reordermap:
		init_reorder_map();
	}

	~Fourier(){
		for(unsigned int dd=0; dd<comms.size(); dd++){
			QMP_comm_free(comms[dd]);
		}
		if(tuned){
			//destroy plan
			for(unsigned int dir=0; dir<(is3dtuned ? Nd-1 : Nd); dir++ ){
				FFT_DESTROY( plans[dir][0] );
#ifndef CUFFT
				FFT_DESTROY( plans[dir][1] );
#endif
			}
			FFT_FREE(allvec);
#ifndef CUFFT
			fftw_cleanup_threads();
#endif
			is3dtuned=false;
			tunedvecsize=0;
		}
		tuned=false;
	}

	void print_flops(const bool& pflops){
		printflops=pflops;
	}

	//allows for creating a fixed plan and using it over and over:
	void tune(const unsigned int& vecsize, const unsigned int& is3d){
		tunedvecsize=vecsize;
		is3dtuned=is3d;

		//determine maximal vector size:
		unsigned int maxsize=1;
		for(unsigned int dir=0; dir<(is3dtuned ? Nd-1 : Nd); dir++ ){
			int dira= (dir+1)%Nd;
			int dirb= (dir+2)%Nd;
			int dirc= (dir+3)%Nd;

			int ntrafo= nsites_dir[dira]*nsites_dir[dirb]*nsites_dir[dirc]*vecsize;
			if ( ntrafo%nnodes_dir[dir] ){
				ntrafo= (ntrafo/nnodes_dir[dir]+1)*nnodes_dir[dir];
			}
			ntrafo/= nnodes_dir[dir];

			maxsize=(ntrafo*ldim[dir]>maxsize ? ntrafo*ldim[dir] : maxsize);
		}

		//allocate memory
		FFT_ALLOC(allvec,maxsize*sizeof(fft_complex));

		//create plans:
#ifndef CUFFT
		fftw_plan_with_nthreads(omp_get_max_threads());
#endif
		for(unsigned int dir=0; dir<(is3dtuned ? Nd-1 : Nd); dir++ ){
			int dira= (dir+1)%Nd;
			int dirb= (dir+2)%Nd;
			int dirc= (dir+3)%Nd;

			int ntrafo= nsites_dir[dira]*nsites_dir[dirb]*nsites_dir[dirc]*vecsize;
			if ( ntrafo%nnodes_dir[dir] ){
				ntrafo= (ntrafo/nnodes_dir[dir]+1)*nnodes_dir[dir];
			}
			ntrafo/= nnodes_dir[dir];

#ifndef CUFFT
			plans[dir][0]=fftw_plan_many_dft( 1, &ldim[dir], ntrafo, allvec, NULL, ntrafo, 1, allvec, NULL, ntrafo, 1, FFTW_FORWARD, FFTW_MEASURE );
			plans[dir][1]=fftw_plan_many_dft( 1, &ldim[dir], ntrafo, allvec, NULL, ntrafo, 1, allvec, NULL, ntrafo, 1, FFTW_BACKWARD, FFTW_MEASURE );
#else
			CUFFT_PLAN_MANY(&plans[dir][0], 1, &ldim[dir], NULL, ntrafo, 1, NULL, ntrafo, 1, CUFFT_Z2Z, 0);
#endif
		}
	}


	//buf has to be a buffer contaning REAL data and 2*vecsize elements (factor 2 comes from complex number) 
	void fourier_array(REAL* buf, const int& vecsize, const int& sign, const bool& is3d){
		unsigned int signid, flops;
		if( tuned && ( (is3d!=is3dtuned) || (vecsize!=tunedvecsize) ) ){
			QDP_error_exit("Fourier: error, you have tuned the FFT for a specific set of parameters but you try to apply it for a different set!");
		}
		else{
			if(sign==FFT_FORWARD) signid=0;
			else signid=1;
		}

#ifdef PROFILE
		StopWatch swatch_local_reord, swatch_global_reord, swatch_fft_core;
		swatch_local_reord.reset();
		swatch_global_reord.reset();
		swatch_global_reord.start();
		swatch_global_reord.stop();
		swatch_fft_core.reset();
#endif

		//OpenMP:
#ifndef CUFFT
		if(!tuned) fftw_plan_with_nthreads(omp_get_max_threads());
#endif
		for(unsigned int dir=0; dir<(is3d ? Nd-1 : Nd); dir++ ){
			int dira= (dir+1)%Nd;
			int dirb= (dir+2)%Nd;
			int dirc= (dir+3)%Nd;

			int ntrafo= nsites_dir[dira]*nsites_dir[dirb]*nsites_dir[dirc]*vecsize;
			if ( ntrafo%nnodes_dir[dir] ){
				ntrafo= (ntrafo/nnodes_dir[dir]+1)*nnodes_dir[dir];
			}
			ntrafo/= nnodes_dir[dir];

			fft_complex *vec, *hostvec;
			FFT_ALLOC(vec,ntrafo*ldim[dir]*sizeof(fft_complex));
#ifndef CUFFT
			hostvec=vec;
#else
			FFT_ALLOC_HOST(hostvec,ntrafo*ldim[dir]*sizeof(fft_complex),cudaHostAllocDefault);
#endif
			memset(reinterpret_cast<char*>(hostvec), 0, ntrafo*ldim[dir]*sizeof(fft_complex));

#ifdef PROFILE
			swatch_local_reord.start();
#endif
			for (unsigned int i=0; i<nsites; i++){
				//get local lexicographic coordinates
				multi1d<int> y=Layout::localLexiCoordFromLinear(i);

				//reorder, such that dir is fastest index after vector index
				int ind3= y[dira]+nsites_dir[dira]*(y[dirb]+nsites_dir[dirb]*y[dirc]);

				for (unsigned int v=0; v<vecsize; v++){
					int indv3= vecsize*ind3+v;
					int itrafo= indv3%ntrafo;
					int block= indv3/ntrafo;

					int ind= ntrafo*(nsites_dir[dir]*block+y[dir])+itrafo;
#ifndef CUFFT
					hostvec[ind][0]= buf[0+2*(v+vecsize*i)];
					hostvec[ind][1]= buf[1+2*(v+vecsize*i)];
#else
					memcpy(reinterpret_cast<char*>(&hostvec[ind]),reinterpret_cast<char*>(&buf[0+2*(v+vecsize*i)]),sizeof(fft_complex));
#endif
				}
			}
#ifdef CUFFT
			CU_MEMCPY(vec,hostvec,ntrafo*ldim[dir]*sizeof(fft_complex),cudaMemcpyHostToDevice);
#endif
#ifdef PROFILE
			swatch_local_reord.stop();
#endif

			/* transpose by alltoall */
			if ( nnodes_dir[dir]>1 ){
				if(!tuned) FFT_ALLOC(allvec,ntrafo*ldim[dir]*sizeof(fft_complex));
#ifdef PROFILE
				swatch_global_reord.start();
#endif
				alltoall_dir(reinterpret_cast<char*>(allvec), reinterpret_cast<char*>(vec), ntrafo*nsites_dir[dir]*sizeof(fft_complex), dir);
#ifdef PROFILE
				swatch_global_reord.stop();
#endif				
			}
			else{
				if(!tuned) allvec=vec;
				else{
#ifndef CUFFT
					memcpy(reinterpret_cast<char*>(allvec),reinterpret_cast<char*>(vec),ntrafo*ldim[dir]*sizeof(fft_complex));
#else
					CU_MEMCPY(allvec,vec,ntrafo*ldim[dir]*sizeof(fft_complex),cudaMemcpyDeviceToDevice);
#endif
				}
			}

			//create plan and execute fft
			if(!tuned){
				fft_plan plan;
#ifndef CUFFT
#ifdef PROFILE
				swatch_fft_core.start();
#endif
				plan=fftw_plan_many_dft( 1, &ldim[dir], ntrafo, allvec, NULL, ntrafo, 1, allvec, NULL, ntrafo, 1, sign, FFTW_ESTIMATE );
				fftw_execute( plan );
#ifdef PROFILE
				swatch_fft_core.stop();
#endif

				if(printflops){
					double add,mul,fma;
					fftw_flops(plan, &add, &mul, &fma);
					QDPIO::cout << "fourier: number of flops= " << add+mul+2*fma << std::endl; 
				}
#else
				int inembed[1]={ldim[dir]};
				int onembed[1]={ldim[dir]};
				CUFFT_PLAN_MANY(&plan, 1, &ldim[dir], inembed, ntrafo, 1, onembed, ntrafo, 1, CUFFT_Z2Z, ntrafo);
				CUFFT_EXEC(plan, allvec, allvec, sign);
				CU_DEVICE_SYNCHRONIZE();
#endif
				FFT_DESTROY( plan );
			}
			else{
#ifndef CUFFT
#ifdef PROFILE
				swatch_fft_core.start();
#endif
				fftw_execute(plans[dir][signid]);
#ifdef PROFILE
				swatch_fft_core.stop();
#endif
				
				if(printflops){
					double add,mul,fma;
					fftw_flops(plans[dir][signid], &add, &mul, &fma);
					QDPIO::cout << "fourier: number of flops= " << add+mul+2*fma << std::endl; 
				}
#else
				CUFFT_EXEC(plans[dir][0], allvec, allvec, sign);
				CU_DEVICE_SYNCHRONIZE();
#endif
			}

			//reorder if necessary by using alltoall
			if ( nnodes_dir[dir]>1 ){
#ifdef PROFILE
				swatch_global_reord.start();
#endif
				alltoall_dir(reinterpret_cast<char*>(vec), reinterpret_cast<char*>(allvec), ntrafo*nsites_dir[dir]*sizeof(fft_complex), dir);
				if(!tuned) FFT_FREE( allvec );
#ifdef PROFILE
				swatch_global_reord.stop();
#endif
			}
#ifndef CUFFT
			else if(tuned){
				memcpy(reinterpret_cast<char*>(vec),reinterpret_cast<char*>(allvec),ntrafo*ldim[dir]*sizeof(fft_complex));
			}
#else
			CU_MEMCPY(hostvec,vec,ntrafo*ldim[dir]*sizeof(fft_complex),cudaMemcpyDeviceToHost);
#endif			

			/* re-reorder with normalize */
#ifdef PROFILE
			swatch_local_reord.start();
#endif
			REAL norm= 1.0/sqrt(static_cast<REAL>(ldim[dir]));
			for (unsigned int i=0; i<nsites; i++){
				//get local lexicographic coordinates
				multi1d<int> y=Layout::localLexiCoordFromLinear(i);

				//reorder, such that dir is fastest index after vector index
				int ind3= y[dira]+nsites_dir[dira]*(y[dirb]+nsites_dir[dirb]*y[dirc]);

				for (unsigned int v=0; v<vecsize; v++){
					int indv3= vecsize*ind3+v;
					int itrafo= indv3%ntrafo;
					int block= indv3/ntrafo;

					int ind= ntrafo*(nsites_dir[dir]*block+y[dir])+itrafo;

#ifndef CUFFT
					buf[0+2*(v+vecsize*i)]=norm*hostvec[ind][0];
					buf[1+2*(v+vecsize*i)]=norm*hostvec[ind][1];
#else
					hostvec[ind]=cuCmul(hostvec[ind],make_cuDoubleComplex(norm,0.));
					memcpy(reinterpret_cast<char*>(&buf[0+2*(v+vecsize*i)]),reinterpret_cast<char*>(&hostvec[ind]),sizeof(fft_complex));
#endif
				}
			}
			FFT_FREE( vec );
#ifdef CUFFT
			CU_FREE_HOST( hostvec );
#endif
#ifdef PROFILE
			swatch_local_reord.stop();
#endif
		} // end dir-loop 
#ifndef CUFFT
		if(!tuned) fftw_cleanup_threads();
#endif
#ifdef PROFILE
		QDPIO::cout << "fourier: local reordering: " << swatch_local_reord.getTimeInSeconds() << std::endl;
		QDPIO::cout << "fourier: globl reordering: " << swatch_global_reord.getTimeInSeconds() << std::endl;
		QDPIO::cout << "fourier: fft core routine: " << swatch_fft_core.getTimeInSeconds() << std::endl;
#endif
	}
	
	

	template<class T>
	OLattice<T> operator()(const OLattice<T>& data, const int& isign, const bool& is3d=true){
		//sanity check and setting sign:
		int sign=0;
		if ( isign==+1 ) sign=FFT_FORWARD;
		else if ( isign==-1 ) sign=FFT_BACKWARD;
		else{
			QDP_error_exit("fourier: error, ivalid sign in fourier transform");
		}

		//first, flatten array of data using linear index, i.e. components,x,y,z,t
		int obj_size=sizeof(T)/(2*float_size); //elements of complex object
		REAL* buf=new REAL[nsites*2*obj_size];
		
		//convert to lexi
#pragma omp parallel for firstprivate(obj_size) shared(buf,data) default(shared)
		for(unsigned int run=0; run<nsites; run++){
			memcpy(reinterpret_cast<char*>(buf+run*2*obj_size),&(data.elem(reordermap[run])),sizeof(T));
		}

		//do fft on reordered array
		fourier_array(buf,obj_size,sign,is3d);

		//write back into lattice field:
		OLattice<T> result=zero;
#pragma omp parallel for firstprivate(obj_size) shared(buf,result) default(shared)
		for(unsigned int run=0; run<nsites; run++){
			memcpy(&(result.elem(reordermap[run])),reinterpret_cast<char*>(buf+run*2*obj_size),sizeof(T));
		}

		//clean up:
		delete [] buf;

		return result;
	} //end operator()
	
	
	//***********************************************************************************************************************************************
	//ATTENTION, THIS IS NOT TESTED		*************************************************************************************************************
	//***********************************************************************************************************************************************
	template<class T>
	void truncFFT(Complex* output, const OLattice<T>& data, const int& isign, multi1d<unsigned int> filter, const bool& is3d=true){
		//sanity check and setting sign:
		int sign=0;
		if ( isign==+1 ) sign=FFT_FORWARD;
		else if ( isign==-1 ) sign=FFT_BACKWARD;
		else{
			QDP_error_exit("fourier: error, ivalid sign in fourier transform!");
		}
		
		//first, flatten array of data using linear index, i.e. components,x,y,z,t
		unsigned int obj_size=filter.size(); //number of elements after filter is applied
		REAL* buf=new REAL[nsites*2*obj_size];
		
#pragma omp parallel firstprivate(obj_size,filter) shared(buf,data) default(shared)
		{
			unsigned int nthreads=MINIMUM(omp_get_num_threads(),nsites);
			unsigned int tid=omp_get_thread_num();
			unsigned int blocksize=static_cast<unsigned int>(ceil(static_cast<double>(nsites)/static_cast<double>(nthreads)));
			REAL* tmpbuf=new REAL[sizeof(T)/float_size];
			
			//convert to lexi
			for(unsigned int run=(tid*blocksize); run<MINIMUM(nsites,(tid+1)*blocksize); run++){
				//fill temporary buffer
				memcpy(reinterpret_cast<char*>(tmpbuf),&(data.elem(reordermap[run])),sizeof(T));
				for(unsigned int f=0; f<obj_size; f++){
					memcpy(reinterpret_cast<char*>(&buf[2*(f+obj_size*run)]),reinterpret_cast<char*>(&tmpbuf[2*filter[f]]),2*float_size);
				}
			}
			
			//clean up:
			delete [] tmpbuf;
		}
		
		//FFT the array
		fourier_array(buf,obj_size,sign,is3d);

		//write back into lattice field:
#pragma omp parallel firstprivate(obj_size,filter) shared(buf,output) default(shared)
		{
			unsigned int nthreads=MINIMUM(omp_get_num_threads(),nsites);
			unsigned int tid=omp_get_thread_num();
			unsigned int blocksize=static_cast<unsigned int>(ceil(static_cast<double>(nsites)/static_cast<double>(nthreads)));
		
			//convert to layout
			for(unsigned int run=(tid*blocksize); run<MINIMUM(nsites,(tid+1)*blocksize); run++){
				memcpy(reinterpret_cast<char*>(&(output[obj_size*reordermap[run]])),reinterpret_cast<char*>(&buf[2*obj_size*run]),obj_size*2*float_size);
				run++;
			}
		}

		//clean up:
		delete [] buf;
	} //end truncFFT

};
#endif
