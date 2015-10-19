#ifndef _NON_INT_RECOVER_GENERAL_H_
#define _NON_INT_RECOVER_GENERAL_H_

/*
 * @file nonint_recover_general.h
 * @brief NonIntRecover class which is used to assist Gradient class to do gradient descent
 * @author Tianwei
 * @date 06/25/2014
 *
 */


#include "projection.h"
#include "common.h"
#include "inv_lib.h"
#include "params.h"
#include <time.h>
#include <math.h>

/**
 * @brief This class NonIntRecover is used to do step 2 in light field reconstruction, which is to get the best possible values for each peak we get from step one to recover original signal.
 *
 * It has a subclass RecoverConfig which is used to represent a system like :
 * Ap * Xp = Xs. where Ap is the exponential matrix with size |S| * |P|, Xp is the coefficients/values we want to get with size |P| * 1, and Xs is the time domain sample values we have with a size of |S| * 1
 *
 * Here we choose Eigen template library to denote matrices and vectors. 
 *
 * There are two noticable changes of this version compared with previous version of permutation or projection implementation:
 * 1) We change the calculation from frequency domain to time domain; i.e. the samples we have will be values of time domain. 
 * 2) As we don't put any assumption on the input, in order to keep the precision, we will always keep all samples in calculation process. i.e. num_equations will always be the same as num_samples;
 *
 *
 * Note that in this specific implementation of algorithm, the time domain refers to the intermediate domain, and the frequency domain refers to the actual 4D frequency domain
 *
 * @author Tianwei
 * @date   06/25/2014
 */
class NonIntRecover {

friend class Gradient;

protected:

/**
 * @brief Configuration to a non-integer recovering
 *
 * @author Tianwei
 * @date   06/25/2014
 */



class RecoverConfig {

    /**
     * The exponential matrix, value in each entry is decided by the equation below
     * Ap((u,v), (Wu, Wv)) = exp(2 * j * pi * (u * Wu + v * Wv) / N) 
     *
     */
    CoeffMat Ap;

    /**
     * The vector that contains useful samples (ie. samples that are currently used in the system)
     */
    BukVec samples;

    /**
     * The current coefficient vector we try to recover.
     */
    CoeffVec coeff;

    /**
     * Each entry is that sample's index to the equations/row 
     */
    //int32_t * sample_index;     //removd as samples will always be in the equations and they will always follow the same order

    /**
     * Each entry is that peak's index to the columns.
     */
    //int32_t * peak_index;       //removed as peaks will always be in the same order

    /**
     * Instead of a peak_index, we need to store a peaklist
     */
    PeaksList peaks;   

    /**
     * The number of peaks in this recovering system
     */
     uint32_t num_peaks;

    /**
     * The total number of equations/samples; i.e.
     * the number of rows in RecoverConfig::Ap
     */
    uint32_t num_samples;

    /**
     * The positions of all samples
     */
    IntPosition * sample_pos;

    /**
     * The value of samples
     */ 
    ComplexPtr sample_values;

    /**
     * The maximum scaling; any scaling larger than that is viewed as singular
     *
     */
    double scaling_threshold;

#if defined(UPDATING_PSEUDOINVERSE)
    InvLib inv_lib;
#endif
	
	/**
     * The residue graph that stores all the residues of sample positions
     * 	
	 */
	 Real ** residue_graph;		

    public: 
	/**
 * Constructor.
 * Initialize scalars to be 0 and pointers to be 0.
	 */
	RecoverConfig() :
		peaks(0), num_peaks(0), num_samples(0), sample_pos(0),  sample_values(0), scaling_threshold(0), residue_graph(0) {}

	/**
	 * Destructor.
	 * Deallocate every dynamically allocated buffer
	 */
	~RecoverConfig() {
	deallocateBuffers();
	}

	/**
	 * setup num_samples, samples
	 */
	void setSamples(Params * params, ComplexPtr x_in) {

		num_samples = params->num_samples;
		samples.resize(num_samples, 1);
		// get sample positions
		if(sample_pos) {
			free(sample_pos);
			sample_pos = NULL;
		}
		sample_pos = (IntPosition *) malloc(sizeof(IntPosition) * num_samples);
		int32_t * cur_pos = params->sample_pos;
		for (uint32_t i = 0; i < num_samples; i ++) {
			sample_pos[i].setXY(cur_pos[2*i], cur_pos[2*i+1]);
		} 
		/*
		for (uint32_t i = 0; i < num_samples; i ++) {
			printf("sample pos %d, X: %d, Y: %d\n", i, sample_pos[i].getX(), sample_pos[i].getY());

		}
		*/
		//get sample values
		ComplexPtr writer = samples.data();
		for (uint32_t i = 0; i < num_samples; i ++) {
			writer[i] = x_in[sample_pos[i].map2Index()];
		}

		// allocate space for residue_graph
		residue_graph = (Real **) malloc(sizeof(Real *) * NonIntSFFT::n_v);
		for (uint32_t i = 0; i < NonIntSFFT::n_v; i ++) {
			residue_graph[i] = (Real *) malloc(sizeof(Real) * NonIntSFFT::n_h);
			for (uint32_t j = 0; j < NonIntSFFT::n_h; j ++) {
				residue_graph[i][j] = 0.0;
			}  
		}
	} 

	/**
     * Used by Gradient to do shadow bucket detection
     */
	void getResidueGraph(in out Real *** residue_graph_ptr) {
		/*
		for (uint32_t i = 0; i < NonIntSFFT::n_v; i ++) {
			for (uint32_t j = 0; j < NonIntSFFT::n_h; j ++) {
				fprintf(stderr, "%f ", residue_graph[i][j]);
			}
			fprintf(stderr, "\n");
		}
		*/
		*residue_graph_ptr = residue_graph;
		if (!(*residue_graph_ptr)) {
			fprintf(stderr, "Error residue graph not set\n");
			exit(-1);
		}
	}



#ifdef UPDATING_PSEUDOINVERSE
        /**
         * Get the InvLib field
         */
        const InvLib & getInvLib() const { return inv_lib; }

        /**
         * Get the InvLib field
         */
        InvLib & getInvLib() { return inv_lib; }
#endif

        /**
         * Replicate the recovery configuration
         *
         * @param[out]  target  the target recovery configuration
         */
        inline void replicate2New(out RecoverConfig * target) forceinline;

        /** 
         * This function resets the peaks in 
         * this class, and recover the values, and set the 
         * values in the peak lists in _gradient. If the rank
         * is deficient, this algorithm returns false.
         *
         * @param[in,out] peaks     the list of peaks
         * @param[in]   _num_peaks  the number of peaks
         * @param[out]  residue     the residue of the system
         *
         * @return  false if this is an under-determined system, or there is 
         *          some mistake happening;
         *          otherwise return true
         *
         */
        inline bool resetPeaks(in out PeaksList * peaks, in const uint32_t _num_peaks, out Real * residue) forceinline;

        /**
         * This function change the position of one peak, 
         * recover and set the values, and calculate the residue.
         *
         * @param[in,out]   peaks           the list of peaks
         * @param[in]       _num_peaks      the number of peaks
         * @param[in]       update_peak     the iterator to the peak that is updated
         * @param[in]       update_index    the index of the peak that is updated
         * @param[out]      residue         the residue of the system
         *
         * @return  false if this is an under-determined system, or there is 
         *          some mistake happening;
         *          otherwise return true
         */
        inline bool updatePeak(in out PeaksList * peaks, in const uint32_t _num_peaks, in PeakIterator update_peak, in uint32_t update_index, out Real * residue) forceinline;

        /**
         * This function remove the position of one peak,
         * recover and set the values, and calculate the residue.
         *
         * @param[in,out]   peaks           the list of peaks
         * @param[in]       new_num_peaks   the new number of peaks
         * @param[in]       remove_index    the index of the peak that is removed
         * @param[out]      residue         the residue of the system
         *
         * @return  false if this is an under-determined system, or there is 
         *          some mistake happening;
         *          otherwise return true
         */
        inline bool removePeak(in out PeaksList * peaks, in const uint32_t new_num_peaks, in uint32_t remove_index, out Real * residue) forceinline;

        /**
         * This function adds one new peak, 
         * recover and set the values, and calculate the residue.
         *
         * @param[in,out]   peaks           the list of peaks
         * @param[in]       new_num_peaks   the new number of peaks
         * @param[in]       new_peak        the iterator to the new peak
         * @param[in]       new_index       the index of the peak that is removed
         * @param[out]      residue         the residue of the system
         *
         * @return  false if this is an under-determined system, or there is 
         *          some mistake happening;
         *          otherwise return true
         */
        inline bool newPeak(in out PeaksList * peaks, in const uint32_t new_num_peaks, in PeakIterator new_peak, in uint32_t new_index, out Real * residue) forceinline;




    protected:

        /**
         * Deallocate the buffers: 
         * RecoverConfig::bucket_indicator,
         * RecoverConfig::bucket_index_all_coeff,
         */
        inline void deallocateBuffers() forceinline;

        /**
         * Write different rows into RecoverConfig::coeff from RecoverConfig::all_coeff
         * according to RecoverConfig::submat_indices
         */
        //inline void writeCoeff() forceinline;     //removed as we don't have similar structure any more

        /**
         * This function assumes that the RecoverConfig::all_coeff, 
         * RecoverConfig::coeff, RecoverConfig::b, RecoverConfig::bucket_b
         * are ready, just solve the system
         *
         * @param[out]  peaks   reset Peaks::v according to the solution of this system
         * @param[out]  residue the residue of the system
         *
         * @return the rank of the system
         */
        inline uint32_t solveSystem(out PeaksList * peaks, out Real * residue) forceinline;

        /**
         * This function update one column in RecoverConfig::all_coeff
         *
         * @param[in]   cur_coeff   pointer to that column
         * @param[in]   update_peak the PeakIterator to that peak
         * 
         * @return      the updated position of cur_coeff. It should be cur_coeff + num_bucket
         */
        inline ComplexPtr updateColumn(in PeakIterator update_peak, in ComplexPtr cur_coeff) forceinline;

        /**
         * Get the value of a column in current coefficient matrix
         *
         * @param[in]   peak_index      the index of the peak (column)
         * @param[out]  column_value    the buffer to store the column vector
         */
		inline void getColumnValue(IN const uint32_t & peak_index, OUT ComplexPtr column_value) forceinline {

			const Complex * matrix_ptr = Ap.data() + num_samples * peak_index;
			for (uint32_t equation_index = 0; equation_index < num_samples; ++ equation_index) {
				*(column_value ++) = matrix_ptr[equation_index];
			} 

		}

	
        /**
         * Get the value of a row in current coefficient matrix
         *
         * @param[in]   equation_index  the index of the equation (row)
         * @param[out]  row_value       the buffer to store the row vector
         
	inline void getRowValue(IN const uint32_t & equation_index, OUT ComplexPtr row_value) forceinline {
	    
	    const Complex * matrix_pre = all_Ap.data() + submat_indices[equation_index];
            for (int32_t peak_index = 0; peak_index < num_peaks; ++ peak_index, matrix_ptr += num_peaks)
                *(row_value ++) = *matrix_ptr;

	} */
	
	/**
	 * Checking if the InvLib has maintaining the right coefficient matrix
	 */
	 inline void checkInvLib() forceinline;

        /**
         * reset inv_lib to a new matrix
         */
	inline void resetInvLib() forceinline;


};


 /** 
     * A circular storage of all currently available RecoverConfig
     */
    RecoverConfig recoveries[NUM_RECOVER_BUFFER];

    /**
     * Beginning of NonIntRecover::recoveries
     */
    int32_t start;

    /**
     * End of NonIntRecover::recoveries
     * It is also the one that is under operation
     */
    int32_t end;

    /**
     * The array of IDs
     */
    int32_t ids[NUM_RECOVER_BUFFER];

    /**
     * The next ID to be allocated
     * The ID starts from zero and will increase
     */
    int32_t next_ID;

    /**
     * The inv_lib's statistics.
     * It is just for debugging. 
     */
    InvOpStatistics statistics;

public: 
    
    /**
     * Constructor that initialize the circular array data
     * structure
     */
    NonIntRecover() : start(0), end(-1), next_ID(0) {}

    /**
     * The deconstructor
	 */
    ~NonIntRecover() {}

    inline void setRecoverConfigSamples(Params * params, ComplexPtr x_in) {
		for(uint32_t i = 0; i < NUM_RECOVER_BUFFER; i ++) {
			recoveries[i].setSamples(params, x_in);
		}
    }


    /**
     * Get the statistic from inv_lib.
     * This function is for debugging
     */
    inline const InvOpStatistics & getOpStatistics() const forceinline {
	return statistics;
    }

    /**
     * Reset the statistics from inv_lib
     * This function is for debugging
     */
    inline void resetOpStatistics() forceinline {
        statistics.reset();
    }

    /**
     * This function either retrieve an old configuration or 
     * initialize a new one, depending on whether the old one
     * is still in the record. This configuration is set to be
     * the configuration at operation. 
     *
     * @param[in]   target_ID   the ID of the old configuration.
     *                          -1 means to initialize a new one
     * @param[out]  retrieve_ID the ID of the retrieved or initialized 
     *                          configuration
     * @param[in,out]   peaks       the list of peaks
     * @param[in]   num_peaks   the number of peaks
     * @param[out]  residue     the residue
     *
     * @return  false if this is an under-determined system, or there is 
     *          some mistake happing;
     *          otherwise return true
     */
    inline bool retrieveConfig(in int32_t target_ID, out int32_t * retrieve_ID, in out PeaksList * peaks, in const uint32_t num_peaks, out Real * residue) forceinline {

        if (target_ID >= 0 && ids[end] == target_ID) {
            *retrieve_ID = target_ID;
            return true;
        }

        if (target_ID >= 0 && (target_ID >= ids[start])) {
            // the reason for me to write this stupid low-efficient loop 
            // is the hope of the smart compiler to unroll it with -funroll-loops
            // option on
            for (uint32_t i = 0; i < NUM_RECOVER_BUFFER; i ++)
                if (ids[i] == target_ID) {
                    *retrieve_ID = target_ID;
                    end = i;


                    return true;
                }
        }
		
        // otherwise, let's just initialize a new one
        end = (end + 1) % NUM_RECOVER_BUFFER;

        if (end == start) 
            start = (start + 1) % NUM_RECOVER_BUFFER;

        *retrieve_ID = next_ID;
        ids[end] = next_ID;

        next_ID ++;

	return recoveries[end].resetPeaks(peaks, num_peaks, residue);
    }

    /**
     * This function replicate the configuration to a new configuration 
     * and set the new one at operation. 
     *
     * This is used to back up the current configuration for future backup
     *
     * @param[out]  old_ID  the ID of the old configuration
     * @param[out]  new_ID  the ID of the new configuration
     */
    inline void replicateConfig(out int32_t * old_ID, out int32_t * new_ID) {
        
        *old_ID = ids[end];
        
        RecoverConfig * old_config = recoveries + end;

        // initialize a new one
        end = (end + 1) % NUM_RECOVER_BUFFER;
        if (end == start)
            start = (start + 1) % NUM_RECOVER_BUFFER;
        *new_ID = next_ID;
        ids[end] = next_ID;
        next_ID ++;
        
        old_config->replicate2New(recoveries + end);
    }

    /**
     * This function update the configuration by providing an 
     * update of the peak
     *
     * @param[in,out]   peaks           the list of peaks
     * @param[in]       _num_peaks      the number of peaks
     * @param[in]       update_peak     the iterator to the peak that is updated
     * @param[in]       update_index    the index of the peak that is updated
     * @param[out]      residue         the residue of the system
     *
     * @return  false if this is an under-determined system, or there is 
     *          some mistake happing;
     *          otherwise return true
     */
    inline bool updateConfig(in out PeaksList * peaks, in const uint32_t _num_peaks, in PeakIterator update_peak, in uint32_t update_index, out Real * residue) forceinline {
        
        return recoveries[end].updatePeak(peaks, _num_peaks, update_peak, update_index, residue);
    }

    /**
     * This function update the configuration by providing a
     * removal of the peak
     *
     * @param[in,out]   peaks           the list of peaks
     * @param[in]       new_num_peaks   the new number of peaks
     * @param[in]       remove_index    the index of the peak that is removed
     * @param[out]      residue         the residue of the system
     *
     * @return  false if this is an under-determined system, or there is 
     *          some mistake happing;
     *          otherwise return true
     */
    inline bool shrinkConfig(in out PeaksList * peaks, in const uint32_t new_num_peaks, in uint32_t remove_index, out Real * residue) forceinline {
        
        return recoveries[end].removePeak(peaks, new_num_peaks, remove_index, residue);
    }

    /**
     * This function update the configuration by providing an 
     * addition of the peak
     *
     * @param[in,out]   peaks           the list of peaks
     * @param[in]       new_num_peaks   the new number of peaks
     * @param[in]       new_peak        the iterator to the new peak
     * @param[in]       new_index       the index of the peak that is removed
     * @param[out]      residue         the residue of the system
     *
     * @return  false if this is an under-determined system, or there is 
     *          some mistake happing;
     *          otherwise return true
     */
    inline bool expandConfig(in out PeaksList * peaks, in const uint32_t new_num_peaks, in PeakIterator new_peak, in uint32_t new_index, out Real * residue) forceinline {
        
        return recoveries[end].newPeak(peaks, new_num_peaks, new_peak, new_index, residue);
    }


};




	/*
	 * set values in Ap and fill in num_peaks
	 * 
	 * param[in] _peaks        pointer to the peaks list 
	 * param[in] _num_peaks    number of peaks and also should be equal to _peaks->size()
	 * param[in] residue       pointer to the residue
	 *
	 * @return      return true if the system is stil solvable
	 */ 
	
	inline bool NonIntRecover::RecoverConfig::resetPeaks(in out PeaksList * _peaks, in const uint32_t _num_peaks, out Real * residue) {

		if(!_peaks || !residue) {
		fprintf(stderr, "resetPeaks: null pointer");
		exit(-1);
		}

		
		if(_num_peaks == 0) {
		fprintf(stderr, "ERROR resetPeaks num_peaks = 0");
		//exit(-1);
		}
	   
		num_peaks = _num_peaks;
		peaks = *_peaks; // STL will help us set the size and everything
		
		//sanity check
		if(num_peaks != peaks.size()) {
		fprintf(stderr, "ERROR stl copy does not work properly");
		}
		
		Ap.resize(num_samples, num_peaks);
		coeff.resize(num_peaks, 1);
		ComplexPtr cur_Ap = Ap.data();
		
		for (PeakIterator it = _peaks->begin(); it != _peaks->end(); ++ it) {
			cur_Ap = updateColumn(it, cur_Ap); // note that here the meaning of updateColumn is to reset the values in cur_Ap
		}
	   
	#ifdef UPDATING_PSEUDOINVERSE
		resetInvLib(); 
		inv_lib.initPinv(); 
		checkInvLib();
	#endif

		return (solveSystem(_peaks, residue) == num_peaks);


	}

	/**
	 * update one column in Ap matrix
	 *
	 * @param[in]   update_peak  iterator to the peak that we want to update
	 * @param[in]   cur_Ap       pointer to the column that we want to update
	 *
	 * @return return the pointer to the next column
	 */
	inline ComplexPtr NonIntRecover::RecoverConfig::updateColumn(in PeakIterator update_peak, in ComplexPtr cur_Ap) {
		
		// fill in the Ap matrix according to the equations that Ap((u, v), (Wu, Wv)) = exp (2*j*pi*(u*Wu + v* Wv)/N)
		ComplexPtr writer = cur_Ap;
		for (uint32_t i = 0; i < num_samples; i ++) {
		*writer =  Utils::unityRootPower(sample_pos[i].getX() * update_peak->getX(), sample_pos[i].getY() * update_peak->getY()) / double(NonIntSFFT::n_v); 
			writer ++;
		}

		return writer;
	}


	/**
	 * Replicate current configuration	
	 *
	 * @param[out]   target   pointer to the new configuration
	 */

	inline void NonIntRecover::RecoverConfig::replicate2New(out RecoverConfig * target) {

		// copy arrays
		target->Ap.resize(num_samples, num_peaks);
		memcpy(target->Ap.data(), Ap.data(), sizeof(Complex) * num_samples * num_peaks);
		target->samples.resize(num_samples, 1);
		memcpy(target->samples.data(), samples.data(), sizeof(Complex) * num_samples);
		target->coeff.resize(num_peaks, 1);
		memcpy(target->coeff.data(), coeff.data(), sizeof(Complex) * num_peaks);
	 
		// copy scalar values
		target->num_samples = num_samples;
		target->num_peaks = num_peaks;
		target->scaling_threshold = scaling_threshold;

	   #ifdef UPDATING_PSEUDOINVERSE
		target->inv_lib.copyFrom(inv_lib);

		checkInvLib();
	   #endif


	}

	inline void NonIntRecover::RecoverConfig::checkInvLib() {

	#if defined(UPDATING_PSEUDOINVERSE) && defined(PSEUDOINVERSE_CHECKING_MODE)
		// actually nothing to do here
	#endif
	}

	/**
	 * Reset the InvLib
	 */
	
	inline void NonIntRecover::RecoverConfig::resetInvLib() {
	#ifdef UPDATING_PSEUDOINVERSE
		Complex * matrix = inv_lib.getMatrix();
		inv_lib.resize(num_samples, num_peaks);
		//printf("num_samples %d, num_peaks %d\n", num_samples, num_peaks);

		Complex * Ap_ptr = NULL;
		for (uint32_t equation_index = 0; equation_index < num_samples; ++ equation_index) {
			//all_coeff_ptr = all_coeff.data() + submat_indices[equation_index];
			Ap_ptr = Ap.data() + equation_index;
			for (uint32_t peak_index = 0; peak_index < num_peaks; ++ peak_index, ++ matrix, Ap_ptr += num_samples) {
				*matrix = *Ap_ptr;   
			}
		}
	#endif
	}

	/**
	 * Solve the system 
	 *
	 * @param[out]  peaks     pointer to the output peaks list
	 * @param[out]  residue   pointer to the output residue
	 *
	 * @return    number of peaks
	 */	


	inline uint32_t NonIntRecover::RecoverConfig::solveSystem(out PeaksList * peaks, out Real * residue) { 

		if (!peaks || !residue) {
			fprintf(stderr, "solveSystem: null pointer");
			exit(-1);
		}

	#if defined(PSEUDOINVERSE_CHECKING_MODE) || !defined(UPDATING_PSEUDOINVERSE)
		// First solve the system
		Decompose qr = Ap.decompose(decomp_args);
		uint32_t rank = qr.rank();
		BukVec solution = qr.solve(samples);
	#endif

	#if defined(UPDATING_PSEUDOINVERSE)
		// the inv lib
		inv_lib.solveSys(num_samples, samples.data(), num_peaks);
	#endif

		// Now set the values
		uint32_t peak_index = 0;


	#ifdef UPDATING_PSEUDOINVERSE
		const Complex * vs = inv_lib.getSolution();
	#else
		const Complex * vs = solution.data();
	#endif
	//
	#if defined(PSEUDOINVERSE_CHECKING_MODE) || !defined(UPDATING_PSEUDOINVERSE)
		// Check the singularity
		Real max_sample = 0.0;
		for(int32_t i = 0; i < num_samples; i ++) {
		if(comAbs(sample_value[i]) > max_sample) 
			max_sample = std::norm(sample_value[i]);
		}    
	 
		for(; peak_index < num_peaks; ++ vs, ++peak_index) {
			if (comAbs(*vs) > scaling_threshold * max_sample)
				rank = 0;
		}
		
	#endif

	#ifdef UPDATING_PSEUDOINVERSE
		vs = inv_lib.getSolution();
	#else
		vs = solution.data();
	#endif

		//vs = solution.data();

		peak_index = 0;
		// Note: why we don't use it != peaks.end() as the condition: 
		// this is just a simple trick (but somewhat ugly because the interface sucks)
		// that cooperates with Gradient::try2Remove(), enabling that in that function, 
		// we are NOT actually removing when evaluating whether we need to remove or not
		for (PeakIterator it = peaks->begin(); peak_index < num_peaks; ++ it, ++ vs, ++ peak_index)
			it->setValue(*vs);

	#ifdef UPDATING_PSEUDOINVERSE

		// time domain version of calculating residues

		ComplexPtr cur_sample;
		cur_sample = (ComplexPtr) malloc(sizeof(Complex) * num_samples);
		memset(cur_sample, 0, sizeof(Complex) * num_samples);
		vs = inv_lib.getSolution();

		Complex temp = 0.0;    
		const Complex * cur_vs = vs;
		ComplexPtr cur_Ap = Ap.data();
		// we need to get Ap * (*vs) and compare it with samples
		for (uint32_t i = 0; i < num_samples; i ++) {
			cur_Ap = Ap.data() + i;
			for(uint32_t j = 0; j < num_peaks; j ++) {
				temp += *(cur_vs) * *(cur_Ap);
				cur_vs ++;
				cur_Ap += num_samples;
			}
			cur_sample[i] = temp;
			temp = 0.0;
			cur_vs = vs;
		}
		
		// calculate residue
		*residue = 0;
		ComplexPtr actual_values = samples.data();
		int32_t pos_x, pos_y;
		for(uint32_t i = 0; i < num_samples; i ++) {
			*residue += std::norm( cur_sample[i]  - actual_values[i]);
			// we also need to fill this value into the residue_graph
			pos_x = sample_pos[i].getX();
			pos_y = sample_pos[i].getY();
			residue_graph[pos_x][pos_y] = std::norm( cur_sample[i] - actual_values[i]);
		}



		
		*residue = sqrt(*residue);

	#else

		// Finally, calculate the residue
		AllBukVec calc_samples;
		calc_samples.data().resize(num_samples, 1);
		calc_samples = Ap * solution;
		*residue = (bucket_cur - samples).norm();
	#endif

		if(cur_sample) {
			free(cur_sample);
			cur_sample = NULL;
		}


	#if defined(PSEUDOINVERSE_CHECKING_MODE) || !defined(UPDATING_PSEUDOINVERSE)
		return rank;
	#else
		return num_peaks;
	#endif


	}

	/**
	 * Update position of one of the peak and re-solve the system
     *
	 * @param[in]  _peaks        pointer to the peaks list
	 * @param[in]  _num_peaks    number of peaks that the system has
	 * @param[in]  _update_peak  iteraotr to the peak that will be updated
	 * @param[out] _residue      pointer to the output residue
	 *
     * @return      if the system is solvable or not
	 */

	inline bool NonIntRecover::RecoverConfig::updatePeak(in out PeaksList * _peaks, in uint32_t const _num_peaks, in PeakIterator update_peak, in uint32_t update_index, out Real * residue) {

		//sanity check 
		if(!_peaks || !residue) {
		fprintf(stderr, "updatePeak: null pointer");
		exit(-1);
		}

		if(_num_peaks != num_peaks) {
		fprintf(stderr, "updatePeak: number of peaks from %u to %d", num_peaks, _num_peaks);
		exit(-1);
		}

		if(update_index >= num_peaks) {
		fprintf(stderr, "updatePeak: num_peaks = %d, update_index = %d", num_peaks, update_index);
		exit(-1);
		}

		peaks = * _peaks;

		// we need to update Ap and re-solve the system to see if it is under-determined or not
		ComplexPtr col = Ap.data();
		col += update_index * num_samples;
		updateColumn(update_peak, col);

	#ifdef UPDATING_PSEUDOINVERSE
		getColumnValue(update_index, inv_lib.getColVec());
		inv_lib.updateCol(update_index, num_samples);
	#endif


	#ifdef PSEUDOINVERSE_CHECKING_MODE
		checkInvLib();
	#endif

		//solve the system
		return (solveSystem(_peaks, residue) == num_peaks);

	}

	/**
	 * Remove one peak from peaks list
	 *
	 * @param[in]   _peaks          pointer to the peaks list we have 
	 * @param[in]   new_num_peaks   number of peaks in the system after we remove one peak
	 * @param[in]   remove_index    index of the peak that we will remove
	 * @param[out]  residue         pointer to the residue of the system
	 *	
	 * @return       if the system is solvable or not
	 */

	inline bool NonIntRecover::RecoverConfig::removePeak(in out PeaksList * _peaks, in const uint32_t new_num_peaks, in uint32_t remove_index, out Real * residue) {

		// sanity check
		if (!_peaks || !residue) {
			fprintf(stderr, "removePeak: null pointer");
			exit(-1);
		}

		if (new_num_peaks + 1 != num_peaks) {
			fprintf(stderr, "removePeak: number of peaks from %u to %d", num_peaks, new_num_peaks);
			exit(-1);
		}

		peaks = * _peaks;
		num_peaks --;
		coeff.resize(num_peaks, Eigen::NoChange);

		ComplexPtr data = Ap.data();
		ComplexPtr Ap_ptr = data + remove_index * num_samples;

		for(; Ap_ptr != data + num_peaks * num_samples; Ap_ptr += num_samples)
			memcpy(Ap_ptr, Ap_ptr + num_samples, sizeof(Complex) * num_samples);

		Ap.resize(Eigen::NoChange, num_peaks);

	#ifdef UPDATING_PSEUDOINVERSE
		inv_lib.removeCol(remove_index);
	#endif

		if(num_peaks == 0) return false;
	 
		return (solveSystem(_peaks, residue) == num_peaks);
	}

	/**
	 * Add one new peak to peaks list
	 *
	 * @param[in]   _peaks          pointer to the peaks list we have 
	 * @param[in]   new_num_peaks   number of peaks in the system after we add one new peak
	 * @param[in]   new_peak        iterator to the peak that we will add to the system
	 * @param[in]   new_index       index of the peak that we will add
	 * @param[out]  residue         pointer to the residue of the system
	 *	
	 * @return       if the system is solvable or not
	 */

	inline bool NonIntRecover::RecoverConfig::newPeak(in out PeaksList * _peaks, in const uint32_t new_num_peaks, in PeakIterator new_peak, in uint32_t new_index, out Real * residue) {

		// sanity check
		if (!_peaks || !residue) {
			fprintf(stderr, "newPeak: null pointer");
			exit(-1);
		}

		if (new_num_peaks - 1 != num_peaks) {
			fprintf(stderr, "newPeak: number of peaks from %u to %d", num_peaks, new_num_peaks);
			exit(-1);
		}

		peaks = * _peaks;
		num_peaks ++;

		// resize coeff, Ap and fill in Ap's new column using updateColumn
		coeff.resize(num_peaks, Eigen::NoChange);
		Ap.resize(Eigen::NoChange, num_peaks);
		ComplexPtr data = Ap.data();
		data += (num_peaks - 1) * num_samples;
		updateColumn(new_peak, data);

	#ifdef UPDATING_PSEUDOINVERSE
		getColumnValue(new_index, inv_lib.getColVec());
		inv_lib.addCol(new_index, num_samples);
	#endif


		return (solveSystem(_peaks, residue) == num_peaks);
	}

	/**
	 * Deallocate buffers in one configuration
	 *
	 *
	 */

	inline void NonIntRecover::RecoverConfig::deallocateBuffers() {

		if(sample_pos) {
		free(sample_pos);
		sample_pos = NULL;
		}

		if(sample_values) {
		free(sample_values);
		sample_values = NULL;
		}

		if (residue_graph) {
			for (uint32_t i = 0; i < NonIntSFFT::n_v; i ++) {
				free(residue_graph[i]);
				residue_graph[i] = NULL;
			}
			free(residue_graph);
			residue_graph = NULL;
		}

	}



#endif
