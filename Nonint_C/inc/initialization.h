#ifndef _INITIALIZATION_H_
#define _INITIALIZATION_H_

/**
 * @file    initialization.h
 * @brief   Implement class Initialization, which contains the process of step one described in light field/SFFT paper
 * @author  Tianwei
 * @date    07/08/2014
 */

#include "common.h"
#include "projection.h"
#include "params.h"
#include <list>
#include "nonint_recover.h"
#include <unistd.h>
#include <stdio.h>
/**
 * @brief   
 *
 * This class implements the step one in our light field reconstruction
 *
 * @author  Tianwei
 * @date    07/08/2014
 */

class Initialization {

    /**
     * NonIntRecover is a friend class of this class
     */
    friend class IntRecover;
    friend class Projection;
	friend class Gradient; 

protected:

    /*
     * The number of projection slices
     */ 
     uint32_t num_projection_lines;

    /**
     * The list of projection slices that are used in this gradient algorithm
     */ 
     Projection * projections;

    /**
     * This records all the parameters that tunes this
     * algorithm
     */
    Params * params;

    /**
     * The time domain signal
     */
    ComplexPtr x_in;

    /**
     * The peaks that are constantly updated by the gradient algorithm. 
     *
     * Since the implemtation of peaks is crucial, 
     * here I would like to say why I want to use C++ STL list for peaks:
     *
     * 1) I would to simply the coding, so I don't want to use compliated 
     * data structure
     *
     * 2) We really don't need to do random indexing for peaks, so we 
     * don't need constant time locating
     *
     * 3) We are going to delete, and insert entries to peaks quite 
     * frequently
     *
     * 4) In order for future optimization, I define the type as a macro
     */
     PeaksList peaks;

    /**
     * The current number of peaks. 
     *
     * The reason why we are maintaining this is that I am not sure how
     * they implement std::list::size().
     */
    uint32_t num_peaks;

    /**
     * The residue of the current system
     */
    Real residue;

    /**
     * The integer recover class that is needed in peak estimation
     */
    IntRecover recovery;
    
    /**
     * This is a pre-allocated position table for voting algorithm
     */
    int32_t * position_table;

    /**
     * The current iteration number
     */
    uint32_t iter;

    /**
     * The fft plan in time-fixing 
     */
    fftw_plan plan;

    /**
     * The input to the time-fixing 2D fft, i.e., 
     * the fixed time domain
     */
    ComplexPtr fftw_input;

    /**
     * This is just an indicator of whether fftw_plan Gradient::plan
     * has been initialized or not.
     */
    bool init_plan;

    /**
     * The time-domain sampling table
     */
    bool *sampling_table;


    /**
	 * The initial peaks we have
	 */	
    PeaksList init_peaks;

	/**
	 * The initial residue calculated with initial peaks we have.
	 */
    Real init_residue;

	// following variables are for shadow bucket detection
	
	// The residue values we have, it has size NonIntSFFT::n_v * NonIntSFFT::n_h
	Real ** residue_graph;

	// The position of all samples we use which is corresponding to non zero positions in residue_graph
	int32_t ** sample_graph;
	
	// The peaks list we have, which contains (mostly) non integer positions and values of current peak list	
	PeaksList * current_peak_list;

	// The number of peaks in current_peak_list
	int32_t current_num_peaks;
	


public:

    /**
     * Constructor. Specifies the construction by 
     * the parameter params. 
     *
     * @param[in]  _params   the parameters for the gradient algorithm
     */
    Initialization(Params * _params) :  num_projection_lines(0), projections(0), params(_params), x_in(0), peaks(PeaksList()), num_peaks(0), residue(0), position_table(0), iter(0),  
        fftw_input(0), init_plan(0), sampling_table(0), init_peaks(PeaksList()), init_residue(0), residue_graph(0), current_peak_list(0), current_num_peaks(0) /*, added_peaks(0) */{

        setParameters(_params);  // we also allocate space for buffers here
        setRecovery();      
    }

    /**
     * Destructor: free the buffers, 
     * note that it will also free the buffers from the projections, 
     * i.e. invoking Filters::deallocateBuffers()
     */
    ~Initialization() {
	
	uint32_t projection_index = 0;
	for (; projection_index < num_projection_lines; ++ projection_index)
	    projections[projection_index].deallocateBuffers();
	
        deallocateBuffers();
        destroyPlan();

    }


    /**
     * This function sets all the parameters used in the algorithm. 
     * This mainly include allocating the Initialization::projections buffer, 
     * Initialization::position_table, Initialization::residue_window, 
     * Initialization::optimal_y, Initialization::fftw_input, Initialization::sampling_table, 
     * Initialization::scores_fake, Initialization::scores_missing, 
     * setting parameters to all of the projections, and set Initialization::num_projection_lines. 
     *
     * @param[in]  _params   the parameters for the gradient algorithm
     */
    void setParameters(Params * _params);


    /**
     * This function serves as the bridge from step one to step two where we pass peaks and residue information to step two. 
     * 
     * @param[in] _num_peaks   pointer to the num_peaks variable that we want to set
	 * @param[in] _peaks       pointer to the peaks variable that we use to store peak list
     * @param[in] _residue     pointer to the residue variable that we want to set
     */
    void getPeaksResidue(uint32_t * _num_peaks, PeaksList * _peaks, Real * _residue);



    /**
     * Set the initial peak positions for the gradient iterating algorithm. 
     * Which method to use depends on Params::InitParams::method.
     *
     * @param[in]   x_in        the time-domain signal
     * @param[out]  recovery_id the id to the recovery configurations
     */
    inline bool setInitState(in ComplexPtr _x_in, out int32_t * recovery_id) forceinline {

        Projection * cur_projection = projections;
        uint32_t projection_index = 0;

        for (; projection_index < num_projection_lines; ++ projection_index, ++ cur_projection) {
	    	cur_projection->projectionSlice(_x_in);
			//cur_projection->printPowers();
		}
        x_in = _x_in;	
		recovery.setBuckets();	
	
        bool success = true;
        num_peaks = 0;
        
        // Based on the projections, let's get the positions
        switch (params->init_param.method) {
        case Params::InitParams::VOTING:
            success = voting(recovery_id);
            break;
        
        default:
            exit(-1);
        }
		if(success) {
			init_peaks = peaks;
			init_residue = residue;
		}
	
        return success;

    }
	
	/** 
	 * Do voting on shadow buckets, i.e. try to find missing buckets
	 *
	 * param[out] added_peak  pointer to the list of peaks added
	 */


	bool shadowBucketVoting(out NonIntPeak * added_peak);


	/*
     * Fill in residue graph and sample graph to prepare for shadow bucket detection
	 *
	 * @param[in]  _residue_graph  pointer to the residue graph of current solution
	 * @param[in]  _sample_graph   pointer to the sample graph of current solution
     */
	
	inline void prepareShadowBucketDetection( Real ** _residue_graph, int32_t ** _sample_graph) {	
		for (uint32_t i = 0; i < NonIntSFFT::n_v; i ++) {
			for (uint32_t j = 0; j < NonIntSFFT::n_h; j ++) {
				residue_graph[i][j] = _residue_graph[i][j];
				sample_graph[i][j] = _sample_graph[i][j];
			}
		}
	}

protected:

    /**
     * Destroy the fftw plan Gradient::plan and set Gradient::init_plan back to 0
     */
    inline void destroyPlan() forceinline {

        if (init_plan) {
            fftw_destroy_plan(plan);
            init_plan = 0;
        }
    }

    /**
     * Set the NonIntRecover object recover
     */
    void setRecovery();

    /**
     * Deallocate the buffers:
     * Gradient::projections,
     * Gradient::position_table,
     * Gradient::residue_window,
     * Gradient::optimal_y, 
     * Gradient::fftw_input, 
     * Gradient::sampling_table
     */
    void deallocateBuffers();


    /**
     * The voting algorithm. 
     *
     * The following parameters tune this voting algorithm:
     * Params::InitParams::soft_power_threshold, 
     * Params::InitParams::hard_power_threshold
     *
     * @param[out]  recovery_id the id to the recovery configurations
     *
     * @return  whether we successfully find the initial state
     */
    inline bool voting(out int32_t * recovery_id) forceinline;

    /**
     * The optimism algorithm: 
     * check the initial values, if they have better residue, then
     * use them.
     *
     * @param[out]  recovery_id the id to the recovery configurations
     * @param[in]   previous    whether the old result is underdetermined or not
     *
     * @return  whether the new result is underdetermined
     */
    inline bool optimismInit(out int32_t * recovery_id, in bool previous) forceinline;

    /**
     * This function tries to vote for all projections, under
     * Projection::soft_power_threshold = cur_ratio, 
     *
     * @param[in]   cur_ratio   the current soft power threshold
     */
    inline void votingProposal(in double cur_ratio) forceinline;

    /**
     * Generate the time domain sampling table 
     * Gradient::sampling_table
     */
   inline void generateSamplingTable() forceinline;

};

inline bool Initialization::voting(out int32_t * recovery_id) {

    Projection * cur_projection = NULL;    
	
	// figure out the maximum ratio we can have 
    double max_ratio = INFINITY;
    Real max_power = 0;
    Real min_power = 0;
    double cur_ratio = 0;
    uint32_t projection_index = 0;
    for (projection_index = 0, cur_projection = projections; projection_index < num_projection_lines; ++ projection_index, ++ cur_projection) {
        max_power = cur_projection->getMaxPower();
        min_power = cur_projection->getMinPower();
        cur_ratio = max_power / min_power;
        if (cur_ratio < max_ratio)
            max_ratio = cur_ratio;
    }

    if (max_ratio > params->init_param.max_soft_power_threshold)
        max_ratio = params->init_param.max_soft_power_threshold;

    cur_ratio = params->init_param.soft_power_threshold;
    double ratio_step = (max_ratio - cur_ratio) / params->init_param.max_voting_iter;
    if (ratio_step < params->init_param.min_soft_threshold_step)
        ratio_step = params->init_param.min_soft_threshold_step;

    bool success = 0;
    uint32_t recover_cap = 0;
    int32_t grad_iter = 0;
    //uint32_t num_samples = params->num_samples; 

    while (!success) {

        if (cur_ratio > max_ratio) {
			// if cur_ratio is already larger than max_ratio, that means we have already tried our best to lower the number of peaks, so we will take the result
			break;
        }
       
        votingProposal(cur_ratio);

        if (num_peaks == 0) {
			// if number of peaks equal to 0 after we set the new ratio, that means the last time ratio is the highest we can have and thus we should take it as final
			// recover ratio and do voting again
			// and then directly take result
			cur_ratio -= ratio_step;
			votingProposal(cur_ratio);
			break; 
        }
          
        recover_cap = 0;
        for (projection_index = 0, cur_projection = projections; projection_index < num_projection_lines; ++ projection_index, ++ cur_projection)
            recover_cap += cur_projection->num_large_buckets; // originally there is another factor (cur_filter->num_shifts - 1) here, but it is always 1

        uint32_t f_v = 0;
        uint32_t f_h = 0;
        int32_t * freq = position_table;

        double pos_threshold = params->init_param.voting_threshold * num_projection_lines;
        
		if (recover_cap >= num_peaks) {
            // Now set peaks
            
            peaks.resize(num_peaks);
            PeakIterator it = peaks.begin();

            for (f_v = 0; f_v < NonIntSFFT::n_v; ++ f_v) {
                for (f_h = 0; f_h < NonIntSFFT::n_h; ++ f_h) {
                    if (*freq >= pos_threshold) { 
                        it->setPosition(NonIntPosition(f_v, f_h));
                        ++ it;
                    }
					freq ++;
				}
			}
 
            // Try to recover the peak
            success = recovery.retrieveConfig(-1, recovery_id, &peaks, num_peaks, &residue);

        } else {
            success = false;
        }
        ++ grad_iter;
        cur_ratio += ratio_step;
	
    }
 
    return true;

}

inline void Initialization::votingProposal(in double cur_ratio) {

    uint32_t projection_index = 0;
    Projection * cur_projection = projections;
    int32_t n = projections->n;
    num_peaks = 0;
    //uint32_t freq_index = 0;
    //double pos_threshold = params->init_param.voting_threshold * num_projection_lines;
#ifdef REAL_INPUT
    pos_threshold *= 2;
#endif
    memset(position_table, 0, sizeof(int32_t) * NonIntSFFT::n);
    int32_t * pos = NULL;
    pos = (int32_t*) malloc(sizeof(int32_t) * 2 * n);
    for (; projection_index < num_projection_lines; ++ projection_index, ++ cur_projection) {

        // first, find the large buckets
        cur_projection->setPowerThreshold(cur_ratio, params->init_param.hard_power_threshold);
        cur_projection->indicateLargeBuckets();
        cur_projection->writeLargeBuckets();
      
    // get the bucket positions that currently map to this bucket and vote them
    uint32_t num_large_buckets = cur_projection->num_large_buckets;
    int32_t bucket_index = -1; 
    int32_t num_pos = 0;
    
    for (uint32_t i = 0; i < num_large_buckets; i ++) {
           bucket_index = cur_projection->large_buckets[i];
           cur_projection->getProj(bucket_index, pos, &num_pos);
		for (int32_t j = 0; j < num_pos; j ++) {
	        position_table[ pos[2*j] * n + pos[2*j+1] ] ++;   // the position of the frequency we want to vote is (pos[2*i], pos[2*i+1])			
	    }
	}
	
    }
    
    // we can count the number of peaks now
    int32_t _num_projection_lines = num_projection_lines;
    for (int32_t i = 0; i < n; i ++) {
		for (int32_t j = 0; j < n; j ++) { 
		  if (position_table[i*n+j] == _num_projection_lines) {
			  num_peaks ++;
		  }		
		}
     }

	// deallocate the space allocated by Projection::getProj
    if(pos) {
		free(pos);
		pos = NULL;
    }   

}





inline bool Initialization::optimismInit(out int32_t * recovery_id, in bool previous) {

    if (params->init_param.num_init_peaks == 0)
        return previous;

    int32_t new_recovery_id = 0;
    Real new_residue = 0;
    bool success = recovery.retrieveConfig(-1, &new_recovery_id, &params->init_param.init_peaks, params->init_param.num_init_peaks, &new_residue);

    if (new_residue < residue) {
        residue = new_residue;
        peaks = params->init_param.init_peaks;
        num_peaks = params->init_param.num_init_peaks;
        * recovery_id = new_residue;
        return success;
    }

    return previous;
}




inline void Initialization::generateSamplingTable() {

    Projection * cur_projection = projections;
    uint32_t projection_index = 0;

    memset(sampling_table, 0, sizeof(bool) * NonIntSFFT::n);

    for (; projection_index < num_projection_lines; ++ cur_projection, ++ projection_index) {
        cur_projection->generateSamplingTable(sampling_table);
    }


}









#endif
