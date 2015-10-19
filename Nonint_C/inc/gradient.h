#ifndef _GRADIENT_H_
#define _GRADIENT_H_

/**
 * @file    gradient.h
 * @brief   Implement gradient search algorithm; i.e. the second part of the gradient algorithm 
 * @author  lixin Tianwei
 * @date    07/06/2014
 */

#include "common.h"
#include "projection.h"
#include "params.h"
#include "initialization.h"
#include <list>
#include "nonint_recover_general.h"
#include <unistd.h>
#include <stdio.h>
#include <cmath>
/**
 * @brief   The gradient algorithm to find the in finding the fractional 
 *          positions of peaks
 *
 * This class implements the gradient algorithm that iteratively find
 * the optimal non-integer peak position. 
 *
 * @author  lixin Tianwei
 * @date    07/06/2014
 */


/**
 * used to do gradient descent to find solution
 */
class Gradient {

    /**
     * NonIntRecover is a friend class of this class
     */
    friend class NonIntRecover;
	friend class Initialization;

protected:

    /**
     * The internal state of the algorithm
     */
    enum STATE {
        INIT,                       /**< initial state */
        GRADIENT_STUCK,             /**< after (joint) gradient search, stucked */
        GRADIENT_UPDATE,            /**< after (joint) gradient search, updated */
        EXHAUSTIVE_STUCK,           /**< after exhaustive search, stucked */
        EXHAUSTIVE_UPDATE,          /**< after exhaustive search, updated */
        SUBTRACT_SHIFT_STUCK,       /**< after subtract_and_shift, stucked */
        SUBTRACT_SHIFT_UPDATE,      /**< after subtract_and_shift, updated */
        CHECK_UPDATED,              /**< after check_esimate, updated */
        FINISHED                    /**< the algorithm finished */
    };


    /**
     * This records all the parameters that tunes this
     * algorithm
     */
    Params * params;

    /**
     * The time domain signal
     */
    //ComplexPtr x_in;

    /**
     * The number of time domain samples we have.
     *
     * Note that in current implementation, we will always keep all time domain samples in the equations
     */
    uint32_t num_samples;

    /**
     * The list contains the sample positions
     */
    IntPosition * sample_pos;

    /**
     * The value of time domain samples
     */
    ComplexPtr sample_values;

    /**
     * The peaks that are constantly updated by the gradient algorithm. 
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
     * Interface to NonIntRecover class that we are using
     */
    NonIntRecover recovery;

    /**
     * The residue of the current system
     */
    Real residue;

    /**
     * The internal state of the algorithm
     */
    STATE state;

    /**
     * The current iteration number
     */
    uint32_t iter;

    /**
     * The current gradient iteration number
     */
    uint32_t grad_iter;

    /**
     * Peak indicator: 
     * An indicator of peaks (0-1). Used in different algorithms as buffering.
     * We put it as a field and with statical allocation
     */
    bool peak_indicator[MAX_PEAKS_NUM];

    /**
     * The residue window, recording the current history of residue
     */
    Real * residue_window;

    /**
     * The current writer to the residue window
     */
    uint32_t residue_writer;

    /**
     * The current optimal residue
     */
    Real optimal_residue;

    /**
     * The current optimal peak
     */
    PeaksList optimal_peaks;

    /**
     * The recovery id for the current optimal peak
     */
    int32_t optimal_recovery_id;

    /**
     * The finally recovered frequency domain
     */
    ComplexPtr optimal_y;

    /**
     * The maximum value in current residue window
     */
    Real maxima_residue;

    /**
     * The minimum value in current residue window
     */
    Real minima_residue;

    /**
     * Whether we have already traced back or not
     */
    bool traced_back;

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
     * Note that compared with Gradient::sampling_table it now has a totally different meaning
     */
    bool *sampling_table;

    /**
     * This stores the initial position of peaks
     */
    PeaksList init_peaks;

    /**
     * This stores the initial residue
     */
    Real init_residue;

	/**
     * Pointer to the first step, as we need it for shadow bucket detection
     */
	Initialization * initialization;

	/**
     * The buffer used to store residue when passing back to step one for shadow bucket detection
     */
	Real ** residue_graph;

	/**
     * The graph used to show that which positions are sampled
     */
	int32_t ** sample_graph;


public:

    /**
     * Constructor. Specifies the construction by 
     * the parameter params. 
     *
     *  THE CORRECT ORDER OF FUNCTION CALL: constructor -> setGradient -> getSampleValues 
     *
     * @param[in]  _params   the parameters for the gradient algorithm
     */
    Gradient(Params * _params, Initialization * _initialization) : params(_params), num_samples(0), sample_pos(0), sample_values(0), peaks(0), num_peaks(0), residue(0), state(INIT), iter(0), 
        grad_iter(0), residue_window(0), residue_writer(0), optimal_residue(INFINITY), optimal_peaks(PeaksList()), optimal_recovery_id(0), optimal_y(0), 
        maxima_residue(0), minima_residue(INFINITY), traced_back(0), fftw_input(0), init_plan(0),  sampling_table(0), init_peaks(PeaksList()), init_residue(0), initialization(_initialization), residue_graph(0), sample_graph(0){

        setParameters(_params);
          
    }

    /**
     * Destructor: free the buffers 
     */
    ~Gradient() {
		initialization = NULL;		
        deallocateBuffers();
		destroyPlan();
    }

    /**
     * This function sets all the parameters used in the algorithm. 
     * This mainly include allocating the Gradient::filters buffer, 
     * Gradient::position_table, Gradient::residue_window, 
     * Gradient::optimal_y, Gradient::sampling_table. 
     *
     * @param[in]  _params   the parameters for the gradient algorithm
     */
    void setParameters(Params * _params);


    /**
     * Used to set up values of num_peaks, peaks
     * This function should be called right after a new gradient is generated 
     * and only after this function is called, we can call setRecovery
     */
    void setGradient(uint32_t _num_peaks, PeaksList _peaks, Real _residue, int32_t * recovery_id);  

    /** 
     * Used to setup NonIntRecover recovery
     * Note that this function can only be called after setGradient
     * and only after this fucntion is called , we can start do gradientSearch
     */
    void setRecovery(Params * params, ComplexPtr x_in);

    /**
     * Set sample values
     * note that this function is based on values of x_in and sample_pos
     */
    void getSampleValues(ComplexPtr x_in);

    /**
     * Run the whole gradient algorithm. 
     */
    void gradientSearch(ComplexPtr x_in, int32_t _recovery_id);

    /**
     * Get the results
     *
     * @param[out]  _optimal_y          the optimal frequency spectrum
     * @param[out]  _optimal_peaks      the optimal peak lists
     * @param[out]  _optimal_residue    the optimal residue
     */
    void getResults(out ComplexPtr * _optimal_y, out PeaksList ** _optimal_peaks, out Real ** _optimal_residue) {

        *_optimal_y = optimal_y;
        *_optimal_peaks = &optimal_peaks;
        *_optimal_residue = &optimal_residue;
    }

    /**
     * Get the initial configuration
     *
     * @param[out]  _init_peaks     the initial peak configuration
     * @param[out]  _init_residue   the initial residue
     */
    void getInitConfig(out PeaksList ** _init_peaks, out Real ** _init_residue) {

        *_init_peaks = &init_peaks;
        *_init_residue = &init_residue;
    }

    /**
     * Get the statistics from inv_lib in recovery
     * This function is for debugging
     */
    inline const InvOpStatistics & getOpStatistics() const forceinline { 
        return recovery.getOpStatistics();
    }

    /**
     * Reset the statistics from inv_lib in recovery
     * This function is for debugging
     */
    inline void resetOpStatistics() forceinline {
        recovery.resetOpStatistics();
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
     * Deallocate buffers:
     * Gradient::position_table,
     * Gradient::residue_window,
     * Gradient::optimal_y, 
     * Gradient::fftw_input, 
     * Gradient::sampling_table.
     */
    void deallocateBuffers();

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
     * This function calculates the gradient decent centered at a given peak
     *
     * @param[in]       cur_peak    the current peak for gradient decent
     * @param[in]       peak_index  the index of that peak
     * @param[in]       cur_residue the current residue for the current position 
     * @param[out]      gradient    the 3x3 gradient matrix
     * @param[out]      residues    the 3x3 residue matrix
     * @param[in,out]   recovery_id the id to the recovery configurations
     * @param[out]      direction_v the optimal vertical direction
     * @param[out]      direction_h the optimal horizontal direction
     *
     * @return  whether the position that achieves smallest gradient is underdetermined
     */
    inline bool calGradient(in PeakIterator cur_peak, in uint32_t peak_index, in Real cur_residue, out Real gradient[][3], 
        out Real residues[][3], in out int32_t * recovery_id, out int32_t * direction_v, out int32_t * direction_h) forceinline;

    /**
     * The gradient algorithm. 
     * It selects between singleGradientSearch() and jointGradientSearch()
     *
     * @param[out]  well_determined whether the system is well-determined or not
     * @param[in,out]  recovery_id     the id to the recovery configurations
     */

    inline void gradientOption(out bool * well_determined, in out int32_t * recovery_id) forceinline {
        
        switch (params->method) {
        case Params::GRADIENT:
            singleGradientSearch(well_determined, recovery_id); 
            break;
        case Params::EXHAUSTIVE:
            exhaustiveSearch(well_determined, recovery_id);
            break;
        default:
			break;
        }
    }

    /**
     * The single gradient search algorithm
     *
     * @param[out]  well_determined whether the system is well-determined or not
     * @param[in,out]  recovery_id     the id to the recovery configurations
     */
    void singleGradientSearch(out bool * well_determined, in out int32_t * recovery_id);

    /**
     * The exhaustive search algorithm. 
     *
     * @param[out]  well_determined whether the system is well-determined or not
     * @param[in,out]  recovery_id     the id to the recovery configurations
     */
    void exhaustiveSearch(out bool * well_determined, in out int32_t * recovery_id);

    /**
     * Merge the pair with smallest distance, if their distance is smaller than 
     * Params::min_peak_distance
     *
     * @param[out]  well_determined whether the system is well-determined or not
     * @param[in,out]  recovery_id     the id to the recovery configurations
     *
     * @return
     *          whether actually merge happens or not
     */
    bool fastMerge(out bool * well_determined, in out int32_t * recovery_id);

    /**
     * Merge peaks that is smaller than Params::min_peak_distance, 
     * and update if the new residue is smaller than old one times Params::amend_peak_tolerance
     *
     * @param[out]  well_determined whether the system is well-determined or not
     * @param[in,out]  recovery_id     the id to the recovery configurations
     *
     * @return
     *          whether actually merge happens or not
     */
    bool merge(out bool * well_determined, in out int32_t * recovery_id);

	/**
	 * Find the nearest peak value
     *
	 * @param[in] cur_peak  iterator that points to the current peak
	 * @param[in] _peaks    pointer to the current peaks list
	 *	
     * @return       return the value of the nearest peak
     */

	Complex nearestPeakValue(PeakIterator cur_peak, PeaksList * _peaks);
	


	/**
     * Check the current estmated shadow bucket against the ground truth bucket, and try to fix the current estimation
     *  	
     * @param[out] well_determined whether the system is well-determined or not 
     * @param[in, out]  recovery_id      the id to the recovery configurations
	 *
	 * @return
	 *            whether there is any change in the peaks list
	 *
	 */
	bool checkEst( in out int32_t * recovery_id);

    /**
     * Add the missing peak found in checkEst()
     *
     * @param[in]   missing_peak_pos    the 1d position of the missing peak
     * @param[out]  well_determined     whether the system is well-determined or not
     * @param[in,out]  recovery_id      the id to the recovery configurations
     *
     * @return
     *      whether this peak is added or not
     */
    inline bool addMissingPeak(in uint32_t missing_peak_pos, out bool * well_determined, in out int32_t * recovery_id) forceinline;

    /**
     * Remove the fake peak found in checkEst()
     *
     * @param[in]   fake_peak_pos       the 1d position of the fake peak
     * @param[out]  well_determined     whether the system is well-determined or not
     * @param[in,out]  recovery_id      the id to the recovery configurations
     *
     * @return
     *      whether this peak is removed or not
     */
    inline bool removeFakePeak(in uint32_t fake_peak_pos, out bool * well_determined, in out int32_t * recovery_id) forceinline;

    /**
     * Merge two peaks together, and update the recovery configurations.
     *
     * @param[in,out]   it1             the peak to be erased
     * @param[in,out]   it2             the peak to be updated (surviving peak)
     * @param[in,out]   new_peaks       the list of peaks
     * @param[in,out]   new_num_peaks   the number of peaks. Note that this number is modfied
     *                                  (decreased by 1) in this function
     * @param[in]       index1          the index to the first peak (it1)
     * @param[in]       index2          the index to the second peak (it2)
     * @param[out]      new_residue     the new residue
     *
     * @return  whether the system is well-determined or not
     */

    inline bool mergeTogether(in out PeakIterator & it1, in out PeakIterator & it2, in out PeaksList * new_peaks, in out uint32_t & new_num_peaks, in uint32_t index1, in uint32_t index2, out Real * new_residue) forceinline {

        it2->setPosition(0.5 * (it1->getX() + it2->getX()), 0.5 * (it1->getY() + it2->getY()));

        recovery.updateConfig(new_peaks, new_num_peaks, it2, index2, new_residue);

        if (NonIntSFFT::error)
            return false;

        new_num_peaks --;

        it1 = new_peaks->erase(it1);
        return recovery.shrinkConfig(new_peaks, new_num_peaks, index1, new_residue);
    }

    /**
     * Update the residue window, and decide whether 
     * the iteration goes into convergence. 
     *
     * @param[out]  periodical      whether any periodical structure is found in the residue window
     * @param[out]  small_change    whether the change of the residue within the window is very 
     *                              small, i.e., smaller than the Params::noticable_change level
     * @param[out]  suboptimal      whether all of the residues in the window is smaller than 
     *                              the current optimal Graident::optimal_residue
     * @param[in]   recovery_id     the current recovery id
     */

    inline void updateResidueWindow(out bool * periodical, out bool * small_change, out bool * suboptimal, in int32_t recovery_id) forceinline {

        bool maxima_update = 0;
        bool minima_update = 0;

        if (residue < optimal_residue) {
            optimal_residue = residue;
            optimal_peaks = peaks;
            optimal_recovery_id = recovery_id;
            traced_back = 0;

        }

        * periodical = 0;
        if (Utils::equal(residue, maxima_residue, params->tolerance)) {
            //printf("aaa\n");
			*periodical = 1; 
            maxima_update = 1;
        } else if (residue > maxima_residue) {
            maxima_residue = residue;
            maxima_update = 1;
        }

        if (Utils::equal(residue, minima_residue, params->tolerance)) {
            //printf("bbb\n");
            * periodical = 1;
            minima_update = 1;
        } else if (residue < minima_residue) {
                minima_residue = residue;
                minima_update = 1;
        }

        Real residue_to_be_removed = residue_window[residue_writer];
        uint32_t residue_index = 0;
        bool window_filled_up = !isnan(residue_to_be_removed);
        if (window_filled_up) {
            if (Utils::equal(residue_to_be_removed, maxima_residue, params->tolerance) && !maxima_update) {
                // need to find new maxima
                maxima_residue = residue;
                for (residue_index = 0; residue_index < params->residue_window_size; ++ residue_index) {
                     if (residue_index != residue_writer && !isnan(residue_window[residue_index]) && residue_window[residue_index] > maxima_residue)
                        maxima_residue = residue_window[residue_index];
                }
            }
            if (Utils::equal(residue_to_be_removed, minima_residue, params->tolerance) && !minima_update) {
                // need to find new minima
                minima_residue = residue;
                for (residue_index = 0; residue_index < params->residue_window_size; ++ residue_index) {
                     if (residue_index != residue_writer && !isnan(residue_window[residue_index]) && residue_window[residue_index] < minima_residue)
                        minima_residue = residue_window[residue_index];
                }
            }
        }

        * suboptimal = (minima_residue > optimal_residue);
        * small_change = (maxima_residue - minima_residue) < params->noticable_change * maxima_residue;

        // if the window is not filled up, we didn't do any check
        // otherwise the check is likely to be frivolous
        * periodical &= window_filled_up;
        * suboptimal &= window_filled_up;
        * small_change &= window_filled_up;

        // update residue window
        residue_window[residue_writer] = residue;
        residue_writer = (residue_writer + 1) % params->residue_window_size;

    }


    /**
     * Try to traceback to the current optimal profile
     *
     * @param[out]  recovery_id     after the trace-back operation, the current recovery id
     *
     * @return  if this is the second time tracing back to the same configuration, we would
     *          say that tracing-back fails and the algorithm is forced to end
     */

    inline bool traceback(out int32_t * recovery_id) forceinline {

        if (traced_back) {
            fprintf(stderr, "SECOND TRACEBACK: Algorithm terminates.");
            return false;
        }

        residue = optimal_residue;
        peaks = optimal_peaks;
        num_peaks = optimal_peaks.size();

        recovery.retrieveConfig(optimal_recovery_id, recovery_id, &peaks, num_peaks, &residue);

        traced_back = 1;

        // update residue window
        minima_residue = residue;
        Real residue_to_be_removed = residue_window[residue_writer];
        uint32_t residue_index = 0;
        if (!isnan(residue_to_be_removed) && Utils::equal(residue_to_be_removed, maxima_residue, params->tolerance)) {
            // need to find new maxima
            maxima_residue = residue;
            for (residue_index = 0; residue_index < params->residue_window_size; ++ residue_index) {
                 if (residue_index != residue_writer && !isnan(residue_window[residue_index]) && residue_window[residue_index] > maxima_residue)
                    maxima_residue = residue_window[residue_index];
            }
        }

        residue_window[residue_writer] = residue;
        residue_writer = (residue_writer + 1) % params->residue_window_size;

        return true;
    }


    /**
     * Try to remove peaks from the current configuration
     *
     * @param[in,out]   recovery_id the id of the current recovery
     *
     * @return  true if one or more peak is actually removed
     */
    bool try2Remove(in out int32_t * recovery_id);


    /**
     * Given optimal_peaks, let's try to recover the frequency domain. 
     * It will write the result into Gradient::optimal_y. 
     *
     * This function also fixes the time domain. 
     *
     * @param[in]   x_in    the time domain signal
     */
    void recoverFreqDomain(in ComplexPtr x_in);

    /**
     * Generate the time domain sampling table 
     * Gradient::sampling_table
     */
   inline void generateSamplingTable() forceinline;

};



inline bool Gradient::optimismInit(out int32_t * recovery_id, in bool previous) {

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


inline bool Gradient::calGradient(in PeakIterator cur_peak, in uint32_t peak_index, in Real cur_residue, 
    out Real gradient[][3], out Real residues[][3], in out int32_t * recovery_id, out int32_t * direction_v, out int32_t * direction_h) {

    double delta = params->delta;
    double shift_v = 0;
    double shift_h = 0;
    Real min_gradient = 0;
    * direction_v = 0;
    * direction_h = 0;
    NonIntPosition pos = cur_peak->getPos();
    Real new_residue = cur_residue;
    Real new_gradient = 0;

    recovery.retrieveConfig(*recovery_id, recovery_id, &peaks, num_peaks, &new_residue);
	
    bool success = true;
    bool cur_success = true;

    for (int32_t i_v = -1; i_v <=1; ++ i_v)
        for (int32_t i_h = -1; i_h <= 1; ++ i_h) {
            shift_v = delta * (double)i_v;
            shift_h = delta * (double)i_h;
            if (i_v != 0 || i_h != 0) {

                if (!isnan(gradient[i_v + 1][i_h + 1])) { 
                    new_residue = gradient[i_v + 1][i_h + 1];
                } else {

                    cur_peak->setPosition(pos.getX() + shift_v, pos.getY() + shift_h);
                    cur_success = recovery.updateConfig(&peaks, num_peaks, cur_peak, peak_index, &new_residue);

                    if (NonIntSFFT::error) return false;
                }

                residues[i_v + 1][i_h + 1] = new_residue;
                new_gradient = (new_residue - cur_residue) / ((i_v == 0 || i_h == 0) ? 1 : sqrt(2)) / delta;
                gradient[i_v + 1][i_h + 1] = new_gradient;

                if (new_gradient < min_gradient) {
                    min_gradient = new_gradient;
                    * direction_v = i_v;
                    * direction_h = i_h;
                    success = cur_success;
                }

            } else {
                residues[i_v + 1][i_h + 1] = cur_residue;
                gradient[i_v + 1][i_h + 1] = 0;
            }
        }
    
	cur_peak->setPosition(pos);
    return success;
}


inline void Gradient::generateSamplingTable() {

    for (uint32_t i = 0; i < num_samples; i ++) {
	sampling_table[sample_pos[i].map2Index()] = 1;	
    }
 }



inline bool Gradient::addMissingPeak(in uint32_t missing_peak_pos, out bool * well_determined, in out int32_t * recovery_id) {

    Real new_residue = 0;
    recovery.retrieveConfig(*recovery_id, recovery_id, &peaks, num_peaks, &new_residue);
    
    double distance = 0;

    IntPosition missing_pos = IntPosition::index2Position(missing_peak_pos);

    NonIntPeak missing_peak(missing_pos.getX(), missing_pos.getY(), 0);

    for (PeakIterator it = peaks.begin(); it != peaks.end(); ++ it) {
        distance = it->gridDistance(missing_peak);
        if (distance < params->checkest_min_distance) {
            return false;
        }
    }

    // add this guy!
    bool new_well_determined = 0;

    peaks.push_back(missing_peak);
    num_peaks ++;
    // uint32_t new_num_peaks = num_peaks + 1;
    new_well_determined = recovery.expandConfig(&peaks, num_peaks, --peaks.end(), num_peaks - 1, &new_residue);

    if (NonIntSFFT::error)
        return false;


    if (new_residue < residue * (1 - params->noticable_change)) {
        
        // accept this change
        residue = new_residue;
        *well_determined = new_well_determined;

        return true;
    }

    // revert back
    peaks.pop_back();
    num_peaks --;
    if (!isinf(new_residue)) {
        *well_determined = recovery.shrinkConfig(&peaks, num_peaks, num_peaks, &new_residue);
    }

    return false;
}


inline bool Gradient::removeFakePeak(in uint32_t fake_peak_pos, out bool * well_determined, in out int32_t * recovery_id) {

    Real new_residue = 0; 
    recovery.retrieveConfig(*recovery_id, recovery_id, &peaks, num_peaks, &new_residue);
    
    NonIntPosition fake_pos = NonIntPosition::index2Position(fake_peak_pos);

    PeakIterator fake_peak = peaks.end();
    int32_t fake_index = -1;
    double distance = 0;
    double min_distance = params->checkest_min_distance;
    uint32_t peak_index = 0;

    for (PeakIterator it = peaks.begin(); it != peaks.end(); ++ it, ++ peak_index) {
        distance = it->getPos().gridDistance(fake_pos);
        if (distance < min_distance) {
            fake_peak = it;
            min_distance = distance;
            fake_index = peak_index;
        }
    }

    // if we don't find peaks to remove
    if (fake_index < 0) {
        return false;
    }

    // try to remove this guy!
    bool new_well_determined = 0;
    NonIntPosition removed_pos(fake_peak->getX(), fake_peak->getY());
    peaks.erase(fake_peak);
    num_peaks --;
    new_well_determined = recovery.shrinkConfig(&peaks, num_peaks, fake_index, &new_residue);

    if (new_residue < residue * (1 - params->noticable_change)) {
        
        // accept this change
        residue = new_residue;
        *well_determined = new_well_determined;

        return true;
    }

    // revert back
    peaks.push_back(NonIntPeak(removed_pos));
    num_peaks ++;
    *well_determined = recovery.expandConfig(&peaks, num_peaks, --peaks.end(), num_peaks - 1, &new_residue);
    
	return false;

}

#endif








