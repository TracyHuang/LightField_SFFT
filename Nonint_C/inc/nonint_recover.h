#ifndef _NON_INT_RECOVER_H_
#define _NON_INT_RECOVER_H_

/**
 * @file    nonint_recover.h
 * @brief   The IntRecover class
 * @author  lixin
 * @date    11/16/2012
 */


/**
 * @modified    The integer recovering class
 * @author      Tianwei
 * @date        07/10/2014
 */

#include "common.h"
#include "projection.h"
#include "inv_lib.h"

#include <time.h>

/**
 * @brief   The integer recovering class. 
 *          This class is heavily optimizied. 
 *
 * This class deals with recovering the values of the 
 * peaks whose positions are non-integer. Simply speaking, 
 * this is a linear-matrix inversion. 
 *
 * This class is heavily optimized for the sparse FFT
 * algorithm, mainly because this is the most time-consuming
 * part in it. 
 *
 * 1) Since most of the time, only one peak in the peaks
 *    list will change each time. Naturally the solution is 
 *    to maintain a data structure and update the data 
 *    structure each time the peak is updated.
 *
 * 2) The code involves matrix inverting. This could be 
 *    optimized because only a constant number of the 
 *    inverted matrix will be changed. TODO: this is a 
 *    future optimization. 
 *
 * 3) But it's not that each time only one peak is updated. 
 *    There are cases when we get fully innovative peaks, but
 *    then we find this is sub-optimal and trace back to the
 *    old scheme. As a result, we would allow to store a 
 *    data structre configuation and have an id to come back 
 *    to this configuation. 
 *
 * 4) A headache problem about matrix operations in C/C++ is
 *    always which library to choose. Here we choose the 
 *    Eigen template library. First, it is free; second, 
 *    benchmark tests show that this is rather fast, 
 *    comparable with the commercial Intel MKL LAPACK library. 
 *
 * @author  lixin
 * @date    11/16/2012
 */
class IntRecover {

    friend class Initialization;


protected:

    /**
     * @brief   Configuation to a integer recovering
     *
     * This is a configuation of a specific recovering. 
     * It includes the coefficient matrix, the vector matrix, 
     * the mapping between current peaks and the buckets
     * under RecoverConfig::coeff and RecoverConfig::b, 
     * and all other information which is necessary. 
     *
     * @author  Tianwei
     * @date    07/10/2014
     */
    class RecoverConfig {
    
    protected:

#if !defined(UPDATING_PSEUDOINVERSE) || defined(PSEUDOINVERSE_CHECKING_MODE)
        /**
         * The coefficient matrix. 
         * The recovery is solving
         * RecoverConfig::coeff \ RecoverConfig::b
         */
        CoeffMat coeff;
#endif
        /**
         * The coefficient matrix that contains all
         * of the buckets as the rows. 
         *
         * Note that the way we arrange the order of the equations, 
         * is first according to which filter it come from, 
         * then which bucket, then which time-domain shift. 
         *
         * This is used to calculate the residue. 
         */
        AllCoeffMat all_coeff;

        /**
         * The vector of the useful buckets
         * Note the recovering is solving
         * RecoverConfig::coeff \ RecoverConfig::b
         */
        BukVec b;
	
		// The list of all projection lines
		Projection * projections;

		// The number of projection lines	
		uint32_t num_projection_lines;

		// The array used to store all the positions that hash to a specific projection_index::bucket_index
        int32_t * hash_pos;

        /**
         * Indicator of whether a bucket is included in this 
         * recovering system. 
         *
         * Each entry is that bucket's index to the equations. 
         * -1 means there is no equation corresponding to that bucket.
         *
         * Further, note that in the way we arrange the order of equations, 
         * same bucket with different time-domain shift will be continous, 
         * so the first entry of this continuous block (We call it A-BLOCK)
         * is recorded in this array.
         */
        int32_t ** bucket_indicator;

        /**
         * This is to record each bucket's index to the rows
         * of all_coeff
         */
        int32_t ** bucket_index_all_coeff;

        /**
         * This is the index of each row of coeff to all_coeff
         */
        int32_t submat_indices[MAX_EQUATIONS];

        /**
         * The number of peaks in this recovering system
         */
        uint32_t num_peaks;

        /**
         * For each of the peak, specifying the index
         * of the bucket this peak hashes into.
         *
         * Note that for each peak, there are num_filters number of 
         * A-BLOCKs, the entries to these A-BLOCKs are stored in a continous 
         * order.
         */
        int32_t peak_buckets[MAX_PEAKS_NUM];

        /**
         * For each of the equation, specifying how many peaks
         * support this equation. 
         *
         * 0 means a large bucket support this entry
         * otherwise it is a positive number of the number of peaks
         */
        int32_t equation_sources[MAX_EQUATIONS];

        /**
         * The total number of equations; i.e.
         * the number of rows in RecoverConfig::coeff
         */
        uint32_t num_equations;

        /**
         * The total number of buckets, i.e.
         * the number of rows in RecoverConfig::all_coeff
         */
        uint32_t num_buckets;

        /**
         * The total bucket coefficient
         */
        AllBukVec bucket_b;

        /**
         * The current reconstructed bucket
         */
        AllBukVec bucket_cur;

        /**
         * The current maximum value in the buckets.
         *
         * We use it to detect singularity
         */
        Real max_bucket;

        /**
         * The maximum scaling; any scaling larger than that is viewed as 
         * singular
         */
        double scaling_threshold;

#if defined(UPDATING_PSEUDOINVERSE)
        InvLib inv_lib;
#endif

    public:

        /**
         * Constructor. 
         * Initialize scalars to be 0 and pointers to be 0. 
         */
        RecoverConfig() : 
            projections(0), num_projection_lines(0), bucket_indicator(0), bucket_index_all_coeff(0), num_peaks(0), num_equations(0), num_buckets(0), max_bucket(0), scaling_threshold(0) {}

        /**
         * Descontructor. 
         * Deallocate every dynamically allocated buffer
         */
        ~RecoverConfig() {
            deallocateBuffers();
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
         * This function reset the gradient algorithm that this recovery
         * is associated with. 
         *
         * This function will allocate space for 
         * RecoverConfig::bucket_indicator, 
         * RecoverConfig::bucket_index_all_coeff, 
         *
         * @param[in]   _projections        the projection list
         * @param[in]   _num_projections    the number of projections
         * @param[in]   _scaling_threshold  the maximum scaling amount that is allowed
         */
        //inline void resetFilters(in const CurFilter * _filters, in const uint32_t _num_filters, in const double _scaling_threshold) forceinline;
         inline void resetProjections(in Projection * _projections, uint32_t _num_projection_lines, double _scaling_threshold) forceinline;	  

        /**
         * The the bucket values
         */
        inline void setBuckets() forceinline;

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
         *          some mistake happing;
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
         *          some mistake happing;
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
         *          some mistake happing;
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
         *          some mistake happing;
         *          otherwise return true
         */
        inline bool newPeak(in out PeaksList * peaks, in const uint32_t new_num_peaks, in PeakIterator new_peak, in uint32_t new_index, out Real * residue) forceinline;

        /**
         * Get the shadow bucket. Each time when we 
         * calculate the residue, we are calculating the 
         * norm of difference between shadow bucket
         * and the real bucket
         *
         * @return  a pointer to the shadow bucket
         */
        /*
        inline ComplexPtr getShadowBucket() forceinline {
            return bucket_cur.data();
        }
        */

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
        inline void writeCoeff() forceinline;

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
         * Replace the same row with another bucket, 
         * by changing RecoverConfig::bucket_indicator, RecoverConfig::submat_indices and RecoverConfig::equation_sources
         *
         * @param[in]   filter_index    the current filter index
         * @param[in]   row_index       the index of the row to be replaced
         * @param[in]   original_bucket the original bucket index
         * @param[in]   new_bucket      the new bucket index
         */
        inline void updateRow(in uint32_t projection_index, in uint32_t row_index, in uint32_t original_bucket, in uint32_t new_bucket) forceinline;

        /**
         * Remove the current row
         * by changing RecoverConfig::bucket_indicator, RecoverConfig::submat_indices and RecoverConfig::equation_sources
         *
         * @param[in]   filter_index    the current filter index
         * @param[in]   row_index       the index of the row to be removed
         * @param[in]   bucket          the bucket index
         */
        inline void removeRow(in uint32_t projection_index, in uint32_t row_index, in uint32_t bucket) forceinline;

        /**
         * Add a new row at the end of the matrix
         * by changing RecoverConfig::bucket_indicator, RecoverConfig::submat_indices and RecoverConfig::equation_sources
         *
         * @param[in]   projection_index    the current projection index
         * @param[in]   bucket          the bucket index
         */
        inline bool addRow(in uint32_t projection_index, in uint32_t bucket) forceinline;

        /**
         * Get the value of a column in current coefficient matrix
         *
         * @param[in]   peak_index      the index of the peak (column)
         * @param[out]  column_value    the buffer to store the column vector
         */
        inline void getColumnValue(IN const uint32_t & peak_index, OUT ComplexPtr column_value) forceinline {


            const Complex * matrix_ptr = all_coeff.data() + num_buckets * peak_index;
            for (uint32_t equation_index = 0; equation_index < num_equations; ++ equation_index) {
                *(column_value ++) = matrix_ptr[submat_indices[equation_index]];
            }

        }

        /**
         * Get the value of a row in current coefficient matrix
         *
         * @param[in]   equation_index  the index of the equation (row)
         * @param[out]  row_value       the buffer to store the row vector
         */
        inline void getRowValue(IN const uint32_t & equation_index, OUT ComplexPtr row_value) forceinline {
            
            const Complex * matrix_ptr = all_coeff.data() + submat_indices[equation_index];
            for (uint32_t peak_index = 0; peak_index < num_peaks; ++ peak_index, matrix_ptr += num_peaks)
                *(row_value ++) = *matrix_ptr;

        }

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
    IntRecover() : start(0), end(-1), next_ID(0){
    }

    /**
     * The deconstructor
     */
    ~IntRecover() {
    }


    /**
     * Set the projections
     */
     inline void setProjections(in Projection * _projections,  in uint32_t _num_projection_lines, in double _scaling_threshold) forceinline {
	
	for (uint32_t i = 0; i < NUM_RECOVER_BUFFER; ++i)
	    recoveries[i].resetProjections(_projections, _num_projection_lines, _scaling_threshold);

	memset(ids, -1, sizeof(int32_t) * NUM_RECOVER_BUFFER);

    }



    /**
     * The the bucket values
     */
    inline void setBuckets() forceinline {

        for (uint32_t i = 0; i < NUM_RECOVER_BUFFER; ++ i)
            recoveries[i].setBuckets();
    }

    /**
     * Get the statistics from inv_lib. 
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


	bool success = recoveries[end].resetPeaks(peaks, num_peaks, residue);

       return success;
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



inline void IntRecover::RecoverConfig::resetProjections(in Projection * _projections, in uint32_t _num_projection_lines, in double _scaling_threshold) {
    
    scaling_threshold = _scaling_threshold;
    
    if(!_projections) {
	fprintf(stderr, "resetProjections: null pointer");
	exit(-1);
    }

    uint32_t projection_index = 0;
    Projection * cur_projection = NULL;

    // the first part is setting parameters
    projections = _projections;
    num_projection_lines = _num_projection_lines;

    num_buckets = 0;
    cur_projection = projections;
    for (projection_index = 0; projection_index < num_projection_lines; ++ projection_index, ++ cur_projection) {
    	num_buckets += cur_projection->n;
    } 

    if (num_buckets > MAX_BUCKET_NUM) {
	fprintf(stderr, "the total number of buckets %d clearly exceeds the allowed matrix size.", num_buckets);

    }

    // we also need to allocate memories as we need
    deallocateBuffers();

    bucket_indicator = (int32_t **)malloc(sizeof(int32_t *) * num_projection_lines);
    bucket_index_all_coeff = (int32_t **)malloc(sizeof(int32_t *) * num_projection_lines);
    cur_projection = projections;
    for (projection_index = 0; projection_index < num_projection_lines; ++ projection_index, ++ cur_projection) {
	bucket_indicator[projection_index] = (int32_t *)malloc(cur_projection->n * sizeof(int32_t));
	bucket_index_all_coeff[projection_index] = (int32_t *)malloc(cur_projection->n * sizeof(int32_t));
    }
    
    hash_pos = (int32_t *) malloc(sizeof(int32_t) * projections->n * 2);          // (x1, y1, x2, y2, x3, y3...)

}




inline void IntRecover::RecoverConfig::setBuckets() {

    bucket_b.resize(num_buckets, 1);
    ComplexPtr cur_bucket = bucket_b.data();

    Projection * cur_projection = projections;
    uint32_t projection_index = 0;
    for (projection_index = 0; projection_index < num_projection_lines; projection_index ++, cur_projection ++) {
	memcpy(cur_bucket, cur_projection->buckets, sizeof(Complex) * cur_projection->n);
        cur_bucket += cur_projection->n;
    }

    max_bucket = 0;
    cur_bucket = bucket_b.data();

    uint32_t bucket_index = 0;
    Real cur_power = 0;
    for (; bucket_index < num_buckets; ++ bucket_index, ++ cur_bucket) {
	cur_power = comAbs(*cur_bucket);
	if (cur_power > max_bucket)
           max_bucket = cur_power;
    }

}


inline bool IntRecover::RecoverConfig::resetPeaks(in out PeaksList * peaks, in const uint32_t _num_peaks, out Real * residue) {

    if (!peaks || !residue) {
        exit(-1);
    }

    num_peaks = _num_peaks;

    Projection * cur_projection = NULL;
    uint32_t projection_index = 0;
    uint32_t cur_bucket = 0;
    //uint32_t shift_index = 0;

    if (num_peaks * num_projection_lines > MAX_PEAKS_NUM) {
        fprintf(stderr, "number of peaks %d is too large\n", num_peaks);
		return false;
    }

    
    // First, let's calculate all_coeff again
    all_coeff.resize(num_buckets, num_peaks);
    bucket_cur.resize(num_buckets, 1);
    ComplexPtr cur_coeff = all_coeff.data();

    
    for (PeakIterator it = peaks->begin(); it != peaks->end(); ++ it) {
		cur_coeff = updateColumn(it, cur_coeff);
    }

    // Then, let's figure out what is a useful bucket: filling out bucket_indicator, peaks_indices and equation_sources
    for (projection_index = 0; projection_index < num_projection_lines; ++ projection_index)
        memset(bucket_indicator[projection_index], -1, sizeof(int32_t) * projections[projection_index].n);
	memset(equation_sources, 0, sizeof(int32_t) * MAX_EQUATIONS);
    
    if (num_peaks == 0) {
        num_equations = 0;
        all_coeff.resize(num_buckets, 0);
        bucket_cur.resize(num_buckets, 0);
#if !defined(UPDATING_PSEUDOINVERSE) || defined(PSEUDOINVERSE_CHECKING_MODE)
        coeff.resize(0, 0);
#endif
        b.resize(0);
       
        return false;
    }

    int32_t * cur_peak_hashing = peak_buckets;
    int32_t block_start = 0;
    int32_t * cur_equation_source = NULL;
    int32_t * cur_submat_index = NULL;
    int32_t index_all_coeff = 0;
    num_equations = 0;
    // go through all peaks
    for (PeakIterator it = peaks->begin(); it != peaks->end(); ++ it) {
        cur_projection = projections;
        for (projection_index = 0; projection_index < num_projection_lines; ++ projection_index, ++ cur_projection) {
            cur_bucket = cur_projection->hashing_table[it->getPos().roundPos().map2Index()];
            
            block_start = bucket_indicator[projection_index][cur_bucket];
            if (block_start < 0) {
                bucket_indicator[projection_index][cur_bucket] = block_start = num_equations;
                num_equations += 1; 
                cur_submat_index = submat_indices + block_start;
                index_all_coeff = bucket_index_all_coeff[projection_index][cur_bucket];
                *cur_submat_index = index_all_coeff; 
				++ cur_submat_index;
            }
            * (cur_peak_hashing ++) = cur_bucket;
            cur_equation_source = equation_sources + block_start;
            (*cur_equation_source) ++;
        }
    }
    
    // go through all large buckets
    cur_projection = projections;
    uint32_t large_bucket_index = 0;
    for (projection_index = 0; projection_index < num_projection_lines; ++ projection_index, ++ cur_projection) {

        for (large_bucket_index = 0; large_bucket_index < cur_projection->num_large_buckets; ++ large_bucket_index) {
            cur_bucket = cur_projection->large_buckets[large_bucket_index];
            block_start = bucket_indicator[projection_index][cur_bucket];
            if (block_start < 0) {
                bucket_indicator[projection_index][cur_bucket] = block_start = num_equations;
                num_equations += 1;

                if (num_equations > MAX_EQUATIONS) {
                    fprintf(stderr, "resetPeaks: the number of equations %d is too large", num_equations);
                    NonIntSFFT::error = 1;
                    return false;
                }

                cur_submat_index = submat_indices + block_start;
                index_all_coeff = bucket_index_all_coeff[projection_index][cur_bucket];
                *cur_submat_index = index_all_coeff;
            } else {
                memset(equation_sources + block_start, 0, sizeof(int32_t));
            }
        }
    }


    // now we are able to construct the coefficient matrix
    writeCoeff();
#ifdef UPDATING_PSEUDOINVERSE
    // update the inv_lib
    //fprintf(stderr, "Can we get here?\n");
    resetInvLib();
    //inv_lib.setMatrix(coeff);
    inv_lib.initPinv();

    checkInvLib();
#endif
    
    return (solveSystem(peaks, residue) == num_peaks);
}

inline void IntRecover::RecoverConfig::replicate2New(out RecoverConfig * target) {

    target->all_coeff.resize(num_buckets, num_peaks);
    target->bucket_cur.resize(num_buckets, 1);
    memcpy(target->all_coeff.data(), all_coeff.data(), sizeof(Complex) * num_buckets * num_peaks);
#if !defined(UPDATING_PSEUDOINVERSE) || defined(PSEUDOINVERSE_CHECKING_MODE)
    target->coeff.resize(num_equations, num_peaks);
    memcpy(target->coeff.data(), coeff.data(), sizeof(Complex) * num_equations * num_peaks);
#endif
    target->b.resize(num_equations, 1);
    memcpy(target->b.data(), b.data(), sizeof(Complex) * num_equations);
    
    memcpy(target->bucket_cur.data(), bucket_cur.data(), sizeof(Complex) * num_buckets);

    // change bucket_index_all and bucket_indicator
    uint32_t projection_index = 0;
    //const CurFilter * cur_filter = filters;
    Projection * cur_projection = projections;

    for (projection_index = 0; projection_index < num_projection_lines; ++ projection_index, ++ cur_projection) {
        memcpy(target->bucket_indicator[projection_index], bucket_indicator[projection_index], sizeof(int32_t) * cur_projection->n);
        memcpy(target->bucket_index_all_coeff[projection_index], bucket_index_all_coeff[projection_index], sizeof(int32_t) * cur_projection->n);
    }

    // change the others
    target->num_peaks = num_peaks;
    target->num_equations = num_equations;
    memcpy(target->submat_indices, submat_indices, sizeof(int32_t) * num_equations);
    memcpy(target->peak_buckets, peak_buckets, sizeof(int32_t) * num_peaks * num_projection_lines);
    memcpy(target->equation_sources, equation_sources, sizeof(int32_t) * num_equations);
#ifdef UPDATING_PSEUDOINVERSE
    target->inv_lib.copyFrom(inv_lib);
    
    checkInvLib();
#endif
    
}

inline void IntRecover::RecoverConfig::checkInvLib() {

#if defined(UPDATING_PSEUDOINVERSE) && defined(PSEUDOINVERSE_CHECKING_MODE)
// actually doing nothing
#endif
}

inline void IntRecover::RecoverConfig::writeCoeff() {
#if defined(PSEUDOINVERSE_CHECKING_MODE) || !defined(UPDATING_PSEUDOINVERSE)
    coeff.resize(num_equations, num_peaks);
#endif
    b.resize(num_equations, 1);

    uint32_t equation_index = 0;
    for (; equation_index < num_equations; ++ equation_index) {
#if defined(PSEUDOINVERSE_CHECKING_MODE) || !defined(UPDATING_PSEUDOINVERSE)
        coeff.row(equation_index) = all_coeff.row(submat_indices[equation_index]);
#endif
        b(equation_index) = bucket_b(submat_indices[equation_index]);
    }
}

inline void IntRecover::RecoverConfig::resetInvLib() {
#ifdef UPDATING_PSEUDOINVERSE
    Complex * matrix = inv_lib.getMatrix();
    inv_lib.resize(num_equations, num_peaks);

    Complex * all_coeff_ptr = NULL;
    for (uint32_t equation_index = 0; equation_index < num_equations; ++ equation_index) {
        all_coeff_ptr = all_coeff.data() + submat_indices[equation_index];
        for (uint32_t peak_index = 0; peak_index < num_peaks; ++ peak_index, ++ matrix, all_coeff_ptr += num_buckets) {
            *matrix = *all_coeff_ptr;
        }
    }
#endif
}

inline uint32_t IntRecover::RecoverConfig::solveSystem(out PeaksList * peaks, out Real * residue) {

    if (!peaks || !residue) {
		fprintf(stderr, "solveSystem: null pointer");
        exit(-1);
    }
#if defined(PSEUDOINVERSE_CHECKING_MODE) || !defined(UPDATING_PSEUDOINVERSE)
    // First solve the system
    Decompose qr = coeff.decompose(decomp_args);
    uint32_t rank = qr.rank();
    BukVec solution = qr.solve(b);
#endif

#if defined(UPDATING_PSEUDOINVERSE)
    // the inv lib
    inv_lib.solveSys(num_equations, b.data(), num_peaks);
#endif

    // Now set the values
    uint32_t peak_index = 0;

#ifdef UPDATING_PSEUDOINVERSE
    const Complex * vs = inv_lib.getSolution();
#else
    const Complex * vs = solution.data();
#endif

#if defined(PSEUDOINVERSE_CHECKING_MODE) || !defined(UPDATING_PSEUDOINVERSE)
    // Check the singularity
    for(; peak_index < num_peaks; ++ vs, ++peak_index) {
        if (comAbs(*vs) > scaling_threshold * max_bucket)
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

    // Get the residue vector
    // bucket_cur.resize(num_buckets, 1);
    Complex * bucket_cur_ptr = bucket_cur.data();
    Complex * all_coeff_ptr = all_coeff.data();
    memset(bucket_cur_ptr, 0, sizeof(Complex) * num_buckets);
    vs = inv_lib.getSolution();

    for (uint32_t i = 0; i < num_peaks; ++ i, ++ vs) {
        bucket_cur_ptr = bucket_cur.data();
        for (uint32_t j = 0; j < num_buckets; ++ j, ++ bucket_cur_ptr, ++ all_coeff_ptr) {
            *bucket_cur_ptr += (*vs) * (*all_coeff_ptr);
        }
    }

    // calculate residue
    *residue = 0;
    const double * dbl_bucket_ptr = (const double *)bucket_cur.data();
    const double * dbl_bucket_b = (const double *)bucket_b.data();
    double diff = 0;
    for (uint32_t i = 0; i < 2 * num_buckets; ++ i, ++ dbl_bucket_ptr, ++ dbl_bucket_b) {
        diff = (*dbl_bucket_b) - (*dbl_bucket_ptr);
        *residue += diff * diff;
    }
    *residue = sqrt(*residue);
#else
    // Finally, calculate the residue
    bucket_cur = all_coeff * solution;
    *residue = (bucket_cur - bucket_b).norm();
#endif

#if defined(PSEUDOINVERSE_CHECKING_MODE) || !defined(UPDATING_PSEUDOINVERSE)
    return rank;
#else
    return num_peaks;
#endif
}

inline bool IntRecover::RecoverConfig::updatePeak(in out PeaksList * peaks, in uint32_t const _num_peaks, in PeakIterator update_peak, in uint32_t update_index, out Real * residue) {

    if (!peaks || !residue) {
        fprintf(stderr, "updatePeak: null pointer");
        exit(-1);
    }

    if (_num_peaks != num_peaks) {
        fprintf(stderr, "updatePeak: number of peaks from %u to %d", num_peaks, _num_peaks);
        exit(-1);
    }

    if (update_index >= num_peaks) {
        fprintf(stderr, "updatePeak: num_peaks = %d, update_index = %d", num_peaks, update_index);
        exit(-1);
    }
    
    //fprintf(stderr, "num_peaks1 = %d, num_peaks2 = %d, update_index = %d\n", peaks->size(), num_peaks, update_index);

    Projection * cur_projection = projections;
    uint32_t projection_index = 0;
    uint32_t hashed_bucket = 0;
    int32_t equ_index = 0;
    int32_t original_bucket = 0;
    int32_t original_index = 0;
    uint32_t index2peak_buckets = 0;
    bool no_change = 1;
    bool pending_add = 0;

    // First, let's update all_coeff
    ComplexPtr data = all_coeff.data();
    // std::cout << "UPDATE ALL" << std::endl;
    // std::cout << all_coeff.col(update_index) << std::endl << std::endl;
    updateColumn(update_peak, data + update_index * num_buckets);
#ifdef UPDATING_PSEUDOINVERSE
    getColumnValue(update_index, inv_lib.getColVec());
    inv_lib.updateCol(update_index, num_equations); //, col_vec_buffer);
#endif
    // Then, let's check if there is any need to update the 
    // bucket that is included in the equations
    for (cur_projection = projections, projection_index = 0; projection_index < num_projection_lines; ++ projection_index, ++ cur_projection) {

        pending_add = 0;

        hashed_bucket = cur_projection->hashing_table[update_peak->getPos().roundPos().map2Index()];
        // hashed_buckets[filter_index] = hashed_bucket;
        equ_index = bucket_indicator[projection_index][hashed_bucket];

        index2peak_buckets = update_index * num_projection_lines + projection_index;
        original_bucket = peak_buckets[index2peak_buckets];
        original_index = bucket_indicator[projection_index][original_bucket];

        peak_buckets[index2peak_buckets] = hashed_bucket;

        // This means that this is a new bucket
        if (equ_index < 0) {
          
            pending_add = 1;

        } else {

            if (equation_sources[equ_index] > 0)
                equation_sources[equ_index] ++;
            fprintf(stderr, "No new equation, equation_sources = %d", equation_sources[equ_index]);
        }

        // If the old one is the only peak in the bucket
        // and that bucket is not large
        if (equation_sources[original_index] == 1) {
            // remove_list[filter_index] = original_index;
            // bucket_indicator[filter_index][original_bucket] = -1;

            if (pending_add) {
                // in this case, we can replace the row
                updateRow(projection_index, original_index, original_bucket, hashed_bucket);

                if (no_change) {
#if !defined(UPDATING_PSEUDOINVERSE) || defined(PSEUDOINVERSE_CHECKING_MODE)
                    coeff.row(original_index) = all_coeff.row(submat_indices[original_index]);
#endif
                    b(original_index) = bucket_b(submat_indices[original_index]);
                }
                pending_add = 0;
            } else {
                removeRow(projection_index, original_index, original_bucket);

                no_change = 0;
            }
        } else {
            if(equation_sources[original_index] > 1)
                equation_sources[original_index] --;
        }

        if (pending_add) {
            if (!addRow(projection_index, hashed_bucket)) {
                NonIntSFFT::error = 1;
                return false;
            }

            no_change = 0;
        }
    }

#if !defined(UPDATING_PSEUDOINVERSE) || defined(PSEUDOINVERSE_CHECKING_MODE)
    // Next step is to conduct any change to the coefficient matrix
    // If there are any change
    uint32_t equation_index = 0;
    ComplexPtr base_coeff = coeff.data() + update_index * num_equations;
    ComplexPtr base_all = all_coeff.data() + update_index * num_buckets;
    if (!no_change) {
        writeCoeff();
    } else {
        // change the column of that peak
        
        for (equation_index = 0; equation_index < num_equations; ++ equation_index) {
            *(base_coeff + equation_index) = *(base_all + submat_indices[equation_index]); 
        }
    }
#else
    if (!no_change)
        writeCoeff();
#endif

#ifdef PSEUDOINVERSE_CHECKING_MODE
    checkInvLib();
#endif

    // Finally, solve the system
    return (solveSystem(peaks, residue) == num_peaks);
}

inline bool IntRecover::RecoverConfig::removePeak(in out PeaksList * peaks, in const uint32_t new_num_peaks, in uint32_t remove_index, out Real * residue) {
    
    if (!peaks || !residue) {
        fprintf(stderr, "removePeak: null pointer");
        exit(-1);
    }

    if (new_num_peaks + 1 != num_peaks) {
        fprintf(stderr, "removePeak: number of peaks from %u to %d", num_peaks, new_num_peaks);
        exit(-1);
	}

    // First, let's update all_coeff
    num_peaks = new_num_peaks;

    ComplexPtr data = all_coeff.data();
    ComplexPtr all_ptr = data + remove_index * num_buckets;
    for (; all_ptr != data + num_peaks * num_buckets; all_ptr += num_buckets)
        memcpy(all_ptr, all_ptr + num_buckets, sizeof(Complex) * num_buckets);

    // if (remove_index != num_peaks)
        // all_coeff.col(remove_index).swap(all_coeff.rightCols<1>());
    all_coeff.resize(Eigen::NoChange, num_peaks);
#ifdef UPDATING_PSEUDOINVERSE
    inv_lib.removeCol(remove_index);
#endif
    // Then, let's check if there is any need to update the 
    // bucket that is included in the equations
    Projection * cur_projection = projections;
    uint32_t projection_index = 0;
    // uint32_t hashed_bucket = 0;
    // int32_t equ_index = 0;
    int32_t original_bucket = 0;
    int32_t original_index = 0;
    uint32_t index2peak_buckets = 0;
    bool no_change = 1;
    // bool pending_add = 0;

    for (projection_index = 0; projection_index < num_projection_lines; ++ projection_index, ++ cur_projection) {
        
        index2peak_buckets = remove_index * num_projection_lines + projection_index;
        original_bucket = peak_buckets[index2peak_buckets];
        original_index = bucket_indicator[projection_index][original_bucket];
        
        // If the old one is the only peak in the bucket
        // and that bucket is not large
        if (equation_sources[original_index] == 1) {

            removeRow(projection_index, original_index, original_bucket);
            no_change = 0;

        } else {
            if(equation_sources[original_index] > 1) 
                equation_sources[original_index] --;
        }
    }

    // also change peak_buckets
    int32_t * cur_peak_bucket = peak_buckets + remove_index * num_projection_lines;
    for (; cur_peak_bucket != peak_buckets + num_projection_lines * num_peaks; cur_peak_bucket += num_projection_lines)
        memcpy(cur_peak_bucket, cur_peak_bucket + num_projection_lines, sizeof(int32_t) * num_projection_lines);

#if !defined(UPDATING_PSEUDOINVERSE) || defined(PSEUDOINVERSE_CHECKING_MODE)
    // Next step is to conduct any change to the coefficient matrix
    // If there are any change
    data = coeff.data();
    ComplexPtr coeff_ptr = data + remove_index * num_equations;
    if (!no_change) {
        writeCoeff();
    } else  {
        //if (remove_index != num_peaks)
            //coeff.col(remove_index).swap(coeff.rightCols<1>());
        
        for (; coeff_ptr != data + num_peaks * num_equations; coeff_ptr += num_equations)
            memcpy(coeff_ptr, coeff_ptr + num_equations, sizeof(Complex) * num_equations);

        coeff.resize(Eigen::NoChange, num_peaks);
    }
#else
    if (!no_change)
        writeCoeff();
#endif
#ifdef PSEUDOINVERSE_CHECKING_MODE
    checkInvLib();
#endif

    if (num_peaks == 0)
        return false;

    // Finally, solve the system
    return (solveSystem(peaks, residue) == num_peaks);

}

inline bool IntRecover::RecoverConfig::newPeak(in out PeaksList * peaks, in const uint32_t new_num_peaks, in PeakIterator new_peak, in uint32_t new_index, out Real * residue) {

    // sanity check
    if (!peaks || !residue) {
        fprintf(stderr, "newPeak: null pointer");
        exit(-1);
    }

    if (new_num_peaks - 1 != num_peaks) {
        fprintf(stderr, "newPeak: number of peaks from %u to %d", num_peaks, new_num_peaks);
        exit(-1);
    }

    if (new_num_peaks * num_projection_lines > MAX_PEAKS_NUM) {
        fprintf(stderr, "newPeak: the size of peak_buckets %d is too large", new_num_peaks * num_projection_lines);
        *residue = INFINITY;
        return false;
    }

    Projection * cur_projection = projections;
    uint32_t projection_index = 0;
    uint32_t hashed_bucket = 0;
    int32_t equ_index = 0;
    uint32_t index2peak_buckets = 0;
    bool no_change = 1;

    // First, let's update all_coeff
    num_peaks = new_num_peaks;
    all_coeff.resize(Eigen::NoChange, num_peaks);
    ComplexPtr data = all_coeff.data();

    // Move the second part forward
    // we cannot use memcpy because of overlapping, and memmove is low-efficient
    ComplexPtr cur_data = data + num_buckets * (num_peaks - 1);
    for (; cur_data != data + num_buckets * new_index; cur_data -= num_buckets) {
        memcpy(cur_data, cur_data - num_buckets, sizeof(Complex) * num_buckets);
    }

    // Update for the new peak
    updateColumn(new_peak, data + new_index * num_buckets);
#ifdef UPDATING_PSEUDOINVERSE
    getColumnValue(new_index, inv_lib.getColVec());
    inv_lib.addCol(new_index, num_equations); //, col_vec_buffer);
#endif
    // Then, let's check if there is any need to update the 
    // bucket that is included in the equations

    // Update peak_buckets to give enough space for this peak
    int32_t * cur_peak_bucket = peak_buckets + num_projection_lines * (num_peaks - 1);
    for(; cur_peak_bucket != peak_buckets + num_projection_lines * new_index; cur_peak_bucket -= num_projection_lines)
        memcpy(cur_peak_bucket, cur_peak_bucket - num_projection_lines, sizeof(int32_t) * num_projection_lines);
    

    for (cur_projection = projections, projection_index = 0; projection_index < num_projection_lines; ++ projection_index, ++ cur_projection) {
        hashed_bucket = cur_projection->hashing_table[new_peak->getPos().roundPos().map2Index()];
        equ_index = bucket_indicator[projection_index][hashed_bucket];
        
        index2peak_buckets = new_index * num_projection_lines + projection_index;
        peak_buckets[index2peak_buckets] = hashed_bucket;

        // This means that this is a new bucket
        if (equ_index < 0) {
          
            if (!addRow(projection_index, hashed_bucket)) {
                NonIntSFFT::error = 1;
                return false;
            }

            no_change = 0;

        } else {

            if (equation_sources[equ_index] > 0)
                equation_sources[equ_index] ++;
        }

    }
#if !defined(UPDATING_PSEUDOINVERSE) || defined(PSEUDOINVERSE_CHECKING_MODE)
    // Next step is to conduct any change to the coefficient matrix
    // If there are any change
    uint32_t equation_index = 0;
    ComplexPtr base_coeff = NULL;
    ComplexPtr base_all = all_coeff.data() + new_index * num_buckets;

    if (!no_change) {
        writeCoeff();
    } else {
        // add that new column
        coeff.resize(Eigen::NoChange, num_peaks);
        base_coeff = coeff.data() + new_index * num_equations;

        // Move the second part forward
        // we cannot use memcpy because of overlapping, and memmove is low-efficient
        cur_data = base_coeff + num_equations * (num_peaks - 1);
        for (; cur_data != base_coeff + num_equations * new_index; cur_data -= num_buckets) {
            memcpy(cur_data, cur_data - num_equations, sizeof(Complex) * num_equations);
        }

        // copy that new column from all_coeffs
        for (equation_index = 0; equation_index < num_equations; ++ equation_index) {
            *(base_coeff + equation_index) = *(base_all + submat_indices[equation_index]); 
        }

    }
#else
    if (!no_change)
        writeCoeff();
#endif

#ifdef PSEUDOINVERSE_CHECKING_MODE
    checkInvLib();
#endif
   
    // Finally, solve the system
    return (solveSystem(peaks, residue) == num_peaks);
}

inline ComplexPtr IntRecover::RecoverConfig::updateColumn(in PeakIterator update_peak, in ComplexPtr cur_coeff) {

    // We cache a sinc template to speed up the calculation
    Utils::generateSincTemplate(update_peak->getX(), update_peak->getY()); 
    Projection * cur_projection = projections;
    uint32_t projection_index = 0;
    //uint32_t shift_index = 0;
    int32_t bucket_index = 0;
    uint32_t num_hashes = 0; // represent how many frequencies ar going to be hasehed into the same bucket
    uint32_t hash_index = 0;
    IntPosition * cur_pos = NULL; // used to store the current bucket position
    int32_t cur_bucket_index = 0;
    
    ComplexPtr writer = cur_coeff;

       
    int32_t _n = 0;
    int32_t * cur_hash_pos = hash_pos; 

    for (projection_index = 0; projection_index < num_projection_lines; ++ projection_index, ++ cur_projection) {
        num_hashes = NonIntSFFT::n / cur_projection->n; // NonIntSFFT::n = cur_projection->n * cur_projection->n
        cur_pos = cur_projection->bucket_table;
        
        for (bucket_index = 0; bucket_index < cur_projection->n; bucket_index ++) {
            *writer = 0;
	    cur_projection->getProj(bucket_index, hash_pos, &_n);         
	    if (_n != cur_projection->n) printf("ERROR: Did not get enough number of positions\n");
            cur_hash_pos = hash_pos; 
            for (hash_index = 0; hash_index < num_hashes; ++ hash_index) {
		
                *writer += Utils::readOffSincTail(*(cur_hash_pos), *(cur_hash_pos + 1)) * cur_projection->coeffs[IntPosition(*(cur_hash_pos), *(cur_hash_pos+1)).map2Index()];
                cur_hash_pos += 2;
            }
            writer ++;
            cur_pos += 1; 
            bucket_index_all_coeff[projection_index][bucket_index] = cur_bucket_index;
            cur_bucket_index += 1; 
        }
    }

    
    return writer;
}

inline void IntRecover::RecoverConfig::updateRow(in uint32_t projection_index, in uint32_t row_index, in uint32_t original_bucket, in uint32_t new_bucket) {

    bucket_indicator[projection_index][original_bucket] = -1;
    bucket_indicator[projection_index][new_bucket] = row_index;

    uint32_t size_block = 1; 
    uint32_t base_index = bucket_index_all_coeff[projection_index][new_bucket];
    uint32_t equation_index = row_index;
    for (; equation_index < row_index + size_block; equation_index ++) {
        submat_indices[equation_index] = base_index ++;
        equation_sources[equation_index] = 1;
#ifdef UPDATING_PSEUDOINVERSE
        getRowValue(equation_index, inv_lib.getRowVec());
        inv_lib.updateRow(equation_index, num_peaks); //, row_vec_buffer);
#endif
    } 

}

inline void IntRecover::RecoverConfig::removeRow(in uint32_t projection_index, in uint32_t row_index, in uint32_t bucket) {

    
    uint32_t size_block = 1; //projections[projection_index].num_shifts - 1;

    // First update all data structures
    bucket_indicator[projection_index][bucket] = -1;
    
    uint32_t cur_projection_index = 0;
    int32_t cur_bucket = 0;
    Projection * cur_projection = projections;
    int32_t * cur_indicator = NULL;
    int32_t row_threshold = row_index + size_block;
    for (cur_projection_index = 0; cur_projection_index < num_projection_lines; ++ cur_projection_index, ++ cur_projection) {
        cur_indicator = bucket_indicator[cur_projection_index];
        for (cur_bucket = 0; cur_bucket < cur_projection->n; ++ cur_bucket, ++ cur_indicator) {
            if (*cur_indicator >= row_threshold) {
                *cur_indicator -= size_block;
            }
        }
    }


    num_equations -= size_block;

    // Just a notice, here we cannot use memcpy, and I don't like memmove
    // so I'd rather move one by one by myself

    int32_t * submat_ptr = submat_indices + row_index;
    int32_t * equ_ptr = equation_sources + row_index;

    for (; submat_ptr != submat_indices + num_equations; ++ submat_ptr)
        *(submat_ptr) = *(submat_ptr + size_block);

    for (; equ_ptr != equation_sources + num_equations; ++ equ_ptr)
        *(equ_ptr) = *(equ_ptr + size_block);

    // memcpy(submat_indices + row_index, submat_indices + row_index + size_block, sizeof(uint32_t) * (num_equations - row_index));
    // memcpy(equation_sources + row_index, equation_sources + row_index + size_block, sizeof(uint32_t) * (num_equations - row_index));
#ifdef UPDATING_PSEUDOINVERSE
    for (uint32_t equation_index = row_index; equation_index < row_index + size_block; ++ equation_index)
        inv_lib.removeRow(equation_index);
#endif
}

inline bool IntRecover::RecoverConfig::addRow(in uint32_t projection_index, in uint32_t bucket) {
    
    bucket_indicator[projection_index][bucket] = num_equations;
    uint32_t size_block = 1; //filters[filter_index].num_shifts - 1;

    uint32_t equation_index = num_equations;
    uint32_t new_num_equations = num_equations + size_block;

    if (new_num_equations > MAX_EQUATIONS) {
        fprintf(stderr, "addRow: the number of equations %d is too large", new_num_equations);
        return false;
    }

    uint32_t base_index = bucket_index_all_coeff[projection_index][bucket];
    for (; equation_index < new_num_equations; equation_index ++) {
        submat_indices[equation_index] = base_index ++;
        equation_sources[equation_index] = 1;
#ifdef UPDATING_PSEUDOINVERSE
        getRowValue(equation_index, inv_lib.getRowVec());
        inv_lib.addRow(equation_index, num_peaks); //, row_vec_buffer);
#endif
    }
    
    num_equations = new_num_equations;

    
    return true;
}


inline void IntRecover::RecoverConfig::deallocateBuffers() {
   
    uint32_t projection_index = 0;
    if (bucket_indicator) {
        for (projection_index = 0; projection_index < num_projection_lines; ++ projection_index)
            free(bucket_indicator[projection_index]);
        free(bucket_indicator);
        bucket_indicator = NULL;
    }

    projection_index = 0;
    if(bucket_index_all_coeff) {
	for(; projection_index < num_projection_lines; ++projection_index) 
	    free(bucket_index_all_coeff[projection_index]);
	free(bucket_index_all_coeff);
	bucket_index_all_coeff = NULL;
    }


    if(hash_pos) {
	free(hash_pos);
	hash_pos = NULL;
    }

}


#endif
