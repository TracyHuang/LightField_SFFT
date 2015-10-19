#ifndef _SFFT_H_
//#define _SFFT_H_


/**
 * @file    sfft.h  
 * @brief   Implementing SFFT class
 * @author  lixin
 * @date    11/20/2012
 */

#include "gradient.h"
#include "initialization.h"
#include "peaks.h"
#include "data_io.h"
#include <time.h>

/**
 * @brief   the sparse FFT class that actually performs sparse FFT
 *
 * This is the class that reads data from some DataSource object
 * perform sparse FFT, and output the result to some DataSink object.
 *
 * @author  lixin
 * @date    11/20/2012
 */
class SFFT {

protected:

    /**
     * The Param class that administrates paramters
     */
    Params params;

    /**
     * The Initialization class that does peak estimation process; i.e. step one
     */
    Initialization initialization;

    /**
     * The Gradient class that does gradient search process; i.e. step twp
     */
    Gradient gradient;

    /**
     * The data source: where we get the input to sparse fft
     */
    DataSource * source;

    /**
     * The data sink: where we write the output of sparse fft
     */
    DataSink * sink;

    /**
     * The pattern for voting. It is based on the dimensionality 
     * gap of 4D light field
     */
    bool * pattern;

    /**
     * The parallel job id that is being run from this instance
     *
     * -1 means it is a job master
     */
    int32_t paral_id;

    /**
     * The total number of parallel job allowed
     */
    uint32_t num_parals;

public: 

    /**
     * The constructor. 
     *
     * It initialize everything.
     *
     * @param[in]   param_file  the file that specifies paramters
     * @param[in]   _paral_id   the parallel id that is used for this job. 
     *                          When between 0 and _num_parals-1, it means a specific job in a job pool;
     *                          -1 means forking into _num_parals jobs and invoking them
     * @param[in]   _num_parals the total number of parallel jobs
     * @param[in]   io_prefix   the folder to which the io is dealing with
     */
    SFFT(in char * param_file, in int32_t _paral_id = 0, in uint32_t _num_parals = 1, in const char * io_prefix = NULL) // , in DataSource * _source, in DataSink * _sink) 
        : params(parseFileName(param_file), io_prefix), initialization(&params), gradient(&params, &initialization), paral_id(_paral_id), num_parals(_num_parals) {

        switch (params.control_params.input_method) {
        case Params::ControlParams::FILE_INPUT:
            source = new FileInput(paral_id, num_parals, parseFileName(params.control_params.input_file_prefix));
            break;
        default:
            fprintf(stderr, "initializing SFFT: wrong input_method");
            exit(-1);
        }

        sink = new FileOutput(paral_id, num_parals, parseFileName(params.control_params.result_file_prefix), parseFileName(params.control_params.init_file_prefix));

        source->init();
        sink->init(source);
       
        uint32_t n1 = 0;
        uint32_t n2 = 0;
        uint32_t n3 = 0;
        uint32_t n4 = 0;
        source->getDims(&n1, &n2, &n3, &n4);
        
        printf("n1 = %u, n2 = %u, n3 = %u, n4 = %u\n", n1, n2, n3, n4);

        if (n1 != NonIntSFFT::n_v || n2 != NonIntSFFT::n_h) {
            fprintf(stderr, "error n1 %d, n2 %d, do not match with NonIntSFFT::n_v %d, NonIntSFFT::n_h %d\n", n1, n2, NonIntSFFT::n_v, NonIntSFFT::n_h);
			exit(-1);
        }

        // allocate space for pattern
        if (params.lightfield_params.impose_dim_gap)
            pattern = (bool *)malloc(sizeof(bool) * NonIntSFFT::n);
    }

    /**
     * Deallocate source and sink
     */
    ~SFFT() {

        delete source;
        delete sink;
        if (pattern)
           free(pattern);
    }


    int32_t getBlockNum(int32_t x, int32_t y, int32_t n) {
	    if( n % 2 == 0) {
            if (y >= n / 2) { // block 2 or 4
				if(x >= n/2) {
				return 4;
				}
				else {
				return 2;
				}
			}
	        else {   // block 1 or 3
				if (x >= n/2) {
				return 3;
				}	
				else {
					return 1;
				}
			}
        } 
        else if ( n % 2 == 1) {
			if ( y >= (n + 1) / 2) {
				if ( x >= (n + 1) / 2) {
				return 4;
				}
				else {
				return 2;
				}
			}
			else {
				if( x >= (n + 1) / 2) {
				return 3;
				}
				else {
				return 1;
				}
			}
        }

		return -1;
    }    



    void getCorrespondPos(int32_t * old_x, int32_t * old_y, int32_t * new_x, int32_t * new_y, int32_t n) {

        //first decide which block it is
        int32_t block_num = getBlockNum(*old_x, *old_y, n);
		if(block_num == -1) {
			fprintf(stderr, "Fatal error, block_num is -1\n");
		}


        if (n % 2 == 0) {
	    int32_t off_x, off_y;
	    switch(block_num) {
	    case 1:
	        off_x = *old_x;
            off_y = *old_y;
	        if (off_x == 0 && off_y == 0) {    // must calc, no corresponding point
		    *new_x = -1;
		    *new_y = -1;
	        }
            else if (off_x == 0) {             // block 2
		    *new_x = 0;
		    *new_y = n - off_y;
	        }
	        else if (off_y == 0) {             // block 3
		    *new_x = n - off_x;
		    *new_y = 0;
                }
	        else {                             // block 4
		    *new_x = n - off_x;
		    *new_y = n - off_y;
	        }
	    break;
            case 2:
	        off_x = *old_x;
	        off_y = *old_y - n / 2;
                if (off_x == 0 && off_y == 0) {     // must calc, no corresponding point
		    *new_x = -1;    
		    *new_y = -1;
	        }
	        else if (off_x == 0) {             // block 1
		    *new_x = 0;
		    *new_y = n / 2 - off_y;
	        }
	        else if (off_y == 0) {             // block 4
		    *new_x = n / 2;
		    *new_y = n - off_x;
	        }
	        else {                             // block 3
		    *new_x = n - off_x;
		    *new_y = n / 2 - off_y;
	        }
	    break;
	    case 3:
	        off_x = *old_x - n / 2;
	        off_y = *old_y;
 	        if (off_x == 0 && off_y == 0) {     // must calc, no corresponding point
		    *new_x = -1;
		    *new_y = -1;
	        }
	        else if (off_x == 0) {              // block 4
		    *new_x = n / 2;
		    *new_y = n - off_y;
	        }
	        else if (off_y == 0) {              // block 1
		    *new_x = n / 2 - off_x;
		    *new_y = 0;
	        }
	        else {                              // block 2
		    *new_x = n / 2 - off_x;
			*new_y = n - off_y;
	        }
            break;
	    case 4:
			off_x = *old_x - n / 2;
			off_y = *old_y - n / 2;
	        if (off_x == 0 && off_y == 0) {     // must calc, no corresponding point
		    *new_x = -1;
		    *new_y = -1;
	        }
	        else if (off_x == 0) {              // block 3
		    *new_x = n / 2;
		    *new_y = n / 2 - off_y;
	        }
	        else if (off_y == 0) {              // block 2
		    *new_x = n / 2 - off_x;
		    *new_y = n / 2;
	        }
	        else {
		    *new_x = n / 2 - off_x;         // block 1
		    *new_y = n / 2 - off_y;
	        }
	    break;
            default: 
	    break;
        }
    }
    else if (n % 2 == 1) {
	int32_t off_x, off_y;
	switch (block_num) {
	    case 1:
		off_x = *old_x;
	    off_y = *old_y;
		if (off_x == 0 && off_y == 0) {     // must calc, no corresponding point
		*new_x = -1;
		*new_y = -1;
	    }
	    else if (off_x == 0) {              // block 2
		*new_x = 0;
		*new_y = n - off_y;
	    }
	    else if (off_y == 0) {              // block 3
		*new_x = n - off_x;
		*new_y = 0;
	    }
	    else {                              // block 4
		*new_x = n - off_x;
		*new_y = n - off_y;
	    }
	    break;
	    case 2:
	    off_x = *old_x;
	    off_y = *old_y - (n + 1) / 2;
		if(off_x == 0) {                    // block 1
		*new_x = 0;
		*new_y = (n - 1) / 2 - off_y;
	    }
	    else {                              // block 3
		*new_x = n - off_x;
		*new_y = (n - 1) / 2 - off_y;
	    }
		break;
	    case 3:
		off_x = *old_x - (n + 1) / 2;
		off_y = *old_y;
		if (off_y == 0) {                   // block 1
		*new_x = (n - 1) / 2 - off_x;
		*new_y = 0;
	    }
	    else {                              // block 2
		*new_x = (n - 1) / 2 - off_x;
		*new_y = n - off_y;
	    }
		break;
	    case 4:                             // block 1
	    off_x = *old_x - (n + 1) / 2;
		off_y = *old_y - (n + 1) / 2;
		*new_x = (n - 1) / 2 - off_x;
		*new_y = (n - 1) / 2 - off_y;
	    break;
		default:
		break;
	}
    }

        return;

}





    /**
     * The runner: 
     *
     * It runs Sparse FFT. 
     *
     */
    void run() {

	
        ComplexPtr input = NULL;
        ComplexPtr output = NULL;
        PeaksList peak_list;
        PeaksList * peaks = &peak_list;
        PeaksList init_peaks_list;
        PeaksList * init_peaks = &init_peaks_list;
        Real init_residue_value;
        Real * init_residue = &init_residue_value;
        Real residue_value;
        Real * residue = &residue_value;

		GO_AWAY_UNUSED_VARIABLE(output);
		GO_AWAY_UNUSED_VARIABLE(init_peaks);
		GO_AWAY_UNUSED_VARIABLE(init_residue);
		GO_AWAY_UNUSED_VARIABLE(residue);

		uint32_t camera_v;
		uint32_t camera_h;
		uint32_t input_v;
		uint32_t input_h;
		source->getDims(&camera_v, &camera_h, &input_v, &input_h);
		
        clock_t start_t, end_t;
        double time;
        
        
        uint32_t * indices = NULL;
        
        PeaksList init_peaks_cpy = params.init_param.init_peaks;
        uint32_t num_init_peaks_cpy = params.init_param.num_init_peaks;

      
		int32_t i_x, i_y;
		GO_AWAY_UNUSED_VARIABLE(i_x);
		GO_AWAY_UNUSED_VARIABLE(i_y); 

		
        while ((input = source->getNextSlice(paral_id))) {

			NonIntSFFT::error = false; // important change 07/09
            start_t = clock();
            
			if (params.lightfield_params.impose_dim_gap) {
				indices = source->getCurIndices(paral_id);
				if (indices) {
					generatePattern(indices);
				} 
			}
			
			if (peaks && (peaks->size() != 0) && params.init_param.use_adjacent /* && ((abs(i_x - indices[0]) % NonIntSFFT::n_v) + (abs(i_y - indices[1]) % NonIntSFFT::n_h) == 1) */) {
				fprintf(stderr, "copy from adjacent--------------------------------------");
				params.init_param.init_peaks = *peaks;
				params.init_param.num_init_peaks = peaks->size();
				for (PeaksList::iterator it = params.init_param.init_peaks.begin(); it != params.init_param.init_peaks.end(); ++ it)
					it->setValue(0);
				fprintf(stderr, "copy finished\n");
			} else {
				params.init_param.init_peaks = init_peaks_cpy;
				params.init_param.num_init_peaks = num_init_peaks_cpy;
			}
			
			if (indices) {
				i_x = indices[0];
				i_y = indices[1];
			}
			
			source->getInitPeaks(params.init_param.init_peaks, params.init_param.num_init_peaks); 
			
			uint32_t num_peaks = 0;
			int32_t recovery_id = 0;
		
			initialization.setInitState(input, &recovery_id);
			initialization.getPeaksResidue(&num_peaks, peaks, residue);
			 
			gradient.setRecovery(&params, input);  
			gradient.setGradient(num_peaks, *peaks, *residue, &recovery_id);
			
			gradient.getSampleValues(input);
			gradient.gradientSearch(input, recovery_id);
			gradient.getResults(&output, &peaks, &residue);
			
			gradient.getInitConfig(&init_peaks, &init_residue);
		
			
			end_t = clock();
			time = (double)(end_t - start_t) / CLOCKS_PER_SEC;
			fprintf(stderr, "TOTAL TIME = %lf\n", time);
			fprintf(stderr, "-------------------------------------------------------------------------------------\n");	
			sink->writeNextSlice(paral_id, output, peaks, residue, &time, &NonIntSFFT::error, init_peaks, init_residue);
			
					
     	}

    }
    
protected:

    /**
     * Given a file name, judge whether it is NULL or not
     * @param[in]   file_name   the input file name
     * @return      NULL if file_name is "NULL" or is empty, 
     *              otherwise file_name itself
     */
    inline char * parseFileName(char * file_name) forceinline {
        
        //fprintf(stderr, "filename = %s\n", file_name);

        if (!file_name)
            return NULL;

        if (strlen(file_name) == 0)
            return NULL;
            
        if (!strcmp(file_name, "NULL"))
            return NULL;
            
        return file_name;
    }

    /**
     * Generate the pattern used for voting given the current position. 
     * The result is written in SFFT::pattern
     *
     * @param   indices     the pointer to the indices. Here we require the length 
     *                      of the indices to be 2
     */
    inline void generatePattern(in const uint32_t * indices) forceinline {

        if (!indices || !pattern)
            return;

        memset(pattern, 0, sizeof(bool) * NonIntSFFT::n);
       
        int32_t w_x = 0; 
        int32_t w_y = 0; 

        if (params.lightfield_params.camera_pixel_consistent) {
            w_x = ((*indices) + params.lightfield_params.n_x / 2) % params.lightfield_params.n_x - params.lightfield_params.n_x / 2;
            w_y = ((*(indices + 1)) + params.lightfield_params.n_y / 2) % params.lightfield_params.n_y - params.lightfield_params.n_y / 2;
        } else {
            w_y = ((*indices) + params.lightfield_params.n_x / 2) % params.lightfield_params.n_x - params.lightfield_params.n_x / 2;
            w_x = ((*(indices + 1)) + params.lightfield_params.n_y / 2) % params.lightfield_params.n_y - params.lightfield_params.n_y / 2;
        }

        fprintf(stderr, "w_x = %d, w_y = %d.\n", w_x, w_y);

        if ((w_x == 0) && (w_y == 0)) {
            pattern[0] = 1;
            return;
        }

        double norm = sqrt(w_x * w_x + w_y * w_y);

        uint32_t w_u = 0; 
        uint32_t w_v = 0;
        int32_t w_u_shifted = 0;
        int32_t w_v_shifted = 0;
        bool * cur_pattern = pattern;
        for (w_u = 0; w_u < NonIntSFFT::n_v; w_u ++) {
            for (w_v = 0; w_v < NonIntSFFT::n_h; w_v ++) {
                w_u_shifted = (w_u + NonIntSFFT::n_v / 2) % NonIntSFFT::n_v - NonIntSFFT::n_v / 2;
                w_v_shifted = (w_v + NonIntSFFT::n_h / 2) % NonIntSFFT::n_h - NonIntSFFT::n_h / 2;
                if (fabs((double)(w_u_shifted * w_y + w_v_shifted * w_x) / norm) < params.lightfield_params.lam_thickness)
                    * cur_pattern = 1;
                cur_pattern ++;
            }
        }


    }

    /**
     * This function is used to print out the results of step one, namely all the peaks in our estimation
     *
     *@in peaks ptr to the result of step one
     */
     inline void printStep1Result(PeaksList * peaks) forceinline {
		PeakIterator cur_peak = peaks->begin();
		int peak_index = 0;
		int x,y;
                int size = peaks->size();
		for (; peak_index < size; peak_index ++, cur_peak ++) {
			x = cur_peak->getX();
			y = cur_peak->getY();
			printf("The pos of the %d peak is X: %d, Y: %d\n", peak_index, x, y); 
		}
		return;
	}


    

};




#endif
