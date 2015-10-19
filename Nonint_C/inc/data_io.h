#ifndef _DATA_IO_H_
#define _DATA_IO_H_

/**
 * @file    data_io.h
 * @brief   Dealing with input/output
 * @author  lixin
 * @date    11/20/2012
 */

#include "common.h"
#include "peaks.h"
#include "gradient.h"

/**
 * @brief   Abstract class for data input
 *
 * DataSource effectively defines the interface for data input
 *
 * @author  lixin
 * @date    11/20/2012
 */
class DataSource {
    
public:
    
    /**
     * The total number of parallel computations
     */
    uint32_t num_parals;

    /**
     * The current parallel id
     *
     * -1 means it is a job master
     */
    int32_t paral_id;

    /**
     * The vertical dimension of the slice
     */
    uint32_t n1;

    /**
     * The horizontal dimension of the slice
     */
    uint32_t n2;

    /**
     * The third dimension
     */
    uint32_t n3;

    /**
     * The forth dimension
     */
    uint32_t n4;

    /**
     * The length of each index
     */
    uint32_t index_length;

    /**
     * Number of data dimensions. 
     */
    uint32_t num_data_dims;

    /**
     * The total number of slices
     */
    uint32_t num_slices;

public:

    /**
     * Default constructor
     * @param[in]   _paral_id   the parallel id that is used for this job. 
     *                          When between 0 and _num_parals-1, it means a specific job in a job pool;
     *                          -1 means forking into _num_parals jobs and invoking them
     * @param[in]   _num_parals the total number of parallel jobs
     */
    DataSource(in int32_t _paral_id, in uint32_t _num_parals) : num_parals(_num_parals), paral_id(_paral_id), n1(0), n2(0), n3(0), n4(0), index_length(0), num_data_dims(0), num_slices(0) {
        
        if (paral_id == -1) {
            exit(-1);
        }
    
    };

    /**
     * Just the destructor
     */
    virtual ~DataSource() {}
    
    /**
     * Initialization of the whole class
     */
    virtual void init() = 0;

    /**
     * Get the next slice for a specific parallel computation
     *
     * @param[in]   paral    the id for the parallel computation, 
     *                       a number between 0 and (DataSource::num_parals - 1)
     * @return  the pointer to the next slice
     *          NULL if there is no more slice
     */
    virtual ComplexPtr getNextSlice(in uint32_t paral) = 0;

    /**
     * Get the current index for a specific parallel computation
     *
     * @param[in]   paral    the id for the parallel computation, 
     *                       a number between 0 and (DataSource::num_parals - 1)
     * @return  the pointer to the indices matrix
     */
    virtual uint32_t * getCurIndices(in uint32_t paral) = 0;

    /**
     * Get the dimensions here
     *
     * @param[out]  _n1     the first dimension
     * @param[out]  _n2     the second dimension
     * @param[out]  _n3     the third dimension
     * @param[out]  _n4     the forth dimension
     */
    void getDims(out uint32_t * _n1, out uint32_t * _n2, out uint32_t * _n3, out uint32_t * _n4) {

        *_n1 = n1;
        *_n2 = n2;
        *_n3 = n3;
        *_n4 = n4;

    }

    /**
     * It will read the next init peaks
     * 
     * @param[out]  peak_list   the PeaksList
     * @param[out]  num_peaks   the number of peaks
     *
     * @return whether we get the peaks or not
     *
     */
    virtual bool getInitPeaks(out PeaksList peak_list, out uint32_t & num_peaks) {
        GO_AWAY_UNUSED_VARIABLE(peak_list);
        GO_AWAY_UNUSED_VARIABLE(num_peaks);
        return false;
    }
};


/**
 * @brief   Abstract class for data output
 *
 * DataSink effectively defines the interface for data output
 *
 * @author  lixin
 * @date    11/20/2012
 */
class DataSink {

protected:
    
    /**
     * The total number of parallel computations
     */
    uint32_t num_parals;

    /**
     * The current parallel id
     *
     * -1 means it is a job master
     */
    int32_t paral_id;

    /**
     * The data source that connects to this sink
     */
    DataSource * source;

public:

    /**
     * Default constructor
     * @param[in]   _paral_id   the parallel id that is used for this job. 
     *                          When between 0 and _num_parals-1, it means a specific job in a job pool;
     *                          -1 means forking into _num_parals jobs and invoking them
     * @param[in]   _num_parals the total number of parallel jobs
     */
    DataSink(in int32_t _paral_id, in uint32_t _num_parals) : num_parals(_num_parals), paral_id(_paral_id), source(0) {
        
        if (paral_id == -1) {
            fprintf(stderr, "DataSource: parallel computation master is not implemented yet");
            exit(-1);
        }
    
    }

    /**
     * Just the destructor
     */
    virtual ~DataSink() {}
    
    /**
     * Initialization of the whole class
     *
     * @param[in]   _source the data source that connects to this sink
     */
    virtual void init(in DataSource * _source) = 0;

    /**
     * Get the next slice for a specific parallel computation
     *
     * @param[in]   paral           a number between 0 and (DataSource::num_parals - 1)
     * @param[in]   x               the result to write
     * @param[in]   peaks           the list of peaks to write
     * @param[in]   residue         the residue to write
     * @param[in]   time            the time to write
     * @param[in]   error           whether there is error or not
     * @param[in]   init_peaks      the initial peak configuration, default is NULL
     * @param[in]   init_residue    the initial residue, default is NULL
     */
    virtual void writeNextSlice( in uint32_t paral, in ComplexPtr x, in PeaksList * peaks, in Real * residue, in double * time, in bool *error, in PeaksList * init_peaks = NULL, in Real * init_residue = NULL) = 0;

};


/**
 * @brief   An instantiation of DataSource, namely reading from a data file
 *
 * @author  lixin   
 * @date    11/20/2012
 */
class FileInput : public DataSource {

protected:

    /**
     * The pointer to the file from which we reads the input
     */
    FILE * file;

    /**
     * The pointer to which the init file will be read
     */
    FILE * init_file;

    /**
     * The current x buffer
     */
    ComplexPtr x;

    /**
     * The current indices buffer
     */
    uint32_t * indices;

    /**
     * Current number of slices
     */
    uint32_t slices_read;

    /**
     * The manual initial position
     *
     * Default : empty
     */
    PeaksList init_peaks;

    /**
     * The number of initial peaks
     *
     * Default: 0
     */
    uint32_t num_init_peaks;

public:
    
    /**
     * The constructor: 
     * It reads the header information of the file
     *
     * @param[in]   _paral_id   the parallel id that is used for this job. 
     *                          When between 0 and _num_parals-1, it means a specific job in a job pool;
     *                          -1 means forking into _num_parals jobs and invoking them
     * @param[in]   _num_parals the total number of parallel jobs
     * @param[in]   file_prefix the file to read
     */
    FileInput(in int32_t _paral_id, in uint32_t _num_parals, in const char * file_prefix) : DataSource(_paral_id, _num_parals), file(NULL), init_file(NULL), x(0), indices(0), slices_read(0), 
        init_peaks(), num_init_peaks(0) {

        char file_name_buffer[MAX_STRING_LENGTH];

        sprintf(file_name_buffer, "%s_%04d.dat", file_prefix, paral_id);

        file = fopen(file_name_buffer, "rb");
        
        if (!file) {
            fprintf(stderr, "Cannot open file %s", file_name_buffer);
            exit(-1);
        } else {
            fprintf(stderr, "Reading file %s\n", file_name_buffer);
        }

        // read the header

        if (fread(&num_slices, sizeof(uint32_t), 1, file) != 1) {
            fprintf(stderr, "IO error: reading header of %s", file_name_buffer);
            exit(-1);
        }

        if (fread(&index_length, sizeof(uint32_t), 1, file) != 1) {
            fprintf(stderr, "IO error: reading header of %s", file_name_buffer);
            exit(-1);
        }
        if (fread(&num_data_dims, sizeof(uint32_t), 1, file) != 1) {
            fprintf(stderr, "IO error: reading header of %s", file_name_buffer);
            exit(-1);
        }
       
        if (num_data_dims < 3 || num_data_dims > 4) {
            fprintf(stderr, "num_data_dims = %d is not supported", num_data_dims);
            exit(-1);
        }

        if (fread(&n1, sizeof(uint32_t), 1, file) != 1) {
            fprintf(stderr, "IO error: reading header of %s", file_name_buffer);
            exit(-1);
        }
        if (fread(&n2, sizeof(uint32_t), 1, file) != 1) {
            fprintf(stderr, "IO error: reading header of %s", file_name_buffer);
            exit(-1);
        }
        if (fread(&n3, sizeof(uint32_t), 1, file) != 1) {
            fprintf(stderr, "IO error: reading header of %s", file_name_buffer);
            exit(-1);
        }

        if (num_data_dims == 3)
            n4 = 1;
        else {
            if (fread(&n4, sizeof(uint32_t), 1, file) != 1) {
                fprintf(stderr, "IO error: reading header of %s", file_name_buffer);
                exit(-1);
            }
        }

        sprintf(file_name_buffer, "%s_init_%04d.dat", file_prefix, paral_id);
        init_file = fopen(file_name_buffer, "rb");

        // initialize the buffer
        x = (ComplexPtr)malloc(sizeof(Complex) * n1 * n2);
        indices = (uint32_t *)malloc(sizeof(uint32_t) * index_length);
    }


    /**
     * Destructor. 
     * It deallocates all buffers:
     * FileInput::x, FileInput::indices
     *
     * Also it closes the file
     */
    ~FileInput() {
        if (x)
            free(x);
        
        if (indices)
            free(indices);

        if (file)
            fclose(file);
    }

    /**
     * It does nothing
     */
    void init() {
        // basically doing nothing
    }

    /**
     * It will read the next init peaks from the file
     * 
     * @param[out]  peak_list   the PeaksList
     * @param[out]  num_peaks   the number of peaks
     *
     * @return whether we get the peaks or not
     *
     */
    bool getInitPeaks(out PeaksList peak_list, out uint32_t & num_peaks) {

        if (!init_file)
            return false;

        uint32_t * init_indices = (uint32_t *)malloc(sizeof(uint32_t) * index_length);
        if (index_length) {
            if (fread(init_indices, sizeof(int32_t), index_length, file) != index_length) {
                fprintf(stderr, "IO error: reading the init file");
                exit(-1);
            }

            if (memcpy(init_indices, indices, sizeof(uint32_t) * index_length) != 0) {
                fprintf(stderr, "IO error: index from init file is different");
                exit(-1);
            }
        }

        if (fread(&num_peaks, sizeof(uint32_t), 1, file) != 1) {
            fprintf(stderr, "IO error: reading the init file");
            exit(-1);
        }

        double * x_indices = (double *)malloc(sizeof(double) * num_peaks);
        double * y_indices = (double *)malloc(sizeof(double) * num_peaks);

        if (fread(x_indices, sizeof(double), num_peaks, file) != 1) {
            fprintf(stderr, "IO error: reading the data");
            exit(-1);
        }

        if (fread(y_indices, sizeof(double), num_peaks, file) != 1) {
            fprintf(stderr, "IO error: reading the data");
            exit(-1);
        }

        peak_list.clear();
        for (uint32_t i = 0; i < num_peaks; i ++) {
            peak_list.push_back(NonIntPeak(x_indices[i], y_indices[i]));
        }

        free(init_indices);
        free(x_indices);
        free(y_indices);

        return true;
    }

    /**
     * Read from data the next slice as well as the indices
     *
     *
     * @param[in]   paral    the id for the parallel computation, 
     *                       a number between 0 and (DataSource::num_parals - 1)
     * @return  the pointer to the next slice; 
     *          NULL if there is no more slice
     */
    ComplexPtr getNextSlice(in uint32_t paral) {

        // fprintf(stderr, "JUST START GET NEXT SLICE\n");
        if (slices_read >= num_slices)
            return NULL;

        if ((paral_id >= 0 && (uint32_t)paral_id != paral) || (paral_id < 0 && paral >= num_parals)) {
            fprintf(stderr, "Invalid parallel computation ID %d", paral); 
            exit(-1);
        }
        
        // read index matrix
        if (index_length) {
            if (fread(indices, sizeof(int32_t), index_length, file) != index_length) {
                fprintf(stderr, "IO error: reading the data");
                exit(-1);
            }
        }

        fprintf(stderr, "%d %d\n", indices[0], indices[1]);

        // read data matrix
        double * inner_data_real = (double *)x;
        double * inner_data_imag = ((double *)x) + 1;
        uint32_t x_index = 0;
        uint32_t slice_size = n1 * n2;

        for (x_index = 0; x_index < slice_size; ++ x_index) {

            if (fread(inner_data_real, sizeof(double), 1, file) != 1) {
                fprintf(stderr, "IO error: reading the data");
                exit(-1);
            }
            inner_data_real += 2;

        }

        for (x_index = 0; x_index < slice_size; ++ x_index) {

            if (fread(inner_data_imag, sizeof(double), 1, file) != 1) {
                fprintf(stderr, "IO error: reading the data");
                exit(-1);
            }
            inner_data_imag += 2;
        }

        slices_read ++;

        return x;
    }

    /**
     * Get the current index for a specific parallel computation
     *
     * @param[in]   paral    the id for the parallel computation, 
     *                       a number between 0 and (DataSource::num_parals - 1)
     * @return  the pointer to the indices matrix
     */
    uint32_t * getCurIndices(in uint32_t paral) {

        if (paral >= num_parals) {
            fprintf(stderr, "Invalid parallel computation ID %d", paral); 
            exit(-1);
        }
        return indices;
    }
};

/**
 * @brief   An instantiation of the DataSink class, writing result to a file
 *
 * @author  lixin   
 * @date    11/20/2012
 */
class FileOutput : public DataSink {

protected:
    
    /**
     * The pointer to the file to which we write the results
     */
    FILE * file;

    /**
     * The pointer to the file to which we write the initial configuraiton information
     */
    FILE * init_config;

public:

    /**
     * The constructor
     * It opens the file necessary for writing.
     *
     * @param[in]   _paral_id   the parallel id that is used for this job. 
     *                          When between 0 and _num_parals-1, it means a specific job in a job pool;
     *                          -1 means forking into _num_parals jobs and invoking them
     * @param[in]   _num_parals the total number of parallel jobs
     * @param[in]   file_prefix         the file to write to
     * @param[in]   debugging_file_prefix the file to write to debugging information
     * @param[in]   init_file_prefix      the file to write to initial configuraiton information
     */
    FileOutput(in int32_t _paral_id, in uint32_t _num_parals, in const char * file_prefix, in const char * init_file_prefix = NULL)
        : DataSink(_paral_id, _num_parals) {

        char file_name_buffer[MAX_STRING_LENGTH];

        if (!file_prefix)
            file = NULL;
        else {

            sprintf(file_name_buffer, "%s_%04d.dat", file_prefix, paral_id);
            file = fopen(file_name_buffer, "wb");
            
            if (!file) {
                fprintf(stderr, "Cannot open file %s", file_name_buffer);
                exit(-1);
            } else {
                fprintf(stderr, "Writing file %s\n", file_name_buffer);
            }
        }

        if (!init_file_prefix)
            init_config = NULL;
        else {
            sprintf(file_name_buffer, "%s_%04d.init", init_file_prefix, paral_id);
            init_config = fopen(file_name_buffer, "wb");

            if (!init_config) {
                fprintf(stderr, "Cannot open file %s", file_name_buffer);
                exit(-1);
            } else {
                fprintf(stderr, "Writing file %s\n", file_name_buffer);
            }
        }
    }

    /**
     * The destructor
     *
     * It closes the file
     */
    ~FileOutput() {
        if (file)
            fclose(file);

        if (init_config)
            fclose(init_config);
    }
    
    /**
     * Initialization of the whole class. 
     *
     * Basically it writes the header.
     *
     * @param[in]   _source the data source that connects to this sink
     */
    void init(in DataSource * _source) {
    
        source = _source;

        if (file) {
            if (fwrite(&source->num_slices, sizeof(uint32_t), 1, file) != 1)
                fprintf(stderr, "IO error: writing header");

            if (fwrite(&source->index_length, sizeof(uint32_t), 1, file) != 1)
                fprintf(stderr, "IO error: writing header");
            if (fwrite(&source->num_data_dims, sizeof(uint32_t), 1, file) != 1)
                fprintf(stderr, "IO error: writing header");
            if (fwrite(&source->n1, sizeof(uint32_t), 1, file) != 1)
                fprintf(stderr, "IO error: writing header");
            if (fwrite(&source->n2, sizeof(uint32_t), 1, file) != 1)
                fprintf(stderr, "IO error: writing header");
            if (fwrite(&source->n3, sizeof(uint32_t), 1, file) != 1)
                fprintf(stderr, "IO error: writing header");
            if (source->num_data_dims > 3) {
                if (fwrite(&source->n4, sizeof(uint32_t), 1, file) != 1)
                    fprintf(stderr, "IO error: writing header");
            }
        }
    }


	/**
	 * Finding the block number of a specific position (x, y) in a square of dimension n*n
     *
	 * @param[in]  x  the x coordinate of the position
	 * @param[in]  y  the y coordinate of the position
	 * @param[in]  n  the width/height of the square
	 *
	 * @return     the block number of the position
	 */

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

	/**
	 * Finding the corresponding position in copying. For more information about copying, please take a look at the detailed documentation
	 *
	 * @param[in] old_x  pointer to the current x position
	 * @param[in] old_y  pointer to the current y position
	 * @param[in] new_x  pointer to the matching x position
	 * @param[in] new_y  pointer to the matching y position
	 * @param[in] n      width/height of the square
	 */

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
		    *new_x = *old_x;
		    *new_y = *old_y;
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
		    *new_x = *old_x;    
		    *new_y = *old_y;
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
		    *new_x = *old_x;
		    *new_y = *old_y;
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
		    *new_x = *old_x;
		    *new_y = *old_y;
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
			*new_x = *old_x;
			*new_y = *old_y;
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

    
    void flipResult(in ComplexPtr x, out ComplexPtr y, int32_t n) {
		
		int32_t new_x = 0;
		int32_t new_y = 0;
		for (int32_t i = 0; i < n; i ++) {
			for (int32_t j = 0; j < n; j ++) {
				getCorrespondPos(&i, &j, &new_x, &new_y, n);
                //fprintf(stderr, "old_pos (%d, %d) new_pos (%d, %d)\n", i , j, new_x, new_y);
				y[i * n + j] = x[new_x * n + new_y];
			}
		}
	}
   





    /**
     * Get the next slice for a specific parallel computation
     *
     * @param[in]   has_peak_info    whether this slice has peaks and residue info or not
     * @param[in]   paral   a number between 0 and (DataSource::num_parals - 1)
     * @param[in]   x       the result to write
     * @param[in]   peaks   the list of peaks to write
     * @param[in]   residue the residue to write
     * @param[in]   time    the time to write
     * @param[in]   error           whether there is error or not
     * @param[in]   init_peaks      the initial peak configuration, default is NULL
     * @param[in]   init_residue    the initial residue, default is NULL
     */
    void writeNextSlice(in uint32_t paral, in ComplexPtr x, in PeaksList * peaks, in Real * residue, in double * time, in bool *error, in PeaksList * init_peaks/* = NULL*/, in Real * init_residue/* = NULL*/) {

        if ((paral_id >= 0 && (uint32_t)paral_id != paral) || (paral_id < 0 && paral >= num_parals)) {
            fprintf(stderr, "Invalid parallel computation ID %d", paral); 
            exit(-1);
        }


        uint32_t slice_size = source->n1 * source->n2;
        double * ptr_real = (double *)x;
        double * ptr_imag = ((double *)x) + 1;
		//double temp = 0.0;
		//double * ptr_imag_temp = &temp;
        uint32_t * ptr_indices = source->getCurIndices(paral);
        
		// new result 
        if (file) {
            // write the index of the slice
            if (source->index_length) {
                if (fwrite(ptr_indices, sizeof(uint32_t), source->index_length, file) != source->index_length)
                    fprintf(stderr, "IO error: writing the data");
            }

            // write the data of that slice
            uint32_t x_index = 0;
            for (x_index = 0; x_index < slice_size; ++ x_index) {

                if (fwrite(ptr_real, sizeof(double), 1, file) != 1)
                    fprintf(stderr, "IO error: writing the data");
                ptr_real += 2;

            }

            for (x_index = 0; x_index < slice_size; ++ x_index) {
                if (fwrite(ptr_imag, sizeof(double), 1, file) != 1)
                    fprintf(stderr, "IO error: writing the data");
                ptr_imag += 2;
            }
			
        }


        if (init_config) {
            uint32_t num_peaks = init_peaks->size();
            if (fwrite(&num_peaks, sizeof(uint32_t), 1, init_config) != 1)
                fprintf(stderr, "IO error: writing the init infomation");

            for (PeakIterator it = init_peaks->begin(); it != init_peaks->end(); ++ it) {
                double dummy;
                dummy = it->getX();
                if (fwrite(&dummy, sizeof(double), 1, init_config) != 1)
                    fprintf(stderr, "IO error: writing the init infomation");
                dummy = it->getY();
                if (fwrite(&dummy, sizeof(double), 1, init_config) != 1)
                    fprintf(stderr, "IO error: writing the init infomation");
                Complex dummy2 = it->getV();
                if (fwrite(&dummy2, sizeof(double), 2, init_config) != 2)
                    fprintf(stderr, "IO error: writing the init infomation");
            }

            if (fwrite(init_residue, sizeof(Real), 1, init_config) != 1)
                fprintf(stderr, "IO error: writing the init infomation");

            double dummy = 0;
            if (fwrite(&dummy, sizeof(double), 1, init_config) != 1)
                fprintf(stderr, "IO error: writing the init infomation");

            bool dummy2 = 0;
            if (fwrite(&dummy2, sizeof(bool), 1, init_config) != 1)
                fprintf(stderr, "IO error: writing the debug infomation");
        }

    }
};




#endif
