#ifndef _FILTER_H_
#define _FILTER_H_

/**
 * @file    projection.h
 * @brief   Defines class Projection
 * @author  Tianwei
 * @date    07/10/2014
 */

#include "common.h"
#include "peaks.h"
#include "fftw3.h"

class NonIntRecover;
class Gradient;

/**
 * @brief This class is used to implement step one in projection slice style.
 * Each projection slice/line is decided by 4 interger a, b, c, d in the following way:
 * u = a * t + b mod N
 * v = c * t + d mod N
 *
 * To use this class, there are two fundamental assumptions: 
 * 1) GCD(a, c) = 1 
 * 2) NonIntSFFT::n_v = NonIntSFFT::n_h
 *
 * @author Tianwei
 * @date 07/10/2014
 */
class Projection {

    friend class NonIntRecover;
    friend class IntRecover;
    friend class Gradient;
    friend class Initialization;

protected:
   // The four parameters in the line equation
   int32_t a;

   int32_t b;

   int32_t c;

   int32_t d;

   // number of buckets also the horizontal/vertical size of the camera location map
   int32_t n;

   /* 
    * This array is used to store the values of invert mod
    */
    int32_t * invert_mod;


   /**
    * This is fftw plan for the fftw library's use.
    * Unless we need to change projection slices, otherwise we don't need to regenerate the fftw plan.
    */
   fftw_plan plan;

   /**
    * This is just an indicator of whether fftw_plan Projection::plan has been initialized or not.
    */
    bool init_plan;

   /**
    * This is a buffer for the execution of fftw (the input)
    */
    ComplexPtr fftw_input;

   /**
    * This is a buffer to store all the positions on the projection line
    */
    IntPosition * bucket_table;

   /**
    * This is used to store the output of the fft for future use, namely voting 
    */
    ComplexPtr  buckets;

   /**
    * The pointer to the coefficient array. It is segementally 2D array but here we just use 1D array to implement it.
    */
    ComplexPtr  coeffs;

   /**
    * The power of all the buckets
    */
    Real * powers;

   
    /**
     * Index to the large buckets. 
     * Note that this array is allocated by findLargeBuckets() function
     */
    uint32_t * large_buckets;


    /**
     * The total number of large buckets, i.e. the length of the array
     * Filter::large_buckets
     */
    uint32_t num_large_buckets;

    /**
     * The indicator array of large buckets, with size w_v * w_h, 
     * where 1 means large bucket and 0 otherwise
     */
    bool * large_bucket_indicator;

    /**
     * The minimum power of all large buckets
     */
    Real min_large_power;

    /**
     * The hard power threshold. Every bucket whose power is larger than this
     * threshold is regarded as a candidate of large bucket.
     *
     * Normally we set hard power threshold to be 0 because hard power threshold
     * is not stable
     */
    double hard_power_threshold;

    /**
     * The soft power threshold. Every bucket whose power is larger than 
     * the smallest bucket mutilplied by the soft power threshold is a 
     * candidate of large bucket
     *
     */
    double soft_power_threshold;

    /**    
     * This is a hashing table: for each frequency, which bucket it hashes into
     */
     uint32_t * hashing_table;

public:

     /**
      * Constructor that initalize the projection slice, set the parameters,
      * and allocate the memory space if necessary
      */
     Projection(int32_t _a = 0, int32_t _b = 0, int32_t _c = 0, int32_t _d = 0)
        :  a(0), b(0), c(0), d(0), n(0),invert_mod(NULL), init_plan(0), fftw_input(0), bucket_table(NULL), buckets(NULL), coeffs(NULL), powers(NULL), large_buckets(NULL), num_large_buckets(0), large_bucket_indicator(NULL), 
	    min_large_power(0.0), hard_power_threshold(0.0), soft_power_threshold(0.0), hashing_table(NULL) {
		a = _a; b = _b; c = _c; d = _d;
	    n = NonIntSFFT::n_v;  // we are pretty sure that according to the paramter equation, we will have NonIntSFFT::n_v|n_h buckets
        allocateBuffers();
    }
	/**
     * Deallocation function.
     */

    ~Projection() {
		destroyPlan();
	}

  
    /**
     * Generate the time domain sampling table
     */
	void generateSamplingTable(bool * sampling_table) {
		IntPosition pos;
		for ( int32_t i = 0; i < n; i ++) {
			pos = bucket_table[i];
			sampling_table[pos.map2Index()] = 1;
		}
    }

    /**
     * Print all related parameters of the projection
     */
    void printProjection() {
	
		printf("----------------------------- PROJECTION INFO -----------------------------\n");
		printf (" a = %d, b = %d, c = %d, d = %d\n", a, b, c, d);
		printf("------------------------------ END INFO -----------------------------------\n"); 

    }

    /**
     * This function is used to calculate and set values in bucket_table, coeffs, and hashing_table
     *
     */
     void calBucketCoeffs() {
	
		int32_t f_v = 0;
		int32_t f_h = 0;

		// calculate coeffs
		ComplexPtr cur_coeff = coeffs;

		for (f_v = 0; f_v < n; f_v ++) {
			for (f_h = 0; f_h < n; f_h ++) {
			* (cur_coeff ++) = Utils::unityRootPower(b * f_v, d * f_h);
			}
		}
	
		// calculate bucket_table
		getPositions();	

        // calculate invert mod
        getInvertMod();        

		// calculating hashing_table
		uint32_t remainder = 0;
		uint32_t * cur_bucket = hashing_table;
        for (f_v = 0; f_v < n; f_v ++) {
			for (f_h = 0; f_h < n; f_h ++) {
				remainder = (f_v * a + f_h * c) % n;
				* (cur_bucket ++ ) = remainder; //why this place should be index instead of position? ------nonint_recover.h:1097
			}
		}

	}

 

    /**
     * Run the projection slice algorithm on the samples. The detailed procedure is as follows:
     * 1) Get all the positions on the line
     * 2) Do fft on the line
     * 3) Calculate power on the line
     * 4) Voting
     * 5) Get result, which is a bunch of peaks.
	 *
	 * @param[in]  x  input data
     */
	void projectionSlice(in ComplexPtr x) {

		if(!init_plan) {
			// For the flags, either FFTW_MEASURE or FFTW_ESTIMATE, 
			// FFTW_MEASURE instructs FFTW to run and measure the execution time of several 
			// FFTs in order to find the best way to compute the transform of size n. This process takes some time (usually a few seconds), 
			// depending on your machine and on the size of the transform. FFTW_ESTIMATE, on the contrary, 
			// does not run any computation and just builds a reasonable plan that is probably sub-optimal.
			plan = fftw_plan_dft_1d(n, (fftw_complex *)((void *)fftw_input), (fftw_complex *)((void *)buckets), FFTW_FORWARD, FFTW_MEASURE);

		}

        // Construct the input for fft
		constructInput(x);

		if(!init_plan) {	
			init_plan = 1;

			//execute fft
			fftw_execute(plan);
		}
		else {
			fftw_execute_dft(plan, (fftw_complex *) ((void *)fftw_input), (fftw_complex * )((void *) buckets));
		}

		// originally we need to scale up by sqrt(n) / w but here this value is 1

		// calculate the powers 
		ComplexPtr cur_bucket = buckets;
		for (int i = 0; i < n; i ++) {
			powers[i] = comAbs(*(cur_bucket ++));

		}
	}


	/**
	 * Print out bucket powers
	 */
    inline void printPowers() forceinline {
		printf("Line with a: %d, b: %d, c: %d, d:%d\n", a, b, c, d);
		for (int32_t i = 0; i < n; i ++) {
			printf("(%d, %d) with power %f\n", bucket_table[i].getX(), bucket_table[i].getY(), powers[i]);
		}
        printf("\n\n\n");
     }


     /**
     * Return the minimal power of all large buckets. 
     *
     * This information is needed by Gradient::checkEst()
     *
     * @return  the minimal power pf all large buckets
     */
    inline Real getMinLargePower() forceinline {
        return min_large_power;
    }

     /*
     * Set both the soft and hard threshold
     */
    inline void setPowerThreshold(in double _soft_power_threshold, in double _hard_power_threshold = 0) forceinline {
        
        soft_power_threshold = _soft_power_threshold;
        hard_power_threshold = _hard_power_threshold;

    }

    /**
     * Get the maximum bucket power in this filter
     *
     * @return  the maximum power
     */
    inline Real getMaxPower() const forceinline {

        Real max_power = 0;
        Real * cur_power = powers;
        int32_t b = 0;

        for (b = 0; b < n; ++b, ++ cur_power) {
            if (*cur_power > max_power)
                max_power = *cur_power;
        }

        return max_power;
    }

    /**
     * Get the minimum bucket power in this filter
     *
     * @return the minimum power
     */
    inline Real getMinPower() const forceinline {

        Real min_power = INFINITY;
        int32_t b = 0;
        Real * cur_power = powers;

        for (b = 0; b < n; ++ b, ++ cur_power) {
            if ((*cur_power) < min_power)
                min_power = *cur_power;
        }

        return min_power;
    }

    /**
     * Find the "large" bucket. 
     *
     * Large bucket means the bucket that might have peaks in it. 
     *
     * This function will write the array large_bucket_indicator, as well as
     * scalar num_large_buckets. 
     *
     * But note that it will <i>not</i> touch the array Projection::large_buckets. The reason is because of
     * optimization, we don't want to allocate and deallocate the Projection::large_buckets again and
     * again. To allocate and write Projection::large_buckets, invoke function Projection::writeLargeBuckets()
     */
    inline void indicateLargeBuckets() forceinline {
       
        Real min_power = getMinPower(); 
		// Let's find the combined power threshold of both soft and hard power threshold
        double power_threshold = Utils::max<double>(hard_power_threshold, soft_power_threshold * min_power);
        min_large_power = INFINITY;
        // now, set everything which is larger than power threshold as large bucket
        bool * cur_indicator = large_bucket_indicator;
        int32_t b = 0;
        Real * cur_power = powers;
        memset(cur_indicator, 0, sizeof(bool) * n);
        num_large_buckets = 0;

        for (b = 0, cur_power = powers; b < n; ++b, ++cur_power, ++cur_indicator) {
            if ((*cur_power) >= power_threshold) {
                *cur_indicator = 1;
                num_large_buckets ++;

			if (*cur_power < min_large_power)	min_large_power = *cur_power;
            }
        }
    }

    /**
     * Given that we have already find the large buckets, i.e.
     * function Projection::indicateLargeBuckets() has already been invoked, we 
     * need to allocate and write to buffer Projection::large_buckets
     */
    inline void writeLargeBuckets() forceinline {

        if (large_buckets)
            free(large_buckets);
            
        large_buckets = (uint32_t *)malloc(sizeof(uint32_t) * num_large_buckets);

        uint32_t * cur_index = large_buckets;

        int32_t b = 0;
        for (; b < n; ++ b)
            if (large_bucket_indicator[b])
                * (cur_index ++) = b;

    }




protected:

  
   /**
    * Use the position array bucket_table to find the value of buckets and set them into buckets
    */
    inline void constructInput(in ComplexPtr x_in) {
		IntPosition cur_pos;
			
		//construct the input to fft
		ComplexPtr cur_input = fftw_input;
		for (int32_t i = 0; i < n; i ++) {
			cur_pos.setPosition(bucket_table[i].getX(), bucket_table[i].getY());
			*(cur_input ++) = x_in[cur_pos.map2Index()];		
		}
		return;
    }


    /**
     * Get all the positions on the line
     * positions are derived from two parameter equations
     */ 
    inline void getPositions() {
		int32_t pos_x, pos_y;
		//int32_t n = NonIntSFFT::n_v;
		for (int32_t i = 0; i < n; i ++) {
			pos_x = (a * i + b) % n;
			pos_y = (c * i + d) % n; 
			bucket_table[i].setPosition(pos_x, pos_y);
		}
		
		return;
     }


     /**
      * Get all invert mod values
      */
      inline void getInvertMod() {
        invert_mod[0] = -1;
        invert_mod[1] = 1;
		for (int32_t i = 2; i < n;  i ++) {
			for (int32_t j = 2; j < n; j ++) {
				if ((i * j) % n == 1) {
					invert_mod[i] = j;
					break;
				}
			}
		}
     }




    /** 
     * Deallocate buffers, namely Projection::fftw_input, bucket_table and buckets
     */
    inline void deallocateBuffers() forceinline {
		if (fftw_input) {
			fftw_free(fftw_input);
			fftw_input = NULL;
		}

		if(bucket_table) {
			free(bucket_table);
			bucket_table = NULL;
		}
		
		if(buckets) {
			free(buckets);
			buckets = NULL;
		}

		if(powers) {
			free(powers);
			powers = NULL;
		}

		if(coeffs) {
			free(coeffs);
			coeffs = NULL;
		}

		if(hashing_table) {
				free(hashing_table);
			hashing_table = NULL;
		}

		if(large_bucket_indicator) {
			free(large_bucket_indicator);
			large_bucket_indicator = NULL;
		}

        if(invert_mod) {
		free(invert_mod);
		invert_mod = NULL;
        }

    }

    /**
     * Allocate the buffers in Projection, basically Projection::fftw_input and pos
     */
    inline void allocateBuffers() forceinline {
		if (fftw_input)
			fftw_free(fftw_input);

		fftw_input = (ComplexPtr)fftw_malloc(sizeof(Complex) * n);

		if(bucket_table) free(bucket_table);
		bucket_table = (IntPosition *) malloc(sizeof(IntPosition) * n); // There are n points
        if(buckets) free(buckets);
		buckets = (ComplexPtr) fftw_malloc(sizeof(Complex) * n); // we have n buckets
        if(powers) free(powers);
		powers = (double *) malloc(sizeof(double) * n); // record the power of each bucket

        if(coeffs) free(coeffs);
		coeffs = (ComplexPtr) malloc(sizeof(Complex) * n * n);

		if(hashing_table) free(hashing_table);
        hashing_table = (uint32_t *) malloc(sizeof(uint32_t) * n * n); // for each frequency position, we need to have one entry to show which bucket it hashes to
	
		if(large_bucket_indicator) free(large_bucket_indicator);
		large_bucket_indicator = (bool *) malloc(sizeof(bool) * n);

		if(invert_mod) free(invert_mod);
        invert_mod = (int32_t *) malloc(sizeof(int32_t) * n);

		calBucketCoeffs();

    }

    /**
     * Destroy the fftw plan Projection::plan and set Projection::init_plan back to 0
     */
    inline void destroyPlan() forceinline {
		if (init_plan) {
			fftw_destroy_plan(plan);
			init_plan = 0;
		}
    }



    public: 

    /*
     * This function is used to calculate the extended gcd of (a, b). We also assume that gcd(a, b) = 1 
     * The result is (s,t) and it satisfy that:
     * a * s + b * t = gcd(a, b) = 1
     */
    inline void extendedGCD(int32_t a, int32_t b, int32_t* s, int32_t* t) {
		if (b == 0) {
			*s = 1;
			*t = 0;
		}
		else {
			int32_t q = a / b;
			int32_t r = a % b;
			int32_t temp1, temp2;
			extendedGCD(b, r, &temp1, &temp2);
			*s = temp2;
			*t = temp1 - q * temp2;
		}
	}

    /*
     * This function is used to find positions that 
     * hash into the same bucket (t = f). We need to assume that n_v = n_h = n.
     * a will be slopeX, c will be slopeY, t will be the corresponding index of the bucket
     */

    inline void getProj(int32_t t, int32_t* _pos, int32_t* num_pos) {
		int32_t initx, inity;
        initx = inity = 0;
		int32_t n = NonIntSFFT::n_v;
		if (a==0 && c == 0) {
			fprintf(stderr, "ERROR Not a projection line\n");	
			return;
    	}
		else {
			if (a != 0 && c != 0) {
				extendedGCD(a, c, &initx, &inity);
				if ((a * initx + c * inity) % n != 1) {
				fprintf(stderr, "ERROR Projection!\n");
				return;
				}
            }
            else if (a == 0) {
				initx = 0;
				inity = invert_mod[c];    	
			}
			else if (c == 0) {
				initx = invert_mod[a];
				inity = 0;
			}
        }
		initx = ( t * initx ) % n;
		inity = ( t * inity ) % n;
		*num_pos = 0;

        int32_t * cur_pos = _pos;
		for( int32_t i = 0; i < n; i ++) {
			*( cur_pos ++) = Utils::mathModular(( - c * i + initx), n);
			*( cur_pos ++) = Utils::mathModular(( a * i + inity), n);
			*num_pos += 1;
        }  
        
	}

      


};

#endif
