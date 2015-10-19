#ifndef _UTILS_H_
#define _UTILS_H_

/**
 * @file    utils.h
 * @brief   Defining class Utils, which serves as a static class
 *          providing different utility functions
 * @author  lixin
 * @date    11/14/2012
 */

#include "common.h"
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>

// #include <list>
// #include "peaks.h"

/**
 * @brief   Utility class
 *
 * This is a static class (all of its members are static) providing
 * different miscellaneous functions
 *
 * @author  lixin
 * @date    11/14/2012
 */
class Utils {

protected:

    /**
     * The buffer for vertical sinc template 
     */
    static ComplexPtr verticalSinc;

    /**
     * The buffer for horizontal sinc template
     */
    static ComplexPtr horizontalSinc;

public:

	



    /**
     * This should be called to initialize the Utils, including
     * allocate space for the sinc template and initialize the random 
     * number generator
     */
    inline static void initUtils() forceinline {
        
        deallocateSincTemplate();
        horizontalSinc = (ComplexPtr)malloc(sizeof(Complex) * NonIntSFFT::n_h);
        verticalSinc = (ComplexPtr)malloc(sizeof(Complex) * NonIntSFFT::n_v);

        srand(time(NULL));
    }

    /**
     * This should be called to deallocate the space allocated by Utils
     * Deallocate space for the sinc template
     */
    inline static void deallocateSincTemplate() forceinline {

        if (horizontalSinc)
            free(horizontalSinc);
        
        if (verticalSinc)
            free(verticalSinc);
        
        horizontalSinc = NULL;
        verticalSinc = NULL;
    }


    /**
     * This is doing the mathematical modular operation
     *
     * Note that difference with the standard C modular (%) is about dealing
     * with negative values. This function will always return a value which 
     * is between 0 and (n - 1), i.e., following the math definition of 
     * modular. 
     *
     * @param[in]   a   the base
     * @param[in]   n   the modulus
     *
     * @return      a mod n, such that the returning value falls within 0 and (n-1);
     */
    inline static int32_t mathModular(in int32_t a, in int32_t n) forceinline {
        // First, test if a is negative
        if ((a & (0x80000000)) != 0)
            return (abs(a * n) + a) % n;
        else
            return a % n;
    }

    /**
     * Caculate the following formula: a^p mod n.
     *
     * Note that we require p is larger than 0, so that inverse modular arithmetic
     * could only be calculate by Fermat's little theorem if n is prime and a is
     * co-prime than a. Otherwise, inverseModular() must be needed to calculate
     * negative power.
     *
     * This algorithm uses binary-division method to calculate power.
     *
     * @param[in]   a   the base
     * @param[in]   p   the exponent. Note that p needs to be larger than 0.
     * @param[in]   n   the modulus
     *
     * @return      the result a^p mod n, which is within 0 and (n-1);
     *              -1 if p is smaller than 0.
     */
    inline static int32_t modularPower(in int32_t a, in int32_t p, in int32_t n) forceinline {

        if (p <= 0) {
            return -1;
        }


        int32_t x = 1;
        int32_t y = mathModular(a, n);
        int32_t b = p;

        while (b > 0) {
            // if odd
            if ((b % 2) == 1) {
                x = x * y;
                if (x > n) 
                    x %= n;
            }

            // just do quadratic of y
            y = y * y;

            if(y > n) 
                y %= n;

            b /= 2;
        }

        return x;
    }


    /**
     * The extended Euclid algorithm. 
     *
     * This algorithm calculates the great common divider gcd of two integers a and b, 
     * as well as coefficents x and y such that ax + by = gcd(a, b). 
     *
     * When some of a and b is negative, the result might be negative or positive. 
     * When both of them are positive, the result is positive.
     *
     * @param[in]   a   the first operand
     * @param[in]   b   the second operand
     * @param[out]  gcd the greatest common divider of a and b
     * @param[out]  x   the first coefficient
     * @param[out]  y   the second coefficient
     */
    inline static void extendedEuclid(in int32_t a, in int32_t b, out int32_t * gcd, out int32_t * x, out int32_t * y) forceinline {

        if (!x || !y || !gcd) {
            return;
        }

        int32_t x_new = 1, y_new = 0;

        *x = 0;
        *y = 1;

        int32_t q, r, m, n;

        while (a != 0) {
            q = b / a;
            r = b % a;
            m = (*x) - q * x_new;
            n = (*y) - q * y_new;
            *x = x_new;
            *y = y_new;
            x_new = m;
            y_new = n;
            b = a;
            a = r;
        }

        // assign the values
        *gcd = b;
    }


    /**
     * This function calculates the greatest common divisor of a and b, 
     * using Euclid algorithm. 
     *
     * @param[in]   a   the first operand
     * @param[in]   b   the second operand
     *
     * @return      the greatest common divisor of a and b
     */
    inline static int32_t gcdEuclid(in int32_t a, in int32_t b) forceinline {

        int32_t gcd = 0;
        int32_t x = 0;
        int32_t y = 0;

        extendedEuclid(a, b, &gcd, &x, &y);

        return gcd;
    }

    /**
     * Calculate the following formula: a^(-1) mod n
     *
     * Note that we exploy extended Euclid algorithm to calculate the invert modulus. 
     * Function extendedEuclid() is invoked. 
     *
     * The results will always with 0 and (n-1) regardless of the sign of a
     *
     * @param[in]   a   the base
     * @param[in]   n   the modulus
     *
     * @return      the inverse modular of a over n;
     *              0 if the inverse modular doesen't exist
     */
    inline static int32_t inverseModular(in int32_t a, in int32_t n) forceinline {

        
        int32_t gcd = 0;
        int32_t x = 0;
        int32_t y = 0;

        extendedEuclid(mathModular(a, n), n, &gcd, &x, &y);

        if (gcd > 1)
            return 0;
        else
            return mathModular(x, n);
    }

    /**
     * Caculate the following formula:
     * 
     * exp(2 * pi * i * (m_v/n_v + m_h/n_h))
     *
     * Note that here we want to hoist out this function as to possible
     * speed-up, for example, look-up table approach.
     *
     * @param[in]   m_v the first operand
     * @param[in]   m_h the second operand
     *
     * @return  the resulting complex value
     */

    inline static Complex unityRootPower(in double m_v, in double m_h) forceinline {
        return comExp((2 * M_PI * (m_v / (double)NonIntSFFT::n_v + m_h / (double)NonIntSFFT::n_h)) * NonIntSFFT::I);
    }



    /**
     * This function scales up an array, i.e., 
     * for any entry to array x, multiple it by a
     *
     * @param[in,out]   x   the array to be scaled up
     * @param[in]       n   the size of the array
     * @param[in]       a   the scaling factor
     */
    template<class BaseType>
    inline static void scaleUp(in out BaseType * x, in int32_t n, in BaseType a) {
        
        int32_t i = 0;
        BaseType * x_ptr = x;
        for (; i < n; ++ x_ptr, ++ i)
            *x_ptr = a * (*x_ptr);

    }

    /**
     * This function returns the maximum of two values
     *
     * @param[in]   x   value one
     * @param[in]   y   value two
     *
     * @return  the maximum of them
     */
    template<class BaseType>
    inline static BaseType max(BaseType x, BaseType y) {
        return (x > y) ? x : y;    
    }

    /**
     * This calculate 1D discrete sinc function:
     * y = exp(-pi * 1i * (n-1) / n * x) .* sin(pi * x) ./ sin(pi * x / n) / n
     *
     * @param[in]   x   the frequency position
     * @param[in]   n   the size of the frequency domain
     *
     * @return  the sinc tail at position x
     */
    inline static Complex sincN(in double x, in uint32_t n) {

        double scalar = sin(M_PI * x) / sin(M_PI * x / n) / n;
        if (std::isnan(scalar))
            return 1;

        if (std::isinf(scalar)) {
            exit(-1);
        }

        Complex phase = comExp(-M_PI * (double)(n - 1) / (double)n * x * NonIntSFFT::I);

        return phase * scalar;
    }

    /**
     * Given non-integer position x and y, 
     * this function generates the sinc template along horizontal and vertical
     * dimension, i.e., fills Utils::horizontalSinc and Utils::verticalSinc. 
     *
     * Moreover, make sure that Utils::allocateSincTemplate is invoked before this.
     *
     * @param[in]   x   the vertical position
     * @param[in]   y   the horizontal position
     */
    inline static void generateSincTemplate(in double x, in double y) forceinline {

        if (!horizontalSinc || !verticalSinc) {
            exit(-1);
        }

        uint32_t m = 0;
        for (m = 0; m < NonIntSFFT::n_v; ++ m)
            verticalSinc[m] = sincN(m - x, NonIntSFFT::n_v);
        for (m = 0; m < NonIntSFFT::n_h; ++ m)
            horizontalSinc[m] = sincN(m - y, NonIntSFFT::n_h);

    }

    /**
     * Given non-integer frequency position x and y, 
     * this function generates the time domain sinusoid tempalte along horizontal and vertical
     * dimensions, and write them into Utils::horizontalSinc and Utils::verticalSinc
     *
     * @param[in]   x   the vetical poision
     * @param[in]   y   the horizontal position
     */
    inline static void generateSinusoidTemplate(in double x, in double y) forceinline {

        if (!horizontalSinc || !verticalSinc) {
            exit(-1);
        }
        
        // return comExp((2 * M_PI * ((double)m_v / (double)NonIntSFFT::n_v + (double)m_h / (double)NonIntSFFT::n_h)) * NonIntSFFT::I);

        uint32_t m = 0;
        for (m = 0; m < NonIntSFFT::n_v; ++ m) 
            verticalSinc[m] = comExp(2 * M_PI * m * x / (double)NonIntSFFT::n_v * NonIntSFFT::I) / sqrt(NonIntSFFT::n_v);
        for (m = 0; m < NonIntSFFT::n_h; ++ m) 
            horizontalSinc[m] = comExp(2 * M_PI * m * y / (double)NonIntSFFT::n_h * NonIntSFFT::I) / sqrt(NonIntSFFT::n_v);
            // horizontalSinc[m] = unityRootPower(0, m * y) / sqrt(NonIntSFFT::n_h);
    }

    /**
     * Given a sinc(sinusoid) template has been generated, read off the values from the template
     *
     * @param[in]   x   the vertical position
     * @param[in]   y   the horizontal position
     */
    inline static Complex readOffSincTail(in int32_t x, in int32_t y) forceinline {

        return horizontalSinc[y] * verticalSinc[x];
    }

    /**
     * The same as readOffSinusoid. 
     * It is logically different because Utils::horizontalSinc and Utils::verticalSinc
     * are re-used for sinc and sinusoid computation
     *
     * @param[in]   x   the vertical position
     * @param[in]   y   the horizontal position
     */
    inline static Complex readOffSinusoid(in int32_t x, in int32_t y) forceinline {

        return horizontalSinc[y] * verticalSinc[x];
    }

    /**
     * Compare two double values, 
     * if their absolute difference is smaller than some threshold, 
     * then say that they are equal. 
     *
     * @param[in]   op1         the first double value
     * @param[in]   op2         the second double value
     * @param[in]   threshold   the threshold
     */
    inline static bool equal(in const double & op1, in const double & op2, in const double & threshold) forceinline {

        return fabs(op1 - op2) < threshold;
    }

    /**
     * Generate a random number unformly distributed between a and b, i.e., [a, b)
     *
     * @param[in]   a   lower bound of the uniform distribution, default is 0
     * @param[in]   b   upper bound of the uniform distribution, default is 1
     *
     * @return
     *          A uniform variable
     */
    inline static double randUniform(double a = 0, double b = 1) forceinline {

        double X = (double)rand() / RAND_MAX;

        return X * (b - a) + a;
    }

    /**
     * Generate Gaussian variable using Marsaglia's method
     *
     * @param[in]   sigma   the standard deviation, default is one
     * @param[in]   mean    the mean, deafult is zero
     *
     * @return  
     *              A gaussian variable
     */
    inline static double randGaussian(double sigma = 1, double mean = 0) forceinline {

        static double V1, V2, S;
        static int phase = 0;

        double X;
        double U1, U2;

        if (phase == 0) {
            do {
                U1 = (double)rand() / RAND_MAX;
                U2 = (double)rand() / RAND_MAX;

                V1 = 2 * U1 - 1;
                V2 = 2 * U2 - 1;
                S = V1 * V1 + V2 * V2;
            } while( S >= 1 || S == 0);

            X = V1 * sqrt(-2 * log(S) / S);
        } else
            X = V2 * sqrt(-2 * log(S) / S);

        phase = 1 - phase;

        X = X * sigma + mean;

        return X;
    }

    /**
     * Given a signal, calculate its power
     *
     * @param[in]   x   the input signal
     * @param[in]   n   the length of that signal
     * @return      the power of x
     */
    inline static Real calPower(in ComplexPtr x, in const uint32_t & n) forceinline {

        ComplexPtr cur_x = x;
        uint32_t i = 0;
        Real signal_power = 0;
        Real mag = 0;

        for (; i < n; ++ i, ++ cur_x) {
            mag = comAbs(*cur_x);
            signal_power += mag * mag;
        }
	
	return signal_power;
    }


    /*
     * @param[in,out]   x   the input signal
     * @param[in]       n   the length of that signal
     * @param[in]       snr the snr of the signal
     *
     * @return
     *          the total power of noise added
     */
    inline static Real addAwgnNoise(in out ComplexPtr x, in out const uint32_t & n, in double snr) forceinline {

        Real signal_power = calPower(x, n);
        double sigma = sqrt(signal_power / n * pow(10, -snr / 10));
        
        ComplexPtr cur_x = x;
        uint32_t i = 0;

        Real noise_power = 0;
        Real noise_magnitude = 0;
        Complex noise;
        for (; i < n; ++ i, ++ cur_x) {
            noise = randGaussian(sigma / 2, 0) + randGaussian(sigma / 2, 0) * NonIntSFFT::I;
            *cur_x += noise;
            noise_magnitude = comAbs(noise);
            noise_power += noise_magnitude * noise_magnitude;
        }

        return noise_power;
    }

	inline static double median(double ** input_array, int32_t n1, int32_t n2) {

		double * array = (double *) malloc(sizeof(double) * n1 * n2);
		for (int32_t i = 0; i < n1; i ++) {
			for (int32_t j = 0; j < n2; j ++) {
				array[i * n2 + j] = input_array[i][j];
			}	
		}
	

		// copy value into vector
		// yes I admit I am lazy 
		std::vector<double> value(array, array + sizeof(array) / sizeof(double)) ;
		/*
		for (int32_t i = 0; i < n1; i ++) {
			for (int32_t j = 0; j < n2; j ++) {
				if (sample_graph[i][j]) {
					value.push_back(input_array[i][j]);
				}
			}
		}
		*/
		std::sort(value.begin(), value.end());

		int size = value.size();
		double result;
		if (size % 2 == 0) {
			result = (value[size / 2] + value[size / 2 - 1]) / 2;
		}
		else {
			result = value[(size - 1) / 2];
		}

		// deallocate space
		value.clear();
		free(array);
		array = NULL;

		return result;

	}

	inline static double max2DArray(double ** input, int32_t n1, int32_t n2, int32_t * pos_x, int32_t * pos_y) {

		double max = 0.0;
		int32_t max_x = 0;
		int32_t max_y = 0;
		for (int32_t i = 0; i < n1; i ++) {
			for (int32_t j = 0; j < n2; j ++) {
				if (input[i][j] > max) {
					max_x = i;
					max_y = j;
					max = input[i][j];
				}
			}		
		}
		
		*pos_x = max_x;
		*pos_y = max_y;
		return max;
	}	


};

#endif
