#ifndef _COMMON_H_
#define _COMMON_H_

/**
 *
 * @file    common.h
 * @brief   Common definitions that is needed
 * @author  lixin
 * @date    11/13/2012
 *
 */

#include "Eigen/Dense"
#include <complex>
#include <stdint.h>

//#define REAL_INPUT

///////////////////////////////
//for debug: 
//#define forceinline

//#define UPDATING_PSEUDOINVERSE
//#define PSEUDOINVERSE_CHECKING_MODE
#if (defined(PSEUDOINVERSE_CHECKING_MODE) && !defined(UPDATING_PSEUDOINVERSE))
    #error "PSEUDOINVERSE_CHECKING_MODE is defined but UPDATING_PSEUDOINVERSE is not"
#endif

/**
 * The maximum length of a C-style string
 */
#define MAX_STRING_LENGTH 256

/**
 * The maximum number of rows of the all_coeff matrix
 *
 * This is an important constant, because any run-time matrix
 * that is larger than this will yield an error. 
 *
 * But this cannot, at the same time, be too large, otherwise
 * we are going to run out of memory. 
 *
 * Note that this constant has some effect in the running time: 
 * more peaks means larger running time, so make it smaller
 * will force the algorithm to use smaller number of peaks
 */
#ifndef MAX_BUCKET_NUM
    #define MAX_BUCKET_NUM 320
#endif

/**
 * The maximum number of columns of the all_coeff matrix
 * and the maximum number of rows of the coeff matrix
 */
#ifndef MAX_PEAKS_NUM
    #define MAX_PEAKS_NUM 100
    #define MAX_EQUATIONS 250
    //TODO: Come back to this later
    //#define MAX_PEAKS_NUM 100
    //#define MAX_EQUATIONS 100
#endif

/**
 * The maximum size of the recovery buffer that is available
 */
#ifndef NUM_RECOVER_BUFFER
    #define NUM_RECOVER_BUFFER 5
#endif

/**
 * The flag for input parameters
 */
#define in
#define IN

/**
 * The flag for outout parameters
 */
#define out
#define OUT

/**
 * Force compiler to inline a function
 * Used in G++
 */
#ifndef forceinline
    #define forceinline __attribute__((always_inline))
#endif

/**
 * I hate unused variable warning...
 * But I don't want to disable all such warnings because
 * sometimes it is useful.
 * So I add this macro.
 */
#define GO_AWAY_UNUSED_VARIABLE(var) ((void)&(var))
#if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 402 
    #define GCC_DIAG_STR(s) #s 
    #define GCC_DIAG_JOINSTR(x,y) GCC_DIAG_STR(x ## y) 
    #define GCC_DIAG_DO_PRAGMA(x) _Pragma (#x) 
    #define GCC_DIAG_PRAGMA(x) GCC_DIAG_DO_PRAGMA(GCC diagnostic x) 
    #if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 406 
        #define GCC_DIAG_OFF(x) GCC_DIAG_PRAGMA(push) \
        GCC_DIAG_PRAGMA(ignored GCC_DIAG_JOINSTR(-W,x)) 
        #define GCC_DIAG_ON(x) GCC_DIAG_PRAGMA(pop) 
    #else 
        #define GCC_DIAG_OFF(x) GCC_DIAG_PRAGMA(ignored GCC_DIAG_JOINSTR(-W,x)) 
        #define GCC_DIAG_ON(x) GCC_DIAG_PRAGMA(warning GCC_DIAG_JOINSTR(-W,x)) 
    #endif 
#else 
    #define GCC_DIAG_OFF(x) 
    #define GCC_DIAG_ON(x) 
#endif

/**
 * The filter type that is being used in this implementation
 */
#define CurFilter MansourFilter

/**
 *
 * @brief The namespace for common variable definitions.
 *
 * This is the major namespace that we are using for 
 * defining enums, variables, memory allocations of
 * Non-Integer Sparse FFT algorithm.
 *
 * @author  lixin
 * @date    11/13/2012
 *
 */
namespace NonIntSFFT {

    /**
     * The error code definition
     */
    enum ERROR {

        NULL_POINTER,           /**< Null pointer error. */
        UNIMPLEMENTED,          /**< Unimplemented function error */
        SMALLER_THAN_ZERO,      /**< Input smaller than zero while expected positive */
        NONSQUARE_MAT,          /**< Non-square domain size while square domain size is expected */
        WRONG_ENUM_NUMBER,      /**< When we come into a wrong enumerator number */
        MATRIX_TOO_LARGE,       /**< The error is when the size of matrix is larger than MAX_MATRIX_SIZE */
        INFINITY_VALUE,         /**< Get an infinity value */
        NUM_PEAKS_NOMATCH,      /**< If the number of peaks doesn't match */
        INTERNAL_ERROR,         /**< Any error that is logically wrong; or indicate a clear bug in the code */
        CANNOT_OPEN_FILE,       /**< Cannot open file for io */
        IO_ERROR,               /**< IO error */
        DIM_UNSUPPORTED,        /**< The dimension is not supported */
        INVALID_PARAL_ID,       /**< Invalid parallel computation ID */
        INVALID_PARAM,          /**< Invalid parameter */
        FORMAT_ERR,             /**< Format error while parsing param file */
        NEED_TO_RECALCULATE,    /**< Need to recalculate the pseudo-inverse */
        TOTAL_NUM_ERRORS        /**< Total number of error codes. */
    };

    /**
     * The names of all the errors
     */
    extern const char ERROR_NAMES[TOTAL_NUM_ERRORS][MAX_STRING_LENGTH];

    /**
     * The corresponding definition for real/imaginary part of 
     * the complex_t definition.
     */
    typedef double real_t;

    /**
     * Definition of the complex class goes here.
     * Note that we use the STL template class complex. 
     * There are two reasons that we use STL instead of GLU complex
     *
     * 1) The eigen library we are using for matrix compuation
     *    doesn't support GLU complex
     *
     * 2) The GLU complex is kind of deprecated
     *
     * 3) In terms of efficiency, these two classes are comparable
     */
    typedef std::complex<real_t> complex_t;

    /**
     * Definition of an array of complex values
     */
    typedef complex_t * complex_ptr;

    /**
     * The vertical size of the 2D spectrum, 
     */
    extern uint32_t n_v;

    /**
     * The horizontal size of the 2D spectrum, 
     */
    extern uint32_t n_h;

    /**
     * The total size; iterally n_v * n_h
     */
    extern uint32_t n;

    /**
     * The imaginary constant
     */
    const complex_t I(0, 1);

    /**
     * The error state
     * 1 means a fatal error occurs
     */
    extern bool error;
}

/** 
 * This definition is just for convenience for denoting complex type
 */
#define Complex NonIntSFFT::complex_t

/** 
 * This definition is just for convenience for denoting real type
 */
#define Real NonIntSFFT::real_t

/**
 * This definition is just for convenience for denoting complex array type
 */
#define ComplexPtr NonIntSFFT::complex_ptr

/**
 * The Position<int32_t> class
 */
#define IntPosition Position<int32_t>

/**
 * The non-integer position class Position<double>
 */
#define NonIntPosition Position<double>

/**
 * The non-integer positioned class Peaks<double>
 */
#define NonIntPeak Peaks<double>

/**
 * The integer-positioned class Peaks<int32_t>
 */
#define IntPeak Peaks<int32_t>

/**
 * The list of non-integer peaks
 */
#define PeaksList std::list<NonIntPeak>

/**
 * The iterator to the peaks list
 */
#define PeakIterator std::list<NonIntPeak>::iterator


/**
 * The Peaks window used to store recent calculation of Peaks
 */
//#define PeaksWindow std::list< std::list<NonIntPeak> >

/**
 * The iterator to the PeaksWindow
 */
//#define PeaksWindowIterator std::list< std::list<NonIntPeak> >::iterator


/**
 * Complex exponential
 */
#define comExp std::exp<Real>

/**
 * Complex magnitude
 */
#define comAbs std::abs<Real>

/**
 * Complex conjugate
 */
#define comConj std::conj<Real>

/**
 * The matrix type that is used
 */
#define AllCoeffMat Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, MAX_BUCKET_NUM, MAX_PEAKS_NUM>

#define CoeffMat Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, MAX_EQUATIONS, MAX_PEAKS_NUM>
#define CoeffMatRowMajor Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor, MAX_EQUATIONS, MAX_PEAKS_NUM>

/**
 * The vector type that is used
 */
#define AllBukVec Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, MAX_BUCKET_NUM, 1>
#define BukVec Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, MAX_EQUATIONS, 1>
#define CoeffVec Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, MAX_PEAKS_NUM, 1>
/**
 * What method we are going to use to decompose the matrix
 * There are two options, SVD and QR decomposition
 * If SVD is defined, then we use SVD decomposition; 
 * otherwise QR is used
 */
//#define SVD

/**
 * The decomposition method that is used to solve linear system. 
 * The alternatives are SVD decomposition and QR decomposition
 */
#ifdef SVD
    #define decompose jacobiSvd
    #define Decompose Eigen::JacobiSVD<CoeffMat>
    #define rank nonzeroSingularValues
    #define decomp_args Eigen::ComputeThinU | Eigen::ComputeThinV
#else
    #define decompose colPivHouseholderQr
    #define Decompose Eigen::ColPivHouseholderQR<CoeffMat>
    #define rank rank
    #define decomp_args 
#endif

#ifndef DEBUG_ON
    GCC_DIAG_OFF(unused-parameter)
#endif

#endif
