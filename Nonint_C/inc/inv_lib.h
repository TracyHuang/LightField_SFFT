/*******************************************************************************
  Copyright [2013] <Zhunping>

  This file implements the rank-1 matrix pseudoinverse update
  algorithm in :
  Carl D. Meyer. Generalized inversion of modified matrices.
  SIAM Journal of Applied Mathematics, 24(3):315–323, May 1973

Author: Zhunping Zhang
 *******************************************************************************/

#ifndef _INT_INC_INVLIB_H_
#define _INT_INC_INVLIB_H_
#include <stdio.h>
//#include <sys/time.h>
#include <Eigen/Dense>
#include <complex>
#include <stdlib.h>
#include "common.h"
#include <stdarg.h>

// *****************************************************************************
// Global constants
#define TRUE                  1
#define FALSE                 0

//#define INV_MAXSIZE           170    // maximal matrix size, can't be too big,
// otherwise segfault
//#define MAXSIZE               INV_MAXSIZE
//#define MAXIMALROWS           MAXSIZE
//#define MAXIMALCOLS           MAXSIZE

#define INV_ROUNDOFF_TOLE     1e-12  // error tolerance for rounding off error
#define TOLE                  1e-6   // error tolerance used for deciding cases
// in Meyer's algorithm
#define EPS_ERROR_CHECKING    1e-3   // error tolerance used for error checking
#define TOLE_PINV             1e-4   // error tolerance used for svd pinv's
// matrix inverting

#define U_THRESHOLD 1e-8
#define V_THRESHOLD 1e-2
#define BETA_THRESHOLD 1e-6

#undef INV_DEBUG
// #define   INV_DEBUG
// #define DEBUG

#define MAX_NUM_ROWS MAX_EQUATIONS
#define MAX_NUM_COLS MAX_PEAKS_NUM
#define MAX_VEC_LENGTH ((MAX_EQUATIONS > MAX_PEAKS_NUM) ? MAX_EQUATIONS : MAX_PEAKS_NUM)


// *****************************************************************************
// Global macros
#define INLINE                inline __attribute__((always_inline))

/**
 * complex number macro
 * real part of (a+bi)*(c+di)
 */
#define cplxmr(a, b, c, d)        ((a)*(c)-(b)*(d))
/**
 * complex number macro
 * imaginary part of (a+bi)*(c+di)
 */
#define cplxmi(a, b, c, d)        ((a)*(d)+(b)*(c))
/**
 * complex number macro
 * squared norm of (a+bi)
 */
#define cplxn2(a, b)              ((a)*(a)+(b)*(b))
/**
 * complex number macro
 * squared distance between a+bi and c+di
 */
#define cplxdis2(a, b, c, d)      (cplxn2((a)-(c), (b)-(d)))
/**
 * complex number macro
 * real part of (a+bi)/(c+di), dn=cc+dd
 */
#define cplxdr(a, b, c, d, dn)    (((a)*(c)+(b)*(d))/(dn))
/**
 * complex number macro
 * imaginary part of (a+bi)/(c+di), dn=cc+dd
 */
#define cplxdi(a, b, c, d, dn)    (((b)*(c)-(a)*(d))/(dn))


// *****************************************************************************
// Inversion library

struct InvCaseStatistics {
    uint32_t num_total_ops;
    uint32_t num_case5;

    inline void reset() forceinline {
        num_total_ops = 0;
        num_case5 = 0;
    }

    inline void newCase(IN const bool & _need2recalculate) forceinline {
        num_total_ops ++;
        num_case5 += (_need2recalculate ? 0 : 1);
    }

    inline void printItself(IN FILE * fout) const forceinline {
        fprintf(fout, "%4d, %4d, %05.2lf%%", num_total_ops, num_case5, 100.0 * (double)num_case5 / (double)num_total_ops);
    }
};

struct InvOpStatistics {

    struct InvCaseStatistics update_cols;
    struct InvCaseStatistics update_rows;
    struct InvCaseStatistics add_cols;
    struct InvCaseStatistics add_rows;
    struct InvCaseStatistics remove_cols;
    struct InvCaseStatistics remove_rows;
    struct InvCaseStatistics overall;

    inline void reset() forceinline {
        update_cols.reset();
        update_rows.reset();
        add_cols.reset();
        add_rows.reset();
        remove_cols.reset();
        remove_rows.reset();
        overall.reset();
    }

    inline void printItself(IN FILE * fout, IN const uint32_t & indent = 0) const forceinline {
        
        fprintf(fout, "%*s%10s: ", indent, "", "Update Col");
        update_cols.printItself(fout);
        fprintf(fout, "\n");
        
        fprintf(fout, "%*s%10s: ", indent, "", "Update Row");
        update_rows.printItself(fout);
        fprintf(fout, "\n");
        
        fprintf(fout, "%*s%10s: ", indent, "", "Add Col");
        add_cols.printItself(fout);
        fprintf(fout, "\n");
        
        fprintf(fout, "%*s%10s: ", indent, "", "Add Row");
        add_rows.printItself(fout);
        fprintf(fout, "\n");
        
        fprintf(fout, "%*s%10s: ", indent, "", "Remove Col");
        remove_cols.printItself(fout);
        fprintf(fout, "\n");
        
        fprintf(fout, "%*s%10s: ", indent, "", "Remove Row");
        remove_rows.printItself(fout);
        fprintf(fout, "\n");
        
        fprintf(fout, "%*s%10s: ", indent, "", "Overall");
        overall.printItself(fout);
        fprintf(fout, "\n");

    }

};

enum MATRIX_OPS {
    UPDATE_COL, 
    UPDATE_ROW, 
    ADD_COL, 
    ADD_ROW, 
    REMOVE_COL, 
    REMOVE_ROW, 
};

class InvLib {

protected:
    /**
     * Number of rows
     */
    uint32_t num_rows;
    /**
     * Number of columns
     */
    uint32_t num_cols;
    /**
     * The matrix A to be maintained
     */

    CoeffMatRowMajor matrix_mat;
    //CoeffMatRowMajor pinv_mat;

    double * matrix;
    double * pinv;                         // inversion of matrix
    double * new_pinv;                        // inversion of new matrix
    //uint32_t seed;  // thread safe random number generator's seed
    //double  m_modelelement;

    uint32_t need2recalculate;

protected:
    // temporary variables
    double *   tmp_matrix;                       // temporary matrix
    double *   c;
    double *   d;

    Complex zero_row_vec[MAX_NUM_ROWS * 2];  // temporary complex vector
    Complex zero_col_vec[MAX_NUM_COLS * 2];  // temporary complex vector

    double singular_vec[MAX_NUM_COLS];

    double k[MAX_NUM_COLS * 2];
    double h[MAX_NUM_ROWS * 2];
    double u[MAX_NUM_ROWS * 2];
    double v[MAX_NUM_COLS * 2];
    Complex beta;

    double k2;
    double h2;
    double u2;
    double v2;
    double beta2;

    double pinv_mul_h[MAX_NUM_COLS * 2];
    double q[MAX_NUM_ROWS * 2]; 

    Complex linear_solution[MAX_NUM_COLS];
    Complex tmp_row_vec[MAX_NUM_ROWS];

    InvOpStatistics * statistics;

    Eigen::JacobiSVD<CoeffMatRowMajor> jacobi_svd;

protected:
    inline double * allocMat() forceinline {
        return (double *)malloc(sizeof(double) * MAX_NUM_COLS * MAX_NUM_ROWS * 2);
    }

    // alloc a vec
    inline double * allocRowVec() forceinline { 
        return (double *)malloc(sizeof(double) * MAX_NUM_COLS * 2);
    }

    inline double * allocColVec() forceinline {
        return (double *)malloc(sizeof(double) * MAX_NUM_ROWS * 2);
    }

    // swap two mat_t* pointers
    inline void swapPointers(IN OUT double ** a, IN OUT double ** b) forceinline {
        double * t;
        t  = *a;
        *a = *b;
        *b = t;
    }
    
    inline void calculatePinv(OUT double * pinv_out) forceinline {

        matrix_mat.resize(num_rows, num_cols); 
        memcpy(matrix_mat.data(), matrix, sizeof(Complex) * num_rows * num_cols);
        jacobi_svd.compute(matrix_mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
        const Eigen::JacobiSVD<CoeffMatRowMajor>::SingularValuesType & singular_values = jacobi_svd.singularValues();
        
		// Get the new singular values
        //uint32_t rank = 0;
        for (uint32_t i = 0; i < num_cols; ++ i) {
            if (singular_values(i) > TOLE_PINV) {
                singular_vec[i] = 1.0 / singular_values(i);
                //rank ++;
            } else
                singular_vec[i] = 0.0;
        }

        // This is O(n^3) complexity matrix multiplication

        ComplexPtr pinv_ptr = (ComplexPtr)pinv_out;
        const Complex * matrix_v = jacobi_svd.matrixV().data();
        const Complex * matrix_u = jacobi_svd.matrixU().data();
        const Complex * matu_ptr = matrix_u;
        const Complex * matv_ptr = matrix_v;
	
        for (uint32_t i = 0; i < num_cols; ++ i) {
            matu_ptr = matrix_u;

            for (uint32_t j = 0; j < num_rows; ++ j) {

                *pinv_ptr = 0;
                matv_ptr = matrix_v;
                for (uint32_t k = 0; k < num_cols; ++ k, ++ matv_ptr, ++ matu_ptr) {
                    *pinv_ptr += (*matv_ptr) * comConj(*matu_ptr) * singular_vec[k];
                }

                pinv_ptr ++;
            }
            matrix_v += num_cols;
        }

    }

    // Computer k, h, u, v, beta in Meyer's paper
    // _k2 is the squared norm of k
    inline static void computeKHUVBeta(IN const uint32_t & num_rows, IN const uint32_t num_cols, 
        IN const double * matrix, IN const double * pinv, IN const double * c, IN const double * d,
        OUT double * k, OUT double * h, OUT double * u, OUT double * v, OUT Complex & beta, 
        OUT double & k2, OUT double & h2, OUT double & u2, OUT double & v2, OUT double & beta2) forceinline {

        const double * c_ptr = NULL;
        const double * d_ptr = NULL;
        const double * pinv_ptr = NULL;
        const double * matrix_ptr = NULL;
        double tmp_real = 0;
        double tmp_imag = 0;

        // compute k(num_cols,1) = pinv(num_cols,num_rows) * c(num_rows,1): k(i) = \sum_j pinv(i, j) * c(j)
        pinv_ptr = pinv;
        double * k_ptr = k;
        k2 = 0;
        for (uint32_t i = 0; i < num_cols; i++) {
            c_ptr = c;
            tmp_real = tmp_imag = 0;
            for (uint32_t j = 0; j < num_rows; j++) {
                tmp_real += cplxmr(*pinv_ptr, *(pinv_ptr+1), *c_ptr, *(c_ptr+1));
                tmp_imag += cplxmi(*pinv_ptr, *(pinv_ptr+1), *c_ptr, *(c_ptr+1));
                pinv_ptr += 2;
                c_ptr += 2;
            }
            *(k_ptr++) = tmp_real;
            *(k_ptr++) = tmp_imag;
            k2 += cplxn2(tmp_real, tmp_imag);
        }

        // compute h(1,num_rows) = d(num_cols,1).adjoint() * pinv(num_cols,num_rows):
        //                                 h(j) = \sum_i conj(d(i)) * pinv(i, j)
        pinv_ptr = pinv;
        d_ptr = d;
        memset(h, 0, sizeof(double) * num_rows * 2);
        double * h_ptr = NULL;
        for (uint32_t i = 0; i < num_cols; i++) {
            h_ptr = h;
            for (uint32_t j = 0; j < num_rows; j++) {
                *(h_ptr++) += cplxmr(*d_ptr, -*(d_ptr+1), *pinv_ptr, *(pinv_ptr+1));
                *(h_ptr++) += cplxmi(*d_ptr, -*(d_ptr+1), *pinv_ptr, *(pinv_ptr+1));
                pinv_ptr += 2;
            }
            d_ptr += 2;
        }
        h_ptr = h;
        h2 = 0;
        for (uint32_t j = 0; j < num_rows; j++) {
            h2 += cplxn2(*h_ptr, *(h_ptr+1));
            h_ptr += 2;
        }

        // compute u(num_rows,1)= c(num_rows,1) - A(num_rows,num_cols) * k(num_cols,1):
        //                                 u(i) = c(i) - \sum_j A(i, j) * k(j)
        c_ptr = c;
        matrix_ptr = matrix;
        double * u_ptr = u;
        u2 = 0;
        for (uint32_t i = 0; i < num_rows; i++) {
            k_ptr = k;
            tmp_real = *(c_ptr++);
            tmp_imag = *(c_ptr++);
            for (uint32_t j = 0; j < num_cols; j++) {
                tmp_real -= cplxmr(*matrix_ptr, *(matrix_ptr+1), *k_ptr, *(k_ptr+1));
                tmp_imag -= cplxmi(*matrix_ptr, *(matrix_ptr+1), *k_ptr, *(k_ptr+1));
                matrix_ptr += 2;
                k_ptr += 2;
            }
            *(u_ptr++) = tmp_real;
            *(u_ptr++) = tmp_imag;
            u2 += cplxn2(tmp_real, tmp_imag);
        }

        // compute v(1,c) = d(num_cols,1).adjoint() - h(1,num_rows) * A(num_rows,num_cols):
        //                              v(j) = conj(d(j)) - \sum_i h(i)*A(i, j)
        d_ptr = d;
        double * v_ptr = v;
        for (uint32_t i = 0; i < num_cols; i++) {
            *(v_ptr++) = *(d_ptr++);
            *(v_ptr++) = -*(d_ptr++);
        }
        matrix_ptr = matrix;
        h_ptr = h;
        for (uint32_t i = 0; i < num_rows; i++) {
            v_ptr = v;
            for (uint32_t j = 0; j < num_cols; j++) {
                *(v_ptr++) -= cplxmr(*h_ptr, *(h_ptr+1), *matrix_ptr, *(matrix_ptr+1));
                *(v_ptr++) -= cplxmi(*h_ptr, *(h_ptr+1), *matrix_ptr, *(matrix_ptr+1));
                matrix_ptr += 2;
            }
            h_ptr += 2;
        }
        v_ptr = v;
        v2 = 0;
        for (uint32_t i = 0; i < num_cols; i++) {
            v2 += cplxn2(*v_ptr, *(v_ptr+1));
            v_ptr += 2;
        }

        // Compute Beta = 1 + h(1,num_rows) * c(num_rows,1): b = 1 + \sum_i h(i)*c(i)
        double & br = ((double *)&beta)[0];
        double & bi = ((double *)&beta)[1];
        br = 1; 
        bi = 0;
        h_ptr = h;
        c_ptr = c;
        for (uint32_t i = 0; i < num_rows; i++) {
            br += cplxmr(*h_ptr, *(h_ptr+1), *c_ptr, *(c_ptr+1));
            bi += cplxmi(*h_ptr, *(h_ptr+1), *c_ptr, *(c_ptr+1));
            h_ptr += 2;
            c_ptr += 2;
        }
        beta2 = cplxn2(br, bi);
    }


    // compute pinv_mul_h(num_cols,1) = pinv(num_cols,num_rows) * h(1,num_rows).adjoint()
    inline static void computePinvMulH(IN const uint32_t & num_rows, IN const uint32_t & num_cols, 
        IN const double * pinv, IN const double * h, OUT double * pinv_mul_h) forceinline {

        double * pinv_mul_h_ptr = pinv_mul_h;
        const double * pinv_ptr = pinv;
        const double * h_ptr = h;
        double tmp_real = 0;
        double tmp_imag = 0;

        for (uint32_t i = 0; i < num_cols; i++) {
            h_ptr = h;
            tmp_real = tmp_imag = 0;
            for (uint32_t j = 0; j < num_rows; j++) {
                tmp_real += cplxmr(*pinv_ptr, *(pinv_ptr+1), *h_ptr, -*(h_ptr+1));
                tmp_imag += cplxmi(*pinv_ptr, *(pinv_ptr+1), *h_ptr, -*(h_ptr+1));
                h_ptr += 2;
                pinv_ptr += 2;
            }
            *(pinv_mul_h_ptr++) = tmp_real;
            *(pinv_mul_h_ptr++) = tmp_imag;
        }
    }

    // update A = A + c * d^T and new_pinv = pinv(A)
    // The "special" means that we are not using the original Meyer method
    // Instead we only use case 5 and give the rest of the case to SVD method
    // We use case 5 only because statistically 99% of the data we use is
    // case 5
    // return if we need to recalculate the pseudo-inverse
    // Note: new_pinv is changed, not pinv
    inline bool updateInverse(/*int num_rows, int num_cols, double * matrix, double * pinv, double * c, double * d,
                double * new_pinv*/) forceinline {

        computeKHUVBeta(num_rows, num_cols, matrix, pinv, c, d, k, h, u, v, beta, k2, h2, u2, v2, beta2);
        bool need2recalculate = (v2 > V_THRESHOLD) || (beta2 <= BETA_THRESHOLD);


        // if it is not case 5
        //if (u2 <= U_THRESHOLD || v2 > V_THRESHOLD || beta2 <= BETA_THRESHOLD) {
        if (need2recalculate) {
            //fprintf(stderr, "u2 = %lf, v2 = %lf, beta2 = %lf\n", u2, v2, beta2);
            return true;
        }

        //m_ncase5s++;
        computePinvMulH(num_rows, num_cols, pinv, h, pinv_mul_h);
        const double & br = ((double *)&beta)[0];
        const double & bi = ((double *)&beta)[1];

        // c0 = 1/conj(beta)
        double c0r = br / beta2;
        double c0i = bi / beta2;
        // c1 = conj(beta) / sigma2, where sigma2 = ||h||^2 ||u||^2 + |beta|^2
        double c1r = h2 * u2 + beta2;
        double c1i = -bi / c1r;
        c1r = br / c1r;
        // c2 = -||u||^2 / conj(beta)
        double c2r = -u2 * c0r;
        double c2i = -u2 * c0i;
        // c3 = -||h||^2 / conj(beta)
        double c3r = -h2 * c0r;
        double c3i = -h2 * c0i;

        const double * u_ptr = u;
        const double * h_ptr = h;
        double * q_ptr = q;
        for (uint32_t i = 0; i < num_rows; i++) {
            *(q_ptr++) = cplxmr(c3r, c3i, *u_ptr, -*(u_ptr+1)) - *(h_ptr++);
            *(q_ptr++) = cplxmi(c3r, c3i, *u_ptr, -*(u_ptr+1)) - *(h_ptr++);
            u_ptr += 2;
        }

        // x = pinv + c0 * pinv_mul_h * u' - c1 * (c2 * pinv_mul_h - k) * q
        double * new_pinv_ptr = new_pinv;
        const double * pinv_ptr = pinv;
        const double * pinv_mul_h_ptr = pinv_mul_h;
        const double * k_ptr = k;

        double c4r = 0;
        double c4i = 0;
        double c5r = 0;
        double c5i = 0;
        double c6r = 0;
        double c6i = 0;
        for (uint32_t i = 0; i < num_cols; i++) {
            // c5 = c0 * pinv_mul_h
            c5r = cplxmr(c0r, c0i, *pinv_mul_h_ptr, *(pinv_mul_h_ptr+1));
            c5i = cplxmi(c0r, c0i, *pinv_mul_h_ptr, *(pinv_mul_h_ptr+1));

            // c6 = c1 *(c2 * pinv_mul_h - k)
            c4r = cplxmr(c2r, c2i, *pinv_mul_h_ptr, *(pinv_mul_h_ptr+1)) - *(k_ptr++);
            c4i = cplxmi(c2r, c2i, *pinv_mul_h_ptr, *(pinv_mul_h_ptr+1)) - *(k_ptr++);
            c6r = cplxmr(c1r, c1i, c4r, c4i);
            c6i = cplxmi(c1r, c1i, c4r, c4i);

            u_ptr = u;
            q_ptr = q;
            for (uint32_t j = 0; j < num_rows; j++) {
                // x = pinv + c5 * u' - c6 * c3ujmh
                *(new_pinv_ptr++) = *(pinv_ptr++)
                    + cplxmr(c5r, c5i, *u_ptr, -*(u_ptr+1))
                    - cplxmr(c6r, c6i, *q_ptr, *(q_ptr+1));
                *(new_pinv_ptr++) = *(pinv_ptr++)
                    + cplxmi(c5r, c5i, *u_ptr, -*(u_ptr+1))
                    - cplxmi(c6r, c6i, *q_ptr, *(q_ptr+1));

                u_ptr += 2;
                q_ptr += 2;
            }
            pinv_mul_h_ptr += 2;
        }

        return false;
    }



public:

    InvLib() : num_rows(0), num_cols(0), statistics(NULL) {

        matrix = allocMat();
        pinv = allocMat();

        new_pinv = allocMat();
        tmp_matrix = allocMat();
        c = allocColVec();
        d = allocRowVec();

        memset(zero_row_vec, 0, sizeof(Complex) * MAX_NUM_COLS);
        memset(zero_col_vec, 0, sizeof(Complex) * MAX_NUM_ROWS);
    }

    ~InvLib() {
        free(new_pinv);
        free(tmp_matrix);
        free(c);
        free(d);
    }

    inline void setStatistics(IN InvOpStatistics * _statisitcs) {
        statistics = _statisitcs;
    }

    inline void copyFrom(IN const InvLib & source) forceinline {
        num_rows = source.num_rows;
        num_cols = source.num_cols;
        memcpy(matrix, source.matrix, sizeof(Complex) * MAX_NUM_ROWS * MAX_NUM_COLS);
        memcpy(pinv, source.pinv, sizeof(Complex) * MAX_NUM_ROWS * MAX_NUM_COLS);
    }

    inline const Complex * getSolution() const forceinline {
        return linear_solution;
    }

    inline const Complex * getMatrix() const forceinline {
        return (Complex *)matrix; 
    }

    inline const uint32_t & getNumRows() const forceinline {
        return num_rows;
    }

    inline const uint32_t & getNumCols() const forceinline {
        return num_cols;
    }

    inline void resize(IN const uint32_t & _num_rows, IN const uint32_t & _num_cols) forceinline {
        num_rows = _num_rows;
        num_cols = _num_cols;
    }

    inline Complex * getMatrix() forceinline {
        return (Complex *)matrix;
    }


    // multiple a matrix A[r][co] with a vector v[co] into a vector re[r], 
    inline static void matVecMul(IN const uint32_t & num_rows, IN const uint32_t & num_cols, 
        IN const double * matrix, IN const double * vec, OUT double * result) forceinline {
        // compute re[r] = A[r][co] * vec[co]    re(i) = \sum_j A(i,j) * vec(j)
        const double * matrix_ptr = matrix;
        double * result_ptr = result;
        const double * vec_ptr = NULL;
        double real = 0;
        double imag = 0;
        //nm = 0;
        for (uint32_t i = 0; i < num_rows; i++) {
            vec_ptr = vec;
            real = imag = 0;
            for (uint32_t j = 0; j < num_cols; j++, vec_ptr += 2, matrix_ptr += 2) {
                real += cplxmr(*matrix_ptr, *(matrix_ptr+1), *vec_ptr, *(vec_ptr+1));
                imag += cplxmi(*matrix_ptr, *(matrix_ptr+1), *vec_ptr, *(vec_ptr+1));
                //matrix_ptr += 2;
                //vec_ptr += 2;
            }
            *(result_ptr++) = real;
            *(result_ptr++) = imag;
            //nm += cplxn2(real, imag);   // compute squared norm
        }
        //return nm;
    }


    // set A
    void setMatrix(IN const CoeffMat & cm) {

        num_rows = cm.rows();
        num_cols = cm.cols();

        double * matrix_ptr = matrix;
        for (uint32_t i = 0; i < num_rows; i++)
            for (uint32_t j = 0; j < num_cols; j++) {
                *(matrix_ptr++) = cm(i, j).real();
                *(matrix_ptr++) = cm(i, j).imag();
            }
    }


    // check if A is the same as the given
    inline bool isSameMatrix(IN const CoeffMat& cm) forceinline {

        if (num_rows != cm.rows() || num_cols != cm.cols()) 
            return false;

        double * matrix_ptr = matrix;
        for (uint32_t i = 0; i < num_rows; i++) {
            for (uint32_t j = 0; j < num_cols; j++) {
                if (fabs(*(matrix_ptr++)-cm(i, j).real()) > INV_ROUNDOFF_TOLE)
                    return false;
                if (fabs(*(matrix_ptr++)-cm(i, j).imag()) > INV_ROUNDOFF_TOLE)
                    return false;
            }
        }

        return true;
    }

    inline bool checkSolution(IN const Complex * their_solution, IN const Complex * b) forceinline {

        double average_error = 0;
        double our_residue = 0;
        double their_residue = 0;

        // calculate the average error
        const Complex * my_result = linear_solution;
        const Complex * their_result = their_solution;
        for (uint32_t i = 0; i < num_cols; ++ i, ++ my_result, ++ their_result)
            average_error += comAbs(*my_result - *their_result);

        average_error /= num_cols;

        // calculate our residue
        memcpy(tmp_row_vec, b, sizeof(Complex) * num_rows);
        Complex * cur_check = tmp_row_vec;
        const Complex * matrix_ptr = (const Complex *)matrix;
        for (uint32_t i = 0; i < num_rows; ++ i, ++ cur_check) {
            my_result = linear_solution;
            for (uint32_t j = 0; j < num_cols; ++ j, ++ my_result, ++ matrix_ptr)
                *cur_check -= (*matrix_ptr) * (*my_result);
        }
        cur_check = tmp_row_vec;
        for (uint32_t i = 0; i < num_rows; ++ i, ++ cur_check)
            our_residue += cur_check->real() * cur_check->real() + cur_check->imag() * cur_check->imag();

        // calculate their residue
        matrix_ptr = (const Complex *)matrix;
        memcpy(tmp_row_vec, b, sizeof(Complex) * num_rows);
        cur_check = tmp_row_vec;
        for (uint32_t i = 0; i < num_rows; ++ i, ++ cur_check) {
            their_result = their_solution;
            for (uint32_t j = 0; j < num_cols; ++ j, ++ their_result, ++ matrix_ptr)
                *cur_check -= (*matrix_ptr) * (*their_result);
        }
        cur_check = tmp_row_vec;
        for (uint32_t i = 0; i < num_rows; ++ i, ++ cur_check)
            their_residue += cur_check->real() * cur_check->real() + cur_check->imag() * cur_check->imag();
        
        // then let's decide if it is the same or better
        if (average_error < EPS_ERROR_CHECKING || our_residue < their_residue) {
            //fprintf(stderr, "different reconstruction: average difference = %g, "
                //"my residue = %g, their residue = %g\n", average_error, our_residue, their_residue);
            return true;
        } else {
            // we have a different result and it has a larger residue
            fprintf(stderr, "different reconstruction: average difference = %g, "
                "my residue = %g, their residue = %g", average_error, our_residue, their_residue);
            return false;
        }
    }


    // initialize Ai to be pinv(A) using SVD
    inline void initPinv() forceinline {
        calculatePinv(pinv);
    }

    inline Complex * getColVec() forceinline {
        return (Complex *)c;
    }

    // update a col[nele] at location upco in A[r][co]
    // side effect: matrix, m_An, pinv, new_pinv are changed
    inline bool updateCol(IN const uint32_t & col_index, IN const uint32_t & num_elements, IN const ComplexPtr & col = NULL) forceinline {

        //const double * cur_col = (const double *)col;
        if (col != NULL && col != (Complex *)c)
            memcpy(c, col, sizeof(Complex) * num_elements);
        // Now A is row-major
        double * matrix_ptr = matrix + col_index * 2;

        double * c_ptr = c;

        for (uint32_t i = 0; i < num_rows; i ++) {
            // c = col - A(:,col_index), so that A(:,col_index) + c gives col
            *(c_ptr ++) -= *matrix_ptr;
            *(c_ptr ++) -= *(matrix_ptr + 1);
            matrix_ptr += num_cols * 2;
        }
        memset(d, 0, sizeof(d[0]) * num_cols * 2);
        d[col_index * 2] = 1;

        need2recalculate = updateInverse(); //num_rows, num_cols, matrix, pinv, c, d, new_pinv);
        // update the matrix
        matrix_ptr = matrix + col_index * 2;
        c_ptr = c;
        for (uint32_t i = 0; i < num_rows; ++ i, c_ptr += 2, matrix_ptr += 2 * num_cols) {
            //*matrix_ptr = *cur_col;
            //*(matrix_ptr + 1) = *(cur_col + 1);
            *matrix_ptr += *c_ptr;
            *(matrix_ptr + 1) += *(c_ptr + 1);
        }
        if (need2recalculate) {
            calculatePinv(new_pinv);
        }
        swapPointers(reinterpret_cast<double **>(&pinv),
                reinterpret_cast<double **>(&new_pinv));
        return need2recalculate;
    }

    inline Complex * getRowVec() forceinline {
        return (Complex *)d;
    }

    // update a row[nele] at location upr in A[r][co]
    // side effect: matrix, m_An, pinv, new_pinv are changed
    inline bool updateRow(IN const uint32_t & row_index, IN const uint32_t & num_elements, IN const ComplexPtr & row = NULL) forceinline {
        // then we compute c and d

        //d_ = d;
        //const double * cur_row = (const double *)row;
        if (row != NULL && row != (Complex *)d)
            memcpy(d, row, sizeof(Complex) * num_elements);
        //else
            //cur_row = d;

        double * d_ptr = d;
        
        double * matrix_ptr = matrix + row_index * num_cols * 2;
        for (uint32_t i = 0; i < num_cols; i ++) {
            // d = (row - A(row_index,:)).adjoint(),
            // so that A(row_index,:) + d.adjoint() gives row
            *(d_ptr ++) -= *(matrix_ptr ++);
            *d_ptr= *(matrix_ptr ++) - *d_ptr;
            d_ptr ++;
        }

        memset(c, 0, sizeof(c[0]) * num_rows * 2);
        c[row_index * 2] = 1;

        need2recalculate = updateInverse(); //num_rows, num_cols, matrix, pinv, c, d, new_pinv);

        // Now update the matrix
        matrix_ptr = matrix + row_index * num_cols * 2;
        d_ptr = d;
        for (uint32_t i = 0; i < num_cols; ++ i, d_ptr += 2, matrix_ptr += 2) {
            *matrix_ptr += *d_ptr;
            *(matrix_ptr + 1) = *(matrix_ptr + 1) - *(d_ptr + 1);
        }

        if (need2recalculate) {
            calculatePinv(new_pinv);
        }

        swapPointers(reinterpret_cast<double **>(&pinv),
                reinterpret_cast<double **>(&new_pinv));
       
        return need2recalculate;
    }



    // add a col[nele] at location adc to A[r][co]
    inline bool addCol(IN const uint32_t & col_index, IN const uint32_t & num_elements, IN const ComplexPtr & col = NULL) forceinline {

        //matrix_mat.resize(Eigen::NoChange, num_cols + 1);
        //pinv_mat.resize(num_cols + 1, Eigen::NoChange);

        // add a col of 0s to matrix
        double * matrix_ptr = matrix;
        double * tmp_matrix_ptr = tmp_matrix;
        for (uint32_t i = 0; i < num_rows; i++) {
            memcpy(tmp_matrix_ptr, matrix_ptr, sizeof(tmp_matrix[0])*col_index*2);
            matrix_ptr += col_index * 2;
            tmp_matrix_ptr += col_index * 2;
            *(tmp_matrix_ptr++) = 0;
            *(tmp_matrix_ptr++) = 0;
            memcpy(tmp_matrix_ptr, matrix_ptr, sizeof(tmp_matrix[0])*(num_cols - col_index)*2);
            matrix_ptr += (num_cols - col_index)*2;
            tmp_matrix_ptr += (num_cols - col_index)*2;
        }
        swapPointers(reinterpret_cast<double **>(&matrix),
                reinterpret_cast<double **>(&tmp_matrix));

        // add a row of 0s to pinv
        double * pinv_ptr = pinv + num_rows * num_cols * 2;
        for (uint32_t i = num_cols; i > col_index; i--) {
            memcpy(pinv_ptr, pinv_ptr - num_rows * 2, sizeof(pinv[0]) * num_rows * 2);
            pinv_ptr -= num_rows * 2;
        }
        pinv_ptr = pinv + col_index * num_rows * 2;
        memset(pinv_ptr, 0, sizeof(pinv[0])*num_rows*2);

        // inc num_cols
        num_cols++;

        // then we call update_col to add in the col
        need2recalculate = updateCol(col_index, num_elements, col);
        return need2recalculate;
    }



    // add a row[nele] to adr location in A[r][co]
    // side effect: matrix, m_An, pinv, new_pinv are changed
    inline bool addRow(IN const uint32_t & row_index, IN const uint32_t & num_elements, IN const ComplexPtr & row = NULL) forceinline {

        // add a row of 0s to matrix
        double * matrix_ptr = matrix + num_rows * num_cols * 2;
        for (uint32_t i = num_rows; i > row_index; i--)      {
            memcpy(matrix_ptr, matrix_ptr - num_cols * 2, sizeof(matrix[0]) * num_cols * 2);
            matrix_ptr -= num_cols * 2;
        }
        matrix_ptr = matrix + row_index * num_cols * 2;
        memset(matrix_ptr, 0, sizeof(matrix[0])*num_cols*2);

        // add a column of 0s to pinv
        double * tmp_matrix_ptr = tmp_matrix;
        double * pinv_ptr = pinv;
        for (uint32_t i = 0; i < num_cols; i++)      {
            memcpy(tmp_matrix_ptr, pinv_ptr, sizeof(tmp_matrix[0])*row_index*2);
            pinv_ptr += row_index*2;
            tmp_matrix_ptr += row_index*2;
            *(tmp_matrix_ptr++) = 0;    // create a column of 0s
            *(tmp_matrix_ptr++) = 0;
            memcpy(tmp_matrix_ptr, pinv_ptr, sizeof(tmp_matrix[0])*(num_rows-row_index)*2);
            // It was A_ here, but it should be Ai_
            pinv_ptr += (num_rows-row_index)*2;
            tmp_matrix_ptr += (num_rows-row_index)*2;
        }
        swapPointers(reinterpret_cast<double **>(&pinv),
                reinterpret_cast<double **>(&tmp_matrix));

        // inc num_rows
        num_rows++;

        // then we call update_row to plug the row in
        need2recalculate = updateRow(row_index, num_elements, row);
        return need2recalculate;
    }



    // remove a column at location col_index from A[r][co]
    inline bool removeCol(IN const uint32_t & col_index) forceinline {

        // clear out the col to be 0s
        //memset0_complex_vector(num_rows, zero_col_vec);
        //memset(zero_col_vec, 0, sizeof(Complex) * num_rows);
        need2recalculate = updateCol(col_index, num_rows, zero_col_vec);

        // check the row in pinv is 0s
        double * pinv_ptr = NULL;

        // remove the col in matrix
        double * matrix_ptr = matrix;                 // remove column
        double * tmp_matrix_ptr = tmp_matrix;
        for (uint32_t i = 0; i < num_rows; i++) {
            memcpy(tmp_matrix_ptr, matrix_ptr, sizeof(tmp_matrix[0])*col_index*2);
            matrix_ptr += (col_index+1)*2;
            tmp_matrix_ptr += col_index*2;
            memcpy(tmp_matrix_ptr, matrix_ptr, sizeof(tmp_matrix[0])*(num_cols-col_index-1)*2);
            matrix_ptr += (num_cols-col_index-1)*2;
            tmp_matrix_ptr += (num_cols-col_index-1)*2;
        }
        swapPointers(reinterpret_cast<double **>(&matrix),
                reinterpret_cast<double **>(&tmp_matrix));

        // remove the row in pinv
        pinv_ptr = pinv+col_index*num_rows*2;
        for (uint32_t i = col_index; i < num_cols-1; i++) {
            memcpy(pinv_ptr, pinv_ptr+num_rows*2, sizeof(pinv[0])*num_rows*2);
            pinv_ptr += num_rows*2;
        }

        // dec num_cols
        num_cols--;
        //matrix_mat.resize(Eigen::NoChange, num_cols);
        //pinv_mat.resize(num_cols, Eigen::NoChange);
        return need2recalculate;
    }



    // remove a row at location row_index from A[r][co]
    // side effect: matrix, m_An, pinv, new_pinv are changed
    inline bool removeRow(IN const uint32_t & row_index) forceinline {
        //static uint32_t i; //, j;
        //static double* a_;
        //assert(row_index < num_rows);

        // clear out the row to be 0s
        //memset0_complex_vector(num_cols, zero_row_vec);
        //memset(zero_row_vec, 0, sizeof(Complex) * num_cols);
        need2recalculate = updateRow(row_index, num_cols, zero_row_vec);

        double * pinv_ptr = NULL;

        // remove the row in matrix
        double * matrix_ptr = matrix + row_index * num_cols * 2;
        for (uint32_t i = row_index; i < num_rows-1; i++) {
            memcpy(matrix_ptr, matrix_ptr + num_cols * 2, sizeof(matrix[0])*num_cols*2);
            matrix_ptr += num_cols * 2;
        }

        // remove a column in pinv
        pinv_ptr = pinv;
        double * tmp_matrix_ptr = tmp_matrix;
        for (uint32_t i = 0; i < num_cols; i++) {
            memcpy(tmp_matrix_ptr, pinv_ptr, sizeof(tmp_matrix[0])*row_index*2);
            pinv_ptr += (row_index+1)*2;
            tmp_matrix_ptr  += row_index*2;
            memcpy(tmp_matrix_ptr, pinv_ptr, sizeof(tmp_matrix[0])*(num_rows-row_index-1)*2);
            pinv_ptr += (num_rows-row_index-1)*2;
            tmp_matrix_ptr += (num_rows-row_index-1)*2;
        }
        swapPointers(reinterpret_cast<double **>(&pinv),
                reinterpret_cast<double **>(&tmp_matrix));

        // dec num_rows
        num_rows--;
        //matrix_mat.resize(num_rows, Eigen::NoChange);
        //pinv_mat.resize(Eigen::NoChange, num_rows);

        return need2recalculate;
    }



    // mul the inverted matrix pinv[co][r] with a given vector v[nele] into the
    // vector re[nreturnele]
    inline void solveSys(IN const uint32_t & num_elements, IN const ComplexPtr & v, IN const uint32_t & num_return_elements, IN ComplexPtr result = NULL) forceinline {

        //assert(num_elements == num_rows);
        //assert(num_return_elements == num_cols);
        if (result)
            matVecMul(num_cols, num_rows, pinv, (double *)v, (double *)result);
        else
            matVecMul(num_cols, num_rows, pinv, (double *)v, (double *)linear_solution);

    }

};


#endif  // INT_INC_INVLIB_H_
