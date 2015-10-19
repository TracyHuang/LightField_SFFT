/**
 * @file    common.cc
 * @brief   The main implementation to the definitions in common.h
 * @author  lixin
 * @date    11/13/2012
 */

#include "common.h"
#include "utils.h"

namespace NonIntSFFT {

    const char ERROR_NAMES[TOTAL_NUM_ERRORS][MAX_STRING_LENGTH] = {
        "Null Pointer",         // NULL_POINTER
        "Unimplemented",        // UNIMPLEMENTED
        "Negative Value",       // SMALLER_THAN_ZERO
        "Non-square matrix",    // NONSQUARE_MAT
        "Wrong enum number",    // WRONG_ENUM_NUMBER
        "Matrix too large",     // MATRIX_TOO_LARGE
        "Infinity value",       // INFINITY_VALUE
        "Peak-num mismatch",    // NUM_PEAKS_NOMATCH
        "Internal error",       // INTERNAL_ERROR
        "Cannot open file",     // CANNOT_OPEN_FILE
        "IO error",             // IO_ERROR
        "Dim unsupported",      // DIM_UNSUPPORTED
        "Invalid paral ID",     // INVALID ID
        "Invalid parameter",    // INVALID_PARAM
        "Format error",         // FORMAT_ERR
        "Need recalculation",   // NEED_TO_RECALCULATE
    };

    uint32_t n_v = 0;
    uint32_t n_h = 0;
    uint32_t n = 0;
    bool error = 0;

}

ComplexPtr Utils::horizontalSinc = NULL;
ComplexPtr Utils::verticalSinc = NULL;


