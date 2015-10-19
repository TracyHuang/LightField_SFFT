/**
 * @file    main.cc
 * @brief   The main entry to the non-integer sparse fft code
 * @author  lixin
 * @date    11/13/2012
 */


#include "sfft.h"
#include "utils.h"

#include <iostream>

int main(int argc, char ** argv) {

    uint32_t num_parals = 1;
    int32_t paral_id = 0;
    char param_file_name[MAX_STRING_LENGTH] = {0};
    char io_prefix[MAX_STRING_LENGTH] = {0};

    switch (argc) {

    default:
    case 5: 
        strcpy(io_prefix, argv[4]); 
        fprintf(stderr, "io_prefix = %s\n", io_prefix);
    case 4:
        num_parals = atoi(argv[3]);
        fprintf(stderr, "num_parals = %d\n", num_parals);
    case 3:
        paral_id = atoi(argv[2]);
        fprintf(stderr, "paral_id = %d\n", paral_id);
    case 2: 
        strcpy(param_file_name, argv[1]);
        fprintf(stderr, "param_file_name = %s\n", param_file_name);
        break;
    case 1:
        param_file_name[0] = 0;
    }

    clock_t start_time, end_time;

    start_time = clock();

    SFFT * sfft = new SFFT(param_file_name, paral_id, num_parals, io_prefix);

    sfft->run();

   // delete sfft;

    end_time = clock();

    double time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    fprintf(stderr, "TOTAL TIME OF WHOLE PROGRAM EXECUTION = %lf\n", time);


    return 0;
}
