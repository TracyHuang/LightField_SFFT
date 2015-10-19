#include "initialization.h"

/**
 * @file    initialization.cc
 * @brief   Implementing the Initialization class
 * @author  Tianwei
 * @date    07/13/2014
 */

using namespace NonIntSFFT;

void Initialization::setParameters(Params * _params) {

    deallocateBuffers();

    params = _params;

    
    // allocate projection array
    num_projection_lines = params->num_projection_lines;
    projections = (Projection *) malloc(sizeof(Projection) * num_projection_lines);
    Projection * cur_projection = projections;

    // set the parameters to the projections 
    uint32_t projection_index = 0;
    int32_t * cur_projection_matrix = params->projection_matrices;
    for (; projection_index < num_projection_lines; projection_index ++, cur_projection ++) {
		*cur_projection = Projection(*(cur_projection_matrix), *(cur_projection_matrix + 1), *(cur_projection_matrix + 2), *(cur_projection_matrix + 3));
		cur_projection->calBucketCoeffs();
		cur_projection->printProjection();
		cur_projection_matrix += 4;
    }
    
    position_table = (int32_t *) malloc (sizeof(int32_t) * NonIntSFFT::n);
    fftw_input = (ComplexPtr)fftw_malloc(sizeof(Complex) * NonIntSFFT::n);
    sampling_table = (bool *)malloc(sizeof(bool) * NonIntSFFT::n);
   
	residue_graph = (Real **) malloc (sizeof(Real *) * NonIntSFFT::n_v);
	for (uint32_t i = 0; i < NonIntSFFT::n_v; i ++) {
		residue_graph[i] = (Real *) malloc(sizeof(Real) * NonIntSFFT::n_h);
		for (uint32_t j = 0; j < NonIntSFFT::n_h; j ++) {
			residue_graph[i][j] = 0.0;
		}
	}
	
	sample_graph = (int32_t **) malloc (sizeof(int32_t *) * NonIntSFFT::n_v);
	for (uint32_t i = 0; i < NonIntSFFT::n_v; i ++) {
		sample_graph[i] = (int32_t *) malloc(sizeof(int32_t) * NonIntSFFT::n_h);
		for (uint32_t j = 0; j < NonIntSFFT::n_h; j ++) {
			sample_graph[i][j] = 0;
		}
	}
 
    generateSamplingTable();

    return;
}

void Initialization::setRecovery() {

    recovery.setProjections(projections, num_projection_lines, params->scaling_threshold);

}

void Initialization::deallocateBuffers() {
 
    if (projections) {
		free(projections);
		projections = NULL;
    }

    if(position_table) {
		free(position_table);
		position_table = NULL;
    }

    if (fftw_input) {
        fftw_free(fftw_input);
        fftw_input = NULL;
    }

    if (sampling_table) {
        fftw_free(sampling_table);
        sampling_table = NULL;
    }
	
	if (residue_graph) {
		for (uint32_t i = 0; i < NonIntSFFT::n_v; i ++) {
			free(residue_graph[i]);
			residue_graph = NULL;
		}
		free(residue_graph);
		residue_graph = NULL;
	}

	if (sample_graph) {
		for (uint32_t i = 0; i < NonIntSFFT::n_v; i ++) {
			free(sample_graph[i]);
			sample_graph[i] = NULL;
		}
		free(sample_graph);
		sample_graph = NULL;
	}
}

void Initialization::getPeaksResidue(uint32_t * _num_peaks, PeaksList * _peaks, Real * _residue) {
    *_num_peaks = num_peaks;
    *_peaks = peaks;
    *_residue = residue;
}




bool Initialization::shadowBucketVoting(out NonIntPeak * added_peak) {

	if (!added_peak) {
		fprintf(stderr, "Fatal error shadowBucketVoting, added_peaks is NULL\n");
		exit(-1);
	}

	double ** voting_graph = (double **) malloc(sizeof(double *) * NonIntSFFT::n_v);
	for (uint32_t i = 0; i < NonIntSFFT::n_v; i ++) {
		voting_graph[i] = (double *) malloc(sizeof(double) * NonIntSFFT::n_h);
		for (uint32_t j = 0; j < NonIntSFFT::n_h; j ++) {
			voting_graph[i][j] = 0.0;
		}
	}

	Projection * cur_projection = projections;
	int32_t * pos = NULL;
	pos = (int32_t *) malloc(sizeof(int32_t) * 2 * cur_projection->n);
	int32_t num_pos = 0;
	IntPosition * bucket_pos;
	for (uint32_t i = 0; i < num_projection_lines; i ++, cur_projection ++) {	
		bucket_pos = cur_projection->bucket_table;
		for (int32_t bucket_index = 0; bucket_index < cur_projection->n; bucket_index ++) {
			cur_projection->getProj(bucket_index, pos, &num_pos);
			for (int32_t j = 0; j < num_pos; j ++) {
				voting_graph[pos[2*j]][pos[2*j+1]] += residue_graph[bucket_pos[bucket_index].getX()][bucket_pos[bucket_index].getY()];
			}
		}
	} 		

	int32_t pos_x = 0;
	int32_t pos_y = 0;
	double max = Utils::max2DArray(voting_graph, NonIntSFFT::n_v, NonIntSFFT::n_h, &pos_x, &pos_y);
	double median = Utils::median(voting_graph, NonIntSFFT::n_v, NonIntSFFT::n_h);

	if (max > median * params->init_param.min_value_ratio_threshold) {
		
		added_peak->setXY(pos_x, pos_y);
		added_peak->setValue(Complex(0.0, 0.0));

		// deallocate space of voting graph
		for (uint32_t i = 0; i < NonIntSFFT::n_v; i ++) {
			free(voting_graph[i]);
			voting_graph[i] = NULL;
		}
		free(voting_graph);
		voting_graph = NULL;
		
		return true;

	}
	else {	
		// deallocate space of voting graph
		for (uint32_t i = 0; i < NonIntSFFT::n_v; i ++) {
			free(voting_graph[i]);
			voting_graph[i] = NULL;
		}
		free(voting_graph);
		voting_graph = NULL;

		return false;
	}		

}



