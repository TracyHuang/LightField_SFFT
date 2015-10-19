#include "gradient.h"

/**
 * @file    gradient.cc
 * @brief   Implementing the Gradient class
 * @author  Tianwei
 * @date    07/11/2014
 */


using namespace NonIntSFFT;

void Gradient::setParameters(Params * _params) {

    deallocateBuffers();

    params = _params;

    // setup samples, sample_values, sample_pos, num_samples
    num_samples = params->num_samples;
    if(sample_pos) {
	free(sample_pos);
	sample_pos = NULL;
    }
    sample_pos = (IntPosition *) malloc(sizeof(IntPosition) * num_samples);

    // @question maybe need to consider whether to input time domain signal or get it from input?
    int32_t * pos = params->sample_pos;
    //Real * values = params->sample_values;
    for (uint32_t i = 0; i < num_samples; i ++) {
	sample_pos[i].setPosition(pos[2 * i], pos[2 * i + 1]);
        //sample_values[i] = values[i];
    }

    residue_window = (Real *)malloc(sizeof(Real) * params->residue_window_size);
    optimal_y = (ComplexPtr)fftw_malloc(sizeof(Complex) * NonIntSFFT::n);
    fftw_input = (ComplexPtr)fftw_malloc(sizeof(Complex) * NonIntSFFT::n);
    sampling_table = (bool *)malloc(sizeof(bool) * NonIntSFFT::n);

	residue_graph = (Real **) malloc(sizeof(Real *) * NonIntSFFT::n_v);
	for (uint32_t i = 0; i < NonIntSFFT::n_v; i ++) {
		residue_graph[i] = (Real *) malloc(sizeof(Real) * NonIntSFFT::n_h);
		for (uint32_t j = 0; j < NonIntSFFT::n_h; j ++) {
			residue_graph[i][j] = 0.0;
		}
	}

	sample_graph = (int32_t **) malloc(sizeof(int32_t *) * NonIntSFFT::n_v);
	for (uint32_t i = 0; i < NonIntSFFT::n_v; i ++) {
		sample_graph[i] = (int32_t *) malloc(sizeof(int32_t) * NonIntSFFT::n_h);
		for (uint32_t j = 0; j < NonIntSFFT::n_h; j ++) {
			sample_graph[i][j] = 0;
		}
	}

	// set up sample_graph

	for (uint32_t i = 0; i < num_samples; i ++) {
		sample_graph[sample_pos[i].getX()][sample_pos[i].getY()] = 1;		
	}

    generateSamplingTable();

    return;
}

void Gradient::getSampleValues(ComplexPtr x_in) {
    if(sample_values) {
	free(sample_values);
	sample_values = NULL;
    } 
    sample_values = (Complex *) malloc(sizeof(Complex) * num_samples);
    for (uint32_t i = 0; i < num_samples; i ++) {
	sample_values[i] = x_in[sample_pos[i].map2Index()];
    }

    return;
}




void Gradient::setGradient(uint32_t _num_peaks, PeaksList _peaks, Real _residue, int32_t * recovery_id) {
    
    num_peaks = _num_peaks;
    peaks = _peaks;
    residue = _residue;
    // we also need to reset recovery
    recovery.retrieveConfig(-1, recovery_id, &peaks, num_peaks, &residue);

}



void Gradient::setRecovery(Params * params, ComplexPtr x_in) {

  recovery.setRecoverConfigSamples(params, x_in); 
}

void Gradient::deallocateBuffers() {

    if (residue_window) {
        free(residue_window);
        residue_window = NULL;
    }

    if (optimal_y) {
        fftw_free(optimal_y);
        optimal_y = NULL;
    }

    if(fftw_input) {
	fftw_free(fftw_input);
        fftw_input = NULL;
    }

    if (sampling_table) {
        fftw_free(sampling_table);
        sampling_table = NULL;
    }

	if (residue_graph) {
		for (uint32_t i = 0; i < NonIntSFFT::n_v; i ++) {
			if (residue_graph[i]) {
				free(residue_graph[i]);
				residue_graph[i] = NULL;
			}
		}
		if (residue_graph) {
			free(residue_graph);
			residue_graph = NULL;
		}
	}

	if (sample_graph) {
		for (uint32_t i = 0; i < NonIntSFFT::n_v; i ++) {
			if (sample_graph[i]) {
				free(sample_graph[i]);
				sample_graph[i] = NULL;
			}
		}
		if (sample_graph) {
			free(sample_graph);
			sample_graph = NULL;
		}
	}


}

void Gradient::gradientSearch(ComplexPtr x_in ,int32_t _recovery_id) {
     
    residue_writer = 0;
    maxima_residue = 0;
    minima_residue = INFINITY;
    optimal_residue = INFINITY;
    memset(residue_window, -1, params->residue_window_size * sizeof(Real));
    traced_back = 0;

    int32_t recovery_id = _recovery_id;

    state = INIT; 
 
    // do gradient decent algorithm
    iter = 0;
    int stucked_num = 0;
    bool merged = 0;
    bool periodical = 0;
    bool small_change = 0;
    bool suboptimal = 0;
    bool try_removed = 0;

    bool well_determined;
    // We also need to update for the initial state
    updateResidueWindow(&periodical, &small_change, &suboptimal, recovery_id);
    for (; (iter < params->max_iter) && (state != FINISHED); ++ iter) {

        switch (state) {
        case INIT: 
        case CHECK_UPDATED:
        case EXHAUSTIVE_UPDATE:
            stucked_num = 0;
            gradientOption(&well_determined, &recovery_id); 
            break;
        case GRADIENT_STUCK:
            stucked_num ++;
            // exhaustive
            exhaustiveSearch(&well_determined, &recovery_id);
            break;
        case GRADIENT_UPDATE:
            stucked_num = 0;
            // (joint) gradient
            gradientOption(&well_determined, &recovery_id);  
            break;
        case EXHAUSTIVE_STUCK:
            stucked_num ++;
            // should consider finish
            state = FINISHED;
            break;
        default:
            ;
        }

        if (NonIntSFFT::error) {
			fprintf(stderr, "error\n" );	
            return;
		}

        // why we need to update here? 
        // because sometimes merger screw up things
        updateResidueWindow(&periodical, &small_change, &suboptimal, recovery_id);
        // first, if they are degenerated, try fast merge
        merged = false;
        merged = fastMerge(&well_determined, &recovery_id);
        
        
        if (NonIntSFFT::error) {
            return;
		}        
        // do merge
        if (!merged || !well_determined) { 
            merged = merge(&well_determined, &recovery_id);
        }
        
        // if merged
        if (merged) {
            stucked_num = 0;
            switch (state) {
            case GRADIENT_STUCK:
                state = GRADIENT_UPDATE;
                break;
            case EXHAUSTIVE_STUCK:
            case FINISHED:
                state = EXHAUSTIVE_UPDATE;
                break;
            case SUBTRACT_SHIFT_STUCK:
                state = SUBTRACT_SHIFT_UPDATE;
                break;
            default: 
                //printf("unknown state %d\n", state);
                ;
            } 
            updateResidueWindow(&periodical, &small_change, &suboptimal, recovery_id);
        }
		else {

			// do check update
			bool check_updated = checkEst( &recovery_id);

			// if success
			if (check_updated) {
				stucked_num = 0;
				state = CHECK_UPDATED;
				updateResidueWindow(&periodical, &small_change, &suboptimal, recovery_id);
			}
		}        


        if (suboptimal) {
            // if this is the second time we trace back the same profile
            if (!traceback(&recovery_id)) {
                break;
			}
        }
        
        // now let's see what to do
        if (periodical || small_change || suboptimal) {
            if (!(try_removed = try2Remove(&recovery_id))) {
                // stop
                state = FINISHED;
            } else {
                // now let's see what to do
                state = CHECK_UPDATED;
            }

        }
    }
    
	// now, good, let's recover the frequency domain based on the values we get
    recoverFreqDomain(x_in);
    
    
}


void Gradient::singleGradientSearch(out bool * well_determined, out int32_t * recovery_id) {

    int32_t new_recovery_id = 0;
    Real new_residue = residue;
   
    recovery.retrieveConfig(*recovery_id, &new_recovery_id, &peaks, num_peaks, &new_residue);

    bool not_stucked = true;
    NonIntPosition old_pos;

    Real gradient[3][3] = {{NAN, NAN, NAN}, {NAN, NAN, NAN}, {NAN, NAN, NAN}};
    Real residues[3][3];
    int32_t direction_v = 0;
    int32_t direction_h = 0;
    int32_t i_v = 0;
    int32_t i_h = 0;
    int32_t new_i_v = 0;
    int32_t new_i_h = 0;
    uint32_t update_index = 0;

    PeaksList ci_peaks;
    PeakIterator ci_it;
    ci_peaks = peaks;
    ci_it = ci_peaks.begin();

    double difference = 0;

    for (PeakIterator it = peaks.begin(); it != peaks.end(); ++ it, ++ update_index) {

        memset(gradient, -1, sizeof(Real) * 9);

        old_pos = it->getPos(); 

        not_stucked = 1;
        while (not_stucked) {

            // calculate the gradient
            *well_determined = calGradient(it, update_index, new_residue, gradient, residues, &new_recovery_id, &direction_v, &direction_h);
            
            if (NonIntSFFT::error) return;

            if (direction_v == 0 && direction_h == 0)
                not_stucked = 0;
            else {
                // update position    
                it->setPosition(it->getX() + params->step * (double)direction_v, it->getY() + params->step * (double)direction_h);

                // calculate the new residue and values
                if (params->delta == params->step) {
                    new_residue = residues[direction_v + 1][direction_h + 1];

                    // more than that: update the residue matrix
                    for (i_v = -1; i_v <= 1; ++ i_v)	{
						for (i_h = -1; i_h <= 1; ++i_h)    {
                            new_i_v = direction_v + i_v;
                            new_i_h = direction_h + i_h;
                            if (new_i_v >= -1 && new_i_v <= 1 && new_i_h >= -1 && new_i_h <= 1)
                                gradient[i_v + 1][i_h + 1] = residues[new_i_v + 1][new_i_h + 1];
                            else
                                gradient[i_v + 1][i_h + 1] = NAN;
                        }
					}

                }
				else {
                    printf("singleGradientSearch: dead end, we don't come here\n");
                    *well_determined = recovery.updateConfig(&peaks, num_peaks, it, update_index, &new_residue);

                    if (NonIntSFFT::error)
                        return;
                }



            }

		}

        difference += it->getPos().gridDistance(old_pos);

        if (!params->commit_step_by_step) {
            // revert everything back

            ci_it->setPosition(it->getPos());
            ++ ci_it;
            
            it->setPosition(old_pos);

            *well_determined = recovery.updateConfig(&peaks, num_peaks, it, update_index, &new_residue);
        }
    }
    
    // commit the change
    if (difference < params->tolerance) {
        *well_determined = true;
        state = GRADIENT_STUCK;
    } else {
        if (!params->commit_step_by_step) {
            peaks = ci_peaks;
            *well_determined = recovery.retrieveConfig(-1, recovery_id, &peaks, num_peaks, &residue);
            if (NonIntSFFT::error)
                return;
        } else {
            residue = new_residue;
        }
        state = GRADIENT_UPDATE;
    }

}


void Gradient::exhaustiveSearch(out bool * well_determined, out int32_t * recovery_id) {
    
    int32_t new_recovery_id = 0;
    Real new_residue = residue;
    
    recovery.retrieveConfig(*recovery_id, &new_recovery_id, &peaks, num_peaks, &new_residue);

    NonIntPosition old_pos;

    Real min_residue = residue;
    double direction_v = 0;
    double direction_h = 0;
    uint32_t update_index = 0;

    PeaksList ci_peaks;
    PeakIterator ci_it;
    if (!params->commit_step_by_step) {
        ci_peaks = peaks;
        ci_it = ci_peaks.begin();
    }

    double difference = 0;
    double shift_v = 0;
    double shift_h = 0;

    for (PeakIterator it = peaks.begin(); it != peaks.end(); ++ it, ++ update_index) {

        old_pos = it->getPos(); 

        min_residue = residue;

        for (shift_v = -params->range; shift_v < params->range; shift_v += params->epsilon) {
            for (shift_h = -params->range; shift_h < params->range; shift_h += params->epsilon) {

                it->setPosition(old_pos.getX() + shift_v, old_pos.getY() + shift_h);
                *well_determined = recovery.updateConfig(&peaks, num_peaks, it, update_index, &new_residue);
                if (NonIntSFFT::error)
                    return;
                
                if (new_residue < min_residue) {
                    min_residue = new_residue;
                    direction_v = shift_v;
                    direction_h = shift_h;
                }

            }
        }

        difference += fabs(direction_v) + fabs(direction_h);

        if (!params->commit_step_by_step) {
            // revert everything back

            ci_it->setPosition(old_pos.getX() + direction_v, old_pos.getY() + direction_h);
            ++ ci_it;
            
            it->setPosition(old_pos);
            recovery.updateConfig(&peaks, num_peaks, it, update_index, &new_residue);

        } else {
            it->setPosition(old_pos.getX() + direction_v, old_pos.getY() + direction_h);
            // reset residue as well
            *well_determined = recovery.updateConfig(&peaks, num_peaks, it, update_index, &residue);
        }
    }

    // commit the change
    if (difference < params->tolerance) {
        *well_determined = true;
        if (params->subtract_and_shift_enabled)
            state = EXHAUSTIVE_STUCK;
        else
            state = FINISHED;
    } else {

        if (!params->commit_step_by_step) {
            peaks = ci_peaks;
            *well_determined = recovery.retrieveConfig(-1, recovery_id, &peaks, num_peaks, &residue);
            if (NonIntSFFT::error)
                return;
        }

        state = EXHAUSTIVE_UPDATE;

    }

}

bool Gradient::fastMerge(out bool * well_determined, in out int32_t * recovery_id) {

    // step1: if two peaks are very near, merge them!
    // In this step, the residue might increase dramatically, just because
    // previously it is a presudo-sufficient system, but actually deficient. 
    // So I would say, ignore the residue, man!!
    PeakIterator it1, it2;
    uint32_t peak_index = 0;
    uint32_t peak_index2 = 0;
    double distance = 0;
    bool merged = false; 
    bool merged_this_time = true;

    for (it1 = peaks.begin(); it1 != peaks.end();) {

        it2 = it1;
        merged_this_time = false;
        for (peak_index2 = peak_index + 1, it2 ++; it2 != peaks.end(); ++ it2, ++ peak_index2) {
            distance = it1->gridDistance(*it2);

            if (distance < params->tolerance) {

                // merge
                *well_determined = mergeTogether(it1, it2, &peaks, num_peaks, peak_index, peak_index2, &residue);
                if (NonIntSFFT::error)
                    return false;

                merged_this_time = true;
                break;

                // update
                merged = true;
                continue;
            } 
        }
        if (!merged_this_time) {
            ++ it1;
            ++ peak_index;
        }
    }

    if (*well_determined)
        return merged;

    // step2: if we still cannot get the system enough rank, let's do our
    // best try: merge the nearest pair of peaks
    double min_distance = params->min_peak_distance;
    int32_t min_peak1 = -1;
    int32_t min_peak2 = -1;
    PeakIterator min_it1, min_it2;

    // PeaksList new_peaks = peaks;

    for (peak_index = 0, it1 = peaks.begin(); it1 != peaks.end(); ++ it1, ++ peak_index) {

        it2 = it1; 
        for (peak_index2 = peak_index + 1, it2 ++; it2 != peaks.end(); ++ it2, ++ peak_index2) {
            distance = it1->gridDistance(*it2);
            if (distance < min_distance) {
                min_distance = distance;
                min_peak1 = peak_index;
                min_peak2 = peak_index2;
                min_it1 = it1;
                min_it2 = it2;
            }
        }

    }

    // If we can merge 
    if (min_peak1 >= 0) {
        recovery.retrieveConfig(*recovery_id, recovery_id, &peaks, num_peaks, &residue);
        * well_determined = mergeTogether(min_it1, min_it2, &peaks, num_peaks, min_peak1, min_peak2, &residue);
        if (NonIntSFFT::error)
            return false;
        
		return true;
    }
    
    return merged;
}

bool Gradient::merge(out bool * well_determined, in out int32_t * recovery_id) {

    PeaksList new_peaks = peaks;
    uint32_t new_num_peaks = num_peaks;

    double min_distance = 0;

    PeakIterator it, it2, min_it;
    uint32_t peak_index = 0;
    uint32_t peak_index2 = 0;
    int32_t min_peak_index = 0;
    double distance = 0;
    Real new_residue = 0;
    bool merged = false;

    int32_t new_recovery_id = -1;
    
    recovery.retrieveConfig(*recovery_id, recovery_id, &peaks, num_peaks, &new_residue);
    
    for (it = new_peaks.begin(); it != new_peaks.end();) {
		
		min_distance = params->min_peak_distance;
        
        min_peak_index = -1;
        it2 = it;
        for (peak_index2 = peak_index + 1, it2 ++; it2 != new_peaks.end(); ++ it2, ++ peak_index2) {
            distance = it->gridDistance(*it2);

            if (distance < min_distance) {
                min_distance = distance;
                min_it = it2;
                min_peak_index = peak_index2;
            }
        }
        if (min_peak_index >= 0) {
            // remove things from it
           
            if (!merged) {
                // replicate the configuration only when necessary
                recovery.replicateConfig(recovery_id, &new_recovery_id);
                merged = 1;
            }

            mergeTogether(it, min_it, &new_peaks, new_num_peaks, peak_index, min_peak_index, &new_residue);
            if (NonIntSFFT::error)
                return false;

        } else {
            ++ it;
            ++ peak_index;
        }
    }

    if (!merged)
        return merged;

    if (new_residue < residue * params->amend_peak_tolerance) {
        
		// in this case, merge!
		peaks = new_peaks;
        num_peaks = new_num_peaks;
        residue = new_residue;
        * well_determined = true;
        * recovery_id = new_recovery_id;
        return true;
    } else {
        // revert back
        recovery.retrieveConfig(*recovery_id, recovery_id, &peaks, num_peaks, &residue);

        return false;
    }

}

bool Gradient::try2Remove(in out int32_t * recovery_id) {

    if (num_peaks == 1) {
        fprintf(stderr, "One peak: not going to remove");
        return false;
    }


    Real removable_score;
    Real max_removable_score = -INFINITY;
    uint32_t best_peak_index = -1;
    
    PeakIterator it, best_peak;
    int32_t new_recovery_id = -1;
    Real new_residue = 0;
    uint32_t peak_index = 0;
    memset(peak_indicator, 0, sizeof(bool) * num_peaks);
    uint32_t num_deletes = 0;

    for (peak_index = 0, it = peaks.begin(); it != peaks.end(); ++ it, ++ peak_index) {

        // we don't actually remove the peak, but rather, just evaluate the removing gain
        // this function will change the Peaks::value field, but who cares
        recovery.retrieveConfig(*recovery_id, recovery_id, &peaks, num_peaks, &residue);
        recovery.replicateConfig(recovery_id, &new_recovery_id);
        recovery.shrinkConfig(&peaks, num_peaks - 1, peak_index, &new_residue);

        removable_score = residue - new_residue;

        if (removable_score > -params->noticable_change * residue) {
            // maybe we really want to remove this guy
            peak_indicator[peak_index] = 1;
            num_deletes ++;
            
            if (removable_score > max_removable_score) {
                best_peak = it;
                max_removable_score = removable_score;
                best_peak_index = peak_index;
            }
        }

    }

    recovery.retrieveConfig(*recovery_id, recovery_id, &peaks, num_peaks, &residue);

    if (num_deletes == 0) {
        return false;
    }

    // now let's see what to do
    recovery.replicateConfig(recovery_id, &new_recovery_id);
    
    // first try: remove everything
    PeaksList new_peaks = peaks;
    uint32_t new_num_peaks = num_peaks;
    uint32_t old_peak_index = 0;
    for (peak_index = 0, it = new_peaks.begin(); it != new_peaks.end(); ++ old_peak_index) {

        if (peak_indicator[old_peak_index]) {
            new_num_peaks --;

            it = new_peaks.erase(it);
            recovery.shrinkConfig(&new_peaks, new_num_peaks, peak_index, &new_residue);

        } else {
            ++ it; 
            ++ peak_index;
        }
    }

    if (new_residue < residue * params->amend_peak_tolerance) {
        // update
        peaks = new_peaks;
        num_peaks = new_num_peaks;
        residue = new_residue;
        * recovery_id = new_recovery_id;
        return true;
    }

    // second try: only remove the most-removable guy

    PeaksList old_peaks = peaks;
    recovery.retrieveConfig(*recovery_id, recovery_id, &peaks, num_peaks, &residue);
    recovery.replicateConfig(recovery_id, &new_recovery_id);
    peaks.erase(best_peak);
    recovery.shrinkConfig(&peaks, --num_peaks, best_peak_index, &new_residue);
    if (new_residue < residue * params->amend_peak_tolerance) {
        // update
        residue = new_residue;
        * recovery_id = new_recovery_id;
        return true;
    } else {

        peaks = old_peaks;
        num_peaks ++;
        return false;
    }

}

void Gradient::recoverFreqDomain(in ComplexPtr x_in) {

    if (!params->fix_time_samples) {
        NonIntPeak::reconstructFrequencyDomain(&optimal_peaks, optimal_y);
        return;
    }

    if (!init_plan) {
        // For the flags, either FFTW_MEASURE or FFTW_ESTIMATE, 
        // FFTW_MEASURE instructs FFTW to run and measure the execution time of several 
        // FFTs in order to find the best way to compute the transform of size n. This process takes some time (usually a few seconds), 
        // depending on your machine and on the size of the transform. FFTW_ESTIMATE, on the contrary, 
        // does not run any computation and just builds a reasonable plan that is probably sub-optimal.
        plan = fftw_plan_dft_2d(NonIntSFFT::n_v, NonIntSFFT::n_h, (fftw_complex *)((void *)fftw_input), (fftw_complex *)((void *)optimal_y), FFTW_FORWARD, FFTW_MEASURE);
    }


    // fixing the time domain samples
    // One import issue is that if we use FFTW_MEASURE, we need to initialize after the plan
    NonIntPeak::reconstructTimeDomain(&optimal_peaks, fftw_input);

    uint32_t sample_index = 0;
    ComplexPtr cur_x = x_in;
    for (; sample_index < NonIntSFFT::n; sample_index ++, cur_x ++)
        if (sampling_table[sample_index]) {
            fftw_input[sample_index] = *cur_x;
        }

    if (!init_plan) {
        // execute fft
        fftw_execute(plan);
        init_plan = 1;
    } else {
        fftw_execute_dft(plan, (fftw_complex *)((void *)fftw_input), (fftw_complex *)((void *)optimal_y));
    }



    // scale up result
    Utils::scaleUp<Complex>(optimal_y, n, 1.0 / sqrt(n));

}


/**
 * Find the value of the nearest peak's value
 */
Complex Gradient::nearestPeakValue(PeakIterator cur_peak, PeaksList * _peaks) {
	double min_dis = INFINITY;
	double cur_dis = 0.0;
	Complex result = 0.0;
	PeakIterator target_peak = _peaks->begin();
	for (; target_peak != _peaks->end(); target_peak ++) { 
		cur_dis = std::sqrt( pow( std::abs(target_peak->getX() - cur_peak->getX()), 2.0) + pow( std::abs(target_peak->getY() - target_peak->getY()), 2.0));			
		if (cur_dis < min_dis) {
			min_dis = cur_dis;
			result = target_peak->getV();
		}
	}

	return result;
}





bool Gradient::checkEst(in out int32_t * recovery_id) {

	if (!params->check_est)
		return false;

	// the way we find missing peak is to try to detect whether there is a large peak in residue graph
	// if there is no missing peak, then the residue graph should be random so that the largest/median ratio should not be too high
	// on the other hand, if there exists missing peak, then there should be a clear peak in the residue graph	
	recovery.recoveries[recovery.end].getResidueGraph(&residue_graph);
	initialization->prepareShadowBucketDetection(residue_graph, sample_graph);


	NonIntPeak added_peak = NonIntPeak();
	bool success = false;
	success = initialization->shadowBucketVoting(&added_peak); 

	if (!success) {
		return success;
	}

	peaks.push_back(added_peak);
	num_peaks ++;

	
	// we actually added shadow bucket/missing peak(s)
	recovery.retrieveConfig(-1, recovery_id, &peaks, num_peaks, &residue);

	return true;

}

