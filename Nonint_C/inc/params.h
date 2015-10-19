#ifndef _PARAMS_H_
#define _PARAMS_H_

/**
 * @file    params.h
 * @brief   Defining class Params, dealing with interfaces with input
 * @author  lixin
 * date     02/15/2013
 */

#include "common.h"
#include <list>
#include "peaks.h"
#include <libconfig.h++>

////////////////////////////////// TEST
/*
extern uint32_t test_ws[10];
extern uint32_t test_mats[20];
extern int32_t test_init_offset[10];
extern int32_t test_freq_offset[10];
extern uint32_t test_shifts[10];
*/

/**
 * @brief Dealing with input parameters that tunes the algorithm
 *
 * This class will read parameters from the configuration file
 * and set them appropriately. 
 *
 * @author  lixin
 * @date    11/13/2012
 *
 */
class Params {

public:
    /**
     * @brief   parameters regarding the initial method of gradient search method
     * @author  lixin
     * @date    11/16/2012
     */
    class InitParams {

    public:

        /**
         * The method to initialize the algorithm
         */
        enum INIT_METHOD {
            VOTING,                 /**< the voting algorithm */
            OFDM_TRICK,             /**< estimating the initial positions by OFDM trick */
            MANUAL,                 /**< manual specify the initial peaks */
            TOTAL_NUM_INIT_METHODS  /**< the total number of initial method */
        };

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

        /**
         * The method to estimate the initial positions. 
         *
         * Default: INIT_METHOD::VOTING
         */
        INIT_METHOD method; 

        /**
         * The soft power threshold. 
         *
         * Default : 1.0
         *
         * @see Filter::soft_power_threshold
         */
        double soft_power_threshold;

        /**
         * The hard power threshold. 
         *
         * Default : 0
         *
         * @see Filter::hard_power_threshold
         */
        double hard_power_threshold;

        /**
         * The maximum soft power threshold we can have. 
         * Since some voting algorithm is actively changing
         * this soft power threshold, we need to set a 
         * maximum
         *
         * Default: 1000
         */
        double max_soft_power_threshold;

        /**
         * The maximum number of iterations we can have 
         * for the voting algorithm
         *
         * Default: 10
         */
        uint32_t max_voting_iter;

        /**
         * The minimum step of increasing soft_power_threshold
         *
         * Default: 0.1
         */
        double min_soft_threshold_step;

        /**
         * Optunism: if the peak is not manually set, and if 
         * there are manually specified peaks,
         * compare the residue with the manually set peaks
         * if they are better, use the manually set peaks
         *
         * Default: 1
         */
        bool optimism;

        /**
         * In voting, what is the threshold: 
         * i.e. to how many percentage of the buckets vote for this frequency
         * it should be between 0 and 1
         */
        double voting_threshold;

        /**
         * Whether to use the adjacent values
         */
        bool use_adjacent;
		
		/**
         * The maximum number of peaks allowed to be added in shadow bucket detection procedure
         * 
		 * Default : 3
		 */
		uint32_t added_peaks_num_threshold;

		/**
         *	The minimum power ratio needed for the shadow buckets to be recognised as "real" missing peaks
         *
		 * Default : 0.1
		 */
		double min_shadow_power_threshold; 

		/**
         *	The minimum ratio of values between new peak and its nearest original peak to be accepted
         *
		 *	Default : 0.1
		 */
		double min_value_ratio_threshold;

        /**
         * The constructor, set all parameters to their default value
         */
        InitParams() : 
            init_peaks(PeaksList()), num_init_peaks(0), method(VOTING), soft_power_threshold(1.0), hard_power_threshold(0), 
            max_soft_power_threshold(1000.0), max_voting_iter(1000), min_soft_threshold_step(0.1), optimism(1), voting_threshold(1), 
            use_adjacent(0), added_peaks_num_threshold(3), min_shadow_power_threshold(0.1), min_value_ratio_threshold(0.1) {}

    };

    /**
     * @brief   parameters regarding controlling of the whole algorithm input and output
     *
     * @author  lixin
     * @date    11/23/2012
     */
    class ControlParams {

    public:
        /**
         * The method to determine input
         */
        enum INPUT_METHOD {
            FILE_INPUT,                 /**< input from file: FileInput */
            RANDOM_GENERATE,            /**< random input: RandomInput */
            TOTAL_NUM_INPUT_METHODS     /**< the total number of input methods */
        };

        /**
         * The result file
         *
         * Default: y_sfft
         */
        char result_file_prefix[MAX_STRING_LENGTH];

        /**
         * The debugging file
         *
         * Default: result
         */
        char debugging_file_prefix[MAX_STRING_LENGTH];

        /**
         * The initial configuration file
         *
         * Default: init
         */
        char init_file_prefix[MAX_STRING_LENGTH];

        /**
         * The logging file appendix
         *
         * Default: log
         */
        char logging_file_prefix[MAX_STRING_LENGTH];

        /**
         * The input method
         *
         * Default: FILE_INPUT
         */
        INPUT_METHOD input_method;

        /**
         * The input file (if needed)
         *
         * Default: x
         */
        char input_file_prefix[MAX_STRING_LENGTH];

        /**
         * The input n3 (if needed)
         *
         * Default: 100
         */
        uint32_t n3;

        /**
         * The input k (if needed)
         *
         * Default: 5
         */
        uint32_t k;

        /**
         * The input snr (if needed)
         *
         * Default: inf
         */
        double snr;

        /**
         * The input min_distance (if needed)
         *
         * Default: 1
         */
        double min_distance;

        /**
         * The full fft file (if needed)
         *
         * Default: y_full.dat
         */
        char full_file[MAX_STRING_LENGTH];

        /**
         * The ground truth debugging file (if needed)
         *
         * Default: ground_truth.debug
         */
        char ground_truth_debugging[MAX_STRING_LENGTH];

        /**
         * Whether the verbos mode is on or not
         * 
         * Default: verbose
         */
        char verbos_file_prefix[MAX_STRING_LENGTH];

        /**
         * Set up control parameters
         * @param[in]   io_prefix   the folder to which the io is dealing with
         */
        ControlParams(in const char * io_prefix = NULL) : input_method(FILE_INPUT), n3(100), k(5), snr(INFINITY), min_distance(1) {

            if (!io_prefix || strlen(io_prefix) == 0) {
                strcpy(result_file_prefix, "y_sfft");
                strcpy(debugging_file_prefix, "result");
                strcpy(init_file_prefix, "init");
                strcpy(logging_file_prefix, "log");
                strcpy(input_file_prefix, "x");
                strcpy(full_file, "y_full.dat");
                strcpy(ground_truth_debugging, "ground_truth.debug");
                strcpy(verbos_file_prefix, "verbose");
            } else {
                sprintf(result_file_prefix, "%s/y_sfft", io_prefix);
                sprintf(debugging_file_prefix, "%s/result", io_prefix);
                sprintf(init_file_prefix, "%s/init", io_prefix);
                sprintf(logging_file_prefix, "%s/log", io_prefix);
                sprintf(input_file_prefix, "%s/x", io_prefix);
                sprintf(full_file, "%s/y_full.dat", io_prefix);
                sprintf(ground_truth_debugging, "%s/ground_truth.debug", io_prefix);
                sprintf(verbos_file_prefix, "%s/verbose", io_prefix);
            }

            //fprintf(stderr, "%s\n", *((char (*)[MAX_STRING_LENGTH])(&(this->result_file_prefix))));
        }


    };

    /**
     * @brief: these parameters are for light field application
     *
     * @author  lixin
     * @date    12/04/2012
     */
    class LightFieldParams {

    public:
        /**
         * The size of pixel dimension one
         *
         * Default: 480
         */
        uint32_t n_x;

        /**
         * The size of pixel dimension two
         *
         * Default: 640
         */
        uint32_t n_y;

        /**
         * The thickness of the lightfield line. 
         * It depends on the perfectness of Lambertion reflection
         *
         * Default: 6
         */
        uint32_t lam_thickness;

        /**
         * This is to indicate that the camera dimension 
         * and the pixel dimension is consistent or not
         *
         * Default: 1
         */
        bool camera_pixel_consistent;

        /**
         * Whether to impose the line-pattern (dimensionality-gap)
         * of lightfield data set
         *
         * Default: 1
         */
        bool impose_dim_gap;

    public:
        LightFieldParams():
            n_x(480), n_y(640), lam_thickness(6), camera_pixel_consistent(1), impose_dim_gap(1) {}

    };

public: 

    /**
     * The gradient decent method
     */
    enum GRAD_METHOD {
        GRADIENT,       /**< the gradient algorithm along each peak separately */
        JOINT_GRADIENT, /**< the gradient algorithm jointly for all peaks */ 
        EXHAUSTIVE,     /**< exhaustive search algorithm for each peak */
        TOTAL_NUM_GRAD_METHOD, /**< the total number of gradient search algorithms */
    };

    /**
     * The vertical dimension size
     *
     * Default: 15
     */
    uint32_t n_v;

    /**
     * The horizontal dimension size
     *
     * Default: 15
     */
    uint32_t n_h;
    
	// The number of projection slices in the program ie. the number of projection lines in the program
    // Default: 0
    uint32_t num_projection_lines;

    /**
     * Specifying the slope and offset of each projection line
     * [a1, b1, c1, d1, a2, b2, c2, d2, .....]
     * Default: NULL
     */
    int32_t * projection_matrices;

    uint32_t num_samples;

    int32_t * sample_pos;

    /**
     * The params regarding the initial positions
     */
    InitParams init_param;

    /**
     * The maximum number of iterations
     *
     * Default: 100
     */
    uint32_t max_iter;

    /**
     * Which method should use
     *
     * Default: GRAD_METHOD::GRADIENT
     */
    GRAD_METHOD method;

    /**
     * Whether the subtract_and_shift process is enabled. 
     *
     * Default: false
     */
    bool subtract_and_shift_enabled;

    /**
     * Should we commit the change step by step, 
     * or commit them all at once at the end
     *
     * Default: false
     */
    bool commit_step_by_step;

    /**
     * The computational tolerance
     *
     * Default: 1e-5
     */
    double tolerance;

    /**
     * When calculating gradient, the delta distance
     *
     * Default: 0.01
     */
    double delta;

    /**
     * When doing single gradient, the update step
     *
     * Default: 0.01
     */
    double step;

    /**
     * In exhaustive search, the range of the search
     *
     * Default: 0.1
     */
    double range;

    /**
     * In exhaustive search, the step within the range
     *
     * Default: 0.01
     */
    double epsilon;

    /**
     * The minimum peak distance that is allowed
     *
     * Default : 0.05
     */
    double min_peak_distance; 

    /**
     * The tolerance of residue increase when merging peaks
     *
     * Default: 2
     */
    double amend_peak_tolerance;

    /**
     * The size of the residue window
     *
     * Default: 10
     */
    uint32_t residue_window_size;

    /**
     * The "noticable change" threshold
     *
     * If change is smaller than that percentage, we think the algorithm
     * converges
     *
     * Default: 0.05
     */
    double noticable_change;

    /**
     * Whether the Graident::checkEst() is enabled or not
     */
    bool check_est;

    /**
     * The threshold of difference in power
     * This is used by Gradient::checkEst()
     *
     * Default: 0.5
     */
    double power_diff_threshold;

    /**
     * The maximal number of changes allowed in Gradient::checkEst()
     *
     * Default: 1
     */
    uint32_t max_change;

    /**
     * The fraction of largest bucket power, such that any bucket
     * whose power is larger than that is a reasonble bucket
     *
     * Default: 0.8
     */
    double reasonable_bucket_power;

    /**
     * The minimal peak distance when trying to insert or remove peak. 
     * Used in Gradient::checkEst()
     *
     * Default: 0.5
     */
    double checkest_min_distance;

    /**
     * The parameters for controlling
     */
    ControlParams control_params;

    /**
     * Whether we fix the time domain samples
     *
     * Default: 0
     */
    bool fix_time_samples;

    /**
     * The parameters for lightfield application
     */
    LightFieldParams lightfield_params;

    /**
     * If the solution is way larger than the maximum bucket, 
     * we say that it is singular. This is the threshold
     *
     * Default: 100
     */
    double scaling_threshold;

protected:

    /**
     * Read the configuration from the configuration file
     * @param[in]   param_file  the file to read parameters
     * @param[in]   io_prefix   the folder to which the io is dealing with
     *
     * @return  success or not 
     */
    bool readAllParams(in const char * param_file, in const char * io_prefix);

public:

    /**
     * The default constructor. 
     * It reads from the configuration file and sets all parameters
     *
     * @param[in]   param_file  the file to read parameters
     * @param[in]   io_prefix   the folder to which the io is dealing with
     */
    Params(in const char * param_file, in const char * io_prefix = NULL) : 
        n_v(17), n_h(17), num_projection_lines(0), projection_matrices(NULL), max_iter(100), 
        method(GRADIENT), subtract_and_shift_enabled(0), commit_step_by_step(0), tolerance(1e-5), delta(0.01), 
        step(0.01), range(0.1), epsilon(0.01), min_peak_distance(0.05), amend_peak_tolerance(2), residue_window_size(10), 
        noticable_change(0.05), check_est(true), 
        power_diff_threshold(0.5), max_change(1), reasonable_bucket_power(0.8), checkest_min_distance(0.5), 
        control_params(io_prefix), fix_time_samples(0), scaling_threshold(100), 
        meta_params(std::list<MetaParam>()) {

            //fprintf(stderr, "param_file = %s\n", param_file);

            // First, let's construct the meta_data data structures
            constructMetaParams();

            // Then, read the parameter file
            readAllParams(param_file, io_prefix);

            // Set the dimension related information in NonIntSFFT
            setDimension(n_v, n_h);

            /*
               setDimension(53, 53);

               num_filters = 5;
               ws = test_ws;
               init_shifts = test_init_offset;
               shifts = test_shifts;
               perm_matrices = test_mats;
               freq_offsets = test_freq_offset;

               lightfield_params.camera_pixel_consistent = 0;

               max_iter = 100;





               setDimension(48, 48);

               num_filters = 4;
               ws = test_ws;
               init_shifts = test_init_offset;
               shifts = test_shifts;
               perm_matrices = test_mats;
               freq_offsets = test_freq_offset;


               max_iter = 100;
               */

        };

    /**
     * Note that it is Params' responsibility to free sinc template memory space
     */
    ~Params() {
        Utils::deallocateSincTemplate();
    }


    /**
     * Set dimension n_v and n_h. 
     * It also allocates the sinc template memory space in Utils
     *
     * @param[in]   _n_v    the size of vertical dimension
     * @param[in]   _n_h    the size of the horizontal dimension
     */
    void setDimension(in double _n_v, in double _n_h) forceinline {
        NonIntSFFT::n_v = _n_v;
        NonIntSFFT::n_h = _n_h;
        NonIntSFFT::n = _n_v * _n_h;

        Utils::initUtils();
    } 

protected: 

    /**
     * @brief   meta data class for each parameter
     * 
     * @author  lixin
     * @date    11/23/2012
     */
    class MetaParam {

    public:
        /**
         * What type this parameter is 
         * Different type means different input format
         */
        enum TYPE {
            INT32_T,        /**< int32_t type */
            UINT32_T,       /**< uint32_t type */
            DOUBLE,         /**< double type */
            STRING,         /**< string type */
            BOOLEAN,        /**< boolean type */
            PEAKSLIST,      /**< peakslist type */
            INT32_ARRAY,    /**< array of int32_t type */
            UINT32_ARRAY,   /**< array of uint32_t type */
            DOUBLE_ARRAY,   /**< array of double type */
            ENUM_INPUT,     /**< enum INPUT_METHOD type */
            ENUM_INIT,      /**< enum INIT_METHOD type */
            ENUM_GRAD,      /**< enum GRAD_METHOD type */
            TOTAL_NUM_TYPES /**< total number of types */
        };

        /**
         * What group this parameter is in
         */
        enum GROUP {
            TOP_GROUP,          /**< the top group; i.e. no group */
            INIT_GROUP,         /**< the InitParams group */
            CONTROL_GROUP,      /**< the ControlParams group */
            LIGHTFIELD_GROUP,   /**< the LightFieldParams group */
            TOTAL_NUM_GROUPS    /**< the total number of groups */
        };

    protected:

        /**
         * The name of the groups
         */
        const static char GROUP_NAMES[TOTAL_NUM_GROUPS][MAX_STRING_LENGTH];

        /**
         * The names of the types 
         */
        const static char TYPE_NAMES[TOTAL_NUM_TYPES][MAX_STRING_LENGTH];

        /**
         *  The names of the types supported by libconfig
         */
        const static char LIBCONFIG_TYPE_NAMES[libconfig::Setting::TypeList + 1][MAX_STRING_LENGTH];

        /**
         * The type of this parameter
         */
        TYPE type;

        /**
         * The group of this parameter
         */
        GROUP group;

        /**
         * The pointer to the field that stores this paramter
         */
        void * storage;

        /**
         * The name of this param.
         */
        char name[MAX_STRING_LENGTH];

        /**
         * The root setting from libconfig
         */
        const libconfig::Setting * root;

        /**
         * The number of elements
         */
        uint32_t num_elements;

        /**
         * The names of the input method enums
         */
        static const char INPUT_ENUM_NAMES[Params::ControlParams::TOTAL_NUM_INPUT_METHODS][MAX_STRING_LENGTH];

        /**
         * The names of the init methed enums
         */
        static const char INIT_ENUM_NAMES[Params::InitParams::TOTAL_NUM_INIT_METHODS][MAX_STRING_LENGTH];

        /**
         * The names of the grad-search methods
         */
        static const char GRAD_ENUM_NAMES[Params::TOTAL_NUM_GRAD_METHOD][MAX_STRING_LENGTH];

        public:

        /**
         * Constructor: assign values to MetaParam::type and MetaParam::storage
         *
         * @param[in]   _type       the MetaParam::type
         * @param[in]   _storage    the MetaParam::storage
         * @param[in]   _name       the MetaParam::name
         * @param[in]   _group      the MetaParam::group
         */
        MetaParam(in TYPE _type, in void * _storage, in const char * _name, in GROUP _group) : type(_type), group(_group), storage(_storage), root(NULL), num_elements(1) {
            strcpy(name, _name); 
        }

        /**
         * Set the root setting
         *
         * @param[in]   _root   the root setting
         */
        inline void setRoot(in const libconfig::Setting * _root) forceinline {
            root = _root;
        }

        /**
         * Actually read the parameter
         * @param[in]   io_prefix   the folder to which the io is dealing with
         *
         * @return  whether reading the parameter is successful or not
         */
        bool readConfig(in const char * io_prefix);

        /**
         * The deconstructor: it frees any pointer dynamically allocated
         */
        ~MetaParam() {
            switch (type) {
                case INT32_ARRAY:
                    if (*((int32_t **)storage))
                        free(*((int32_t **)storage));
                    break;
                case UINT32_ARRAY:
                    if (*((uint32_t **)storage))
                        free(*((uint32_t **)storage));
                    break;
                default:
                    ;
            }

        }
    };

    /**
     * The array of meta parameters
     */
    std::list<MetaParam> meta_params;

    /**
     * Constructs the array of meta parameters meta_params
     */
    void constructMetaParams();

};


#endif
