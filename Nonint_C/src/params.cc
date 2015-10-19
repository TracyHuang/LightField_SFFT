#include "params.h"
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/stringize.hpp>

/**
 * @file    params.cc
 * @brief   Implementation of class Params
 * @author  lixin
 * @date    02/15/2013
 */

using namespace libconfig;

/**
 * The macro to generate one meta parameter data from top group
 */
#define META_TOP_ENTRY(unused, type, name)       \
    meta_params.push_back(MetaParam(MetaParam::type, &(name), BOOST_PP_STRINGIZE(name), MetaParam::TOP_GROUP));

/**
 * The macro to generate one meta parameter data from init_params group
 */
#define META_INIT_ENTRY(unused, type, name)       \
    meta_params.push_back(MetaParam(MetaParam::type, &(init_param.name), BOOST_PP_STRINGIZE(name), MetaParam::INIT_GROUP));

/**
 * The macro to generate one meta parameter data from control_params group
 */
#define META_CONTORL_ENTRY(unused, type, name)       \
    meta_params.push_back(MetaParam(MetaParam::type, &(control_params.name), BOOST_PP_STRINGIZE(name), MetaParam::CONTROL_GROUP));

/**
 * The macro to generate one meta parameter data from lightfield_params group
 */
#define META_LIGHTFIELD_ENTRY(unused, type, name)       \
    meta_params.push_back(MetaParam(MetaParam::type, &(lightfield_params.name), BOOST_PP_STRINGIZE(name), MetaParam::LIGHTFIELD_GROUP));

void Params::constructMetaParams() {

    // init_params
    BOOST_PP_SEQ_FOR_EACH(META_INIT_ENTRY, PEAKSLIST, (init_peaks))
    BOOST_PP_SEQ_FOR_EACH(META_INIT_ENTRY, UINT32_T, (num_init_peaks)(max_voting_iter))
    BOOST_PP_SEQ_FOR_EACH(META_INIT_ENTRY, ENUM_INIT, (method))
    BOOST_PP_SEQ_FOR_EACH(META_INIT_ENTRY, DOUBLE, (soft_power_threshold)(hard_power_threshold)(max_soft_power_threshold)(min_soft_threshold_step)(voting_threshold))
    BOOST_PP_SEQ_FOR_EACH(META_INIT_ENTRY, BOOLEAN, (optimism)(use_adjacent))
	BOOST_PP_SEQ_FOR_EACH(META_INIT_ENTRY, UINT32_T, (added_peaks_num_threshold))
	BOOST_PP_SEQ_FOR_EACH(META_INIT_ENTRY, DOUBLE, (min_shadow_power_threshold) (min_value_ratio_threshold))
    // control params
    BOOST_PP_SEQ_FOR_EACH(META_CONTORL_ENTRY, STRING, (result_file_prefix)(debugging_file_prefix)(init_file_prefix)(logging_file_prefix)(input_file_prefix)(full_file)(ground_truth_debugging)(verbos_file_prefix))
    BOOST_PP_SEQ_FOR_EACH(META_CONTORL_ENTRY, ENUM_INPUT, (input_method))
    BOOST_PP_SEQ_FOR_EACH(META_CONTORL_ENTRY, UINT32_T, (n3)(k))
    BOOST_PP_SEQ_FOR_EACH(META_CONTORL_ENTRY, DOUBLE, (snr)(min_distance))
    
    // lightfield_params
    BOOST_PP_SEQ_FOR_EACH(META_LIGHTFIELD_ENTRY, UINT32_T, (n_x)(n_y)(lam_thickness))
    BOOST_PP_SEQ_FOR_EACH(META_LIGHTFIELD_ENTRY, BOOLEAN, (camera_pixel_consistent)(impose_dim_gap))

    //BOOST_PP_SEQ_FOR_EACH(META_TOP_ENTRY, UINT32_T, (n_v)(n_h)(num_filters)(max_iter)(residue_window_size)(max_change))
    //BOOST_PP_SEQ_FOR_EACH(META_TOP_ENTRY, UINT32_ARRAY, (ws)(shifts)(perm_matrices))
    //BOOST_PP_SEQ_FOR_EACH(META_TOP_ENTRY, INT32_ARRAY, (init_shifts)(freq_offsets))
    BOOST_PP_SEQ_FOR_EACH(META_TOP_ENTRY, UINT32_T, (num_projection_lines))
    BOOST_PP_SEQ_FOR_EACH(META_TOP_ENTRY, INT32_ARRAY, (projection_matrices))
    BOOST_PP_SEQ_FOR_EACH(META_TOP_ENTRY, INT32_ARRAY, (sample_pos))
    //BOOST_PP_SEQ_FOR_EACH(META_TOP_ENTRY, DOUBLE_ARRAY, (sample_values))
    BOOST_PP_SEQ_FOR_EACH(META_TOP_ENTRY, UINT32_T, (num_samples))
    BOOST_PP_SEQ_FOR_EACH(META_TOP_ENTRY, ENUM_GRAD, (method))
    BOOST_PP_SEQ_FOR_EACH(META_TOP_ENTRY, BOOLEAN, (subtract_and_shift_enabled)(commit_step_by_step)(fix_time_samples)(check_est))
    BOOST_PP_SEQ_FOR_EACH(META_TOP_ENTRY, DOUBLE, 
        (tolerance)(delta)(step)(range)(epsilon)(min_peak_distance)(amend_peak_tolerance)(noticable_change)(power_diff_threshold)(reasonable_bucket_power)(checkest_min_distance)(scaling_threshold))
}

const char Params::MetaParam::INPUT_ENUM_NAMES[Params::ControlParams::TOTAL_NUM_INPUT_METHODS][MAX_STRING_LENGTH] = {

    "file_input",       // FILE_INPUT
    "random_generate"   // RANDOM_GENERATE
};

const char Params::MetaParam::INIT_ENUM_NAMES[Params::InitParams::TOTAL_NUM_INIT_METHODS][MAX_STRING_LENGTH] = {
    
    "voting",       // VOTING
    "ofdm_trick",   // OFDM_TRICK
    "manual"        // MANUAL
};

const char Params::MetaParam::GRAD_ENUM_NAMES[TOTAL_NUM_GRAD_METHOD][MAX_STRING_LENGTH] = {

    "gradient",         // GRADIENT
    "joint_gradient",   // JOINT_GRADIENT
    "exhaustive"        // EXHAUSTIVE
};

const char Params::MetaParam::GROUP_NAMES[Params::MetaParam::TOTAL_NUM_GROUPS][MAX_STRING_LENGTH] = {

    "",                     // TOP_GROUP
    "init_params",          // INIT_GROUP
    "control_params",       // CONTROL_GROUP
    "lightfield_params"     // LIGHTFIELD_GROUP
};

const char Params::MetaParam::TYPE_NAMES[Params::MetaParam::TOTAL_NUM_TYPES][MAX_STRING_LENGTH] = {

    "integer",              // INT32_T
    "unsigned int",         // UINT32_T
    "double",               // DOUBLE
    "string",               // STRING
    "boolean",              // BOOLEAN
    "peakslist",            // PEAKSLIST
    "integer array",        // INT32_ARRAY
    "unsigned int array",   // UINT32_ARRAY
    "enum INPUT_METHOD",    // ENUM_INPUT
    "enum INIT_METHOD",     // ENUM_INIT
    "enum GRAD_METHOD"      // ENUM_GRAD
};

const char Params::MetaParam::LIBCONFIG_TYPE_NAMES[Setting::TypeList + 1][MAX_STRING_LENGTH] = {

    "none",         // TypeNone
    "Int",          // TypeInt
    "Int64",        // TypeInt64
    "Float",        // TypeFloat
    "String",       // TypeString
    "Boolean",      // TypeBoolean
    "Group",        // TypeGroup
    "Array",        // TypeArray
    "List"          // TypeList
};


bool Params::MetaParam::readConfig(in const char * io_prefix) {

    const Setting * setting = NULL;

    const Setting & root_set = *root;

    try {
        switch (group) {

        case INIT_GROUP:
            //setting = &root_set[GROUP_NAMES[INIT_GROUP]][(const char *)name];
            //break;
        case CONTROL_GROUP:
            //setting = &root_set[GROUP_NAMES[CONTROL_GROUP]][(const char *)name];
            //break;
        case LIGHTFIELD_GROUP:
            //setting = &root_set[GROUP_NAMES[LIGHTFIELD_GROUP]][(const char *)name];
            setting = &root_set[GROUP_NAMES[group]][(const char *)name];
            break;
        case TOP_GROUP:
            setting = &root_set[(const char *)name];
            break;
        default:
            return false;
        }
    } catch(const SettingNotFoundException) { }

    // Now convert the content back to what it should be
    double x, y = 0;
    int int_x, int_y = 0;
    uint32_t actual_size = 0;
    const char * enum_name = NULL;
    bool matched = false;
    bool getXY = true;
    PeaksList::iterator it;
    try {
        switch (type) {
        case INT32_T:
            *((int32_t *) storage) = (int32_t)*setting;
            break;
        case UINT32_T:
            *((uint32_t *) storage) = (uint32_t)*setting;
            break;
        case DOUBLE: 
            if (setting->getType() == Setting::TypeString) {
                if (!strncmp((const char *)(*setting), "inf", MAX_STRING_LENGTH))
                    *((double *)storage) = INFINITY;
                else if (!strncmp((const char*)(*setting), "-inf", MAX_STRING_LENGTH))
                    *((double *)storage) = -INFINITY;
            } else if (setting->getType() == Setting::TypeInt)
                *((double *)storage) = (double)((int)*setting);
            else
                *((double *) storage) = (double)*setting;
            break;
        case STRING: 
            // Note that we need to cast as char(*)[MAX_STRING_LENGTH] instead of char **
            // otherwise there might be problems with aliasing
            if (!io_prefix || strlen(io_prefix) == 0)
                strncpy(*((char (*)[MAX_STRING_LENGTH])storage), (const char *)(*setting), MAX_STRING_LENGTH);
            else
                snprintf(*((char (*)[MAX_STRING_LENGTH])storage), MAX_STRING_LENGTH, "%s/%s", io_prefix, (const char *)(*setting));
            break;
        case BOOLEAN:
            *((bool *) storage) = (bool)*setting;
            break;
        case PEAKSLIST:
            num_elements = setting->getLength();
            ((PeaksList *)storage)->resize(num_elements);
            it = ((PeaksList *)storage)->begin();
            for (uint32_t i = 0; i < num_elements; ++ i) {
                if (!(*setting)[i].lookupValue("x", x)) {
                    if (!(*setting)[i].lookupValue("x", int_x)) {
                        getXY = false;
                    } else 
                        x = int_x;
                }
                if( !(*setting)[i].lookupValue("y", y)) {
                    if (!(*setting)[i].lookupValue("y", int_y)) {
                        getXY = false;
                    } else
                        y = int_y;
                } 
                if (getXY) {
                    it->setXY(x, y);
                    ++ it;
                    ++ actual_size;
                }
            }
            num_elements = actual_size;
            ((PeaksList *)storage)->resize(actual_size);
            //*((PeaksList *)storage) = peaks_list_ptr;
            break;
        case INT32_ARRAY:
            num_elements = setting->getLength();
            *((int32_t **)storage) = (int32_t *)malloc(sizeof(int32_t) * num_elements); 
            for (uint32_t i = 0; i < num_elements; ++ i) {
                *((*(int32_t **)storage) + i) = (*setting)[i];
            }
            break;
        case UINT32_ARRAY:
            num_elements = setting->getLength();
            *((uint32_t **)storage) = (uint32_t *)malloc(sizeof(uint32_t) * num_elements); 
            for (uint32_t i = 0; i < num_elements; ++ i) {
                *((*(uint32_t **)storage) + i) = (*setting)[i];
            }
            break;
        case ENUM_INPUT:
            enum_name = (const char *)(*setting);
             for (uint32_t i = 0; i < Params::ControlParams::TOTAL_NUM_INPUT_METHODS; ++ i) {
                if (!strncmp(enum_name, INPUT_ENUM_NAMES[i], MAX_STRING_LENGTH)) {
                    *((uint32_t *)storage) = i;
                    matched = true;
                    break;
                }
             }
             if (!matched) {
                return false;
             }
            break;
        case ENUM_INIT:
            enum_name = (const char *)(*setting);
             for (uint32_t i = 0; i < Params::InitParams::TOTAL_NUM_INIT_METHODS; ++ i) {
                if (!strncmp(enum_name, INIT_ENUM_NAMES[i], MAX_STRING_LENGTH)) {
                    *((uint32_t *)storage) = i;
                    matched = true;
                    break;
                }
             }
             if (!matched) {
                return false;
             }
            break;
        case ENUM_GRAD:
            enum_name = (const char *)(*setting);
             for (uint32_t i = 0; i < Params::TOTAL_NUM_GRAD_METHOD; ++ i) {
                if (!strncmp(enum_name, GRAD_ENUM_NAMES[i], MAX_STRING_LENGTH)) {
                    *((uint32_t *)storage) = i;
                    matched = true;
                    break;
                }
             }
             if (!matched) {
                return false;
             }
            break;
        default:
            ;
        }
    } catch(SettingTypeException type_err) {
        return false;
    }

    return true;
}

bool Params::readAllParams(in const char * param_file, in const char * io_prefix) {


    Config cfg;

    // Read the file. If there is an error, report it and exit.
    try {
        cfg.readFile(param_file);
    } catch(FileIOException &fioex) {
        return false;
    } catch(ParseException &pex) {
        return false;
    } 

    const Setting& root = cfg.getRoot();

    // Now we can process all meta data 
    for (std::list<MetaParam>::iterator it = meta_params.begin(); it != meta_params.end(); ++ it) {
        it->setRoot(&root);
        it->readConfig(io_prefix);
    }

    // This is sanity check, mainly check the number of elements is equal to the length of array
    bool passed = true;
    if (init_param.num_init_peaks != init_param.init_peaks.size()) {
        passed = false;
    }
    return passed;
}

