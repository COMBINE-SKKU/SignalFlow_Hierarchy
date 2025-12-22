/*
===============================================================================
BEI MODEL - Brain E-to-E Inhibition Model
===============================================================================

Simplified brain network model with:
1. E-to-E connections only (no FFI)
2. FIC control with burn-in phase
3. Variable Balloon-Windkessel parameters across regions
4. Convergence-driven burn-in followed by clean data collection

Key Components:
- Neural dynamics: DMF equations with E/I populations
- BOLD generation: Region-specific Balloon-Windkessel model
- FIC: Vogels plasticity rule for local inhibition (3Hz target)
- Burn-in: Automatic convergence detection

Usage: ./bei_model <param_file> <subject_id> <random_seed> <results_dir> [bw_file] [bw_variability]

===============================================================================
*/

#include <stdio.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <omp.h>

// ============================================================================
// DATA STRUCTURES
// ============================================================================

struct ConnectionStrengths {
    float *weights;
}; 

struct ConnectionSources {
    int *source_regions;
};

struct BWParams {
    float *rho;
    float *alpha;
    float *tau;
    float *kappa;
    float *gamma;
    float *V0;
    float *k1;
    float *k2;
    float *k3;
};


// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================


static inline void generate_gaussian_noise(float *random_numbers) {
    static double V1 = 0.0, V2 = 0.0, S = 0.0, U1 = 0.0, U2 = 0.0;
    
    S = 0.0;
    do {
        U1 = (double)rand() / RAND_MAX;
        U2 = (double)rand() / RAND_MAX;
        V1 = 2 * U1 - 1;
        V2 = 2 * U2 - 1;
        S = V1 * V1 + V2 * V2;
    } while(S >= 1 || S == 0);
    
    random_numbers[0] = (float)(V1 * sqrt(-2 * log(S) / S));
    random_numbers[1] = (float)(V2 * sqrt(-2 * log(S) / S));
    
    S = 0.0;
    do {
        U1 = (double)rand() / RAND_MAX;
        U2 = (double)rand() / RAND_MAX;
        V1 = 2 * U1 - 1;
        V2 = 2 * U2 - 1;
        S = V1 * V1 + V2 * V2;
    } while(S >= 1 || S == 0);
    
    random_numbers[2] = (float)(V1 * sqrt(-2 * log(S) / S));
    random_numbers[3] = (float)(V2 * sqrt(-2 * log(S) / S));
}

// ============================================================================
// BW PARAMETER LOADING
// ============================================================================

int load_bw_parameters(char *bw_file, int num_regions, struct BWParams *bw_params, float variability) {
    // Try to load from file first
    FILE *file = NULL;
    if (bw_file != NULL) {
        file = fopen(bw_file, "r");
    }
    
    // Default values
    float default_rho = 0.34, default_alpha = 0.32, default_tau = 0.98;
    float default_kappa = 1.0/0.65, default_gamma = 1.0/0.41, default_V0 = 0.02;
    
    if (file) {
        printf("Loading BW parameters from file: %s\n", bw_file);
        int file_regions;
        if (fscanf(file, "%d", &file_regions) != 1 || file_regions != num_regions) {
            printf("ERROR: BW file region count mismatch\n");
            fclose(file);
            return -1;
        }
        
        for (int i = 0; i < num_regions; i++) {
            int region_id;
            if (fscanf(file, "%d %f %f %f %f %f %f", &region_id,
                      &bw_params->rho[i], &bw_params->alpha[i], &bw_params->tau[i],
                      &bw_params->kappa[i], &bw_params->gamma[i], &bw_params->V0[i]) != 7) {
                printf("ERROR: Could not read BW parameters for region %d\n", i);
                fclose(file);
                return -1;
            }
        }
        fclose(file);
        printf("BW parameters loaded from file\n");
    } else {
        // Generate with variability
        printf("Generating BW parameters with %.1f%% variability\n", variability * 100);
        
        for (int i = 0; i < num_regions; i++) {
            bw_params->rho[i] = default_rho * (1.0f + variability * ((float)rand()/RAND_MAX - 0.5f) * 2.0f);
            bw_params->alpha[i] = default_alpha * (1.0f + variability * ((float)rand()/RAND_MAX - 0.5f) * 2.0f);
            bw_params->tau[i] = default_tau * (1.0f + variability * ((float)rand()/RAND_MAX - 0.5f) * 2.0f);
            bw_params->kappa[i] = default_kappa * (1.0f + variability * ((float)rand()/RAND_MAX - 0.5f) * 2.0f);
            bw_params->gamma[i] = default_gamma * (1.0f + variability * ((float)rand()/RAND_MAX - 0.5f) * 2.0f);
            bw_params->V0[i] = default_V0 * (1.0f + variability * ((float)rand()/RAND_MAX - 0.5f) * 2.0f);
            
            // Ensure positive values
            bw_params->rho[i] = fmaxf(0.1f, bw_params->rho[i]);
            bw_params->alpha[i] = fmaxf(0.1f, bw_params->alpha[i]);
            bw_params->tau[i] = fmaxf(0.1f, bw_params->tau[i]);
            bw_params->kappa[i] = fmaxf(0.1f, bw_params->kappa[i]);
            bw_params->gamma[i] = fmaxf(0.1f, bw_params->gamma[i]);
            bw_params->V0[i] = fmaxf(0.001f, bw_params->V0[i]);
        }
    }
    
    // Calculate derived parameters
    for (int i = 0; i < num_regions; i++) {
        bw_params->k1[i] = 7.0f * bw_params->rho[i];
        bw_params->k2[i] = 2.0f;
        bw_params->k3[i] = 2.0f * bw_params->rho[i] - 0.2f;
    }
    
    return 0;
}

// ============================================================================
// STRUCTURAL DATA LOADING
// ============================================================================

void load_structural_data(char *structural_connectivity_file, char *region_file, int num_regions,
                         int **connection_count_per_region, float coupling_strength, int **region_index_lookup,
                         struct ConnectionStrengths **current_weights,
                         struct ConnectionSources **connection_sources) {
    
    FILE *sc_file = fopen(structural_connectivity_file, "r");
    FILE *region_file_handle = fopen(region_file, "r");
    
    if (!sc_file || !region_file_handle) {
        printf("ERROR: Could not open structural connectivity files\n");
        exit(1);
    }
    
    int header_sc, header_region;
    if (fscanf(sc_file, "%d", &header_sc) == EOF || 
        fscanf(region_file_handle, "%d", &header_region) == EOF) {
        printf("ERROR: Could not read file headers\n");
        exit(1);
    }
    
    if (header_sc != num_regions || header_region != num_regions) {
        printf("ERROR: Inconsistent number of regions across files\n");
        exit(1);
    }
    
    // Allocate memory
    *connection_count_per_region = (int *)malloc(num_regions * sizeof(int));
    *region_index_lookup = (int *)malloc(num_regions * sizeof(int));
    
    *current_weights = (struct ConnectionStrengths *)malloc(num_regions * sizeof(struct ConnectionStrengths));
    *connection_sources = (struct ConnectionSources *)malloc(num_regions * sizeof(struct ConnectionSources));
    
    if (!*connection_count_per_region) {
        printf("ERROR: Memory allocation failed\n");
        exit(1);
    }
    
    // Load structural connectivity data
    for (int i = 0; i < num_regions; i++) {
        int region_idx, connection_count;
        if (fscanf(region_file_handle, "%d %d", &region_idx, &connection_count) != 2) {
            printf("ERROR: Could not read region %d header\n", i);
            exit(1);
        }
        
        if (region_idx != i) {
            printf("ERROR: Region index mismatch at region %d\n", i);
            exit(1);
        }
        
        int sc_region_idx, sc_connection_count;
        if (fscanf(sc_file, "%d %d", &sc_region_idx, &sc_connection_count) != 2) {
            printf("ERROR: Could not read SC header for region %d\n", i);
            exit(1);
        }
        
        if (sc_region_idx != i || sc_connection_count != connection_count) {
            printf("ERROR: Connection count mismatch for region %d\n", i);
            exit(1);
        }
        
        (*connection_count_per_region)[i] = connection_count;
        (*region_index_lookup)[i] = i;
        
        // Allocate connection arrays
        (*current_weights)[i].weights = (float *)malloc(connection_count * sizeof(float));
        (*connection_sources)[i].source_regions = (int *)malloc(connection_count * sizeof(int));
        
        if (!(*current_weights)[i].weights || !(*connection_sources)[i].source_regions) {
            printf("ERROR: Memory allocation failed for region %d\n", i);
            exit(1);
        }
        
        // Read connection data
        for (int j = 0; j < connection_count; j++) {
            int source_region;
            float connection_strength;
            
            if (fscanf(region_file_handle, "%d", &source_region) != 1 ||
                fscanf(sc_file, "%f", &connection_strength) != 1) {
                printf("ERROR: Could not read connection data for region %d\n", i);
                exit(1);
            }
            
            (*connection_sources)[i].source_regions[j] = source_region;
            (*current_weights)[i].weights[j] = coupling_strength * connection_strength;
        }
    }
    
    fclose(sc_file);
    fclose(region_file_handle);
}

// ============================================================================
// MEMORY CLEANUP
// ============================================================================

void cleanup_data(int num_regions, int *connections_per_region,
                  int *region_lookup,
                  struct ConnectionStrengths *current_weights,
                  struct ConnectionSources *connection_sources, struct BWParams *bw_params,
                  float *bold_timeseries, float *excitatory_synaptic_activity,
                  float *inhibitory_synaptic_activity, float *excitatory_firing_rate,
                  float *inhibitory_firing_rate, float *long_range_excitatory_input,
                  float *feedback_inhibition_strength, float *mean_excitatory_firing_rate,
                  float *mean_inhibitory_firing_rate, float *bw_vasodilatory_signal,
                  float *bw_blood_flow, float *bw_blood_volume, float *bw_deoxyhemoglobin,
                  float *bold_timeseries_ctrl, float *bw_vasodilatory_signal_ctrl,
                  float *bw_blood_flow_ctrl, float *bw_blood_volume_ctrl, float *bw_deoxyhemoglobin_ctrl) {
    
    // Free per-region arrays
    if (current_weights && connection_sources) {
        for (int i = 0; i < num_regions; i++) {
            if (connections_per_region && connections_per_region[i] > 0) {
                free(current_weights[i].weights);
                free(connection_sources[i].source_regions);
            }
        }
    }
    
    // Free structural arrays
    if (connections_per_region) free(connections_per_region);
    if (region_lookup) free(region_lookup);
    if (current_weights) free(current_weights);
    if (connection_sources) free(connection_sources);
    
    // Free BW parameter arrays
    if (bw_params) {
        free(bw_params->rho);
        free(bw_params->alpha);
        free(bw_params->tau);
        free(bw_params->kappa);
        free(bw_params->gamma);
        free(bw_params->V0);
        free(bw_params->k1);
        free(bw_params->k2);
        free(bw_params->k3);
    }
    
    // Free neural state arrays
    if (excitatory_synaptic_activity) _mm_free(excitatory_synaptic_activity);
    if (inhibitory_synaptic_activity) _mm_free(inhibitory_synaptic_activity);
    if (excitatory_firing_rate) _mm_free(excitatory_firing_rate);
    if (inhibitory_firing_rate) _mm_free(inhibitory_firing_rate);
    if (long_range_excitatory_input) _mm_free(long_range_excitatory_input);
    if (feedback_inhibition_strength) _mm_free(feedback_inhibition_strength);
    if (mean_excitatory_firing_rate) _mm_free(mean_excitatory_firing_rate);
    if (mean_inhibitory_firing_rate) _mm_free(mean_inhibitory_firing_rate);
    
    // Free BW state arrays
    if (bw_vasodilatory_signal) _mm_free(bw_vasodilatory_signal);
    if (bw_blood_flow) _mm_free(bw_blood_flow);
    if (bw_blood_volume) _mm_free(bw_blood_volume);
    if (bw_deoxyhemoglobin) _mm_free(bw_deoxyhemoglobin);
    
    // Free control BW state arrays
    if (bw_vasodilatory_signal_ctrl) _mm_free(bw_vasodilatory_signal_ctrl);
    if (bw_blood_flow_ctrl) _mm_free(bw_blood_flow_ctrl);
    if (bw_blood_volume_ctrl) _mm_free(bw_blood_volume_ctrl);
    if (bw_deoxyhemoglobin_ctrl) _mm_free(bw_deoxyhemoglobin_ctrl);
    
    // Free BOLD timeseries
    if (bold_timeseries) _mm_free(bold_timeseries);
    if (bold_timeseries_ctrl) _mm_free(bold_timeseries_ctrl);
}

// ============================================================================
// MAIN SIMULATION
// ============================================================================

int main(int argc, char *argv[]) {
    time_t simulation_start_time = time(NULL);
    
    // ========================================================================
    // COMMAND LINE ARGUMENTS
    // ========================================================================
    if (argc < 5 || argc > 7) {
        printf("Usage: %s <param_file> <subject_id> <random_seed> <results_dir> [bw_file] [bw_variability]\n", argv[0]);
        exit(1);
    }
    
    char *parameter_file = argv[1];
    char *subject_id = argv[2];
    int random_seed = atoi(argv[3]);
    char *results_dir = argv[4];
    char *bw_file = (argc > 5) ? argv[5] : NULL;
    float bw_variability = (argc > 6) ? atof(argv[6]) : 0.15; // Default 15% variability
    
    printf("===============================================================================\n");
    printf("BEI MODEL - Brain E-to-E Inhibition Model\n");
    printf("===============================================================================\n");
    printf("Subject: %s\n", subject_id);
    printf("Random seed: %d\n", random_seed);
    printf("BW variability: %.1f%%\n", bw_variability * 100);
    if (bw_file) {
        printf("BW parameter file: %s\n", bw_file);
    }
    
    // ========================================================================
    // MODEL PARAMETERS
    // ========================================================================
    
    const float integration_timestep = 1.0;
    const float model_sampling_rate = 0.001;
    const int vectorization_width = 4;
    
    int total_simulation_timesteps;
    int bold_repetition_time = 1000;
    int actual_regions = 40;
    int padded_regions = 40;
    float global_coupling_strength = 0.5;
    float excitatory_coupling_strength = 0.15;
    
    // Neural parameters
    const float exc_gain_factor = 310;
    const float exc_bias_current = 125;
    const float exc_membrane_time = 0.16;
    const float inh_gain_factor = 615;
    const float inh_bias_current = 177;
    const float inh_membrane_time = 0.087;
    const float excitatory_time_constant = 100;
    const float inhibitory_time_constant = 10;
    const float exc_kinetic_factor = 0.641/1000.0;
    const float inh_kinetic_factor = 1.0/1000.0;
    float noise_amplitude = 0.00316228;
    const float external_input_baseline = 0.382;
    const float exc_external_scaling = 1.0;
    const float inh_external_scaling = 0.7;
    float initial_inhibitory_strength = 0.0;
    
    // FIC parameters
    const float fic_learning_rate = 0.01;
    const float target_excitatory_rate = 3.0; // 3Hz target
    const float convergence_threshold = 0.5; // 0.5 Hz tolerance
    const float fic_change_threshold = 0.005; // 0.5% change
    const int required_stable_steps = 10;
    const int max_burnin_time = 1000000; // 5 minutes max
    
    // Calculate simulation time for exactly 490 BOLD timepoints
    int target_bold_points = 490;
    total_simulation_timesteps = target_bold_points * bold_repetition_time;
    
    // Control (homogeneous) BW parameters
    const float default_rho = 0.34f;
    const float default_alpha = 0.32f;
    const float default_tau = 0.98f;
    const float default_kappa = 1.0f/0.65f;
    const float default_gamma = 1.0f/0.41f;
    const float default_V0 = 0.02f;
    
    // Derived constants for control model
    const float default_k1 = 7.0f * default_rho;
    const float default_k2 = 2.0f * default_rho;
    const float default_k3 = 2.0f * default_rho - 0.2f;
    
    // ========================================================================
    // LOAD PARAMETERS FROM FILE
    // ========================================================================
    
    FILE *param_file = fopen(parameter_file, "r");
    if (!param_file) {
        printf("ERROR: Could not open parameter file %s\n", parameter_file);
        exit(1);
    }
    
    float local_excitatory_recurrence = 1.5f;
    int params_read = fscanf(param_file, "%d %f %f %f %f %d %d", 
               &actual_regions, &global_coupling_strength, &excitatory_coupling_strength,
               &initial_inhibitory_strength, &noise_amplitude,
               &total_simulation_timesteps, &bold_repetition_time);
    
    if (params_read < 7) {
        printf("ERROR: Could not read parameters (read %d/7)\n", params_read);
        exit(1);
    }
    fclose(param_file);
    
    // Override simulation time to get exactly 490 BOLD points
    total_simulation_timesteps = target_bold_points * bold_repetition_time;
    
    srand((unsigned)random_seed);
    
    printf("Regions: %d\n", actual_regions);
    printf("Target BOLD timepoints: %d\n", target_bold_points);
    printf("Simulation time: %.1f minutes\n", total_simulation_timesteps/60000.0);
    
    // ========================================================================
    // ALLOCATE BW PARAMETERS
    // ========================================================================
    
    struct BWParams bw_params;
    bw_params.rho = (float *)malloc(actual_regions * sizeof(float));
    bw_params.alpha = (float *)malloc(actual_regions * sizeof(float));
    bw_params.tau = (float *)malloc(actual_regions * sizeof(float));
    bw_params.kappa = (float *)malloc(actual_regions * sizeof(float));
    bw_params.gamma = (float *)malloc(actual_regions * sizeof(float));
    bw_params.V0 = (float *)malloc(actual_regions * sizeof(float));
    bw_params.k1 = (float *)malloc(actual_regions * sizeof(float));
    bw_params.k2 = (float *)malloc(actual_regions * sizeof(float));
    bw_params.k3 = (float *)malloc(actual_regions * sizeof(float));
    
    if (load_bw_parameters(bw_file, actual_regions, &bw_params, bw_variability) != 0) {
        printf("ERROR: Failed to load BW parameters\n");
        exit(1);
    }
    
    // ========================================================================
    // ALLOCATE NEURAL VARIABLES
    // ========================================================================
    
    float *excitatory_synaptic_activity = (float *)_mm_malloc(padded_regions * sizeof(float), 16);
    float *inhibitory_synaptic_activity = (float *)_mm_malloc(padded_regions * sizeof(float), 16);
    float *excitatory_firing_rate = (float *)_mm_malloc(padded_regions * sizeof(float), 16);
    float *inhibitory_firing_rate = (float *)_mm_malloc(padded_regions * sizeof(float), 16);
    float *long_range_excitatory_input = (float *)_mm_malloc(padded_regions * sizeof(float), 16);
    float *feedback_inhibition_strength = (float *)_mm_malloc(padded_regions * sizeof(float), 16);
    float *mean_excitatory_firing_rate = (float *)_mm_malloc(padded_regions * sizeof(float), 16);
    float *mean_inhibitory_firing_rate = (float *)_mm_malloc(padded_regions * sizeof(float), 16);
    
    // BOLD variables
    int bold_timeseries_length = target_bold_points;
    float *bw_vasodilatory_signal = (float *)_mm_malloc(padded_regions * sizeof(float), 16);
    float *bw_blood_flow = (float *)_mm_malloc(padded_regions * sizeof(float), 16);
    float *bw_blood_volume = (float *)_mm_malloc(padded_regions * sizeof(float), 16);
    float *bw_deoxyhemoglobin = (float *)_mm_malloc(padded_regions * sizeof(float), 16);
    float *bold_timeseries = (float *)_mm_malloc(actual_regions * bold_timeseries_length * sizeof(float), 16);
    
    // Control BOLD variables (homogeneous parameters)
    float *bw_vasodilatory_signal_ctrl = (float *)_mm_malloc(padded_regions * sizeof(float), 16);
    float *bw_blood_flow_ctrl = (float *)_mm_malloc(padded_regions * sizeof(float), 16);
    float *bw_blood_volume_ctrl = (float *)_mm_malloc(padded_regions * sizeof(float), 16);
    float *bw_deoxyhemoglobin_ctrl = (float *)_mm_malloc(padded_regions * sizeof(float), 16);
    float *bold_timeseries_ctrl = (float *)_mm_malloc(actual_regions * bold_timeseries_length * sizeof(float), 16);
    
    if (!excitatory_synaptic_activity || !inhibitory_synaptic_activity || !bold_timeseries || !bold_timeseries_ctrl) {
        printf("ERROR: Memory allocation failed\n");
        exit(2);
    }
    
    // ========================================================================
    // LOAD STRUCTURAL CONNECTIVITY
    // ========================================================================
    
    char sc_strengths_file[256], sc_regions_file[256];
    snprintf(sc_strengths_file, sizeof(sc_strengths_file), "../input/FLN40_SC_strengths.txt");
    snprintf(sc_regions_file, sizeof(sc_regions_file), "../input/FLN40_SC_regionids.txt");
    
    int *connections_per_region, *region_lookup_table;
    struct ConnectionStrengths *current_weights;
    struct ConnectionSources *connection_sources;
    
    load_structural_data(sc_strengths_file, sc_regions_file, actual_regions,
                         &connections_per_region, global_coupling_strength * excitatory_coupling_strength, 
                         &region_lookup_table, &current_weights, 
                         &connection_sources);
    
    float local_recurrent_strength = local_excitatory_recurrence * excitatory_coupling_strength;
    
    // ========================================================================
    // INITIALIZE NEURAL STATES
    // ========================================================================
    
    for (int j = 0; j < padded_regions; j++) {
        excitatory_synaptic_activity[j] = 0.001;
        inhibitory_synaptic_activity[j] = 0.001;
        long_range_excitatory_input[j] = 0.001;
        mean_excitatory_firing_rate[j] = 0.0f;
        mean_inhibitory_firing_rate[j] = 0.0f;
        feedback_inhibition_strength[j] = initial_inhibitory_strength;
        
        bw_vasodilatory_signal[j] = 0.0;
        bw_blood_flow[j] = 1.0;
        bw_blood_volume[j] = 1.0;
        bw_deoxyhemoglobin[j] = 1.0;
        
        // Initialize control BOLD states
        bw_vasodilatory_signal_ctrl[j] = 0.0;
        bw_blood_flow_ctrl[j] = 1.0;
        bw_blood_volume_ctrl[j] = 1.0;
        bw_deoxyhemoglobin_ctrl[j] = 1.0;
    }
    
    // ========================================================================
    // VECTORIZATION SETUP
    // ========================================================================
    
    const int vectorized_region_count = padded_regions / vectorization_width;
    const float negative_exc_membrane_time = -1.0 * exc_membrane_time;
    const float negative_inh_membrane_time = -1.0 * inh_membrane_time;
    const float negative_inv_exc_time_const = -1.0 / excitatory_time_constant;
    const float negative_inv_inh_time_const = -1.0 / inhibitory_time_constant;
    const float scaled_external_exc_input = exc_external_scaling * external_input_baseline;
    const float scaled_external_inh_input = inh_external_scaling * external_input_baseline;
    const float scaled_noise = noise_amplitude * integration_timestep;
    const float one = 1.0;
    
    // SIMD vectors
    const __m128 vec_dt = _mm_load1_ps(&integration_timestep);
    const __m128 vec_local_recurrent = _mm_load1_ps(&local_recurrent_strength);
    const __m128 vec_exc_gain = _mm_load1_ps(&exc_gain_factor);
    const __m128 vec_exc_bias = _mm_load1_ps(&exc_bias_current);
    const __m128 vec_neg_exc_membrane = _mm_load1_ps(&negative_exc_membrane_time);
    const __m128 vec_inh_gain = _mm_load1_ps(&inh_gain_factor);
    const __m128 vec_inh_bias = _mm_load1_ps(&inh_bias_current);
    const __m128 vec_neg_inh_membrane = _mm_load1_ps(&negative_inh_membrane_time);
    const __m128 vec_exc_kinetic = _mm_load1_ps(&exc_kinetic_factor);
    const __m128 vec_inh_kinetic = _mm_load1_ps(&inh_kinetic_factor);
    const __m128 vec_neg_inv_exc_tau = _mm_load1_ps(&negative_inv_exc_time_const);
    const __m128 vec_neg_inv_inh_tau = _mm_load1_ps(&negative_inv_inh_time_const);
    const __m128 vec_ext_exc_input = _mm_load1_ps(&scaled_external_exc_input);
    const __m128 vec_ext_inh_input = _mm_load1_ps(&scaled_external_inh_input);
    const __m128 vec_noise_scale = _mm_load1_ps(&scaled_noise);
    const __m128 vec_one = _mm_load1_ps(&one);
    const __m128 vec_nmda_coupling = _mm_load1_ps(&excitatory_coupling_strength);
    
    __m128 *vec_exc_synaptic = (__m128*)excitatory_synaptic_activity;
    __m128 *vec_inh_synaptic = (__m128*)inhibitory_synaptic_activity;
    __m128 *vec_exc_firing = (__m128*)excitatory_firing_rate;
    __m128 *vec_inh_firing = (__m128*)inhibitory_firing_rate;
    __m128 *vec_lre_input = (__m128*)long_range_excitatory_input;
    __m128 *vec_feedback_inhibition = (__m128*)feedback_inhibition_strength;
    __m128 *vec_mean_exc_fr = (__m128*)mean_excitatory_firing_rate;
    __m128 *vec_mean_inh_fr = (__m128*)mean_inhibitory_firing_rate;
    
    float temp_exponential_exc[4] __attribute__((aligned(16)));
    float temp_exponential_inh[4] __attribute__((aligned(16)));
    float random_noise[4] __attribute__((aligned(16)));
    __m128 *vec_temp_exp_exc = (__m128*)temp_exponential_exc;
    __m128 *vec_temp_exp_inh = (__m128*)temp_exponential_inh;
    __m128 *vec_random_noise = (__m128*)random_noise;
    __m128 temp_current_exc, temp_current_inh, temp_transfer_exc, temp_transfer_inh;
    
    // ========================================================================
    // BURN-IN PHASE VARIABLES
    // ========================================================================
    
    int firing_rate_sample_count = 0;
    int bold_timestep_counter = -1;
    int bold_length = -1;
    int consecutive_stable_steps = 0;
    int burnin_complete = 0;
    int burnin_timesteps = 0;
    float *previous_J_I = (float *)malloc(actual_regions * sizeof(float));
    
    printf("\n===============================================================================\n");
    printf("STARTING BURN-IN PHASE\n");
    printf("Target firing rate: %.1f Hz\n", target_excitatory_rate);
    printf("Convergence threshold: +/- %.1f Hz\n", convergence_threshold);
    printf("===============================================================================\n");
    
    // ========================================================================
    // MAIN SIMULATION LOOP
    // ========================================================================
    
    int timestep = 0;
    while (1) {
        
        // STEP 1: COMPUTE LONG-RANGE E-TO-E INPUTS ONLY (NO FFI)
        for (int target_region = 0; target_region < actual_regions; target_region++) {
            float total_lre_input = 0;
            
            for (int connection_idx = 0; connection_idx < connections_per_region[target_region]; connection_idx++) {
                int source_region = connection_sources[target_region].source_regions[connection_idx];
                float source_activity = excitatory_synaptic_activity[source_region];
                total_lre_input += source_activity * current_weights[target_region].weights[connection_idx];
            }
            
            long_range_excitatory_input[target_region] = total_lre_input;
        }
        
        // STEP 2: NEURAL DYNAMICS (VECTORIZED)
        for (int vector_idx = 0; vector_idx < vectorized_region_count; vector_idx++) {
            
            // Excitatory population
            temp_current_exc = _mm_sub_ps(
                _mm_mul_ps(vec_exc_gain, 
                    _mm_add_ps(_mm_add_ps(vec_ext_exc_input, 
                                         _mm_mul_ps(vec_local_recurrent, vec_exc_synaptic[vector_idx])),
                              _mm_sub_ps(vec_lre_input[vector_idx], 
                                        _mm_mul_ps(vec_feedback_inhibition[vector_idx], vec_inh_synaptic[vector_idx])))),
                vec_exc_bias);
            
            *vec_temp_exp_exc = _mm_mul_ps(vec_neg_exc_membrane, temp_current_exc);
            temp_exponential_exc[0] = expf(temp_exponential_exc[0]);
            temp_exponential_exc[1] = expf(temp_exponential_exc[1]);
            temp_exponential_exc[2] = expf(temp_exponential_exc[2]);
            temp_exponential_exc[3] = expf(temp_exponential_exc[3]);
            
            *vec_temp_exp_exc = _mm_load_ps(temp_exponential_exc);
            temp_transfer_exc = _mm_div_ps(temp_current_exc, _mm_sub_ps(vec_one, *vec_temp_exp_exc));
            
            vec_exc_firing[vector_idx] = temp_transfer_exc;
            vec_mean_exc_fr[vector_idx] = _mm_add_ps(vec_mean_exc_fr[vector_idx], temp_transfer_exc);
            
            // Inhibitory population (NO FFI INPUT)
            temp_current_inh = _mm_sub_ps(
                _mm_mul_ps(vec_inh_gain,
                    _mm_sub_ps(_mm_add_ps(vec_ext_inh_input,
                                         _mm_mul_ps(vec_nmda_coupling, vec_exc_synaptic[vector_idx])),
                              vec_inh_synaptic[vector_idx])),
                vec_inh_bias);
            
            *vec_temp_exp_inh = _mm_mul_ps(vec_neg_inh_membrane, temp_current_inh);
            temp_exponential_inh[0] = expf(temp_exponential_inh[0]);
            temp_exponential_inh[1] = expf(temp_exponential_inh[1]);
            temp_exponential_inh[2] = expf(temp_exponential_inh[2]);
            temp_exponential_inh[3] = expf(temp_exponential_inh[3]);
            
            *vec_temp_exp_inh = _mm_load_ps(temp_exponential_inh);
            temp_transfer_inh = _mm_div_ps(temp_current_inh, _mm_sub_ps(vec_one, *vec_temp_exp_inh));
            
            vec_inh_firing[vector_idx] = temp_transfer_inh;
            vec_mean_inh_fr[vector_idx] = _mm_add_ps(vec_mean_inh_fr[vector_idx], temp_transfer_inh);
            
            // Update synaptic activities
            generate_gaussian_noise(random_noise);
            vec_inh_synaptic[vector_idx] = _mm_add_ps(
                _mm_add_ps(_mm_mul_ps(vec_noise_scale, *vec_random_noise), vec_inh_synaptic[vector_idx]),
                _mm_mul_ps(vec_dt, _mm_add_ps(_mm_mul_ps(vec_neg_inv_inh_tau, vec_inh_synaptic[vector_idx]),
                                             _mm_mul_ps(temp_transfer_inh, vec_inh_kinetic))));
            
            generate_gaussian_noise(random_noise);
            vec_exc_synaptic[vector_idx] = _mm_add_ps(
                _mm_add_ps(_mm_mul_ps(vec_noise_scale, *vec_random_noise), vec_exc_synaptic[vector_idx]),
                _mm_mul_ps(vec_dt, _mm_add_ps(_mm_mul_ps(vec_neg_inv_exc_tau, vec_exc_synaptic[vector_idx]),
                                             _mm_mul_ps(_mm_mul_ps(_mm_sub_ps(vec_one, vec_exc_synaptic[vector_idx]), vec_exc_kinetic),
                                                       temp_transfer_exc))));
        }
        
        // Apply bounds
        for (int j = 0; j < padded_regions; j++) {
            excitatory_synaptic_activity[j] = fmaxf(0.0f, fminf(1.0f, excitatory_synaptic_activity[j]));
            inhibitory_synaptic_activity[j] = fmaxf(0.0f, fminf(1.0f, inhibitory_synaptic_activity[j]));
        }
        
        firing_rate_sample_count++;
        
        // STEP 3: BALLOON-WINDKESSEL DYNAMICS
        // Model 3a (Variable parameters)
        for (int j = 0; j < actual_regions; j++) {
            bw_vasodilatory_signal[j] += model_sampling_rate * 
                (excitatory_synaptic_activity[j] - bw_params.kappa[j] * bw_vasodilatory_signal[j] - 
                 bw_params.gamma[j] * (bw_blood_flow[j] - 1.0));
            
            float new_flow = bw_blood_flow[j] + model_sampling_rate * bw_vasodilatory_signal[j];
            
            bw_blood_volume[j] += model_sampling_rate / bw_params.tau[j] * 
                (bw_blood_flow[j] - powf(bw_blood_volume[j], 1.0f/bw_params.alpha[j]));
            
            bw_deoxyhemoglobin[j] += model_sampling_rate / bw_params.tau[j] * 
                (bw_blood_flow[j] * (1.0 - powf(1.0f - bw_params.rho[j], 1.0f/bw_blood_flow[j])) / bw_params.rho[j] - 
                 powf(bw_blood_volume[j], 1.0f/bw_params.alpha[j]) * bw_deoxyhemoglobin[j] / bw_blood_volume[j]);
            
            bw_blood_flow[j] = new_flow;
        }
        
        // Model 3b (Control - homogeneous parameters)
        for (int j = 0; j < actual_regions; j++) {
            bw_vasodilatory_signal_ctrl[j] += model_sampling_rate * 
                (excitatory_synaptic_activity[j] - default_kappa * bw_vasodilatory_signal_ctrl[j] - 
                 default_gamma * (bw_blood_flow_ctrl[j] - 1.0));
            
            float new_flow_ctrl = bw_blood_flow_ctrl[j] + model_sampling_rate * bw_vasodilatory_signal_ctrl[j];
            
            bw_blood_volume_ctrl[j] += model_sampling_rate / default_tau * 
                (bw_blood_flow_ctrl[j] - powf(bw_blood_volume_ctrl[j], 1.0f/default_alpha));
            
            bw_deoxyhemoglobin_ctrl[j] += model_sampling_rate / default_tau * 
                (bw_blood_flow_ctrl[j] * (1.0 - powf(1.0f - default_rho, 1.0f/bw_blood_flow_ctrl[j])) / default_rho - 
                 powf(bw_blood_volume_ctrl[j], 1.0f/default_alpha) * bw_deoxyhemoglobin_ctrl[j] / bw_blood_volume_ctrl[j]);
            
            bw_blood_flow_ctrl[j] = new_flow_ctrl;
        }
        
        // STEP 4: BOLD TIMEPOINT AND FIC UPDATE
        bold_timestep_counter++;
        
        if (bold_timestep_counter % bold_repetition_time == 0) {
            bold_length++;
            
            // Calculate mean firing rates
            float overall_mean_exc_rate = 0;
            for (int j = 0; j < actual_regions; j++) {
                mean_excitatory_firing_rate[j] /= firing_rate_sample_count;
                mean_inhibitory_firing_rate[j] /= firing_rate_sample_count;
                overall_mean_exc_rate += mean_excitatory_firing_rate[j];
            }
            overall_mean_exc_rate /= actual_regions;
            
            // FIC UPDATE (only during burn-in or if not converged)
            if (!burnin_complete) {
                // Store previous J_I values for convergence check
                for (int j = 0; j < actual_regions; j++) {
                    previous_J_I[j] = feedback_inhibition_strength[j];
                }
                
                // FIC plasticity update targeting 3Hz excitatory firing
                for (int j = 0; j < actual_regions; j++) {
                    float firing_rate_error = mean_excitatory_firing_rate[j] - target_excitatory_rate;
                    feedback_inhibition_strength[j] += fic_learning_rate * firing_rate_error;
                    feedback_inhibition_strength[j] = fmaxf(0.0f, feedback_inhibition_strength[j]);
                }
                
                // Check convergence criteria
                int firing_rate_converged = 1;
                float max_J_change = 0;
                
                for (int j = 0; j < actual_regions; j++) {
                    if (fabsf(mean_excitatory_firing_rate[j] - target_excitatory_rate) > convergence_threshold) {
                        firing_rate_converged = 0;
                    }
                    float j_change = fabsf(feedback_inhibition_strength[j] - previous_J_I[j]);
                    if (previous_J_I[j] > 0) {
                        j_change /= previous_J_I[j];
                    }
                    if (j_change > max_J_change) {
                        max_J_change = j_change;
                    }
                }
                
                int fic_stable = (max_J_change < fic_change_threshold);
                
                if (firing_rate_converged && fic_stable) {
                    consecutive_stable_steps++;
                } else {
                    consecutive_stable_steps = 0;
                }
                
                // Check if burn-in is complete
                if (consecutive_stable_steps >= required_stable_steps || timestep >= max_burnin_time) {
                    burnin_complete = 1;
                    burnin_timesteps = timestep;
                    printf("\n===============================================================================\n");
                    if (consecutive_stable_steps >= required_stable_steps) {
                        printf("BURN-IN CONVERGED after %.2f minutes\n", timestep/60000.0);
                    } else {
                        printf("BURN-IN TIMEOUT after %.2f minutes\n", timestep/60000.0);
                    }
                    printf("Final mean firing rate: %.2f Hz\n", overall_mean_exc_rate);
                    printf("Max FIC change: %.4f%%\n", max_J_change * 100);
                    printf("Starting main simulation phase...\n");
                    printf("===============================================================================\n");
                    
                    // Reset BOLD counter and length for main phase
                    bold_length = -1;
                }
                
                // Progress during burn-in
                if (timestep % 10000 == 0) {
                    printf("Burn-in: t=%.1fmin, FR=%.2fHz, stable=%d/%d\n", 
                           timestep/60000.0, overall_mean_exc_rate, consecutive_stable_steps, required_stable_steps);
                }
            }
            
            // RECORD BOLD DATA (only during main phase)
            if (burnin_complete && bold_length >= 0) {
                if (bold_length < bold_timeseries_length) {
                    for (int j = 0; j < actual_regions; j++) {
                        // Variable BOLD signal
                        bold_timeseries[j * bold_timeseries_length + bold_length] = 
                            100.0f / bw_params.rho[j] * bw_params.V0[j] * 
                            (bw_params.k1[j] * (1 - bw_deoxyhemoglobin[j]) + 
                             bw_params.k2[j] * (1 - bw_deoxyhemoglobin[j] / bw_blood_volume[j]) + 
                             bw_params.k3[j] * (1 - bw_blood_volume[j]));
                        
                        // Control BOLD signal (homogeneous parameters)
                        bold_timeseries_ctrl[j * bold_timeseries_length + bold_length] = 
                            100.0f / default_rho * default_V0 * 
                            (default_k1 * (1 - bw_deoxyhemoglobin_ctrl[j]) + 
                             default_k2 * (1 - bw_deoxyhemoglobin_ctrl[j] / bw_blood_volume_ctrl[j]) + 
                             default_k3 * (1 - bw_blood_volume_ctrl[j]));
                    }
                }
                
                // Progress during main phase
                if (bold_length % 50 == 0) {
                    printf("Main phase: BOLD=%d/%d, FR=%.2fHz\n", 
                           bold_length + 1, target_bold_points, overall_mean_exc_rate);
                }
                
                // Stop when we have enough BOLD data
                if (bold_length >= target_bold_points - 1) {
                    break;
                }
            }
            
            // Reset firing rate accumulators
            for (int j = 0; j < padded_regions; j++) {
                mean_excitatory_firing_rate[j] = 0.0f;
                mean_inhibitory_firing_rate[j] = 0.0f;
            }
            firing_rate_sample_count = 0;
        }
        
        timestep++;
    }
    
    // ========================================================================
    // SAVE RESULTS
    // ========================================================================
    
    char output_filename[256];
    snprintf(output_filename, sizeof(output_filename), "%s/%s_BOLD_timeseries.txt", results_dir, subject_id);
    
    FILE *bold_file = fopen(output_filename, "w");
    if (bold_file) {
        fprintf(bold_file, "%d %d\n", actual_regions, bold_length + 1);
        for (int j = 0; j < actual_regions; j++) {
            for (int t = 0; t <= bold_length; t++) {
                fprintf(bold_file, "%.6f ", bold_timeseries[j * bold_timeseries_length + t]);
            }
            fprintf(bold_file, "\n");
        }
        fclose(bold_file);
        printf("BOLD timeseries saved to: %s\n", output_filename);
    }
    
    // Save control BOLD timeseries
    char control_output_filename[256];
    snprintf(control_output_filename, sizeof(control_output_filename), "%s/%s_BOLD_control.txt", results_dir, subject_id);
    
    FILE *bold_file_ctrl = fopen(control_output_filename, "w");
    if (bold_file_ctrl) {
        fprintf(bold_file_ctrl, "%d %d\n", actual_regions, bold_length + 1);
        for (int j = 0; j < actual_regions; j++) {
            for (int t = 0; t <= bold_length; t++) {
                fprintf(bold_file_ctrl, "%.6f ", bold_timeseries_ctrl[j * bold_timeseries_length + t]);
            }
            fprintf(bold_file_ctrl, "\n");
        }
        fclose(bold_file_ctrl);
        printf("Control BOLD timeseries saved to: %s\n", control_output_filename);
    }
    
    // ========================================================================
    // CLEANUP AND RESULTS
    // ========================================================================
    
    float total_runtime = (float)(time(NULL) - simulation_start_time);
    printf("\n===============================================================================\n");
    printf("SIMULATION COMPLETE\n");
    printf("Total runtime: %.2f seconds\n", total_runtime);
    printf("Burn-in time: %.2f minutes\n", burnin_timesteps/60000.0);
    printf("Main simulation: %.2f minutes\n", (timestep - burnin_timesteps)/60000.0);
    printf("BOLD output: %d regions x %d timepoints\n", actual_regions, target_bold_points);
    printf("===============================================================================\n");
    
    // Cleanup
    cleanup_data(actual_regions, connections_per_region, region_lookup_table, 
                current_weights, connection_sources,
                &bw_params, bold_timeseries, excitatory_synaptic_activity, inhibitory_synaptic_activity,
                excitatory_firing_rate, inhibitory_firing_rate, long_range_excitatory_input,
                feedback_inhibition_strength, mean_excitatory_firing_rate, mean_inhibitory_firing_rate,
                bw_vasodilatory_signal, bw_blood_flow, bw_blood_volume, bw_deoxyhemoglobin,
                bold_timeseries_ctrl, bw_vasodilatory_signal_ctrl, bw_blood_flow_ctrl, 
                bw_blood_volume_ctrl, bw_deoxyhemoglobin_ctrl);
    
    free(previous_J_I);
    
    return 0;
}