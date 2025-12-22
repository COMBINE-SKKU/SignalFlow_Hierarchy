/*
===============================================================================
SINGLE REGION VARIABILITY STUDY - HRF Parameter Effect Analysis
===============================================================================

Simplified single-region brain model for HRF parameter variability analysis:
1. Single region neural dynamics (no network connectivity)
2. Same neural signal applied to multiple HRF models
3. Variability sweep: 10%, 28%, 46%, 64%, 82%, 100% (6 linear steps)
4. Outputs: Neural signal (raw & downsampled), BOLD signals for each variability level

Key Components:
- Neural dynamics: DMF equations for single region with E/I populations
- BOLD generation: Multiple Balloon-Windkessel models with varying HRF parameters
- FIC: Vogels plasticity rule for local inhibition (3Hz target)
- Burn-in: Automatic convergence detection

Usage: ./single_region_variability <random_seed> <results_dir>

===============================================================================
*/

#include <stdio.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

// ============================================================================
// DATA STRUCTURES
// ============================================================================

struct BWParams {
    float rho;
    float alpha;
    float tau;
    float kappa;
    float gamma;
    float V0;
    float k1;
    float k2;
    float k3;
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
// BW PARAMETER GENERATION WITH VARIABILITY
// ============================================================================

void generate_bw_parameters(struct BWParams *bw_params, float variability) {
    // Default values
    float default_rho = 0.34f;
    float default_alpha = 0.32f;
    float default_tau = 0.98f;
    float default_kappa = 1.0f/0.65f;
    float default_gamma = 1.0f/0.41f;
    float default_V0 = 0.02f;
    
    if (variability == 0.0f) {
        // No variability - use defaults
        bw_params->rho = default_rho;
        bw_params->alpha = default_alpha;
        bw_params->tau = default_tau;
        bw_params->kappa = default_kappa;
        bw_params->gamma = default_gamma;
        bw_params->V0 = default_V0;
    } else {
        // Apply variability
        bw_params->rho = default_rho * (1.0f + variability * ((float)rand()/RAND_MAX - 0.5f) * 2.0f);
        bw_params->alpha = default_alpha * (1.0f + variability * ((float)rand()/RAND_MAX - 0.5f) * 2.0f);
        bw_params->tau = default_tau * (1.0f + variability * ((float)rand()/RAND_MAX - 0.5f) * 2.0f);
        bw_params->kappa = default_kappa * (1.0f + variability * ((float)rand()/RAND_MAX - 0.5f) * 2.0f);
        bw_params->gamma = default_gamma * (1.0f + variability * ((float)rand()/RAND_MAX - 0.5f) * 2.0f);
        bw_params->V0 = default_V0 * (1.0f + variability * ((float)rand()/RAND_MAX - 0.5f) * 2.0f);
        
        // Ensure positive values
        bw_params->rho = fmaxf(0.1f, bw_params->rho);
        bw_params->alpha = fmaxf(0.1f, bw_params->alpha);
        bw_params->tau = fmaxf(0.1f, bw_params->tau);
        bw_params->kappa = fmaxf(0.1f, bw_params->kappa);
        bw_params->gamma = fmaxf(0.1f, bw_params->gamma);
        bw_params->V0 = fmaxf(0.001f, bw_params->V0);
    }
    
    // Calculate derived parameters
    bw_params->k1 = 7.0f * bw_params->rho;
    bw_params->k2 = 2.0f;
    bw_params->k3 = 2.0f * bw_params->rho - 0.2f;
}

// ============================================================================
// MAIN SIMULATION
// ============================================================================

int main(int argc, char *argv[]) {
    time_t simulation_start_time = time(NULL);
    
    // ========================================================================
    // COMMAND LINE ARGUMENTS
    // ========================================================================
    if (argc != 3) {
        printf("Usage: %s <random_seed> <results_dir>\n", argv[0]);
        exit(1);
    }
    
    int random_seed = atoi(argv[1]);
    char *results_dir = argv[2];
    
    printf("===============================================================================\n");
    printf("SINGLE REGION VARIABILITY STUDY - HRF Parameter Effect Analysis\n");
    printf("===============================================================================\n");
    printf("Random seed: %d\n", random_seed);
    printf("Results directory: %s\n", results_dir);
    
    srand((unsigned)random_seed);
    
    // ========================================================================
    // MODEL PARAMETERS
    // ========================================================================
    
    const float integration_timestep = 1.0;
    const float model_sampling_rate = 0.001;
    int total_simulation_timesteps;
    int bold_repetition_time = 1000;
    
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
    const float excitatory_coupling_strength = 0.15;
    float local_recurrent_strength = 1.5 * excitatory_coupling_strength;
    
    // FIC parameters
    const float fic_learning_rate = 0.01;
    const float target_excitatory_rate = 3.0; // 3Hz target
    const float convergence_threshold = 0.5; // 0.5 Hz tolerance
    const float fic_change_threshold = 0.005; // 0.5% change
    const int required_stable_steps = 10;
    const int max_burnin_time = 1000000; // Max burn-in time
    
    // Calculate simulation time for exactly 490 BOLD timepoints
    int target_bold_points = 490;
    total_simulation_timesteps = target_bold_points * bold_repetition_time;
    
    // Variability levels: 0% (default), 30%, 60%, 90% (4 models total)
    const int num_variability_levels = 4;
    float variability_levels[4] = {0.00, 0.30, 0.60, 0.90};
    
    printf("Simulation settings:\n");
    printf("  Target BOLD timepoints: %d\n", target_bold_points);
    printf("  Simulation time: %.1f minutes\n", total_simulation_timesteps/60000.0);
    printf("  Variability levels: ");
    for (int i = 0; i < num_variability_levels; i++) {
        printf("%.0f%% ", variability_levels[i] * 100);
    }
    printf("\n");
    
    // ========================================================================
    // ALLOCATE NEURAL VARIABLES (SINGLE REGION)
    // ========================================================================
    
    float excitatory_synaptic_activity = 0.001f;
    float inhibitory_synaptic_activity = 0.001f;
    float excitatory_firing_rate = 0.0f;
    float inhibitory_firing_rate = 0.0f;
    float mean_excitatory_firing_rate = 0.0f;
    float mean_inhibitory_firing_rate = 0.0f;
    float feedback_inhibition_strength = initial_inhibitory_strength;
    
    // Neural signal recording (raw and downsampled)
    int neural_raw_length = total_simulation_timesteps + max_burnin_time;
    int neural_downsampled_length = target_bold_points;
    
    float *neural_signal_raw = (float *)malloc(neural_raw_length * sizeof(float));
    float *neural_signal_downsampled = (float *)malloc(neural_downsampled_length * sizeof(float));
    
    if (!neural_signal_raw || !neural_signal_downsampled) {
        printf("ERROR: Memory allocation failed for neural signals\n");
        exit(2);
    }
    
    // BOLD variables for each variability level
    struct BWParams bw_models[4]; // One for each variability level (including default)
    float *bw_vasodilatory_signal[4];
    float *bw_blood_flow[4]; 
    float *bw_blood_volume[4];
    float *bw_deoxyhemoglobin[4];
    float **bold_timeseries = (float **)malloc(num_variability_levels * sizeof(float *));
    
    if (!bold_timeseries) {
        printf("ERROR: Memory allocation failed for BOLD arrays\n");
        exit(2);
    }
    
    // Initialize BW parameters and state variables for each model
    for (int model = 0; model < num_variability_levels; model++) {
        // Generate BW parameters with specific variability
        generate_bw_parameters(&bw_models[model], variability_levels[model]);
        
        // Allocate BW state variables
        bw_vasodilatory_signal[model] = (float *)malloc(sizeof(float));
        bw_blood_flow[model] = (float *)malloc(sizeof(float));
        bw_blood_volume[model] = (float *)malloc(sizeof(float));
        bw_deoxyhemoglobin[model] = (float *)malloc(sizeof(float));
        bold_timeseries[model] = (float *)malloc(target_bold_points * sizeof(float));
        
        if (!bw_vasodilatory_signal[model] || !bw_blood_flow[model] || 
            !bw_blood_volume[model] || !bw_deoxyhemoglobin[model] || !bold_timeseries[model]) {
            printf("ERROR: Memory allocation failed for model %d\n", model);
            exit(2);
        }
        
        // Initialize BW state variables
        *bw_vasodilatory_signal[model] = 0.0f;
        *bw_blood_flow[model] = 1.0f;
        *bw_blood_volume[model] = 1.0f;
        *bw_deoxyhemoglobin[model] = 1.0f;
        
        printf("Model %d (%.0f%% variability): rho=%.3f, alpha=%.3f, tau=%.3f\n", 
               model + 1, variability_levels[model] * 100,
               bw_models[model].rho, bw_models[model].alpha, bw_models[model].tau);
    }
    
    // ========================================================================
    // VECTORIZATION CONSTANTS
    // ========================================================================
    
    const float negative_exc_membrane_time = -1.0 * exc_membrane_time;
    const float negative_inh_membrane_time = -1.0 * inh_membrane_time;
    const float negative_inv_exc_time_const = -1.0 / excitatory_time_constant;
    const float negative_inv_inh_time_const = -1.0 / inhibitory_time_constant;
    const float scaled_external_exc_input = exc_external_scaling * external_input_baseline;
    const float scaled_external_inh_input = inh_external_scaling * external_input_baseline;
    const float scaled_noise = noise_amplitude * integration_timestep;
    
    // ========================================================================
    // BURN-IN PHASE VARIABLES
    // ========================================================================
    
    int firing_rate_sample_count = 0;
    int bold_timestep_counter = -1;
    int bold_length = -1;
    int consecutive_stable_steps = 0;
    int burnin_complete = 0;
    int burnin_timesteps = 0;
    float previous_J_I = feedback_inhibition_strength;
    int neural_raw_index = 0;
    
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
        
        // STEP 1: NEURAL DYNAMICS (SINGLE REGION)
        float random_noise[4];
        
        // Excitatory population
        float current_exc = exc_gain_factor * 
            (scaled_external_exc_input + local_recurrent_strength * excitatory_synaptic_activity - 
             feedback_inhibition_strength * inhibitory_synaptic_activity) - exc_bias_current;
        
        float exp_exc = expf(negative_exc_membrane_time * current_exc);
        float transfer_exc = current_exc / (1.0f - exp_exc);
        
        excitatory_firing_rate = transfer_exc;
        mean_excitatory_firing_rate += transfer_exc;
        
        // Inhibitory population 
        float current_inh = inh_gain_factor * 
            (scaled_external_inh_input + excitatory_coupling_strength * excitatory_synaptic_activity - 
             inhibitory_synaptic_activity) - inh_bias_current;
        
        float exp_inh = expf(negative_inh_membrane_time * current_inh);
        float transfer_inh = current_inh / (1.0f - exp_inh);
        
        inhibitory_firing_rate = transfer_inh;
        mean_inhibitory_firing_rate += transfer_inh;
        
        // Update synaptic activities
        generate_gaussian_noise(random_noise);
        inhibitory_synaptic_activity += scaled_noise * random_noise[0] + 
            integration_timestep * (negative_inv_inh_time_const * inhibitory_synaptic_activity + 
                                   transfer_inh * inh_kinetic_factor);
        
        generate_gaussian_noise(random_noise);
        excitatory_synaptic_activity += scaled_noise * random_noise[0] + 
            integration_timestep * (negative_inv_exc_time_const * excitatory_synaptic_activity + 
                                   (1.0f - excitatory_synaptic_activity) * exc_kinetic_factor * transfer_exc);
        
        // Apply bounds
        excitatory_synaptic_activity = fmaxf(0.0f, fminf(1.0f, excitatory_synaptic_activity));
        inhibitory_synaptic_activity = fmaxf(0.0f, fminf(1.0f, inhibitory_synaptic_activity));
        
        // Record neural signal (raw)
        if (neural_raw_index < neural_raw_length) {
            neural_signal_raw[neural_raw_index] = excitatory_synaptic_activity;
            neural_raw_index++;
        }
        
        firing_rate_sample_count++;
        
        // STEP 2: BALLOON-WINDKESSEL DYNAMICS (ALL MODELS)
        for (int model = 0; model < num_variability_levels; model++) {
            struct BWParams *bw = &bw_models[model];
            
            *bw_vasodilatory_signal[model] += model_sampling_rate * 
                (excitatory_synaptic_activity - bw->kappa * (*bw_vasodilatory_signal[model]) - 
                 bw->gamma * ((*bw_blood_flow[model]) - 1.0));
            
            float new_flow = (*bw_blood_flow[model]) + model_sampling_rate * (*bw_vasodilatory_signal[model]);
            
            *bw_blood_volume[model] += model_sampling_rate / bw->tau * 
                ((*bw_blood_flow[model]) - powf(*bw_blood_volume[model], 1.0f/bw->alpha));
            
            *bw_deoxyhemoglobin[model] += model_sampling_rate / bw->tau * 
                ((*bw_blood_flow[model]) * (1.0 - powf(1.0f - bw->rho, 1.0f/(*bw_blood_flow[model]))) / bw->rho - 
                 powf(*bw_blood_volume[model], 1.0f/bw->alpha) * (*bw_deoxyhemoglobin[model]) / (*bw_blood_volume[model]));
            
            *bw_blood_flow[model] = new_flow;
        }
        
        // STEP 3: BOLD TIMEPOINT AND FIC UPDATE
        bold_timestep_counter++;
        
        if (bold_timestep_counter % bold_repetition_time == 0) {
            bold_length++;
            
            // Calculate mean firing rate
            mean_excitatory_firing_rate /= firing_rate_sample_count;
            mean_inhibitory_firing_rate /= firing_rate_sample_count;
            
            // FIC UPDATE (only during burn-in or if not converged)
            if (!burnin_complete) {
                // Store previous J_I value for convergence check
                previous_J_I = feedback_inhibition_strength;
                
                // FIC plasticity update targeting 3Hz excitatory firing
                float firing_rate_error = mean_excitatory_firing_rate - target_excitatory_rate;
                feedback_inhibition_strength += fic_learning_rate * firing_rate_error;
                feedback_inhibition_strength = fmaxf(0.0f, feedback_inhibition_strength);
                
                // Check convergence criteria
                int firing_rate_converged = (fabsf(mean_excitatory_firing_rate - target_excitatory_rate) <= convergence_threshold);
                
                float j_change = fabsf(feedback_inhibition_strength - previous_J_I);
                if (previous_J_I > 0) {
                    j_change /= previous_J_I;
                }
                int fic_stable = (j_change < fic_change_threshold);
                
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
                    printf("Final mean firing rate: %.2f Hz\n", mean_excitatory_firing_rate);
                    printf("Final FIC strength: %.4f\n", feedback_inhibition_strength);
                    printf("Starting main simulation phase...\n");
                    printf("===============================================================================\n");
                    
                    // Reset BOLD counter and length for main phase
                    bold_length = -1;
                }
                
                // Progress during burn-in
                if (timestep % 10000 == 0) {
                    printf("Burn-in: t=%.1fmin, FR=%.2fHz, J_I=%.4f, stable=%d/%d\n", 
                           timestep/60000.0, mean_excitatory_firing_rate, feedback_inhibition_strength,
                           consecutive_stable_steps, required_stable_steps);
                }
            }
            
            // RECORD BOLD DATA (only during main phase)
            if (burnin_complete && bold_length >= 0) {
                if (bold_length < target_bold_points) {
                    for (int model = 0; model < num_variability_levels; model++) {
                        struct BWParams *bw = &bw_models[model];
                        
                        bold_timeseries[model][bold_length] = 
                            100.0f / bw->rho * bw->V0 * 
                            (bw->k1 * (1 - (*bw_deoxyhemoglobin[model])) + 
                             bw->k2 * (1 - (*bw_deoxyhemoglobin[model]) / (*bw_blood_volume[model])) + 
                             bw->k3 * (1 - (*bw_blood_volume[model])));
                    }
                    
                    // Record downsampled neural signal
                    if (bold_length < neural_downsampled_length) {
                        neural_signal_downsampled[bold_length] = excitatory_synaptic_activity;
                    }
                }
                
                // Progress during main phase
                if (bold_length % 50 == 0) {
                    printf("Main phase: BOLD=%d/%d, FR=%.2fHz, E_syn=%.4f\n", 
                           bold_length + 1, target_bold_points, mean_excitatory_firing_rate, excitatory_synaptic_activity);
                }
                
                // Stop when we have enough BOLD data
                if (bold_length >= target_bold_points - 1) {
                    break;
                }
            }
            
            // Reset firing rate accumulators
            mean_excitatory_firing_rate = 0.0f;
            mean_inhibitory_firing_rate = 0.0f;
            firing_rate_sample_count = 0;
        }
        
        timestep++;
    }
    
    // ========================================================================
    // SAVE RESULTS
    // ========================================================================
    
    printf("\n===============================================================================\n");
    printf("SAVING RESULTS\n");
    printf("===============================================================================\n");
    
    // Save neural signal (raw)
    char neural_raw_filename[512];
    snprintf(neural_raw_filename, sizeof(neural_raw_filename), "%s/neural_signal_raw.txt", results_dir);
    FILE *neural_raw_file = fopen(neural_raw_filename, "w");
    if (neural_raw_file) {
        fprintf(neural_raw_file, "# Raw neural signal (excitatory synaptic gating variable)\n");
        fprintf(neural_raw_file, "# Time points: %d, Sampling rate: %.6f ms\n", neural_raw_index, integration_timestep);
        for (int t = 0; t < neural_raw_index; t++) {
            fprintf(neural_raw_file, "%.6f\n", neural_signal_raw[t]);
        }
        fclose(neural_raw_file);
        printf("Neural signal (raw) saved to: %s\n", neural_raw_filename);
    }
    
    // Save neural signal (downsampled to BOLD resolution)
    char neural_downsampled_filename[512];
    snprintf(neural_downsampled_filename, sizeof(neural_downsampled_filename), "%s/neural_signal_downsampled.txt", results_dir);
    FILE *neural_downsampled_file = fopen(neural_downsampled_filename, "w");
    if (neural_downsampled_file) {
        fprintf(neural_downsampled_file, "# Downsampled neural signal (excitatory synaptic gating variable)\n");
        fprintf(neural_downsampled_file, "# Time points: %d, BOLD resolution (TR=1s)\n", bold_length + 1);
        for (int t = 0; t <= bold_length; t++) {
            fprintf(neural_downsampled_file, "%.6f\n", neural_signal_downsampled[t]);
        }
        fclose(neural_downsampled_file);
        printf("Neural signal (downsampled) saved to: %s\n", neural_downsampled_filename);
    }
    
    // Save BOLD timeseries for each variability level
    for (int model = 0; model < num_variability_levels; model++) {
        char bold_filename[512];
        
        if (variability_levels[model] == 0.0f) {
            // Save default signal separately
            snprintf(bold_filename, sizeof(bold_filename), "%s/bold_default.txt", results_dir);
        } else {
            snprintf(bold_filename, sizeof(bold_filename), "%s/bold_variability_%02d_percent.txt", 
                    results_dir, (int)(variability_levels[model] * 100));
        }
        
        FILE *bold_file = fopen(bold_filename, "w");
        if (bold_file) {
            if (variability_levels[model] == 0.0f) {
                fprintf(bold_file, "# Default BOLD signal with 0%% HRF parameter variability\n");
            } else {
                fprintf(bold_file, "# BOLD signal with %.0f%% HRF parameter variability\n", variability_levels[model] * 100);
            }
            fprintf(bold_file, "# BW parameters: rho=%.4f, alpha=%.4f, tau=%.4f, kappa=%.4f, gamma=%.4f, V0=%.6f\n",
                   bw_models[model].rho, bw_models[model].alpha, bw_models[model].tau,
                   bw_models[model].kappa, bw_models[model].gamma, bw_models[model].V0);
            fprintf(bold_file, "# Time points: %d\n", bold_length + 1);
            
            for (int t = 0; t <= bold_length; t++) {
                fprintf(bold_file, "%.6f\n", bold_timeseries[model][t]);
            }
            fclose(bold_file);
            
            if (variability_levels[model] == 0.0f) {
                printf("Default BOLD signal (0%% variability) saved to: %s\n", bold_filename);
            } else {
                printf("BOLD signal (%02d%% variability) saved to: %s\n", 
                       (int)(variability_levels[model] * 100), bold_filename);
            }
        }
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
    printf("Neural signal length: %d timepoints (raw), %d timepoints (downsampled)\n", 
           neural_raw_index, bold_length + 1);
    printf("BOLD outputs: %d signals with variability levels 0%%, 30%%, 60%%, 90%%\n", num_variability_levels);
    printf("===============================================================================\n");
    
    // Cleanup
    free(neural_signal_raw);
    free(neural_signal_downsampled);
    
    for (int model = 0; model < num_variability_levels; model++) {
        free(bw_vasodilatory_signal[model]);
        free(bw_blood_flow[model]);
        free(bw_blood_volume[model]);
        free(bw_deoxyhemoglobin[model]);
        free(bold_timeseries[model]);
    }
    free(bold_timeseries);
    
    return 0;
}