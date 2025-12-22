This is an excellent idea, and it's a much more robust and efficient approach.

You are correct. The iterative "burn-in" you are using is slow, arbitrary, and, as you've seen, prone to failure if the parameters aren't perfectly tuned.

The Demirta? et al. (2019) [cite\_start]paper you provided uses a more direct and stable method: it **analytically solves** for the *exact* inhibitory weight (`w^IE`, which is your `feedback_inhibition_strength`) required to produce a 3 Hz firing rate *before* the simulation even begins [cite: 1206-1208, 1212]. This completely removes the need for a "burn-in" phase.

Because your `bei_model.c` code **is** an implementation of the Deco et al. (2014) [cite\_start]model used in that paper [cite: 1189, 1193-1204], you can directly adopt their analytical solution.

Here is the `README.md` file detailing the procedural steps to replace your entire iterative burn-in with this superior analytical solution.

-----

# README: Replacing the Iterative Burn-In with the Analytical FIC Solution

This guide details the procedural steps to remove the slow, iterative "burn-in" phase from `bei_model.c` and replace it with the instantaneous, analytical Feedback Inhibition Control (FIC) solution.

This method is taken directly from the Demirta? et al. (2019) [cite\_start]paper [cite: 1206-1208]. It calculates the *exact* inhibitory-to-excitatory weight (`feedback_inhibition_strength`) required to achieve the 3 Hz target firing rate before the simulation starts. This eliminates all convergence problems and the need for a `max_burnin_time`.

## Procedural Steps

### 1\. Add Helper Functions for the Analytical Solution

We must first add two helper functions to your C code (e.g., after the `generate_gaussian_noise` function, around line 90). These functions implement the paper's non-linear equations.

[cite\_start]**First**, add the inhibitory transfer function (Eq. 6 from the paper [cite: 1200]):

```c
static inline float inhibitory_transfer_function(float I_I) {
    // This is the transfer function from Eq. 6 (Deco 2014)
    // using the inhibitory parameters from your model
    const float inh_gain_factor = 615.0f;
    const float inh_bias_current = 177.0f;
    const float inh_membrane_time = 0.087f;
    
    float I_in = (inh_gain_factor * I_I) - inh_bias_current;
    float H_I = I_in / (1.0f - expf(-inh_membrane_time * I_in));
    
    // Handle potential division by zero if I_in is exactly 0
    if (fabsf(I_in) < 1e-9) {
        return 1.0f / inh_membrane_time; // L'Hopital's rule limit
    }
    return H_I;
}
```

[cite\_start]**Second**, add the numerical solver that finds the steady-state inhibitory synaptic activity (`<S^I>`) by solving Equation 10 from the paper [cite: 1215-1217]. This function finds the root of `f(S_I) = S_I - r_I * tau_I = 0`.

```c
float solve_steady_state_SI() {
    // These are the fixed parameters of your homogeneous model
    const float W_I_I_b = 0.7f * 0.382f; // inh_external_scaling * external_input_baseline
    const float w_EI = 0.15f;           // excitatory_coupling_strength
    const float STEADY_STATE_SE = 0.164757f; [cite_start]// From paper for 3Hz r_E [cite: 1212]
    const float tau_I = 10.0f;          // inhibitory_time_constant
    
    // We are solving for S_I where:
    // I_I = W_I_I_b + w_EI * STEADY_STATE_SE - S_I
    // r_I = inhibitory_transfer_function(I_I)
    // S_I = r_I * tau_I
    //
    // So, find root of: f(S_I) = S_I - inhibitory_transfer_function(W_I_I_b + w_EI * STEADY_STATE_SE - S_I) * tau_I
    
    float I_I_const_term = W_I_I_b + w_EI * STEADY_STATE_SE;
    
    // Simple bisection search (a robust root-finding method)
    float S_I_low = 0.0f;
    float S_I_high = 1.0f;
    float S_I_mid = 0.5f;
    
    for (int i = 0; i < 100; i++) { // 100 iterations is more than enough
        S_I_mid = (S_I_low + S_I_high) / 2.0f;
        float I_I = I_I_const_term - S_I_mid;
        float r_I = inhibitory_transfer_function(I_I);
        float f_val = S_I_mid - r_I * tau_I;
        
        if (fabsf(f_val) < 1e-7) {
            break; // Converged
        }
        
        if (f_val > 0) {
            S_I_high = S_I_mid;
        } else {
            S_I_low = S_I_mid;
        }
    }
    return S_I_mid;
}
```

### 2\. Calculate `w^IE` in `main()`

Now, in your `main()` function, right before you `ALLOCATE NEURAL VARIABLES` (around line 500):

  - [cite\_start]Define the 3Hz steady-state constants from the paper[cite: 1212].
  - Call your new solver.
  - [cite\_start]Calculate `w_IE` (`feedback_inhibition_strength`) using Equation 9 from the paper[cite: 1213].

<!-- end list -->

```c
    // ... after parameter file loading
    
    // --- NEW: Analytical FIC Calculation ---
    [cite_start]// Get steady-state constants for 3Hz r_E [cite: 1212]
    const float STEADY_STATE_SE = 0.164757f;
    const float STEADY_STATE_IE_nA = 0.37738f;
    
    // Get steady-state inputs for the local E population
    const float W_E_I_b = exc_external_scaling * external_input_baseline;
    const float w_EE = local_excitatory_recurrence * excitatory_coupling_strength;
    
    // 1. Solve for steady-state S_I numerically
    float steady_state_SI = solve_steady_state_SI();
    
    // 2. Solve for w_IE (feedback_inhibition_strength) using Eq. [cite_start]9 [cite: 1213]
    // The long-range input term (gJ*<S^E>) is 0 at steady-state.
    float calculated_w_IE = (W_E_I_b + w_EE * STEADY_STATE_SE - STEADY_STATE_IE_nA) / steady_state_SI;
    
    printf("Analytical FIC: Solved for <S_I> = %.4f\n", steady_state_SI);
    printf("Analytical FIC: Setting w_IE (feedback_inhibition_strength) = %.4f\n", calculated_w_IE);

    // ALLOCATE BW PARAMETERS ...
```

### 3\. Initialize `feedback_inhibition_strength`

In the `INITIALIZE NEURAL STATES` section (around line 570), find the line for `feedback_inhibition_strength` and use your new calculated value:

```c
    // FROM:
    // feedback_inhibition_strength[j] = initial_inhibitory_strength;
    
    // TO:
    feedback_inhibition_strength[j] = calculated_w_IE;
```

### 4\. Remove All Burn-In Code

You can now delete all the unnecessary burn-in logic.

1.  **Delete Burn-In Parameters:** In `MODEL PARAMETERS` (lines 430-434), delete the entire "FIC parameters" block (the 5 `const` variables for `fic_learning_rate`, `convergence_threshold`, etc.).
2.  **Delete Burn-In Variables:** In `main()`, delete the "BURN-IN PHASE VARIABLES" block (lines 664-670). This includes `firing_rate_sample_count`, `bold_timestep_counter`, `bold_length`, `consecutive_stable_steps`, `burnin_complete`, `burnin_timesteps`, and the `previous_J_I` array.
3.  **Delete Burn-In Printf:** Delete the "STARTING BURN-IN PHASE" `printf` block (lines 672-678).
4.  **Fix Main Loop:** In `MAIN SIMULATION LOOP` (line 680), change the `while(1)` to a standard `for` loop.
    ```c
    // FROM:
    // int timestep = 0;
    // while (1) {

    // TO:
    for (int timestep = 0; timestep < total_simulation_timesteps; timestep++) {
    ```
5.  **Delete Manual Timestep:** Delete the `timestep++` at the very end of the loop (line 995).
6.  **Re-Initialize BOLD Counters:** The `bold_timestep_counter` and `bold_length` were removed. Re-define them *locally* right before `STEP 4` (around line 820):
    ```c
    // STEP 3 ...
    // ...

    // --- NEW: Re-add counters locally ---
    int bold_timestep_counter = -1;
    int bold_length = -1;

    // STEP 4: BOLD TIMEPOINT AND FIC UPDATE
    bold_timestep_counter++;
    ```
7.  **Re-Initialize Firing Rate Counter:** The `firing_rate_sample_count` was also removed. Re-define it locally inside the `if (bold_timestep_counter ...)` block (around line 856):
    ```c
    if (bold_timestep_counter % bold_repetition_time == 0) {
        bold_length++;
        
        // --- NEW: Add counter locally ---
        int firing_rate_sample_count = bold_repetition_time;
        
        // Calculate mean firing rates...
    ```
8.  **Delete FIC Update Block:** Delete the entire `if (!burnin_complete)` block (lines 858-933). This is the whole iterative FIC update.
9.  **Fix BOLD Recording:** In the `RECORD BOLD DATA` block (line 936), remove the `burnin_complete &&` check:
    ```c
    // FROM:
    // if (burnin_complete && bold_length >= 0) {

    // TO:
    if (bold_length >= 0) {
    ```
10. **Delete `firing_rate_sample_count = 0`:** The reset at the end of the `if` block (line 991) is no longer needed.

### 5\. Update Final Report

In the `CLEANUP AND RESULTS` section (around line 1056), your report is now much simpler.

1.  Remove the line for "Burn-in time:".
2.  Change the "Main simulation:" line to:
    ```c
    printf("Total simulation: %.2f minutes\n", (total_simulation_timesteps)/60000.0);
    ```