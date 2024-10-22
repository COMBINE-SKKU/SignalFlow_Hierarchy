## References

Oh YH, Ann YJ, Lee JJ, Ito T, Froudist-Walsh S, Paquola C, Milham M, Spreng RN, Margulies D, Bernhardt B, Woo CW, Hong SK. [**In vivo cartography of state-dependent signal flow hierarchy in the human cerebral cortex **], Arxiv 2024.

----

## Background (Abstract)

Understanding the principle of information flow across distributed brain networks is of paramount importance in neuroscience. Here, we introduce a novel neuroimaging framework leveraging integrated effective connectivity (iEC) and unconstrained signal flow mapping for data-driven discovery of the human cerebral functional hierarchy, one of the principal brain scaffolds for macroscale neural dynamics. Simulation and empirical validation demonstrated the high fidelity of iEC in recovering connectome directionality and its potential relationship with histologically defined feedforward and feedback pathways. Notably, the iEC-based hierarchy revealed a functional axis composed of sequentially ordered sensorimotor-association-paralimbic areas, a pattern supported by the Structural Model of laminar connectivity. This hierarchy was further demonstrated to flexibly reorganize according to brain state, elevating exteroceptive regions during external focus and interoceptive regions during an internally oriented state. Our study highlights the unique role of macroscale directed functional connections in uncovering a neurobiologically grounded, state-dependent signal flow hierarchy.

## Code

This foldder contains all matlab codes for iEC framework and the replication of the results in the paper. We use Bayesian optimization to estimate the iEC and Hopf model to simulate the brain signals (i.e., BOLD).

In our current study, we chose three effective connectivity algorithms to estimate the iEC, namely rDCM, VAR, and FASK. The iEC estimated from these algorithms are used for the simulation of BOLD signals and the subsequent analysis. (Note that users can choose algorithms other than the ones used in the current study, offering a flexibility to the framework.)

The codes are organized as follows:

- 'figures' folder: codes for reproducing the results of the entire study.
- 'data' folder: codes for data used in the analysis and simulation.
- 'utils' folder: codes for the iEC estimation and signal flow mapping.

Please read 'README.md' in each sub-folder for the details of the codes.

## Dependencies
This code repository was tested and ran successfully on Linux.

- MATLAB R2020a or later
- [cifti-matlab toolbox](https://github.com/Washington-University/cifti-matlab) (add path to the cifti-matlab toolbox when running the scripts)
- Mac users need the appropriate version of XCode to run the C file in rDCM.
- All other necessary toolboxes are included in the utils folder.
