## mDLAG (multi-population DLAG)

[![][license-img]][license-url]

[license-img]: https://img.shields.io/github/license/mashape/apistatus.svg
[license-url]: https://github.com/egokcen/mDLAG/blob/main/LICENSE

This repository accompanies the following [paper](https://nips.cc/virtual/2023/poster/70171):
- Gokcen, E., Jasper, A. I., Xu, A., Kohn, A., Machens, C. K., & Yu, B. M. Uncovering motifs of concurrent signaling across multiple neuronal populations.
  _Advances in Neural Information Processing Systems_, 36 (2023).

Please read it carefully before using the code, as it describes all of the
terminology and usage modes. Please cite the above reference if using any
portion of this code for your own purposes.

This codepack includes an implementation of the primary contribution of the paper, mDLAG.
It also includes a custom implementation of group factor analysis (GFA), and—to keep the codepack self-contained—a 
minimal implementation of [DLAG](https://www.nature.com/articles/s43588-022-00282-5.epdf?sharing_token=hyFsoNFZfDyh-EVI76hoddRgN0jAjWel9jnR3ZoTv0Nm1Wps5WZ1Nog-dORHPLUG97YnGS0JZBkvhpO7c5pblBICIHRXMKZ04hmmro2Tn12HIbx2e2LrperSJc6bwzqptnPIaVOrqvl8DcloXzDaOBhLlAqzUvwM4uMyl96KvTE%3D), full repository [here](https://github.com/egokcen/DLAG).

The codepack directory structure is as follows:
- `DLAG`
  - Minimal implementation of DLAG
- `gfa`
  - Implementation of GFA. `gfa/demo/demo_gfa.m` is a script that provides a
    generic demonstration of GFA.
- `mDLAG`
  - Implementation of mDLAG. `mDLAG/demo/demo_mdlag.m` is a script
    that provides a generic demonstration of mDLAG.
- `sim1`
  - Includes the data, results, and code to reproduce 
    "Simulation 1: Uncovering directed interactions across 
    multiple populations" (Section 3, Fig. 3)
- `sim2`
  - Includes the data, results, and code to reproduce
    "Simulation 2: Disentangling concurrent, bidirectional 
    signaling" (Section 3, Fig. 4)
- `V1V2array`
  - Includes part of the data, results, and code to reproduce
    Supplementary Figs. S1 & S2a. Full data is available [here](https://doi.org/10.6080/K0B27SHN).
- `V1V2V3npx`
  - Includes part of the data, results, and code to reproduce
    Fig. 5 and Supplementary Figs. S2bc, S3, S4, and S5. Full data
    cannot be made available at this time.

## System requirements

This codepack was written in Matlab (The MathWorks, Inc.), and must be run in 
Matlab.

It has been tested on Matlab 2019a, 2019b, 2020a, 2022b, and 2023a, on Linux 
(Pop!_OS 22.04 LTS, Red Hat Enterprise Linux release 7.9) and Windows (Windows
10) operating systems.

Note that all runtimes provided below are estimates based on these tests. Your
runtimes may vary.

Some functions rely on the C/MEX Matlab interface for speedup, whereby C code 
can be compiled and used from within Matlab. A native C compiler is necessary
to take advantage of this functionality. Windows users may require extra 
installation of, for example, Microsoft Visual C++ or MinGW. The code is written
to default to a native Matlab code implementation if mex files cannot be 
properly executed, and will still work correctly.

## Installation guide

Install [Matlab](https://www.mathworks.com/products/matlab.html).

For one script, `sim1/crossval_noARD_sim1.m`, you may need to specifically 
install the Matlab [Bioinformatics Toolbox](https://www.mathworks.com/help/bioinfo/index.html), 
if it's not already installed in your Matlab build.

C/MEX Compilation:
1. Enter DLAG directory. Run startup.m.
2. Enter mDLAG directory. Run startup.m.

Assuming Matlab is already installed on your machine, setup should not take 
more than a few minutes.

## Instructions for use

### GFA Generic Demo

1. Enter `gfa` directory.
2. Run `startup.m` to add all necessary dependencies to the Matlab path.
3. Open `gfa/demo/demo_gfa.m` (remain in `gfa` directory).
4. Run `demo_gfa.m` cell-by-cell. The script contains directions and descriptions
   of relevant user-defined parameters.
   
With the current data and settings, `demo_gfa.m` should take less than 1 min to
run.

### mDLAG Generic Demo

1. Enter `mDLAG` directory.
2. Run `startup.m` to add all necessary dependencies to the Matlab path.
3. Open `mDLAG/demo/demo_mdlag.m` (remain in `mDLAG` directory).
4. Run `demo_mdlag.m` cell-by-cell. The script contains directions and 
   descriptions of relevant user-defined parameters.
   
With the current data and settings, `demo_mdlag.m` should take less than 10 min to
run.

### Simulation 1: Uncovering directed interactions across multiple populations

1. Enter `sim1` directory.
2. Run `startup.m` to set up common directories and other constants.
3. mDLAG experiments (with ARD)
   1. To see pre-saved results, open `results_summary_mdlag_sim1.m` and run
      it cell-by-cell.
   2. To fit mDLAG from scratch, run `fit_mdlag_sim1.m`. The script should 
      take ~45 min to run.
4. mDLAG experiments (no ARD)
    1. To see pre-saved results, open `result_summary_noARD_sim1.m` and run
       it cell-by-cell.
    2. To perform the mDLAG (no ARD) experiments from scratch:
         1. To use cross-validation to find the total latent
            dimensionality, run `crossval_noARD_sim1.m`. The script should
            take ~75 min to run without parallelization.
         2. To fit the model with the correct number of latent dimensions,
            run `fit_noARD_sim1.m`. The script should take ~2 hours to run.

### Simulation 2: Disentangling concurrent, bidirectional signaling

1. Enter `sim2` directory.
2. Run `startup.m` to set up common directories and other constants.
3. mDLAG experiments
   1. To see pre-saved results, open `results_summary_mdlag_sim2.m` and run
      it cell-by-cell.
   2. To fit mDLAG from scratch, run `fit_mdlag_sim2.m`. The script should
      take less than 5 min to run.
4. GFA experiments
    1. To see pre-saved results, open `results_summary_gfa_sim2.m` and run
       it cell-by-cell.
    2. To fit GFA from scratch, run `fit_gfa_sim2.m`. The script should take
       less than 1 min to run.
          
### Validating mDLAG on recordings from V1 and V2

1. Enter `V1V2array` directory.
2. Run `startup.m` to set up common directories and other constants.
3. mDLAG experiments
    1. To see pre-saved results, open 
       `./mdlag/mdlag_exampledataset_v1v2array.m` (remain in `V1V2array` 
          directory) and run it cell-by-cell. 
    2. To fit mDLAG from scratch, run `./mdlag/fit_mdlag_v1v2array.m`
       (remain in `V1V2array` directory). The script should take ~16 hours to
       run.
4. DLAG experiments
    1. To see pre-saved results, open `./dlag/dlag_exampledataset_v1v2array.m`
       (remain in `V1V2array` directory) and run it cell-by-cell.
    2. To fit DLAG from scratch (here we're providing the optimal
       dimensionalities, selected through cross-validation), run
       `./dlag/fit_dlag_v1v2array.m` (remain in `V1V2array` directory). 
       The script should take ~5 hours to run. 
5. mDLAG vs DLAG performance comparison: Run `./mdlag_vs_dlag_v1v2array.m` to 
   produce Supplementary Fig. S2a.
   
### Interactions across laminar compartments of V1, V2, and V3d

1. Enter `V1V2V3npx` directory.
2. Run `startup.m` to set up common directories and other constants.
3. Open `mdlag_exampledataset_npx.m`. Run it cell-by-cell to see mDLAG results
   on the example Neuropixels dataset (Fig. 5).
4. Run `performance_summary.m` to compare performance of mDLAG to alternative
   methods (Supplementary Fig. S2bc).
5. Run `data_efficiency_summary.m` to show mDLAG test performance as a function
   of the number of available training trials (Supplementary Fig. S3).
6. Run `runtime_summary.m` to show runtimes of mDLAG and GFA (Supplementary
   Fig. S4).
7. Run `init_sensitivity_summary.m` to show sensitivity of mDLAG Neuropixels 
   results to random initialization (Supplementary Fig. S5).
   
