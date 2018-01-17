# Sync-deSync-model
The Sync/deSync (SDS) model: How a synchronized hippocampus and a de-synchronized neocortex code memories. This repository houses the matlab files that were used for this project.

To use this code, copy all files to the same folder directory. Using Matlab, run the make_figures file, specifying the desired name of the directory and a modifier p (0<p<=1) for variable sampling and number of trials. e.g. make_figures('TEST', 1);

This will run many simulations to evaluate the variability of the model and sample relevant variables. Note: if running all trials (i.e. p==1) simulation time on a single machine might take up to 40 hours. For simple testing, recomended p < 0.2. Recommended disk space for maximum sampling is ~50GB.

The data that was simulated from this code that is used to make the Figures in the repository is available upon request.
