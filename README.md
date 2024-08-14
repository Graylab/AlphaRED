# AlphaRED
AlphaFold-initiated replica exchange protein docking is a pipeline to predict protein complex structures from sequences. This pipeline follows a sequence to structure to complex paradigm while employing AlphaFold for structure prediction followed by ReplicaDock2.0 for protein-protein docking. 

For more details, please check out our paper [here](https://doi.org/10.1101/2023.07.28.551063):
Harmalkar A, Lyskov S, Gray JJ, "Reliable protein-protein docking with AlphaFold, Rosetta and replica-exchange", bioRxiv, July 2023. 

## Pipeline

Currently, this approach is distributed as a step-wise pipeline involving sequence-to-structure prediction, analysis of prediction accuracy, and docking. 
For structure prediction, we have equipped AlphaFold, however, ESMFold and OmegaFold structures could also be used as all predictive methods deposit structure confidence values (e.g., pLDDT) in the B-factor column of the generated models.
For docking, we use ReplicaDock 2.0 with in-built upgrades to select flexible ("mobile") residues based on pLDDT values.

### How to use AlphaRED?

1. Go to `AlphaRED/pipeline/1-structure_prediction/` to set-up local installation of ColabFold. This would be used for structure prediction. Alternative structure prediction tools such as ESMFold, OmegaFold, etc could also be used in place of AlphaFold based on the preference of the user.
2. Next, go to `AlphaRED/pipeline/2-prediction_analysis/` to analyze the state of the prediction, i.e. is the binding site correctly identified or not? If not, perform global docking. Otherwise, perform local docking and refinement.
3. Finally, in `AlphaRED/pipeline/3-docking/`, based on the binding mode identified, perform the respective docking analysis. If you perform global docking first, the top decoys from global docking are down-selected for local docking. If you perform local docking directly, select the top-scoring decoys (based on Interface scores) for side-chian refinement and relaxation.

To run docking simulation, use command:
##   mpirun $ROSETTA_DIR/source/bin/rosetta_scripts.mpi.linuxgccrelease @flags_replica_dock
where $ROSETTA_DIR is your Rosetta's directory

## Installation

The installation instructions for each stage of the pipeline are detailed in its respective directory in `pipeline` folder. 
To navigate, the folder hierarchy is as follows:

* benchmark
    * difficult_targets
    * medium_targets
    * rigid_targets
* pipeline
    * 1-structure_prediction
    * 2-prediction_analysis
    * 3-docking
* utilities


## Web Server

We are working to set up a web-server for structure prediction and docking. This service would be shortly available on ROSIE. Updates coming soon!

## Bug reports

If you run into any problems while using AlphaRED, please create a [Github issue](https://github.com/Graylab/AlphaRED/issues) with a description of the problem and the steps to reproduce it.

## References

Please use the following references to cite our work and corresponding literature:

1. Structure prediction

* [AlphaFold](https://github.com/Graylab/AlphaRED/issues) and [AlphaFold-multimer](https://www.biorxiv.org/content/10.1101/2021.10.04.463034v1)

* [ColabFold](https://www.nature.com/articles/s41592-022-01488-1)

2. Docking

* [AlphaRED](https://doi.org/10.1101/2023.07.28.551063)

* [ReplicaDock 2.0](https://doi.org/10.1371/journal.pcbi.1010124)
