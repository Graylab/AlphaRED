# AlphaRED
AlphaFold-initiated replica exchange protein docking 

## Pipeline

Currently, this approach is distributed as a step-wise pipeline involving sequence-to-structure prediction, analysis of prediction accuracy, and docking. 
For structure prediction, we have equipped AlphaFold, however, ESMFold and OmegaFold structures could also be used as all predictive methods deposit structure confidence values (e.g., pLDDT) in the B-factor column of the generated models.
For docking, we use ReplicaDock 2.0 with in-built upgrades to select flexible ("mobile") residues based on pLDDT values.


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