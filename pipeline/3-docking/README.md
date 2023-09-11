# Docking

For the final step of the pipeline, we initiate a docking protocol. Based on the docking accuracy, you could initiate either global docking or local docking. 

## Global Docking

Rigid global docking protocol for a search through the protein energy landscape to determine putative binding modes. Rigid search is performed with ReplicaDock 2.0. The input files are deposited in the `global_docking` directory with an example. 


## Local Docking

Flexible backbone local docking protocol for local refinement and backbone sampling. An update to the global docking protocol described prior, here we utilize backbone samping using Rosetta Backrub and KIC protocols. 


## Refinement

Side-chain packing and refinement of generated decoys to obtain high-resolution scores with Rosetta. Top models can be refined with Rosetta and obtained for comparisons. 