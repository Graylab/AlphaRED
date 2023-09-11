# Structure prediction using AlphaFold

For structure prediction, we default to AlphaFold. For easier (and faster) access, ColabFold is recommended, particularly if you would like to predict a small number of naturally occuring proteins.
In case you would like to predict large numer of sequences of complexes, would like batch processing or any of the other advanced applications, LocalColabFold is a useful resource.

## Resources

AlphaFold: https://github.com/google-deepmind/alphafold

ColabFold: https://github.com/sokrypton/ColabFold

LocalColabFold: https://github.com/YoshitakaMo/localcolabfold


## Examples

Here, we provide a SLURM script to use local colabfold for a bunch of fastas ( `structure_prediction.slurm`). The results alongwith inputs and outputs are deposited in the examples folder. Alternatively, once you have loaded the conda modules, you can essentially run the following commands to generate a structure:

```
OUTPUT_DIR=results/
FASTA_DIR=fastas/
mkdir -p $OUTPUT_DIR

# We run colabfold for generating 5 predictions using the latest version of AlphaFold_v3
colabfold_batch $FASTA_DIR $OUTPUT_DIR --num-models 1 --model-type alphafold2_multimer_v3
```

Once you have a structure from any of the predicted services (AlphaFold, ColabFold, LocalColabFold), you can proceed to the next step of the pipeline. Alternatively, if you have structure from ESMFold or OmegaFold that records the structural confidence in the b-factor column, please feel free to test our following docking pipeline with those structures. 
