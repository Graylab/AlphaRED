# Prediction Analysis

Code to evaluate prediction accuracy. Here, we provide a simple python script to evaluate the confidence of docked complex prediction.


## Requirements

This script would require Biopython package (https://biopython.org/)
```pip install biopython```


## Evaluating model confidence

```
python interface_plddt.py prediction.pdb -partners A_B -cutoff 10
```

`prediction.pdb` is the PDB file generated from any of the structural prediction methods. The output of this would be a prompt indicating which docking method to use, global docking or local docking.

An example prompt is below:
```
Avg pLDDT 90.329
STATUS: AlphaFold has HIGH overall confidence in the structure
Interface pLDDT 64.472
STATUS: Docking site has lower confidence. Use GLOBAL Docking and refine using LOCAL docking.
```