# Create a file with the tags of the structures you would like to extract
# here for example, I am sorting based on the Interface score (I_sc) column and 
# extracting the top 5 structures and depositing them to tag_low file. 
# Uncomment the next line if you would prefer to do that or change column header accordingly.
# The next steps will fail if there is either no tag_file or if the tag_file is empty
# sort -n -k14 decoys_1.fsc | awk '{print $NF}' | head -n 5 > tag_low  # [Uncomment to extract from low res files]
# Top files based on Irms
# sort -n -k11 high_res.sc | awk '{print $NF}' | head -n 3 >> tag_low  # [Uncomment to extract from high res files based on Irms]
# Top files based on scores
# sort -n -k10 high_res.sc | awk '{print $NF}' | head -n 3 >> tag_low  # [Uncomment to extract from high res files based on Isc]

mkdir out_pdbs; 
# Extract the PDB files specified in the tag_low file. Note: Please change path to Rosetta.
for i in $(ls decoys_P_000?_traj.out); do echo $i; ~/Rosetta/demos/protocol_capture/replica_docking/scripts/extract_tagged_decoys.py $i tag_low > low$i; done; 
for i in $(ls lowdecoys_P_000?_traj.out); do echo $i; ~/Rosetta/main/source/bin/extract_pdbs.mpi.linuxgccrelease -in:file:silent_struct_type binary -out:path:pdb out_pdbs -out:pdb_gz -in:file:silent $i ; done
mv P*pdb out_pdbs;