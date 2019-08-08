# revised_epistasis_pipeline

make scripts folder in home directory. This will contain (atlas.tar python.tar and node.py).


make epistasis_results folder in home directory (empty folder)
make epistasis_condor_out folder in home directory (empty folder)

make data folder in home directory. This will contain the plink files to run. (also with xxx.pheno.txt file)

python epistasis_submit_DAGman_v7.py xxx -g 600 -m 1

*xxx is the plink file name
