# revised_epistasis_pipeline

make scripts folder in home directory. This will contain (atlas.tar python.tar and two of node.py files).


make epistasis_results folder in home directory (empty folder)
make epistasis_condor_out folder in home directory (empty folder)

make data folder in home directory. This will contain the plink files to run. (also with xxx.pheno.txt file)

python epistasis_submit_DAGman_v7.py xxx -g 600 -m 1

run    ./release.sh 2048 1024 600 &  (if want to reset request memory to 2g disk to 1g and waiting time as 10 min)

*xxx is the plink file name
