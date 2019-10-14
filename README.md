# revised_epistasis_pipeline


epistasis_submit.py  is used for re-run failed filed  tasks.



make scripts folder in home directory. This will contain (atlas.tar python.tar and two of node.py files).


make epistasis_results folder in home directory (empty folder)
make epistasis_condor_out folder in home directory (empty folder)

make data folder in home directory. This will contain the plink files to run. (also with xxx.pheno.txt file)

python epistasis_submit_DAGman_v7.py xxx -g 600 -m 1

 ./release.sh 2048 1024 600 &  (if want to reset request memory to 2g disk to 1g and waiting time as 10 min)


* python epistasis_submit_DAGman_v7_osg.py xxx -g 600 -m 1  (if run script in osg, but need to upload file to squid folder first)
* In order to run in osg, run file in chtc first and remove job and then run it in osg (the same file name, but using _osg.py to run).
* need to change the squid/username to your username in chtc .
* And to run jobs in osg, set connect project to  osg.UWMadison_Parks  before running python script.
* xxx is the plink file name
