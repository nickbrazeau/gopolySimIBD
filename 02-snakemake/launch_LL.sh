#! /bin/bash

ROOT=/proj/ideel/meshnick/users/NickB/Projects/gopolySimIBD # root directory for project (non-scratch)
WD=/work/users/n/f/nfb/Projects/gopolySimIBD/ # scratch directory
NODES=1028 # max number of cluster nodes
WAIT=30 # number of seconds to wait for files to appear, absorbing some file system latency

snakemake \
	--snakefile $ROOT/02-snakemake/run_polySimIBD.snake \
	--configfile $ROOT/02-snakemake/config_batch_LL.yaml \
	--printshellcmds \
	--directory $WD \
	--cluster $ROOT/02-snakemake/launch_lite.py \
	-j $NODES \
	--rerun-incomplete \
	--keep-going \
	--latency-wait $WAIT \
	--dryrun -p
