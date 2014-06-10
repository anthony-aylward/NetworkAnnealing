#!/bin/bash
#PBS -l nodes=30:ppn=8

# Run the network annealing process. 

# During the annealing process, each parameter will either be free to move or be held 
# fixed. This choice is specified by how the parameters are entered as arguments. For
# example, --size 100 indicates a network of fixed size 100, while --freesize 100 
# indicates an initial network size of 100 that will change during the annealing process.

# This script assumes that network_annealing.py and seed_seqs.fas are contained in the
# directory ~/NetworkAnnealing/, while ~/HIVClustering/ is a clone of the eponymous git 
# repository.

/opt/python-3.3.1/bin/python3 ~/NetworkAnnealing/network_annealing.py \
--sequences ~/NetworkAnnealing/seed_seqs.fas \
--fasta ~/NetworkAnnealing/sampled_seqs.fas \
--tn93 ~/NetworkAnnealing/inferred_network.csv \
--size 500 \
--freedays 7 \
--freerate 0.0007 \
--freelineages 10 \
--freesplit 0.02 \
--freerandom 0.2 \
--subset 65 \
--freebias 0.1 \
--freesampling 30 \
--freeburst 2

exit
