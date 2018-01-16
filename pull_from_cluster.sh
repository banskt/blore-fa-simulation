#!/bin/bash
rsync -av --exclude-from 'rsync-exclude.txt' cluster1:quasi_laplace_gwas/simulated_phenotype/ simulated_phenotype/
