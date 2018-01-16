#!/bin/bash
rsync -av --exclude-from 'rsync-exclude.txt' simulated_phenotype cluster1:quasi_laplace_gwas/
