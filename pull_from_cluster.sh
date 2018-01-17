#!/bin/bash
rsync -av --exclude-from 'rsync-exclude.txt' cluster1:quasi_laplace_gwas/simulation/ simulation/
rsync -av --exclude-from 'rsync-exclude.txt' cluster1:quasi_laplace_gwas/analysis/   analysis/
