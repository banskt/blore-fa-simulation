#!/bin/bash
rsync -av --exclude-from 'rsync-exclude.txt' simulation cluster1:quasi_laplace_gwas/
rsync -av --exclude-from 'rsync-exclude.txt' analysis   cluster1:quasi_laplace_gwas/
