# G2G analysis on EBV

This repository contains all scripts relevant for G2G analysis on Epstein–Barr virus (EBV) in HIV patients. For EBV data preparation, please see [gitlab.com/ezlab/vir_var_calling](https://gitlab.com/ezlab/vir_var_calling). Results are available here: https://doi.org/10.5281/zenodo.4011995. 

## Usage

In this order: 

1. `./make_1_computation.sh`: runs the G2G analysis (= many GWASs).
1. `./make_2_results.sh`: extracts top results. 
1. `./make_3_genomic_lambda.sh`: computes the genomic lambdas for all GWASs.
1. `./make_4_figures_tables.sh`: makes figures + tables for report.

_Note:_ All bash files have some parameters to choose in the header.

## Dependencies

- [FINEMAP v1.3.1](http://www.christianbenner.com/)
- [LDStore v1.1](http://www.christianbenner.com/)
- [PLINK2](https://www.cog-genomics.org/plink/2.0/)
- [GCTA v1.91.7](https://cnsgenomics.com/software/gcta/)
