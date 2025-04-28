# Codes for the project *Crystallography, Group Cohomology, and Lieb–Schultz–Mattis Constraints* 
# 

[![arXiv](https://img.shields.io/badge/arXiv-2410.03607-b31b1b.svg)](https://arxiv.org/abs/2410.03607)

This repository contains the following three parts

### GAP codes: 
#
Read both SpaceGroupCohomologyData.gi and SpaceGroupCohomologyFunctions.gi
Main command that outputs the mod-2 cohomology and LSM anomaly classes:
SpaceGroupCohomologyRingGapInterface(IT);
where IT is the numbering of the space group (1<=IT<=230).

### Mathematica codes:
#
Space_Group_Cohomology_Data.nb contains the standard inhomogeneous functions for the 1-,2-, and 3-cocycles (except for No. 225, 227, and 229)
Execute the "Preparations" section first.

### Sage codes:
#
Codes to obtain the explicit cochain representatives for degree-3 generators of several complicated space groups


## Reference

If this repository is useful for your research, please consider citing the [arXiv preprint](https://arxiv.org/abs/2410.03607):

```bibtex
@misc{liu2024crystallography,
    title={Crystallography, Group Cohomology, and Lieb–Schultz–Mattis Constraints},
    author={Chunxiao Liu and Weicheng Ye},
    year={2024},
    eprint={2410.03607},
    archivePrefix={arXiv},
    primaryClass={cond-mat.str-el}
}
```
