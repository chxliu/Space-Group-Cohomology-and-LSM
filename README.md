# Codes for the project *Crystallography, Group Cohomology, and Lieb–Schultz–Mattis Constraints* 
# 

[![arXiv](https://img.shields.io/badge/arXiv-2410.03607-b31b1b.svg)](https://arxiv.org/abs/2410.03607)

## Content

This repository contains the following three parts

### GAP codes

Read both SpaceGroupCohomologyData.gi and SpaceGroupCohomologyFunctions.gi
Main command that outputs the mod-2 cohomology and LSM anomaly classes:
SpaceGroupCohomologyRingGapInterface(IT);
where IT is the numbering of the space group (1<=IT<=230).

### Mathematica codes

Space_Group_Cohomology_Data.nb contains the standard inhomogeneous functions for the 1-,2-, and 3-cocycles (except for No. 225, 227, and 229)
Execute the "Preparations" section first.

### Sage codes

Codes to obtain the explicit cochain representatives for degree-3 generators of several complicated space groups


## Reference

If this repository is useful for your research, please consider citing the [Scipost](https://www.scipost.org/SciPostPhys.18.5.161) article:

```bibtex
@Article{10.21468/SciPostPhys.18.5.161,
	title={{Crystallography, group cohomology, and Lieb–Schultz–Mattis constraints}},
	author={Chunxiao Liu and Weicheng Ye},
	journal={SciPost Phys.},
	volume={18},
	pages={161},
	year={2025},
	publisher={SciPost},
	doi={10.21468/SciPostPhys.18.5.161},
	url={https://scipost.org/10.21468/SciPostPhys.18.5.161},
}
```
