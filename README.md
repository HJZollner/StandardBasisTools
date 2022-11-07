# StandardBasisTools
 Tools to generate standardized basis spectra and LCM fitting results representations
 
We have released initial guidelines for a standardized data format to store basis spectra and LCM results. This standard will: 
- provide a common format for data sharing
- reduce workload to support file formats specific to a single density matrix simulation software
- increase transparency and reproducibility by adding relevant metadata
- increase between-software compatibility; and e) maximize the integration of MRS into existing software- and data-sharing efforts.

The development is spearheaded by the Code and Data Sharing Committee of the ISMRM MRS study group, encouraging community feedback and contributions via the community forum (forum.mrshub.org).

## Preliminary Example data
This repository currently contains an example analysis of an [HERCULES MRS dataset]( 10.1016/j.neuroimage.2018.10.002
). The data is provided in NIfTI-MRS format ([MRSHub forum thread](https://forum.mrshub.org/t/nifti-mrs-discussion-thread/443) or full [manuscript]( https://doi.org/10.1002/mrm.29418) for reference) and anlyzed using [Osprey](https://github.com/schorschinho/osprey).

The example script generates the following outputs:
- Processed data in NIfTI-MRS format including model provenance
- Linear-combination model results in NIfTI-MRS format following a preliminiary proposal to store LCM fitting results
- low-level HERCULES basis set in NIfTI-MRS format 

## Planned Features
- Generalize guidelines to leverage NIfTI-MRS to store low-level basis sets and LCM fitting results
- Introduce a data format to store high-level basis sets
- Write conversion scripts for LCM algorithms with discontinued support (LCModel)

## Contact & Feedback
If you are interested to contribute please contact [Helge J. Zöllner](mailto:hzoelln2@jhu.edu).

## Contributors
These efforts are spearheaded by the Code and Data Sharing Committee of the ISMRM MRS study group:

- Helge J. Zöllner (Johns Hopkins University, Baltimore, MD)
- Kelley M. Swanberg (Columbia University, New York, NY)
- John LaMaster (Technische Universität München, München, Germany)
- Antonia Kaiser (University of Amsterdam, Amsterdam, The Netherlands)
- Jamie Near (University of Toronto, Toronto, Canada)
- Candace Fleischer (Georgia Institute of Technology and Emory University, Atlanta, GA)
- Georg Oeltzschner (Johns Hopkins University, Baltimore, MD)

And by other interested MRS software developers:
- Brian J. Soher (University of Oxford, Oxford, United Kingdom)
- William T. Clarke (Duke University Medical Center, Durham, NC)




