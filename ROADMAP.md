# What is next in QMRITools

The main aim of the current roadmap is to optimize automatization of the processing and analysis pipelins.
For this the BIDS dat structure is adoped and optimized for muscle also known as Muscle-Bids.

Currently the toobox allows processing from Dicom to merged fully processed datasets for T2-mapping, Dixon based imaging, and DTI/IVIM.
Next steps will be muscle segentation and per muscle analysis of quantitative values and tractogrphy.

## Upcomming code projects

### Muscle faciculation analsyis

- Optimize faciculation detection
- Optimize faciculation analysis

### Create common tensor space

- Tensor align
  - create methods for aligning between and within subject tensors
  - registration based on distance maps
  - registration based on tensor information
- Muscle atlas cration

### Tractography

- create tract quantification and analysis tools
  - muscle force line definitions
    - pennation anlge definitions
    - muscle to model files
    - tract deformations

### Mask alignment

- tools for mask homogenization
  - Mask templating between subjects/cohorts

### Muscle Bids

- add MuscleBidsQualityReport
- add MuscleBidsAnalysisReport
- add facicualtion analysis in MuscleBidsProcess

### Muscle labeling

- Related to muscle atlas project
  - create muscle codes that lead to names and regions.
  - add function to translate segmentation to muscle codes
    - will allow automated labeling for analysis
    - will allow grouping of segmentation in groups

### Histology based simualations

- If histology is availible allow diffusion simulations of cell segmentation
  - based on the work of David Berry
  - Allow parameterized generation of realistic muscle cells
  - potentially add internal muscle organells
