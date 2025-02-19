# phot-d

**phot-d** is a tool for estimating stellar spectral types and distances. 
It was originally developed for [Tang et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...862..106T/abstract) 
and has since been adopted in other studies, including [Hunter et al. 2024](https://ui.adsabs.harvard.edu/abs/2024AJ....168..211B/abstract). 

The core algorithm for spectral type estimation is based on `photo-type` from 
[Skrzypek et al. 2015](https://ui.adsabs.harvard.edu/abs/2015A%26A...574A..78S/abstract). 
For further details, please refer to [Tang et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...862..106T/abstract).

## Installation and Usage

### 1. Set up the environment
To install the required dependencies, create a conda environment using the provided `environment.yml` file:

```bash
conda env create -f environment.yml
```

### 2. Activate the environment
Activate the phot-d environment:
```bash
conda activate phot-d
```

### 3. Run phot-d
Execute phot-d.py with the desired input file:
```bash
python phot-d.py <input_file>
```
To view available input options, use:
```bash
python phot-d.py -h
```

### 4. Example usage
You can run an example using:
```bash
python phot-d.py auto_query_example.csv -q
```

If you encounter any issues or have questions, feel free to open an issue in this repository or contact me at sytang@rice.edu.

