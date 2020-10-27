# Installation:
Regardless of whether you are an ordinary user or developer, using a conda environment (preventing dependency conflicts) is a very good idea. Install Miniconda from [here](https://docs.conda.io/en/latest/miniconda.html), or install Anaconda if you prefer.

### For Users
* First, create and activate a conda environment with dependencies (named "gem" in this example). This environment must always be activated when running `gemconvert`.
```
conda create -n gem python=3.8
```

* Then, activate that environment and install gemlog into it using pip:
```
conda activate gem
pip install gemlog
```
Note that this will install all of gemlog's dependencies. These include the ubiquitous numpy and pandas, along with the super-useful and strongly recommended obspy package for seismic/infrasound data processing.

### For Developers
* First, set up the conda environment with dependencies.
```
conda config --add channels conda-forge
conda create -n gem python=3.8 numpy obspy pandas matplotlib scipy cython pytest
```

* Next, from the right project folder, clone and install the gemlog python package from github. 
```
conda activate gem
git clone https://github.com/ajakef/gemlog/
pip install --upgrade ./gemlog
```

