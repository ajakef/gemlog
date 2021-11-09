# Installation:
Regardless of whether you are an ordinary user or developer, using a conda environment (preventing dependency conflicts) is a very good idea. Install Miniconda from [here](https://docs.conda.io/en/latest/miniconda.html), or install Anaconda if you prefer.

### For Users
* First, create and activate a conda environment with dependencies (named "gem" in this example). This environment must always be activated when running `gemconvert`.
```
conda create -y -n gem python=3.10
```

* Then, activate that environment and install gemlog into it using pip:
```
conda activate gem
pip install cython
pip install gemlog
```

Note that this will install all of gemlog's dependencies. These include the ubiquitous numpy and pandas, along with the super-useful and strongly recommended obspy package for seismic/infrasound data processing.

* If you will be doing any scientific computing in python besides converting data, consider installing a user-friendly computing environment; for example, ```pip install spyder ipython```.

* To update gemlog to a newer version, run ```pip install --upgrade gemlog```. This will not upgrade your python version. 

### For Developers
* First, set up the conda environment with dependencies.
```
conda config --add channels conda-forge
conda create -y -n gem python=3.10 numpy obspy pandas matplotlib scipy cython pytest
```

* Next, from the right project folder, clone and install the gemlog python package from github. 
```
conda activate gem
git clone https://github.com/ajakef/gemlog/
pip install --upgrade ./gemlog
```

If you are new to git, checkout the [git cheat sheet](https://github.com/ajakef/gemlog/tree/main/git_instructions.md).

## Supported Python versions:
Minimum 3.7. Python 2.7 (obsolete) is not supported and will not be supported in the future.

Expect a lag of at least several weeks between when a new Python version is released and when gemlog becomes available for it. This is because we need to wait for dependencies to start supporting the new Python version. 