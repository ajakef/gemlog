# Installation:

* First, create and activate a conda environment with dependencies (named "gem" in this example). This environment must always be activated when running gem2ms.
```
conda config --add channels conda-forge
conda create -n gem python=3.7.6 numpy obspy pandas matplotlib cython
conda activate gem
```

* Next, install the gemlog python package from github. This will create an executable gem2ms in a folder where it will be recognized.
```
pip install --upgrade https://github.com/ajakef/gemlog/archive/main.zip
```

* Get the syntax to run gem2ms.
```
gemconvert -h # print the help page
```

* Finally, run the demo [here](https://github.com/ajakef/gemlog/tree/main/demo), both to ensure that everything works on your system, and as an example of a typical workflow.

* If you've found a bug, please raise an issue on github. Please be sure you've installed all the dependencies correctly with conda, that the environment is activated on your computer, and that you're including information about your system and enough info for your issue to be reproduced. If you have an idea for an improvement, you can raise an issue for that too--especially if you're willing to implement it! gemlog has no dedicated employees and discretionary development time is limited in general.

# Contributing
User contributions are welcome! If you want to improve performance or address outstanding issues, please follow the guidelines below. Like any group that benefits from volunteers, we very much appreciate contributors but require that contributions follow this framework in order to be manageable.
##### Code requirements:
* New functions must have docstrings (follow the format used in, for example, read_gem).
* Variable names should have meaningful names. Single-letter exceptions: t is allowed for 'time', and p ('pressure') for infrasound time series.
* 100-character limit per line. This can be relaxed if it's burdensome for you because you're new to this sort of thing.

##### Workflow requirements for contributing:
* Set up the conda environment just like above, but add the testing module 'pytest' at the end: `pip install pytest`.
* Contributing code via github pull request is *required* unless, for a good reason, you have a complete and significant contribution and you are totally unable to share using github. Code contributed via other channels is less transparent and more burdensome on maintainers.
* Begin your contribution by raising an issue to get quick feedback from maintainers and users.
* Be sure you are editing the most current gemlog version! 
* After finishing your changes, commit them and push them to the gemlog fork on your github account. A suite of tests will run after the push. Your code must pass all tests before it can be considered for merging.
* Send a pull request when your code is ready for review. This is often at the point when you think you're done, but it can be earlier if you need to share your code to get specific feedback.

