# [Installation](https://github.com/ajakef/gemlog/tree/main/Installation.md) 
* Follow the directions [here](https://github.com/ajakef/gemlog/tree/main/Installation.md).

# Getting started after installing gemlog
* Click [here](https://sites.google.com/view/geminfrasound) for more info on the Gem Infrasound Logger, including information for loans and purchases.
* Most users will only need the `gemconvert` tool for converting data to standard formats. `gemconvert` is run from the terminal. This call will give you the syntax, options, and version number.
```
gemconvert -h # print the help page
```

* Run `pip install --upgrade gemlog` on the command line to update gemlog to the newest version.

* Run [this demo](https://github.com/ajakef/gemlog/tree/main/demo) to ensure that everything works on your system and as an example of a typical workflow. On a typical laptop, conversion may take on the order of 10 wall-clock seconds per day of data for one station (or, a unitless ratio of approximately 10^-4 between conversion time and data duration).

* If you have data files that lack GPS information (e.g., because they were recorded on a high-altitude balloon, indoors, or underground), AND you're ok with having imprecise sample timing, run [this demo](https://github.com/ajakef/gemlog/tree/main/demo_missing_gps) too.

* If you've found a bug, please raise an issue on github. Please be sure that the conda environment is activated on your computer, that the dependencies are installed correctly, and that you're including information about your system and enough info for your issue to be reproduced. If you have an idea for an improvement, you can raise an issue for that too--especially if you're willing to implement it! 

# Contributing
User contributions are welcome! If you want to improve performance or address outstanding issues, please follow the guidelines below. Like any group that benefits from volunteers, we very much appreciate contributors but require that contributions follow this framework in order to be manageable.
##### Code requirements:
* New functions must have docstrings (follow the format used in, for example, read_gem).
* Variable names should have meaningful names. Single-letter exceptions: t is allowed for 'time', and p ('pressure') for infrasound time series.
* 100-character limit per line.
* Following normal python convention, variable names use snake_case and class names use CamelCase.

##### Workflow requirements for contributing:
* Follow the installation instructions for developers in the link above.
* If you are new to `git`, check out the [git cheat sheet](https://github.com/ajakef/gemlog/tree/main/git_instructions.md).
* Contributing code via github pull request is *required* unless, for a good reason, you have a complete and significant contribution and you are totally unable to share using github. Code contributed via other channels (e.g., email) is less transparent and more burdensome on maintainers.
* Begin your contribution by raising an issue to get quick feedback from maintainers and users.
* Be sure you are editing the most current gemlog version! 
* After finishing your changes, commit them and push them to a gemlog fork on your github account. A suite of tests and style checks will run after the push. Your code must pass all automated checks before it can be considered for merging.
* Send a pull request when your code is ready for review. This is often at the point when you think you're done, but it can be earlier if you need to share your code to get specific feedback.

# Need help?
* Check the demonstrations mentioned above to see if it addresses your issue.
* Look through the issues in this github site to see if anyone else is having the same trouble you are.
* If those don't work, email project lead Jake Anderson at jacobanderson152@boisestate.edu.

# Statement of need
The Gem Infrasound Logger [@Anderson2018] is an approach to infrasound recording where the sensor and data logger are built into a single cable-free package that is easy to conceal and permits arbitrary sensor network geometries. Additionally, it has a low cost and weight, runs for months on alkaline batteries, and has a simple, fast installation procedure. These characteristics make it a good choice for temporary infrasound campaigns. By contrast, most campaigns that do not use the Gem use analog infrasound sensors that, via long cables, connect to multichannel data loggers built to record seismometers. This approach yields high-quality data but has several disadvantages: seismic data loggers are expensive, and sensor cables constrain the sensor network's geometry, make a station prone to animal damage and vandalism due to their visibility and exposure, and account for a large share of a recording site's budget for weight, bulk, and setup time. These disadvantages are especially acute for temporary recording campaigns (as opposed to permanent installations), which account for a large share of infrasound research.

Like many geophysical data loggers, the Gem writes data in a non-standard raw format intended to balance firmware simplicity and performance, human readability, and compact file size. Although it is a human-readable text format that an expert can read and understand on a line-by-line basis, data files consist of hundreds of thousands of lines with complicated formatting, meaning that reading it as a spreadsheet or data frame is impractical. Further, operations like clock drift corrections, data decompression, and conversion to standard file formats or classes must be performed to make the data accessible in standard visualization and analysis software. Users often need to convert data from 10 or more infrasound loggers spanning several weeks, meaning that thousands of files including billions of data points must be processed efficiently. The `gemlog` Python package is a cross-platform tool to facilitate data conversion, and is essential for all Gem infrasound logger users.

Some other geophysical data loggers used to record infrasound, for example the [DiGOS DataCube-3](https://digos.eu/seismology/) and Reftek [RT-130](https://www.passcal.nmt.edu/content/instrumentation/dataloggers/3-channel-dataloggers/reftek-rt-130-datalogger), have non-standard raw data formats that must be converted to standard formats by software distributed by the manufacturer. Like `gemlog`, they conduct data conversion as a simple command line operation. However, to the authors' knowledge, none of those software packages (or the raw data formats they convert) are open-source.