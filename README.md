# [Installation](https://github.com/ajakef/gemlog/tree/main/Installation.md) 
* Follow the directions [here](https://github.com/ajakef/gemlog/tree/main/Installation.md).

# Getting started after installing gemlog
* Most users will only need the `gemconvert` tool for converting data to standard formats. `gemconvert` is run from the terminal. This call will give you the syntax, options, and version number.
```
gemconvert -h # print the help page
```

* Run `pip install --upgrade gemlog` on the command line to update gemlog to the newest version.

* Run [this demo](https://github.com/ajakef/gemlog/tree/main/demo) to ensure that everything works on your system and as an example of a typical workflow.

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