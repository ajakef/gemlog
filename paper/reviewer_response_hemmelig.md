I’ve gone through the JOSS checklist and, bar where indicated below, everything looks in order to me. The manuscript set out the motivations for—and capabilities of—the gemlog package in a clear manner.

On a purely stylish note, would it be possible to highlight the package name and command line utility names with inline code formatting in the manuscript?
*Done.*

The reference list seems thorough, though my knowledge of the infrasound field is somewhat limited. The Data Report (line 115 of the original proof, Dannemann Dugick and Bowman) should probably read “TurboWave I and II…” as opposed to “TurboWave i and II”.
*Done.*

You might indicate somewhere (probably in the documentation, rather than setting it in stone in the JOSS paper) a rough definition of ‘low-cost’, based on the prices of the components used.
*I agree that setting it in stone is unwise as it tends to change. I added a link to the Gem website where info on purchasing can be found. Like most (all?) geophysical instrument manufacturers that I know of, my prices and availability are dynamic enough that I currently only provide prices by request via quotes, though I do intend to begin posting estimated prices publicly on the website this year.*

General comments
Who is TLSatterwhite? They have made some contributions to the code in the repository. Besides that, the contribution and authorship looks good to me.
*TLSatterwhite is the third author.*

You might consider releasing the code and assigning it its own DOI using something like Zenodo.


Per the GitHub Community Standards tab (under Insights), the repository is missing a few items:

    Code of conduct
    Contributing
    Security policy
    Issue templates
    Pull request template

The repository readme does include (good + clear) information about contribution guidelines, so you might consider putting these into the corresponding file as well just for the sake of completeness. I don't think you need all of the above, it was just a [sic]

*Added a CODE_OF_CONDUCT.md file (the standard Contributor's Covenant) and, as suggested, moved the contributing info from the README to a dedicated CONTRIBUTING.md. Given the small scale of the project, the limited ability/intent to issue patches for old versions, and the small number of issues and PRs, we think that including a security policy, issue templates, and PR template is more formal than is needed or helpful in this case.*


You might consider incorporating the documentation (e.g. for installation, for the demos, etc) into the existing readthedocs. You could then add some of the information from the paper regarding the problem the software is designed for/target audience to the landing page of the readthedocs site.
*This is a good idea that I will implement in a future version.*


How long does it take to convert, say, a 1 month campaign, per station? Just a ballpark figure would be useful (minutes, hours, days?). I get a rough estimate on the order of 10 minutes per station.
*We added this info to the README.md and to the paper (on the order of 10 seconds per station-day).*

Installation
The installation process was smooth and clearly documented. I was able to install gemlog both directly from the Python Package Index (using pip install gemlog) and also downloading the source code and installing from there (pip install . from the source directory). The authors uses a GitHub Action to build and release the package to the Python Package Index.

Tests
Again, the authors use a GitHub Actions to automatically test and lint any changes made and pushed to either a branch of the main repo, or a fork. However, in order to run these tests locally on my own machine, I had to install pytest (and also had to clone the source repo). You might consider adding some information to the documentation detailing what someone might need to do to test any changes made locally, without having to push changes to a remote repo.
*Done; we added more explicit info on pytest to Installation.md and CONTRIBUTING.md.*

Tests all ran locally without any issues, though there were a number of warnings that seemed to be about the implementation of the tests themselves.


Demos
In general, you might consider converting these markdown files to a series of Jupyter notebooks, which can be hosted and run in the browser using a (free) service like mybinder. This would also make it straightforward to include these notebooks in the readthedocs site. Besides that, they cover the basic applications of the package well. I have some minor, less code-related, comments that you might consider things to consider in the future.
*Done*

Demo 1
How might a user determine whether a file is bad?
*I added some explanation.*

Is there a reason why you can’t specify a station_info.txt file as an input to the gemconvert tool that will populate the SEED headers at this stage? Rather than
*That's a good idea and it's on the to-do list.*

summarize_gps indexes from 0 rather than 1—presumably an i rather than i + 1 in the print statement.
*Fixed.*

Demo 2
Great concept, love to see it. Something we’ve been meaning to do for some time with a number of the QC scripts we have developed over the years for seismic network data validation/huddle tests/pre-archive QC! The explainers at the end of the report are excellent, but I think should ultimately have a dedicated section in the online docs?
*I added an explanation section to the workflow notebook.*

Point the reader to the data convert tutorial as a reminder.
*Done.*

061 BATTERY WARNING – this doesn’t seem to be correct? Claiming 2.81 V is within 0.5 V of limit (1.7 V).
*That was a bug; it's fixed now.*

Is there a way you suggest validating whether data agree in timing? Visually assessing the waveform similarity is fine, but validating down to 0.01 s is a little more difficult.
*It's actually not hard if you have a high-frequency event like a door slam and filter the data above 20 Hz. I added this info to the demo text.*

Comments on the code
Below I’ve listed some general comments that I consider important enough to warrant addressing prior to approving this submission. There is a subsequent section of comments that are more just things to consider as tasks to be tackled in the future, time permitting.


I think it would be worth providing a set of aliases for the command line tools that are all prefixed by gem. This would make things more consistent and avoid any potential namespace clashes (though unlikely). Further to this, you might consider adding something like gemlog as an entry point that lists all of the entry point options?
*Good idea; done.*

See previous comment about providing the demos as Jupyter notebooks.
*Done.*

The default library package pathlib is preferable over os.path. I’ve not tested the package on Windows, but I did note a number of instances of hardcoded path delimiters (/) that might fail if someone tries to use the package on Windows. No requirement to support Windows, but could be useful to be aware of this.
*I refactored the code to completely replace os with pathlib without changing functionality.*

--------------------------
Python thoughts—entirely optional
The comments below go beyond what I think would be required to accept this in JOSS:

*Thank you for these notes. I have applied some of these, noted as follows. Others are intended to be added in future versions.*

Apply some form of autoformatting, e.g. black (can be installed with pip install black). I see there is linting in the GH Actions, but passing everything through black often helps me organise my source code. This can be configured to, for example, allow 100 character lines (per the README) by adding a pyproject.toml file to the repo, containing:

[tool.black]
line-length = 100
*Thank you for this suggestion; I've run black on the code and set up the config file you provided.*

The import section in many of the files is a bit disorganised, often with a number of unnecessary imports. I don’t think there’s a definitive PEP standard for how to format this part of a Python source file, but often something like the following is used:

import standard_lib1
from standard_lib2 import foo

import third_party_lib1
import third_party_lib2

import local_package

Many(/all) of the top-level functions are documented, but many of the underlying/internal functions are not. There seems also to be an inconsistent level of documentation across the package (presumably just as habits changed through time/different individual contributions).

Many bare exceptions – might make things easier in short term, but might create latent issues that are hard to pick up on/track down.

Mixed use of getopt and argparse—I think the latter is slightly nicer to work with/idiomatic.
*I agree and am using argparse for new code. I am leaving the getopt code in place because it's already written, works well, and doesn't cause maintenance problems.*

fstrings make for cleaner print statements, and also implicitly handle conversion of variables to their string representation.
*I changed most old formatted strings with % to fstrings.*

Some function names don’t follow PEP conventions (snake_case)—a quick example is in gem_cat.py with the AppendFile function.

*I found that the functions AppendFile, PlotAmp, CheckDiscontinuity are not the standard snake_case and renamed them accordingly.*

I decided to run through some of the above suggestions for one of the command line tools (gem_cat, as this seemed to only have a few uses). I validated the outputs remain the same using the refactored version using the example data provided in demo 3. I can either attach it as a comment below, or open an Issue—no expectation to merge it (since it’s not been tested against any edge cases), more included to highlight some of the above suggestions.