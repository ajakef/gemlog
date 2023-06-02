# Gem Quality Control Workflow ("huddle test")
###
This notebook is an example workflow for testing Gems. It is *not* intended for routine field data conversion (for that, use the `demo_conversion` workflow).

Be sure to follow the gemlog installation procedure [here](https://github.com/ajakef/gemlog/blob/main/README.md). Note that the resulting conda environment includes relevant packages like obspy and pandas.

Geophysicists test their sensors using "huddle tests", where multiple instruments are placed in a location, subjected to similar conditions, and checked to ensure they record similar signals. The terminal command `verify_huddle_test` automates the process of comparing waveform, GPS, and state-of-health data among the instruments, resulting in a more accurate, comprehensive, and consistent evaluation of test results.
