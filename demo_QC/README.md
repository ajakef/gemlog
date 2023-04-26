# Gem Quality Control Workflow ("huddle test")
###
Geophysicists test their sensors using "huddle tests", where multiple instruments are placed in a location, subjected to similar conditions, and checked to ensure they record similar signals. The terminal command `verify_huddle_test` automates the process of comparing waveform, GPS, and state-of-health data among the instruments, resulting in a more accurate, comprehensive, and consistent evaluation of test results. 

### Installing the gemlog software
Follow the gemlog installation procedure [here](https://github.com/ajakef/gemlog/blob/main/README.md). Note that the resulting conda environment includes relevant packages like obspy and pandas.

### Getting Started
First, download the zip file in the list of files above; this contains all the inputs needed for this demonstration. Move the file to some convenient folder, unzip it, and cd into the project folder. Then, set your conda environment to whatever you set up during the installation using a terminal command like `conda activate gem`.

This demo includes data files from four different Gems run in a huddle test. Raw data files are not included because they are not inspected by this process; rather, only the derived data files are inspected (gps, state-of-health/metadata, and waveform/mseed). When running this procedure on your own test data, raw data must be converted before running `verify_huddle_test`.
```
Project_Folder/
|__gps/		# gps
|__metadata/	# state-of-health data
|__mseed/	# waveform data
```

### Run the automated test evaluation
Use this command in the terminal. The `-s` and `-x` options may be used to only process certain Gem serial numbers.

```verify_huddle_test```

If successful, a pdf report will be generated with several tables and figures showing performance of the Gems in the test and highlighting any problems.

```
Project_Folder/
|__gps/		# gps
|__metadata/	# state-of-health data
|__mseed/	# waveform data
|__figures/	# figures created for the report
|__reports/	# pdf reports showing test results
```

In this example, the tables highlight potential problems:

- A low battery warning for Gem #061. The battery voltage is low but within the Gem's operating range (for now). This is probably not a cause for concern for the Gem's battery sensor accuracy, but is intended to warn the user that unexpected behavior could occur on that Gem if the batteries are not changed.

- A high battery error for Gem #077. The battery voltage is too high, outside the Gem's operating range (1.7-15 V). It would be very unlikely for the Gem to operate at such a high battery voltage; rather, this should be interpreted that the battery sensor has a problem that must be fixed.

- A2 and A3 errors for Gems #061 and #065. A2 and A3 are auxiliary analog inputs, which can be connected to analog sensors of the user's choice (e.g., external thermistors, snow depth sensors, absolute pressure sensors, etc.). The A2 and A3 results for these Gems are considered suspicious because of long periods with exactly zero change, which could be due to firmware bugs (as in this case) or short circuits.

Explanations for all failure conditions and recommended troubleshooting tips are provided at the end of the report.

### How to conduct an effective huddle test
- Test 3-8 Gems at a time; fewer than 3 makes it hard to tell which Gem is wrong in a disagreement, and more than 8 is cumbersome. Ensure that all batteries have enough life that Gems will not drop out during the test. Start with empty memory cards.

- Place all Gems in a place with GPS reception, in as small of a circle as possible, with barbs facing inwards. This helps ensure that they all record similar signals, even for short-wavelength signals due to turbulence/wind. 

- If needed, give the Gems time to equilibrate to the ambient temperature with lids open. Avoid testing them in sunlight, which could affect internal temperatures differently and trigger false errors.

- Start all Gems at about the same time, ideally within one minute. Before leaving them, check the LED codes to ensure they all get a GPS signal.

- Run the Gems for around a couple hours to a couple days. Less than a couple hours can be insufficient to characterize their GPS performance. More than a couple days is cumbersome. Some infrasound event should occur during this period (doors closing work very well).

- Copy the data to your computer and convert it, resulting in a folder structure as shown above.

- Visually inspect some waveform data (this task is difficult to automate) using a viewer like PQL or the `Stream.plot()` method in obspy. A high-pass filter of at least 1 Hz is generally needed to find events efficiently. Find an event that is highly coherent among the different Gems' traces. Data should agree in timing (down to 0.01 sec); this can be checked visually by zooming in on a high-frequency event (e.g., door slam) and filtering above 20 Hz. Additionally, data should agree in waveform shape and spectrum shape.

- Finally, run `verify_huddle_test` as above, and check the report.