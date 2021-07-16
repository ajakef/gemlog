# Setting Up PQL
PQL (PASSCAL Quick Look) is a great tool for quickly perusing seismic/infrasound data. It is not a great tool for data analysis--use a scientific programming language for that instead. However, to my knowledge, PQL is the fastest way to scan data in order to look for events and obtain basic spectral information.

### PQL Installation (Ubuntu):
This procedure downloads an rpm file and processes it. It might be best to make a new folder in ~/Downloads and run the code from there. You can delete the download after running everything.
```
wget https://www.passcal.nmt.edu/ftp/software/pql/linux/x86_64/PQL-2010-246.x86_64.rpm
sudo apt install alien
sudo alien --to-deb *rpm
sudo dpkg -i --force-overwrite *deb
```

Then, open ~/.bashrc with the text editor of your choice, and add these lines to the end. This tells your system where to look to find the right PQL command.
```
PATH=$PATH:/opt/passcal/bin/
alias pql='pql -l' # uses setting that makes it resizeable on small screens
```
then run `source ~/.bashrc` to implement your changes.

### PQL Installation (Mac)
Download the .pkg file from the latest release here: https://www.passcal.nmt.edu/ftp/software/passoft/osx/

Find the downloaded file and double-click to install. If the "unidentified source" error comes up, right click on the installation package, and choose "open".... Another warning will pop up, but you can hit "open" again, and it will start the installation process.  

Find out what your default shell is; newer macs typically use zsh (config file .zshrc), and older macs typically use bash (config file ~/.bash_profile). Open the hidden file `~/.zshrc` or `~/.bash_profile` using your favorite text editor, and add the following lines to the very end. This tells your system where to look to find the PQL command.
```
PATH=$PATH:/opt/passcal/bin/
alias pql='pql -l' # uses setting that makes it resizeable on small screens
```
Then, run this command in the terminal to implement your changes: `source ~/.zshrc` or `source ~/.bash_profile`.

Finally, run the command `pql` in the terminal. If you get a display error, you might need to run it through XTerm instead (which may require installing XQuartz).


### Configuring PQL
PQL's default configuration is not very convenient in my typical applications. I recommend changing the settings as follows.

First, change the settings on the left side of the main screen for each of the following tabs.
* "Trace" tab: set window scale to "window", set time axis to "absolute"
* "Magnify" tab: set window scale to "window"
* "Spectra" tab: set window scale to "window"

Then, click the "Controls" button near the bottom left, and make the following changes:
* "General" tab: click "Edit" next to "Sort Definitions" to 1 absolute start time, 2 network, 3 station, 4 location, 5 channel
* "Overlay" tab: change color #3 (originally yellow) to anything that's more visible against pink background
* "Spectra" tab: set both amplitude and frequency to log scale

Finally, click "Set Defaults" in the top left to save all the settings you just changed, then "Continue" to return to the main screen.

### Using PQL
#### Starting PQL
Start PQL from the terminal using a command like `pql data/*`, where `data/*` is a list of files or a wildcard expresion that tells pql what files to open. If your file names include date, station, and channel info, you can even run commands like this:
```
# open files starting at 15:00, March 15 2020, station ANMO, channel HDF (infrasound)
pql data/2020-03-14T15*ANMO*HDF*
```

#### Trace Tab
You should see several traces on the screen. If you have many traces, you might need to click "Next" or "Previous" in the top right to navigate through them. To look at a set of traces in more detail, click their trace info on the left side of the plot to select them, and then click twice on the plot to set the plot beginning and end times (this will draw vertical lines). Then, click Magnify on the left, which takes you to the Magnify tab.

#### Magnify Tab
You can zoom in again by clicking two times and then "Magnify". You can also zoom in and out with the scroll wheel on your mouse.

Filtering may be the most important thing you do in the Magnify tab. Infrasound especially is affected by low-frequency noise that needs to be filtered out. To create a filter, click the "Filters" menu on the left, select "Manage", click the "NEW" button, and enter a name for your filter. For a basic high-pas filter that works for many applications, set "High Pass Filter" poles to 4 and "Cutoff Freq" to 1; leave the low-pass filter parameters alone. Click "SAVE", then look for your filter name in the "Filters" menu. You can toggle the filter on and off with the radio buttons below the "Filter menu.

You can also measure time intervals on the Magnify screen. If you left-click two points on a trace while holding Ctrl, PQL will tell you the exact times you clicked, as well as the interval between them.


#### Spectra Tab
If you click the Spectra tab from the Magnify screen, it'll show you the spectra of whatever data are currently shown on the Magnify screen. If you have a filter implemented, then your spectra will also be filtered. Note that clicking two times on the Magnify screen without actually zooming in has no effect on the window used for calculating spectra.

You can use Ctrl-left click to find precise frequencies in the spectra. Consider changing the "log/linear" settings in the spectra to make the plot more useful to you. Finally, set "Overlay" to "On" in order to plot the spectra on top of each other, making comparison easier.
