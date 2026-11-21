# Aspen settings file (CONFIG.TXT)

The Aspen logger reads its settings from a text file named `CONFIG.TXT` in the top folder of its microSD card. This
page explains every setting, what the logger does with mistakes, and how to check a card before a deployment. It
applies to Aspen firmware 1.0.0-34 and later; see [Older firmware](#older-firmware) for earlier versions.

## Quick start

A card without `CONFIG.TXT` records all four channels at 200 samples per second, gain 128, with the GPS cycled every
30 minutes. To change anything, put a `CONFIG.TXT` on the card. The simplest useful one is a single line, for example
infrasound only at 100 samples per second:

```
gains = 128               # channel 1 only
sample_rate = 100
```

The standard file, with every setting listed and explained, is
[aspen_CONFIG.TXT](https://github.com/ajakef/gemlog/blob/main/gemlog/data/aspen_CONFIG.TXT) (it is also in your
gemlog installation's `gemlog/data` folder). Save it on the card as `CONFIG.TXT` and edit the lines you need.

Before a deployment, check the card (see [Checking a card](#checking-a-card-before-a-deployment)):

```
check_aspen_config E:/
```

## How the file is written

- One setting per line: `name = value`. Spaces around `=` and around commas don't matter.
- Anything after `#` is a comment, and the logger ignores it.
- Order doesn't matter, and every line is optional. A setting that isn't in the file takes its default (listed below).
  To switch a line off without deleting it, put a `#` in front of it.
- A value of `-` also means "the default".
- Names are not case sensitive (`Sample_Rate` works). Blank lines, Windows line endings, and the invisible marker some
  editors put at the start of a file are all fine.
- Some settings take words as well as numbers: `gps_mode = cycled`, `on`, `off` or `adaptive`, and any on/off setting
  takes `on`/`off`, `yes`/`no` or `true`/`false`.
- Save the file as plain text. In Windows Notepad, choose "UTF-8" or "ANSI", not "Unicode" (UTF-16), which the logger
  can't read.
- A line may be up to 199 characters long, not counting its comment.

## Settings

### What to record

| Setting | Values | Default | What it does |
|---|---|---|---|
| `gains` | 1 to 4 of 0, 1, 2, 4, 8, 16, 32, 64, 128, separated by commas | `128, 128, 128, 128` | Gain of channels 1-4 (channel 1 is the infrasound sensor; 2-4 the geophone's Z, N, E). 0 turns a channel off. Listing fewer than four turns the rest off, so `gains = 128` records infrasound only. At least one channel must be on. |
| `sample_rate` | 100, 200, 400, 800, 1000, 2000, 4000 | 200 | Samples per second. At 4000, only channel 1 records (and Bluetooth stays off). |

### GPS

| Setting | Values | Default | What it does |
|---|---|---|---|
| `gps_mode` | 1 or `cycled`, 2 or `on`, 3 or `off`, 4 or `adaptive` | 1 | `cycled`: the GPS wakes up every `gps_cycle_minutes` to correct the clock, then sleeps (saves power). `on`: always on. `off`: no GPS (no absolute time). `adaptive`: cycled, with the cycle length adjusted to keep timing errors small. |
| `gps_cycle_minutes` | 1-255 | 30 | Minutes between GPS wake-ups (modes 1 and 4). |
| `gps_quota` | 1-255 | 20 | How many GPS time pulses the clock collects at each wake-up before the GPS sleeps again. |
| `gps_settle_s` | 0-120 | 20 | Seconds after each wake-up before GPS time pulses are used (the first ones can be off by a millisecond). |

### Names written in the data files (SEED codes)

The data files are miniSEED files named by these codes. Codes are letters and digits only.

| Setting | Values | Default | What it does |
|---|---|---|---|
| `network` | up to 2 letters/digits | blank | Network code. |
| `station` | up to 5 letters/digits | the logger's serial number | Station code. |
| `location` | up to 2 letters/digits | blank | Location code. |
| `channels` | 1 to 4 codes of exactly 3 letters/digits, separated by commas; `-` for one channel's default | depends on the sample rate (below) | Channel codes for channels 1-4. Channels not listed keep their codes. |
| `channel1` ... `channel4` | one code of 3 letters/digits | as above | The same, one channel at a time. |

Default channel codes follow the sample rate, as the SEED standard's band codes do:

| Sample rate | Channel 1 | Channels 2-4 |
|---|---|---|
| 100-200 | HDF | EHZ, EHN, EHE |
| 400-800 | CDF | DHZ, DHN, DHE |
| 1000-4000 | FDF | GHZ, GHN, GHE |

Two channels that are both on can't have the same code (the logger changes the second one to `XX2`, `XX3` or `XX4`).

### Files, serial port, Bluetooth

| Setting | Values | Default | What it does |
|---|---|---|---|
| `file_minutes` | 1-1440 | 1440 | Length of each file in minutes. Files start at whole multiples of it in UTC, so 1440 gives one set of files per UTC day. |
| `serial_data` | 0/1 (`off`/`on`) | 0 | Also send the data over the USB cable (for the real-time plotters, which normally turn this on themselves). |
| `serial_gps` | 0/1 | 0 | Also send the GPS lines over the USB cable. |
| `serial_health` | 0/1 | 1 | Send the health lines (battery, temperature, timing) over the USB cable. |
| `health_s` | 1-60 | 5 | Seconds between health lines. |
| `bluetooth` | 0/1 | 1 | 1: a phone with the Aspen app can connect and show the data. 0: the radio stays off. |

### Cell accessory (daily check-in)

For the optional cell/satellite check-in accessory. Leave these out if you don't have one.

| Setting | Values | Default | What it does |
|---|---|---|---|
| `checkin_minutes` | 0-30 | 0 | Longest daily check-in, in minutes. 0 means no accessory. |
| `checkin_hour` | 0-23 | 12 | UTC hour of the daily check-in. |

## Mistakes

The logger never refuses to record because of a mistake in `CONFIG.TXT`. Instead, for each line it can't use, it:

1. writes a message in its log file (the `.log` file next to the data), with the line number, for example
   `Error: CONFIG.TXT line 9: sample_rate = 300 isn't allowed (allowed: 100, 200, 400, 800, 1000, 2000 or 4000; default 200); ignored`
2. ignores that line, so that setting keeps its default (or its value from an earlier line), and
3. flashes the red LED 6 times, over and over (error code 6).

The log also lists every setting line it read (`CONFIG.TXT line 8: gains = 128`) and, after startup, one `Config:`
line with the settings it ended up using.

What counts as a mistake:

- a line without `=` (for example `sample_rate 400` or `sample_rate: 400`);
- a name that isn't a setting (the message suggests the closest one: `did you mean sample_rate?`);
- a value that isn't allowed, or two values where one belongs (`sample_rate = 400 800`);
- a setting given twice (the later line is used);
- a line longer than 199 characters;
- `run_minutes` or a `test_` setting in `CONFIG.TXT` (they only work in config sequence files, below);
- channels 2-4 turned on at 4000 samples per second;
- two channels that are on with the same code;
- a file saved as UTF-16 ("Unicode"): the logger uses every default.

Infrasound-only loggers (a wire inside grounds the logger's type pin) always leave channels 2-4 off, whatever the file
says; that isn't reported as a mistake.

## Checking a card before a deployment

A mistake shows up in the field only as red LED flashes, so check the card before you go. Both checkers use the
logger's own rules and give the same messages as the log would.

**With gemlog** (versions after 1.11.0), from a terminal:

```
check_aspen_config E:/CONFIG.TXT
check_aspen_config E:/            # the whole card: CONFIG.TXT and any config sequence files
```

or from Python:

```python
import gemlog
result = gemlog.check_aspen_config('E:/CONFIG.TXT')  # prints a report
result.ok          # False if the logger would report a mistake
result.settings    # every setting as the logger would use it
```

The report lists the logger's error messages, warnings about things the logger accepts silently but that are probably
not what you meant (for example a line in an old positional file that isn't a number), and the settings the logger
would use. The terminal command ends with exit status 1 if any file has a mistake.

**With the Aspen phone app:** connect the logger by USB, open the Check tab, and look at the "CONFIG.TXT on the card"
section. "Fix CONFIG.TXT" opens a form that saves corrected settings to the card.

## Saving settings from the app (the SVC command)

The app's "Save to logger" sends the new settings to the logger, which writes them into `CONFIG.TXT` itself. The old
file is kept on the card as `CONFIG.BAK`. The logger edits the old file line by line: your own comments and lines
stay where they were, each setting is changed in place, a setting whose new value is its default is switched off with
a `#` (so it keeps following the default), and settings the file didn't have are added at the end. A line the logger
would reject (a misspelled name, for example) is switched off with a `#` and a note saying why. An old positional file
is replaced by the standard file with your settings.

## Config sequence files: several settings one after another (testing)

For unattended tests, such as noise or power at several sample rates, put files named `CONFIG01.TXT`, `CONFIG02.TXT`,
... on the card (`CONFIG` followed by 1-3 digits, `.TXT`, any case). Each is an ordinary settings file plus a line
saying how long to run it:

```
# CONFIG01.TXT
run_minutes = 60          # 1-10080
sample_rate = 100
gains = 128
gps_mode = off
```

At startup the logger uses the lowest-numbered one **instead of** `CONFIG.TXT` (so repeat anything that should apply
to every run, such as the station code). After `run_minutes` it renames the file `CONFIG01_DONE.TXT`, closes its
files, and restarts with the next one. With none left it goes back to `CONFIG.TXT`. A restart in the middle of a run
(power cut, unplugging) runs the same file again from the start. To run a sequence again, rename the `_DONE` files back.

A file with such a name but no valid `run_minutes` line is ignored (the log says so), so a copy of an old
`CONFIG.TXT` saved as `CONFIG1.TXT` doesn't take over a deployment. Sequence files may also hold the logger's
`test_...` fault-injection settings, which are for bench testing the firmware; they aren't allowed in `CONFIG.TXT`.

## Older firmware

The logger's startup log names its firmware version (`Firmware: 1.0.0-34_...`).

- **Before 1.0.0-30**: only the original positional format (below). A named file in such a logger is silently ignored
  and every setting takes its default. Update the firmware or use a positional file.
- **1.0.0-30 to 1.0.0-33**: named files work, with a few differences: words (`on`, `cycled`, ...) aren't accepted; a
  mistake in the first setting line (`gains: 128`) makes the logger read the whole file as positional and use every
  default without a message; a sequence file without `run_minutes` runs for 60 minutes; `-` after an earlier value of
  the same setting doesn't restore the default; and the default gains are 128, 32, 32, 32.
- **1.0.0-34 and later**: as described on this page.

### The original positional format

Older cards hold one value per line, in a fixed order, with comments as labels. The logger still reads them: a file
whose first setting line starts with a number (or `-`) is positional. The order is: channel 1-4 gains, sample rate,
gps_mode, gps_cycle_minutes, gps_quota, use_compression (no longer used, but the line must be there), network code,
station code, location code, channel 1-4 codes (`-` for the default), checkin_minutes, checkin_hour. Lines may be left
off the end. Because the position decides the meaning, deleting or adding one line shifts every value after it, often
without any message; `check_aspen_config` warns when a line's comment names a different setting than the one it is read
as. For example:

```
128 # channel 1 gain
64 # channel 2 gain
64 # channel 3 gain
64 # channel 4 gain
400 # sample rate
1 # gps_mode
30 # gps_cycle_minutes
20 # gps_quota
1 # use_compression (not used, but keep this line so the lines below are read in order)
- # network code
- # station code
- # location code
```
