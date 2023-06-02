The gemlog package is something that I've been working with for a few years now so I'm happy to see it progress to this state. I'm coming at this review from a user's perspective since I don't have the programming chops of the other reviewer. Overall the package works as expected and with the help of the checklist and I think the submission is mostly ready to publish. I've just included a few nit picky things below that might be useful below. Nice work!

The "project folder" in installation is ill-defined (especially if they don't read the "Pre-processing Workflow" first)
*Clarified "project folder" in the installation instructions and also in the data conversion workflow notebook.

This warning is thrown with gemconvert:

/opt/miniconda3/envs/gem/lib/python3.10/site-packages/gemlog/core.py:612: FutureWarning: Calling int on a single element Series is deprecated and will raise a TypeError in the future. Use int(ser.iloc[0]) instead
config = {key:int(line[key]) for key in list(line.keys())[1:]}

*Weirdly, I can't reproduce this warning. But I looked up the cause and re-wrote the offending line as recommended, and it still works.*

Elevation API:
https://developers.google.com/maps/documentation/elevation/overview
https://open-elevation.com/

*Thanks for these recommendations. open-elevation.com seems to be exactly what I need.*

Glad to see that the development on the huddle test has come along. This is a really useful function. When I ran it against an old huddle test, things worked as expected, however the report contained 25 pages of errors (the same two) repeated each minute:

SN ['158'] recorded temperatures greater than 1 on either side of the temperature median 21.8 on 2022-06-08 02:36:00. Total temperature range = 2.4 (alpha)

SN ['105', '158'] recorded temperatures greater than 1 on either side of the temperature median 18.8 on 2022-06-07 07:53:00. Total temperature range = 4.8 (alpha)

This is a little strange because all 6 instruments in the huddle test were within 1.5 degrees the whole time (according to the graph). In any case, it might be worth splitting the error outputs to a separate file from the plots and tables. There is information at the end of the report that is otherwise buried when there is an output with a bunch of errors in it.

*If you send me the old huddle test data, I can try to debug this.*

There is no obvious statement of need in the documentation, only in the paper.

*I reproduced the statement of need in the main README.md file.*