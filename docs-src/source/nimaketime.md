# Running nimaketime in a batch

The `nimaketime_command_create` function creates a run.sh containing the command for running the `nimaketime` task. To learn about this task, visit https://heasarc.gsfc.nasa.gov/lheasoft/ftools/headas/nimaketime.html and https://heasarc.gsfc.nasa.gov/lheasoft/ftools/headas/maketime.html. The need to run this task may arise when one wants to exclude/include a set of Good Time Intervals (GTIs) while processing the nicer data. The file created by this task can be supplied to nicerl2 using the input gtifiles. Note that, the gtifiles input will always correspond to the GTIs that one wants to include. However,`nimaketime_command_create` creates the output gti file in accordance with the user's choice to include/exclude a set of specified time intervals, which can then be supplied to nicerl2.

First, create a file that has all of the names of the ObsID for which you want to run nimaketime, for example:

ObsID.dat:

    3627010101
    3627010201
    
This can be the same file as for `nicerl2_command_create`, or `xselect_command_create`. Let us call the directory that contains the directories with names given by the Observation IDs as folder. Inside each folder/obsID directory (where obsID=observation IDs), create a file named as `obsID_time.txt`. This is the file that contains the time intervals that the user wants to include/exclude while processing the data. An example of the content of such a file is given below:

2.25215685e+08 2.25216562e+08
2.25221269e+08 2.25222122e+08

These times are in the same format as that for the GTIs in the event file. For simplicity, let us assume that the file ObsID.dat is inside the folder directory (though this is not necessary). Also, suppose you want to remove the above time intervals.

Then in python:

    >>> import nicer.nicer as ni
    >>> ni.nimaketime_command_create(folder, filename='ObsID.dat', include=False, run_f='run.sh')

This will create a file called `run.sh` in the folder directory. In addition, this function also writes down the expression used in nimaketime inside each folder/obsID directory under the name of `obsID_nimaketime_expr.txt`. To run nimaketime, open a terminal and go to folder. Now type:

    > bash run.sh

In case you want to make a gtifile than includes the time intervals instead of excluding it, set include=True when you call the `nimaketime_command_create' function.


