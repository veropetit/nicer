# Running nicerl2 in a batch

## Creating a nicerl2_command.sh bash script

The `nicerl2_command_create` function will create an all purpose bash script that can be used in combination with an file containing enviroment variables to run `nicerl2` for all the ObsIDs with specific options. 

First, create a file that has all of the names of the ObsID for example:

ObsID.dat:

    3627010101
    3627010201

Then in a python shell:

    >>> import nicer.nicer as n
    >>> n.nicerl2_command_create(filename='ObsID.dat', output='nicerl2_command.sh')

will create a `nicerl2_command.sh` file that looks like this:

    nicerl2 indir=3627010101 clobber=YES ${FLAGS}
    mv 3627010101/xti/event_cl/ni3627010101_0mpu7_cl.evt 3627010101/xti/event_cl/ni3627010101_0mpu7_cl_${TYPE}.evt

    nicerl2 indir=3627010201 clobber=YES ${FLAGS}
    mv 3627010201/xti/event_cl/ni3627010201_0mpu7_cl.evt 3627010201/xti/event_cl/ni3627010201_0mpu7_cl_${TYPE}.evt

As you can see, for each ObsID in the list, `nicerl2` is executed, with options contained in the environement variable `$FLAG` (see below), and then the cleaned event list is renamed `ni3627010201_0mpu7_cl_${TYPE}.evt` where the environement variable `$TYPE` contains the suffix to apply. 

The reason being renaming the files is because `nicerl2` clobbers the cleaned event list. 

## Using the nicerl2_command.sh

To process data, you first need to define the two environment variables. I do this by creating a file `nicerl2_source.sh` that could look like this:

    
    # For default processing. 
    export TYPE="default"
    export FLAGS=""
    
    # For changing the SAA region
    #export TYPE="SAA"
    #export FLAGS="nicersaafilt=NO saafilt=YES"

Then to process the data:

    > source nicerl2_source.sh
    > bash nicerl2_command.sh

The end