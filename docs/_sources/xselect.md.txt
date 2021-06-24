# Running xselect in a batch

## Creating xselect scripts from a template

The `xselect_command_create` function will create a set of xselect scripts from a template, for all the specified ObsIDs.  

First, create a file that has all of the names of the ObsID for example:

ObsID.dat:

    3627010101
    3627010201
    
This can be the same file as for `nicerl2_command_create`.

Create a folder that will contain the template and the created scripts. 

In this folder, place a `template.dat` file that has the xselect commands with the ObsID replaced by {}. For example, to create light curves with a pha cut, it could look like this:

    {}
    read events
    ../../{}/xti/event_cl
    ni{}_0mpu7_ufa.evt
    yes
    set bin 25
    set phaname PI
    filter pha_cut
    40
    200
    extract curve
    save curve
    ni{}
    exit
    no

Then in python:

    >>> import nicer.nicer as ni
    >>> ni.xselect_command_create('folder', filename='ObsID.dat', template_f='Template.dat', run_f='run.sh')

or simply
    
    >>> import nicer.nicer as ni
    >>> ni.xselect_command_create('folder') 
    
if using the default values and file locations. 

This will create a set of .xco scripts in the specified folder, one for each ObsID. A single script can be executed like so:

    > xselect @ni3627010101.xco

There is also a `run.sh` file that gets created in the folder, that contains all of the commands, that is useful to run everything at once:

    > bash run.sh

