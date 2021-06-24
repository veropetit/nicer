import numpy as np

def nicerl2_command_create(filename='ObsID.dat', output='nicerl2_command.sh'):
    """
    Creates a bash script that contains command to run Nicerl2 for each ObsID in the list (read from filename='ObsID.dat'.
    
    The options for the Nicerl2 are listed as ${FLAGS} in the script. By defining this environement variable, the same script can be used for different processing options.
    
    The cleaned event list are renamed to cl_${TYPE}, where the environement variable ${TYPE} contains the desired suffix. This is so that the Nicerl2 can be used with different options, for comparison (the default behavior of Nicerl2 is to clobber).
    
    See the nicerl2_source.sh template to see how to set the environement variables.
    
    :param filename: ('ObsID.dat'), the filename that contains the list of ObsID.
    :param output: ('nicerl2_command.sh'), the filename of the output bash script
    :rtype: None. Creates a bash script called output='nicerl2_command.sh'
    
    """
    
    obsid = np.loadtxt(filename, dtype=str)
    
    f = open( output  ,'w')
    
    for obs in obsid:

        f.write( 'nicerl2 indir={} clobber=YES ${{FLAGS}}\n'.format(obs)  )
        f.write( 'mv {}/xti/event_cl/ni{}_0mpu7_cl.evt {}/xti/event_cl/ni{}_0mpu7_cl_${{TYPE}}.evt\n\n'.format(obs, obs, obs, obs) )

            
    f.close()



def xselect_command_create(folder, filename='ObsID.dat', template_f='template.dat', run_f='run.sh'):
    '''
    This function reads in a xselect template and substitute the ObsID numbers
    
    In the template, the obsID is replaced with {}.
    The template is expected to be located in the specified folder.
    A xselect .xco script will be created for each ObsID.
    Also, a run file (default 'run.sh') will be created, to execute all of the .xco script at once, with a `bash run.sh`
    
    :param folder: the folder containing the template, and that will contain the created scripts
    :param filename: ('ObsID.dat') name of the file containing the list of ObsID
    :param template_f: ('template.dat') name of the file containing the template to use (expected location in `folder`)
    :param run_f: ('run.sh') name of the file containing all of the run commands (created in `folder`)
    '''

    obsid = np.loadtxt(filename, dtype=str)


    template_file = open(folder + '/' + template_f, 'r')
    template = template_file.readlines()
    template_file.close()

    run = open( folder + '/'+ run_f, 'w')

    for obs in obsid:

        f = open( '{}/ni{}.xco'.format(folder,obs)  ,'w')
        run.write('xselect @ni{}.xco\n'.format(obs))
        
        for line in template:
        
            # This will place the obs in any {} in the template.
            # It there isn't any, it just returns the string.
            # Also here, I do not need a /n to break the line,
            # as it is already included in the string from reading the template file
            f.write( line.format(obs, obs, obs, obs, obs, obs, obs, obs) )
            # MAKE SURE that there are more "obs" in there than the mak # of {} in one line.
            
        f.close()
        
    run.close()


