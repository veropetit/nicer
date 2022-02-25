import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
import os

###########################
###########################
# Section on running NICER tools in a batch

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

###########################
###########################
# General data information

def gti_info(gti, verbose=False):
    """
    Returns the duration of each GTI and the gap in between each GTI
    
    :param gti: a GTI object read from the GTI fits extension.
    :param verbose: (False) flag to print out the GTI information
    :rtype duration: array with the duration of each GTI (in sec)
    :rtype gap: array with the gap (in sec) between each GTI
    
    """

    n = gti.shape[0]
    duration = gti[:,1] - gti[:,0]
    gap = gti[1:,0] - gti[:-1,1]
    gap = np.append(0,gap)
    
    if verbose:
        print('{:10} {:10} {:10} {:10}'.format('start', 'stop', 'duration', 'gap'))
        for i in range(0,n):
            print('{:10.1f} {:10.1f} {:10.1f} {:10.1f}'.format(gti[i,0], gti[i,1], duration[i], gap[i] ))

    return(duration, gap)


###########################
###########################
# Spectra

def add_pha(filename, ax=None, **kargs):
    '''
    Add a Xselect PHA spectrum from a file to a graph, using the plt.step function
    
    :param filename: name of the xselect pha file
    :param ax: (None) the ax to which the spectrum should be added
    :param kargs: any additional parameters to pass to the step function.
    
    '''
    
    hdul = fits.open(filename)
    head = hdul[1].header
    data = hdul[1].data
    hdul.close()
    if ax is None:
        # If no ax is passed, use the last ax used or create one
        ax = plt.gca()
    spec_type=head['TTYPE2']

    if spec_type=='COUNTS':
        x   =data[spec_type]/head['EXPOSURE']
    if spec_type=='RATE':
        x   =data[spec_type]
    
    ax.step(data['CHANNEL']/100, x, **kargs )
    return(ax)



def load_pha(filename):
    '''
    Add a Xselect PHA spectrum from a file to a graph, using the plt.step function
    
    :param filename: name of the xselect pha file
    :param ax: (None) the ax to which the spectrum should be added
    :param kargs: any additional parameters to pass to the step function.
    
    '''
    
    hdul = fits.open(filename)
    head = hdul[1].header
    data = hdul[1].data
    hdul.close()
    spec_type=head['TTYPE2']

    if spec_type=='COUNTS':
        x   =data[spec_type]/head['EXPOSURE']
    if spec_type=='RATE':
        x   =data[spec_type]
    
    return(data['CHANNEL']/100,x)


###########################
###########################
# Light curves


def load_curve(filename, T0=0):
    """
    Load a light curve from Xselect.
    
    The light curve starts from the first GTI.
    
    :param filename: The name of the file containing the Xselect LC.
    :param T0: (0) Difference between the start of the GTI in the Xselect LC, and that of the desired "zero" time stamp.
    :rtype time_LC: = Xselect Time stamp + start of Xselect first GTI - T0.
    :rtype LC: The light curve
    
    """
    # Open the light curve
    hdul = fits.open(filename)
    LC = hdul[1].data
    gti_LC = hdul[2].data # get the gti
    hdul.close()
    # The light curve starts from the first GTI.
    # Need to add it back then reset to our T0
    time_LC = LC['TIME']-T0+gti_LC[0][0]
    
    return(time_LC, LC)

def add_lc(filename, T_ref=0, ax=None, **kargs):
    # T_ref is a reference 'NICER' time which we want as the zero of the time-axis
    # Open the light curve

    '''
    Add a Xselect lightcurve from a file to a graph using the plt.errorbar
    
    :param filename: name of the xselect .lc file
    :param T_ref: a reference 'NICER' time which we want as the zero of the time-axis
    :param ax: (None) the ax to which the lightcurve should be added
    :param kargs: any additional parameters to pass to the plt.errorbar function.
    
    '''

    hdul = fits.open(filename)
    LC = hdul[1].data
    LC_head=hdul[1].header
    hdul.close()
 
    lc_type=LC_head['TTYPE2']
    if lc_type=='COUNTS':
        x   =LC[lc_type]/LC_head['EXPOSURE']
    if lc_type=='RATE':
        x   =LC[lc_type]
    
    T0=LC_head['TIMEZERO']
    time_LC =LC['TIME']+T0-T_ref
    
    if ax is None:
        # If no ax is passed, use the last ax used or create one
        ax = plt.gca()

    ax.errorbar(time_LC, x, yerr=LC['ERROR'], **kargs )      
    return(ax)






###########################
###########################
# Background diagnostics


def plot_mk2(filename='ObsID.dat', obspath='.', clean_type='', LC_path='.',LC_clean_type='', PDF_path='.', merge_gap=300):
    """
    Make a graph of the enviromental condition plus Xselect light curves.
    
    Needed:
    
    * The ufa files: obspath/obsIDfolder/xti/event_cl/niObsID_0mpu7_ufa.evt -- they are the base for the global GTI information. The start of the first GTI defined t=0 (T0 in the code)
    * Some cleaned evt list: obspath/obsIDfolder/xti/event_cl/niObsID_0mpu7_cl(_type).evt. This will be use to mark the good GTIs.
    * The mask 3 files: obspath/obsIDfolder/auxil/ni'+obs+'.mkf3
    * A set of LC made from the UFA (LC-ufa-full/, LC-ufa-highcut/, LC-ufa-midcut/). Named ObsID.lc
    * A set of LC made from the clean evt list (LC-LC_CLEAN_TYPE-full/, LC-LC_CLEAN_TYPE-highcut/, LC-LC_CLEAN_TYPE-midcut/). Named ObsID.lc
    
    :param filename: ('ObsID.dat') File containing the list of ObsID
    :param obspath: ('.') Path where the ObsID folders are located.
    :param clean_type: ('') the suffix to use for the clean event files (niObsID_0mpu7_cl(CLEAN_TYPE).evt
    :param merge_gap: (300) If two raw GTIs are separated by less than merge_gap (in sec), they will be printed on the same page.
    :param LC_path: ('.') the path where the LC folders are located.
    :param LC_clean_type: ('') the suffix to use for the folder containing the cleaned LC (e.g. LC-LC_CLEAN_TYPE-full/)
    :param PDF_path: ('.') the path where the PDFs will be saved (e.g. pdf_path/ObsID.pdf)
    
    
    """

    obsid = np.loadtxt(filename, dtype=str)
    
    for obs in obsid:
            
        # Open the ufa event file to get the original GTI information
        filename = obspath+obs+'/xti/event_cl/ni'+obs+'_0mpu7_ufa.evt'
        hdul = fits.open(filename)
        gti_ufa = hdul['GTI'].data
        hdul.close()
        
        # get the GTI into a (6, 2) array to facilitate manipulation later on.
        n_gti_ufa = len(gti_ufa)
        gti_ufa = np.array([np.array(x) for x in gti_ufa])
        # defining a T0 as the start of the first ufa GTI.
        # The graph will use that T0 as the zero point of the time axis.
        T0 = gti_ufa[0,0]
        gti_ufa = gti_ufa - T0

        print()
        print('UFA GTI')
        print()
        
        # get the duration of each GTI, and the gap in between.
        gti_ufa_duration, gti_ufa_gap = gti_info(gti_ufa)

        # Open a cleaned event file to get the GTI information
        filename = '../'+obs+'/xti/event_cl/ni'+obs+'_0mpu7_cl'+clean_type+'.evt'
        hdul = fits.open(filename)
        gti_clean = hdul['GTI'].data
        hdul.close()

        print()
        print('CL GTI')
        print()
        
        # get the duration of each clean GTI, and the gap in between.
        gti_clean_duration, gti_clean_gap = gti_info(gti_clean)

        # Get the mk data
        filename = obspath+obs+'/auxil/ni'+obs+'.mkf3'
        hdul = fits.open(filename)
        head_mk = hdul[1].header
        data = hdul[1].data
        hdul.close()
        # Make t=0 the start of the first UFA GTI.
        MKtime = data['TIME']-T0

        # This piece of code merge together some GTIs that have very small gap
        # so that we can make a better use of a page size.
        # If the gap is rather large, then the next GTI will be printed on a new page.
        gti_plot = np.array( (1,2) )
        gti_plot.shape = (1,2)
        gti_plot[0,0] = gti_ufa[0,0] # set the start of the plit gti to that of the first GTI.
        k=0 # iteration variable for changes
        for i in range(0,n_gti_ufa):
            if gti_ufa_gap[i] > merge_gap: # if the gap with the last seg is more than 5 minutes
                gti_plot[k,1] = gti_ufa[i-1,1] #Set the stop of the k GTI to that of the previous segment
                gti_plot = np.vstack( [gti_plot, np.array( [gti_ufa[i,0], 0] ) ] ) # append a new plot_GTI
                                # and set the star time to the current GTI.
                k = k+1
        gti_plot[-1,1]=gti_ufa[-1,1] # set the stop of last gti
        n_gti_plot = k+1
              
        print()
        print('plot GTI')
        print()
        
        gti_plot_duration, gti_plot_gap = gti_info(gti_plot)

        # Open the light curve
        filename = LC_path+'/LC-ufa-full/ni'+obs+'.lc'
        time_ufa_full, LC_ufa_full = load_curve(filename, T0)
        
        filename = LC_path+'/LC-ufa-highcut/ni'+obs+'.lc'
        time_ufa_hi, LC_ufa_hi = load_curve(filename, T0)

        filename = LC_path+'/LC-ufa-midcut/ni'+obs+'.lc'
        time_ufa_mid, LC_ufa_mid = load_curve(filename, T0)

        filename = LC_path+'/LC'+LC_clean_type+'-full/ni'+obs+'.lc'
        time_full, LC_full = load_curve(filename, T0)

        filename = LC_path+'/LC'+LC_clean_type+'-midcut/ni'+obs+'.lc'
        time_mid, LC_mid = load_curve(filename, T0)

        filename = LC_path+'/LC'+LC_clean_type+'-highcut/ni'+obs+'.lc'
        time_hi, LC_hi = load_curve(filename, T0)


        with PdfPages('{}/ni{}_mk.pdf'.format(PDF_path, obs)) as pdf:
            for i in range(0,n_gti_plot):
            #for i in range(1,2):
                tmin = gti_plot[i,0]
                tmax = gti_plot[i,1]

                n = np.where( np.logical_and( MKtime >= tmin, MKtime <= tmax ) )
                fig, ax = plt.subplots(4,1, figsize=(8,10))
                #-----------------------
                # Overonly
                ax[0].scatter(MKtime[n], data['FPM_OVERONLY_COUNT'][n], s=3, zorder=500)
                ax[0].axhline(y=1.0, c='k', ls='--', label='abs cutoff')
                ax[0].plot( MKtime[n], 1.52*data['COR_SAX'][n]**(-0.633), c='orchid', lw=2, zorder=1000, label='overonly_exp' )
                ax[0].set_ylim(0,3)
                ax[0].set_ylabel('overonly')
                ax[0].set_xlabel('Time since T0 (s)')
                ax[0].set_title('{}, segment {}'.format(obs, i))
                ax[0].legend(loc=0)
                #-----------------------
                # Underonly
                ax[1].scatter(MKtime[n], data['FPM_UNDERONLY_COUNT'][n], s=3, zorder=500)
                ax[1].set_ylim(0,300)
                ax[1].set_ylabel('underonly')
                ax[1].axhline(y=200, c='k', ls='--', label='def cutoff')
                ax[1].set_xlabel('Time since T0 (s)')
                #ax[0].set_title('{}, segment {}'.format(obs, i))
                ax[1].legend(loc=0)
                #-----------------------
                # The light curves from the ufa
                #    full energy range
                n2 = np.where( np.logical_and(  time_ufa_full >= tmin, time_ufa_full <= tmax ) )
                ax[2].errorbar(time_ufa_full[n2], LC_ufa_full['RATE'][n2], yerr=LC_ufa_full['ERROR'][n2], fmt='.', ms=2, color="0.5", label='all keV', zorder=500 )
                #    mid energy range
                n2 = np.where( np.logical_and(  time_ufa_mid >= tmin, time_ufa_mid <= tmax ) )
                ax[2].errorbar(time_ufa_mid[n2], LC_ufa_mid['RATE'][n2], yerr=LC_ufa_mid['ERROR'][n2], fmt='.', ms=2, color="red", label='0.4-2.0 keV', zorder=500 )
                #    high energy range
                n2 = np.where( np.logical_and(  time_ufa_hi >= tmin, time_ufa_hi <= tmax ) )
                ax[2].errorbar(time_ufa_hi[n2], LC_ufa_hi['RATE'][n2], yerr=LC_ufa_hi['ERROR'][n2], fmt='.', ms=2, color="blue", label='12-15 keV', zorder=500 )

                ax[2].plot(MKtime[n], data['SAA'][n], c='k', label='SAA (1:yes, 0:no)' )

                ax[2].set_ylabel('UFA Count Rate')
                ax[2].set_ylim(0,3.5)
                ax[2].axhline(y=1.0, ls='--', c='0.5')
                ax[2].axhline(y=2.0, ls='--', c='0.5')
                ax[2].legend(loc=0)
                #-----------------------
     
                # The light curves from the current filtering
                #    full energy range
                n2 = np.where( np.logical_and(  time_full >= tmin, time_full <= tmax ) )
                ax[3].errorbar(time_full[n2], LC_full['RATE'][n2], yerr=LC_full['ERROR'][n2], fmt='.', ms=2, color="0.5", label='all keV', zorder=500 )
                #    mid energy range
                n2 = np.where( np.logical_and(  time_mid >= tmin, time_mid <= tmax ) )
                ax[3].errorbar(time_mid[n2], LC_mid['RATE'][n2], yerr=LC_mid['ERROR'][n2], fmt='.', ms=2, color="red", label='0.4-2.0 keV', zorder=500 )
                #    high energy range
                n2 = np.where( np.logical_and(  time_hi >= tmin, time_hi <= tmax ) )
                ax[3].errorbar(time_hi[n2], LC_hi['RATE'][n2], yerr=LC_hi['ERROR'][n2], fmt='.', ms=2, color="blue", label='12-15 keV', zorder=500 )

                ax[3].set_ylabel('Filtered Count Rate')
                ax[3].set_ylim(0,3.5)
                ax[3].axhline(y=1.0, ls='--', c='0.5')
                ax[3].axhline(y=2.0, ls='--', c='0.5')
                ax[3].legend(loc=0)


                #-----------------------
                # Overplotting the clean GTIs
                for j in range(0,n_gti_clean):
                    ttmin = gti_clean[j,0]
                    ttmax = gti_clean[j,1]
                    if np.logical_and( ttmin >= tmin, ttmin <= tmax ):
                        ax[0].axvspan(ttmin, ttmax, alpha=0.1, color='green', zorder=100)
                        ax[1].axvspan(ttmin, ttmax, alpha=0.1, color='green', zorder=100)
                        ax[2].axvspan(ttmin, ttmax, alpha=0.1, color='green', zorder=100)
                        ax[3].axvspan(ttmin, ttmax, alpha=0.1, color='green', zorder=100)

                for item in ax:
                    item.set_xlim(tmin, tmax)
                plt.tight_layout()
                pdf.savefig(fig)  # or you can pass a Figure object to pdf.savefig
                plt.close()
                
        return

#####################

def plot_mkf(filename='ObsID.dat', obspath='.',mkf_ext='.mkf', clean_type='', LC_ufa_path='.',LC_clean_path='.', PDF_path='.', merge_gap=300):
    """
    Make a graph of the enviromental condition plus Xselect light curves.
    
    Needed:
    
    * The ufa files: obspath/obsIDfolder/xti/event_cl/niObsID_0mpu7_ufa.evt -- they are the base for the global GTI information. The start of the first GTI defined t=0 (T0 in the code)
    * Some cleaned evt list: obspath/obsIDfolder/xti/event_cl/niObsID_0mpu7_cl(_type).evt. This will be use to mark the good GTIs.
    * The mask 3 files: obspath/obsIDfolder/auxil/ni'+obs+'.mkf3
    * A set of LC made from the UFA (LC-ufa-full/, LC-ufa-highcut/, LC-ufa-midcut/). Named ObsID.lc
    * A set of LC made from the clean evt list (LC-LC_CLEAN_TYPE-full/, LC-LC_CLEAN_TYPE-highcut/, LC-LC_CLEAN_TYPE-midcut/). Named ObsID.lc
    
    :param filename: ('ObsID.dat') File containing the list of ObsID
    :param obspath: ('.') Path where the ObsID folders are located.
    :param mkf_ext: The extension of the filter file, it can also be an array of size equal to the number of Observation IDs in the ObsID.dat file.
    :param clean_type: ('') the suffix to use for the clean event files (niObsID_0mpu7_cl(CLEAN_TYPE).evt
    :param LC_ufa_path: ('.') the path where the file containing the ufa LC filenames and labels, is located. It can be None (no lc will be plotted), a single string, or a list/array.
    :param LC_clean_path: ('.') the path where the file containing the clean LC filenames and labels, is located. It can be None (no lc will be plotted), a single string, or a list/array.
    :param LC_clean_type: ('') the suffix to use for the folder containing the cleaned LC (e.g. LC-LC_CLEAN_TYPE-full/)
    :param PDF_path: ('.') the path where the PDFs will be saved (e.g. pdf_path/ObsID.pdf)
    :param merge_gap: (300) If two raw GTIs are separated by less than merge_gap (in sec), they will be printed on the same page.
    
    """

    obsid = np.loadtxt(filename, dtype=str)
    if len(np.shape(obsid))==0: #if there is just one observation ID
        obsid   =np.array([obsid])
    
    mkf_ext_check=0 #a variable that tells whether the mkf_ext is a string or is a list/array
    if isinstance(mkf_ext, str)==True:
        mkf_ext_check   =1

    LC_ufa_path_check=1 #a variable that tells whether the LC_ufa_path is a string or is a list/array    
    if isinstance(LC_ufa_path, np.ndarray)==True or isinstance(LC_ufa_path, list)==True:
        LC_ufa_path_check=0

    LC_clean_path_check=1 #a variable that tells whether the LC_clean_path is a string or is a list/array
    if isinstance(LC_clean_path, np.ndarray)==True or isinstance(LC_clean_path, list)==True:
        LC_clean_path_check=0    

    clean_type_check=0
    if isinstance(clean_type, str)==True:
        mkf_ext_check   =1

    pdf_path_check=0
    if isinstance(PDF_path, str)==True:
        pdf_path_check   =1

    
    for N in range(len(obsid)): 
        obs=obsid[N]
            
        # Open the ufa event file to get the original GTI information
        os.chdir(obspath)
        filename = obs+'/xti/event_cl/ni'+obs+'_0mpu7_ufa.evt'
        hdul = fits.open(filename)
        gti_ufa = hdul['GTI'].data
        hdul.close()
        
        # get the GTI into a (x, 2) array to facilitate manipulation later on.
        n_gti_ufa = len(gti_ufa)
        gti_ufa = np.array([np.array(x) for x in gti_ufa])
        # defining a T0 as the start of the first ufa GTI.
        # The graph will use that T0 as the zero point of the time axis.
        T0 = gti_ufa[0,0]
        gti_ufa = gti_ufa - T0
        print('Obs ID and T0 are respectively: ',obs,T0)
        print()
        print('UFA GTI')
        print()
        
        # get the duration of each GTI, and the gap in between.
        gti_ufa_duration, gti_ufa_gap = gti_info(gti_ufa)

        # Open a cleaned event file to get the GTI information
        if clean_type_check==0:
            my_clean_type=clean_type[N]
        else:
            my_clean_type=clean_type

        filename = obs+'/xti/event_cl/ni'+obs+'_0mpu7_cl'+my_clean_type+'.evt'
        hdul = fits.open(filename)
        gti_clean = hdul['GTI'].data
        hdul.close()

        # get the GTI into a (x, 2) array to facilitate manipulation later on.
        n_gti_clean = len(gti_clean)
        gti_clean = np.array([np.array(x) for x in gti_clean])
        gti_clean = gti_clean - T0 # The zero point is the start of the UFA ****

        print()
        print('CL GTI')
        print()
        
        # get the duration of each clean GTI, and the gap in between.
        gti_clean_duration, gti_clean_gap = gti_info(gti_clean)

        # Get the mk data        
        if mkf_ext_check==0:
            mkf_extension   =mkf_ext[N]
        else:
            mkf_extension   =mkf_ext
        
        
        filename = obs+'/auxil/ni'+obs+mkf_extension
        hdul = fits.open(filename)
        head_mk = hdul[1].header
        data = hdul[1].data
        hdul.close()
        # Make t=0 the start of the first UFA GTI.
        MKtime = data['TIME']-T0

        # This piece of code merge together some GTIs that have very small gap
        # so that we can make a better use of a page size.
        # If the gap is rather large, then the next GTI will be printed on a new page.
        gti_plot = np.array( (1,2) )
        gti_plot.shape = (1,2)
        gti_plot[0,0] = gti_ufa[0,0] # set the start of the plit gti to that of the first GTI.
        k=0 # iteration variable for changes
        for i in range(0,n_gti_ufa):
            if gti_ufa_gap[i] > merge_gap: # if the gap with the last seg is more than 5 minutes
                gti_plot[k,1] = gti_ufa[i-1,1] #Set the stop of the k GTI to that of the previous segment
                gti_plot = np.vstack( [gti_plot, np.array( [gti_ufa[i,0], 0] ) ] ) # append a new plot_GTI
                                # and set the start time to the current GTI.
                k = k+1
        gti_plot[-1,1]=gti_ufa[-1,1] # set the stop of last gti
        n_gti_plot = k+1
              
        print()
        print('plot GTI')
        print()
        
        gti_plot_duration, gti_plot_gap = gti_info(gti_plot)

        # Open the light curve
        if LC_ufa_path_check==0:
            ufa_path    =LC_ufa_path[N]
        else:
            ufa_path    =LC_ufa_path

        if LC_clean_path_check==0:
            clean_path    =LC_clean_path[N]
        else:
            clean_path    =LC_clean_path  
        
        num_lc_panel  =0
        count_ufa=0
        if ufa_path!=None: 
            num_lc_panel+=1
            os.chdir(ufa_path)
            ufa_lc_filenames_labels=np.loadtxt('ni'+obs+'_ufa_lc_info.txt',dtype=str) #niobsIDlc_file_info.txt contains the filenames for the lcs to be overplotted and the corresponding labels

            if len(np.shape(ufa_lc_filenames_labels[0]))==0:
                ufa_lc_filenames,ufa_lc_labels  =np.array([ufa_lc_filenames_labels[0]]),np.array([ufa_lc_filenames_labels[1]])
            else:
                ufa_lc_filenames,ufa_lc_labels  =ufa_lc_filenames_labels[:,0],ufa_lc_filenames_labels[:,1]

            time_ufa,LC_ufa=[],[]
            for j in range(len(ufa_lc_filenames)):
                filename = ufa_lc_filenames[j]
                time, LC = load_curve(filename, T0)
                time_ufa.append(time)
                LC_ufa.append(LC)
                count_ufa+=1
                del(time)
                del(LC)
                del(filename)
            
            
        count_clean=0
        if clean_path!=None: 
            num_lc_panel+=1
            os.chdir(clean_path)
            clean_lc_filenames_labels=np.loadtxt('ni'+obs+'_clean_lc_info.txt',dtype=str) #niobsIDlc_file_info.txt contains the filenames for the lcs to be overplotted

            if len(np.shape(clean_lc_filenames_labels[0]))==0:
                clean_lc_filenames,clean_lc_labels  =np.array([clean_lc_filenames_labels[0]]),np.array([clean_lc_filenames_labels[1]])
            else:
                clean_lc_filenames,clean_lc_labels  =clean_lc_filenames_labels[:,0],clean_lc_filenames_labels[:,1]

            time_clean,LC_clean=[],[]
            for j in range(len(clean_lc_filenames)):
                filename = clean_lc_filenames[j]
                time, LC = load_curve(filename, T0)
                time_clean.append(time)
                LC_clean.append(LC)
                count_clean+=1
                del(time)
                del(LC)
                del(filename)

        if pdf_path_check==0:
            my_PDF_path =PDF_path[N]
        else:
            my_PDF_path =PDF_path

        with PdfPages('{}/ni{}_mk.pdf'.format(my_PDF_path, obs)) as pdf:
            for i in range(0,n_gti_plot):
            #for i in range(1,2):
                tmin = gti_plot[i,0]
                tmax = gti_plot[i,1]

                n = np.where( np.logical_and( MKtime >= tmin, MKtime <= tmax ) )
                fig, ax = plt.subplots(2+num_lc_panel,1, figsize=(8,10))
                #-----------------------
                # Overonly
                ax[0].scatter(MKtime[n], data['FPM_OVERONLY_COUNT'][n], s=3, zorder=500)
                ax[0].axhline(y=1.0, c='k', ls='--', label='abs cutoff')
                ax[0].plot( MKtime[n], 1.52*data['COR_SAX'][n]**(-0.633), c='orchid', lw=2, zorder=1000, label='overonly_exp' )
                ax[0].set_ylim(0,3)
                ax[0].set_ylabel('overonly')
                ax[0].set_xlabel('Time since T0 (s)')
                ax[0].set_title('{}, segment {}'.format(obs, i))
                ax[0].legend(loc=0)
                #-----------------------
                # Underonly
                ax[1].scatter(MKtime[n], data['FPM_UNDERONLY_COUNT'][n], s=3, zorder=500)
                ax[1].set_ylim(0,300)
                ax[1].set_ylabel('underonly')
                ax[1].axhline(y=200, c='k', ls='--', label='def cutoff')
                ax[1].set_xlabel('Time since T0 (s)')
                #ax[0].set_title('{}, segment {}'.format(obs, i))
                ax[1].legend(loc=0)
                #-----------------------
                # The light curves from the ufa
                #    full energy range
                cmap=plt.get_cmap('rainbow')
                if count_ufa>0:
                    ymax=0
                    col_arr=np.linspace(0,0.9999,count_ufa)
                    for count in range(count_ufa):
                        my_label=ufa_lc_labels[count].replace('_',' ')
                        n2 = np.where( np.logical_and(  time_ufa[count] >= tmin, time_ufa[count] <= tmax ) )
                        ax[2].errorbar(time_ufa[count][n2], LC_ufa[count]['RATE'][n2], yerr=LC_ufa[count]['ERROR'][n2], fmt='.', ms=2, color=cmap(col_arr[count]), label=my_label, zorder=500 )
                        if max(LC_ufa[count]['RATE'][n2])>ymax:
                            ymax=max(LC_ufa[count]['RATE'][n2])

                    ax[2].plot(MKtime[n], data['SAA'][n], c='k', label='SAA (1:yes, 0:no)' )
                    ax[2].set_ylabel('UFA Count Rate')
                    ax[2].set_ylim(0,ymax+0.5)
                    ax[2].axhline(y=1.0, ls='--', c='0.5')
                    ax[2].axhline(y=2.0, ls='--', c='0.5')
                    ax[2].legend(loc=0)
                #-----------------------
     
                # The light curves from the current filtering
                if count_clean>0:
                    ymax=0
                    col_arr=np.linspace(0,0.9999,count_clean)
                    for count in range(count_clean):
                        my_label=clean_lc_labels[count].replace('_',' ')
                        n2 = np.where( np.logical_and(  time_clean[count] >= tmin, time_clean[count] <= tmax ) )
                        ax[num_lc_panel+1].errorbar(time_clean[count][n2], LC_clean[count]['RATE'][n2], yerr=LC_clean[count]['ERROR'][n2], fmt='.', ms=2, color=cmap(col_arr[count]), label=my_label, zorder=500 )
                        if max(LC_clean[count]['RATE'][n2])>ymax:
                            ymax=max(LC_clean[count]['RATE'][n2])

                    ax[num_lc_panel+1].set_ylabel('Filtered Count Rate')
                    ax[num_lc_panel+1].set_ylim(0,ymax+0.5)
                    ax[num_lc_panel+1].axhline(y=1.0, ls='--', c='0.5')
                    ax[num_lc_panel+1].axhline(y=2.0, ls='--', c='0.5')
                    ax[num_lc_panel+1].legend(loc=0)


                #-----------------------
                # Overplotting the clean GTIs
                for j in range(0,n_gti_clean):
                    ttmin = gti_clean[j,0]
                    ttmax = gti_clean[j,1]
                    if np.logical_and( ttmin >= tmin, ttmin <= tmax ):
                        ax[0].axvspan(ttmin, ttmax, alpha=0.1, color='green', zorder=100)
                        ax[1].axvspan(ttmin, ttmax, alpha=0.1, color='green', zorder=100)
                        for num in range(num_lc_panel):
                            ax[2+num].axvspan(ttmin, ttmax, alpha=0.1, color='green', zorder=100)
                        
                for item in ax:
                    item.set_xlim(tmin, tmax)
                plt.tight_layout()
                pdf.savefig(fig)  # or you can pass a Figure object to pdf.savefig
                plt.close()
                
    return


##################
##################
## The following function is useful if one wants to use the nimaketime command to create a gti
## It constructs a string of expression, written to a file outfile, that can be provided to the nimaketime expr input
## The string is written in a file outfile, and also returned directly
## The infile is expected to have a list of tiem intervals, an example is given below:

##################

def nimaketime_expr(infile,outfile,include=False):
    # include determines whether the time intervals given in infile are to be included or excluded in the gtifile
   
    data=np.loadtxt(infile)
 
    if len(np.shape(data[0]))>0:
        t_start,t_end=data[:,0],data[:,1]
    else:
        t_start,t_end   =np.array([data[0]]),np.array([data[1]])

    f=open(outfile,'w')
    n   =len(t_start)
    x=''
    for i in range(n):
        if i==0:
            if include==True:
                f.write('(time>='+str(t_start[i])+' && time<='+str(t_end[i])+')')
                x+='(time>='+str(t_start[i])+' && time<='+str(t_end[i])+')'
            else:
                f.write('!((time>='+str(t_start[i])+' && time<='+str(t_end[i])+')')
                x+='!((time>='+str(t_start[i])+' && time<='+str(t_end[i])+')'
        else:
            f.write(' || (time>='+str(t_start[i])+' && time<='+str(t_end[i])+')')
            x+=' || (time>='+str(t_start[i])+' && time<='+str(t_end[i])+')'

    if include==False:
        f.write(')')
        x+=')'
            
    f.close()
    return x


def nimaketime_command_create(folder,filename='ObsID.dat', include=False, run_f='run.sh'):

    '''
    This function creates the command for the nimaketime task.
    The commands are written in run.sh inside the folder (see below) directory.
    After running bash run.sh from the terminal, the gti file with name niObsID.gti will be created inside each of the ObsID directory.
    Inside folder/ObsID/, there should be text file with name obsid_time.txt, 
    the obsid_time.txt is expected to have a list of time intervals which are to be included (include=True) or excluded (include=False) while running nicerl2.

    In addition, this function also writes down the expression used in nimaketime inside the ObsID directory under the name of ObsID_nimaketime_expr.txt.
   
    :param folder: folder is the directory containing directories with names given in filename. 
    :param filename: ('ObsID.dat') name of the file containing the list of ObsID
    :param include: whether to include the gtis, or exclude them
    :param run_f: ('run.sh') name of the file containing all of the run commands (created in `folder`)

    '''

    os.chdir(folder)
    obsid = np.loadtxt(filename, dtype=str)
    run = open(run_f, 'w')

    for obs in obsid:
        infile  =obs+'/'+obs+'_time.txt'
        outfile =obs+'/'+obs+'_nimaketime_expr.txt'

        expr    =nimaketime_expr(infile,outfile,include)
        run.write('nimaketime infile=\'{}/{}/auxil/ni{}.mkf\' outfile=\'{}/{}/ni{}.gti\' expr=\'{}\' clobber=YES\n'.format(folder,obs,obs,folder,obs,obs,expr))
        
                
    run.close()


