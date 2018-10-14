#!/usr/bin/env python

'''
PyASB FTP module.

Get images iteratively from FTP server.
____________________________

This module is part of the PyASB project, 
created and maintained by Mireia Nievas [UCM].
____________________________
'''


__author__ = "Mireia Nievas"
__copyright__ = "Copyright 2012, PyASB project"
__credits__ = ["Mireia Nievas"]
__license__ = "GNU GPL v3"
__shortname__ = "PyASB"
__longname__ = "Python All-Sky Brightness pipeline"
__version__ = "1.99.0"
__maintainer__ = "Mireia Nievas"
__email__ = "mirph4k[at]gmail[dot]com"
__status__ = "Prototype" # "Prototype", "Development", or "Production"


try:
    import sys,os,inspect
    import datetime
    import time
    import signal
    import urllib
    import gzip
    import ftputil
    import pyfits
except:
    print(str(inspect.stack()[0][2:4][::-1])+': One or more modules missing')
    raise SystemExit


''' User defined options '''

base_dir = "/usr/users/mnievas/PyASB"
temporary_path = "/tmp/"
skymap_path = base_dir+"/skymaps/"
photometry_table_path = base_dir+"/starphotometry/"
bouguerfit_path = base_dir+"/bouguerfits/"
skybrightness_map_path = base_dir+"/skybrightnessmaps/"
skybrightness_table_path = base_dir+"/skybrightnesstable/"
cloudmap_path = base_dir+"/cloudmaps/"
clouddata_table_path = base_dir+"/clouddata/"
summary_path = base_dir+"/summary/"
register_analyzed_files = base_dir+"/register.txt"

ftp_images_since_reconnect = 0
ftp_max_images_until_reconnect = 20


'''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~ Halt handler ~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''

def handler(signum, frame):
        print 'Signal handler called with signal', signum
        print "CTRL-C pressed"
        sys.exit(0)

signal.signal(signal.SIGTERM, handler)
signal.signal(signal.SIGINT, handler)


class FtpSession():
    def __init__(self,ftp_server,ftp_user,ftp_pass,ftp_basedir,analysis_basedir):
        self.ftp_server  = ftp_server
        self.ftp_user    = ftp_user
        self.ftp_pass    = ftp_pass
        self.ftp_basedir = ftp_basedir
        
        self.ftp_connect()
        self.fitslist = []
        self.fits_search_remotepath(basedir=ftp_basedir)
        self.ftp_disconnect()
        
        self.update_dirs_path(analysis_basedir)
        self.create_analysis_paths()
        #self.perform_analysis()
    
    def ftp_connect(self):
        ''' Establish FTP connection. Returns ftp session '''
        self.ftp = ftputil.FTPHost(ftp_server,ftp_user,ftp_pass)
    
    def ftp_disconnect(self):
        ''' End FTP connection '''
        self.ftp.close()
    
    def test_is_directory(self,path):
        try:
            self.ftp.chdir(path)
            self.ftp.chdir('..')
            return(True)
        except:
            return(False)
    
    def list_files(self):
        ''' List files in FTP remote dir '''
        try: listfiles = self.ftp.listdir()
        except: return []
        else: return listfiles
    
    def update_dirs_path(self,analysis_basedir):
        try:
            global base_dir,temporary_path,skymap_path,photometry_table_path,\
              bouguerfit_path,skybrightness_map_path,skybrightness_table_path,\
              cloudmap_path,clouddata_table_path,summary_path,register_analyzed_files
            
            base_dir = analysis_basedir
            temporary_path = "/tmp/"
            skymap_path = base_dir+"/skymaps/"
            photometry_table_path = base_dir+"/starphotometry/"
            bouguerfit_path = base_dir+"/bouguerfits/"
            skybrightness_map_path = base_dir+"/skybrightnessmaps/"
            skybrightness_table_path = base_dir+"/skybrightnesstable/"
            cloudmap_path = base_dir+"/cloudmaps/"
            clouddata_table_path = base_dir+"/clouddata/"
            summary_path = base_dir+"/summary/"
            register_analyzed_files = base_dir+"/register.txt"
        except:
            print(str(inspect.stack()[0][2:4][::-1])+' Dont update base_dir')

    def download_file(self,remote_filename,local_filename):
        ''' Download remote file to local temporary copy '''
        if not os.path.exists(temporary_path): os.makedirs(temporary_path)
        #local_filename = temporary_path + tempfilename
        #temporary_file = open(tempfilename,'wb')
        if os.path.isfile(local_filename):
            os.remove(local_filename)
        
        # Binary files are easily corrupted if downloaded from Win->Linux through ftp.
        # Download instead using urllib.    
        urllib.urlretrieve(\
            'ftp://'+str(self.ftp_user)+':'+str(self.ftp_pass)+\
            '@'+str(self.ftp_server)+str(remote_filename),str(local_filename))
        
        #self.ftp.download(remote_filename, local_filename,'w')
        #self.ftp.retrbinary('RETR '+filename,temporary_file.write)
        #temporary_file.close()
    
    def create_analysis_paths(self):
        ''' Create neccesary file path. Must be defined at the beginning of the script'''
        for directory in [\
         base_dir, skymap_path, photometry_table_path, \
         bouguerfit_path, skybrightness_map_path, skybrightness_table_path, \
         summary_path, temporary_path, cloudmap_path, clouddata_table_path]:
            if not os.path.exists(directory): os.makedirs(directory)
    
    def fits_search_remotepath(self,basedir=''):
        ''' Recursively search fits files in the remote server. Return a self.list of fits files'''
        recursive = self.ftp.walk(basedir,topdown=True,onerror=None)
        for root,dirs,files in recursive:
            for name in files:
                thefile = os.path.join(root, name)
                #FileName, FileExt = os.path.splitext(thefile)
                valid_ext = ['.fts','.fit','.fits','.FTS','.FIT','.FITS']
                invalid_prefixes = ['SBJo']
                
                use_file = False
                for known_ext in valid_ext:
                    if known_ext in thefile:
                        use_file = True
                
                if use_file == True:
                    for bad_prefix in invalid_prefixes:
                        #if FileExt == known_ext:
                        if bad_prefix in thefile:
                            use_file = False
                
                if use_file == True:
                    self.fitslist.append("/"+thefile)
        
    
    def perform_analysis(self,pyasb_fullanalysis,pyasb_overwrite):
        ''' 
        For each file in the remote list:
         1.- Download the image to a tmp file.
         2.- Perform the analysis.
         3.- write results and append the file to the list of analyzed files.
        '''
        
        ftp_images_since_reconnect = 0;
        
        for each_fitsfile in self.fitslist:
            ftp_images_since_reconnect+=1;
            if "Johnson_U" in each_fitsfile:
                print('Johnson U file detected, SNR will be too low, discarding')
                continue
            
            if not os.path.isfile(register_analyzed_files):
                register_analyzed = open(register_analyzed_files,"w+")
                register_analyzed.close()
            
            register_analyzed = open(register_analyzed_files,"r+").readlines()
            
            if os.path.split(each_fitsfile)[1] in str(register_analyzed):
                if pyasb_overwrite==True: 
                    print 'Previous analysis results detected, Overwrite mode is ON'
                elif pyasb_overwrite==False: 
                    print 'Previous analysis results detected, Overwrite mode is OFF'
                    continue
            
            self.ftp_connect()
            print('---> Downloading file '+str(each_fitsfile))
            actual_localdir = os.getcwd()
            # Download the image
            temp_filename = temporary_path+'/pyasb_'+\
             str(datetime.datetime.now())\
             .replace("-","").replace(":","").\
             replace(" ","_").replace(".","_")+\
             '.fits'
            
            try:
                self.download_file(each_fitsfile,temp_filename)
            except:
                print(inspect.stack()[0][2:4][::-1])
                raise
                print 'File cannot be downloaded, continue with next one';
                self.ftp_disconnect()
            
            print('---> Launching PyASB to analyze '+str(each_fitsfile))
            # Analize the image
            complete_analysis_exec = " "
            if pyasb_fullanalysis==True:
                complete_analysis_exec = \
                 " -os " + skybrightness_map_path +\
                 " -om " + skymap_path +\
                 " -ost " + skybrightness_table_path +\
                 " -ocm " + cloudmap_path +\
                 " -oct " + clouddata_table_path
            try:
                os.system("python pyasb_launcher.py -i "+\
                 temp_filename+\
                 " -ob " + bouguerfit_path + \
                 " -ot " + photometry_table_path +\
                 " -or " + summary_path + " -c "+ pyasb_config + \
                 complete_analysis_exec)
                os.remove(temp_filename)
            except:
                print('Error performing pyasb analysis, please check the file.')
                #pass
            
            register_analyzed = open(register_analyzed_files,"a+")
            register_analyzed.write(each_fitsfile+"\r\n");
            self.ftp_disconnect()
            time.sleep(2)
            

def show_help():
    ''' Print Help message on std output '''
    print(\
     'PyASB ftp module. This program tries to recursively analize '+\
     'located on a remote FTP repository'+\
     'Parameters:'+\
     ' -h : print this help message and exit\n'+\
     ' -s server_ip: the FTP server IP\n'+\
     ' -u ftp_username: the FTP username\n'+\
     ' -p ftp_userpass: the FTP user passwd\n'+\
     ' -b base_dir: base dir to start the iterative search\n'+\
     ' -c config_file: pyasb config file to use\n'+\
     ' -d analysis_basedir: pyasb analysis base directory\n'+\
     ' --overwrite: overwrite any previous analysis data in dir.\n'+\
     '              Default is to keep the old data\n'+\
     ' --full: perform full analysis (generates sky brightness map, takes more time).\n'+\
     '         Default is to do a simple analysis (no SB map generation)\n')
    exit(0)
    

if __name__ == '__main__':
    ''' Read input options and start the analysis 
        Results will be stored in the same directory'''
    try: 
        input_options = sys.argv[1:]
        print('Input Options: '+str(input_options))
    except: 
        print('ERR. No INPUT parameters detected\n')
        show_help()
    else:
        pyasb_fullanalysis = False
        pyasb_overwrite = False
        ftp_user = 'anonymous'
        ftp_pass = ''
        ftp_basedir = ''
        pyasb_config = 'pyasb_config.cfg'
        
        while len(input_options)>0:
            input_option = input_options[0]
            if input_option == '-h':
                show_help()
            elif input_option == '-s':
                ftp_server = input_options[1]
                input_options.pop(1)
            elif input_option == '-u':
                ftp_user = input_options[1]
                input_options.pop(1)
            elif input_option == '-p':
                ftp_pass = input_options[1]
                input_options.pop(1)
            elif input_option == '-b':
                ftp_basedir = input_options[1]
                input_options.pop(1)
            elif input_option == '-c':
                pyasb_config = input_options[1]
                input_options.pop(1)
            elif input_option == '-d':
                                analysis_basedir = input_options[1]
                                input_options.pop(1)
            elif input_option == '--full':
                pyasb_fullanalysis = True
            elif input_option == '--overwrite':
                pyasb_overwrite = True
            input_options.pop(0)
    
    try: ftp_server
    except:
        print(inspect.stack()[0][2:4][::-1])
        print('ERR. No FTP server specified')
        show_help()
    
    try:
        FtpRemoteSession = FtpSession(\
          ftp_server,\
          ftp_user,\
          ftp_pass,\
          ftp_basedir,\
          analysis_basedir)
        
        print('Number of files found: '+str(len(FtpRemoteSession.fitslist)))
        FtpRemoteSession.perform_analysis(pyasb_fullanalysis,pyasb_overwrite)
    except:
        print(inspect.stack()[0][2:4][::-1])
        raise
