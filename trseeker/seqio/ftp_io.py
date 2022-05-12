#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 05.06.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
"""
Classes:
    
- AbstractFtpIO(object)
    
"""
import ftplib
import os
import gzip

class AbstractFtpIO(object):
    """ Class for working with data via FTP.
        
    Public methods:
    
    - __init__(self, ftp_address=None)
    - connect(self)
    - cd(self, list_of_dirs)
    - ls(self)
    - get(self, file, output_file)
    - unzip(self, file_name)
        
    """

    def __init__(self, ftp_address=None):
        """ Set ftp parameters.
        
        Keyword arguments:
        
        - ftp_address   -- address of ftp server
        """
        self.ftp_address = ftp_address
        self.ftp = None

    def connect(self):
        """ Connect anonymously to server."""
        print "Connect to", self.ftp_address
        self.ftp = ftplib.FTP(self.ftp_address)
        self.ftp.login()

    def cd(self, list_of_dirs):
        """ Go to dir by given path (list of dirs)."""
        for dirname in list_of_dirs:
            self.ftp.cwd(dirname)

    def ls(self):
        """ Return list of files filtered with filter_func in folder.
        """
        return self.ftp.nlst()

    def get(self, file, output_file):
        """ Save file in output file."""
        with open(output_file, "ab") as fh:
            print("Downloading: %s" % file)
            self.ftp.retrbinary("RETR %s" % file, lambda x: fh.write(x))

    def get_aspera_ncbi(self, source, destination):
        """ Download with aspera."""
        download_with_aspera_from_ncbi(source, destination)

    def size(self, file_name):
        """ Return size of file.
        """
        self.ftp.sendcmd("TYPE i")
        return self.ftp.size(file_name)

    def unzip(self, file_name):
        ''' Unzip file.'''
        if not file_name.endswith(".gz"):
            return
        fh = gzip.GzipFile(fileobj=open(file_name, 'rb'))
        fw = file(file_name.replace('.gz', ''), 'wb')
        for line in fh:
            fw.write(line)
        fh.close()
        fw.close()

def download_with_aspera_from_ncbi(source, destination):
    '''
    anonftp@ftp-private.ncbi.nlm.nih.gov:%(source)s
    where source for example /1Gb
    '''
    data = {
        'source': source,
        'dest': destination,

    }
    command = "~/.aspera/connect/bin/ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -k1 -Tr -l100m anonftp@ftp-private.ncbi.nlm.nih.gov:%(source)s %(dest)s" % data
    print command
    os.system(command)



