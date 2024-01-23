#!/bin/sh

#Download from https://www.hdfgroup.org/downloads/hdf5/
#todo Check again for a new download
#actualçy this is wrong firs configure and then say make install
echo "Installing and Moving HDF5 files"


#todo if it already exists do not configure or make install but direcly access the file !
#
# Default for this version 1.12.0
#if [ -d /HDF5Files/HDF5 ]; then cd /HDF5Files/HDF5/hdf5-1.12.0 && \
#this part is for installation but change it
	#./configure  && make install &&\
#	cp -r /HDF5Files/HDF5/hdf5-1.12.0/hdf5 /usr/local/; fi


# Default for this version 1.12.0
if [ -d /HDF5Files/HDF5/hdf5-1.12.0/hdf5 ]; then cd /HDF5Files/HDF5/hdf5-1.12.0 && \
	cp -r /HDF5Files/HDF5/hdf5-1.12.0/hdf5 /usr/local/ ;

elif [ -d /HDF5Files/HDF5/hdf5-1.12.0 ] ; then cd /HDF5Files/HDF5/hdf5-1.12.0 && \
    ./configure  && make install && \
    cp -r /HDF5Files/HDF5/hdf5-1.12.0/hdf5 /usr/local/;

else
    echo "Please download hdf5 library"

fi;


echo "Installed HDF5"