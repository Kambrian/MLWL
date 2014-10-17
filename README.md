Maximum Likelihood code accompanying this paper:
Galaxy and Mass Assembly (GAMA): The halo mass of galaxy groups from maximum-likelihood weak lensing
Jiaxin Han et al., 2014
http://arxiv.org/abs/1404.6828

Prerequisites:
GSL c library
ctypesgsl python library
hdf5 c library

Get started:

Prepare the datafile in hdf5 format, like Null.hdf5. Only the x and v are needed, other datasets are optional. x and v need to be in kpc and km/s, with center position and velocity already subtracted.

in Makefile, modify the rootdir to be your data directory, i.e., change this line
    ROOTDIR='"/cosma/home/jvbq85/data/DynDistr"'
to the directory, so that your datafile is found at $ROOTDIR/data/your_data_file.hdf5.

Then
make lib

Then fit with
fit.py
under py/ folder.
