#The arguments of this function
#First argument: Example file of the scan
#Second argument: save directory for the cxi file
#Third argument: positions dat-file
#Third argument (Optional): Flatfield

import numpy as np
import os
import sys
import h5py as h5
import scipy.constants as sc
#Additional parameters
dbeam_cut_range=10
roi_size=10


#input arguments
print("Arguments should be 1: scan number, 2: target (Mo/Cu/Rh), 3: detector distance in meters")
scanno=sys.argv[1]
targetnam=str(sys.argv[2])
detector_dist=float(sys.argv[3])
print("All arguments are set.")

#Process input arguments
target_dict={
"Mo": 17.48,
"Cu": 8,
"Rh": 20.2,
}
energy=target_dict[targetnam]*1000*sc.e
scan_name="Scan_{}".format(scanno)

#Paths
logdir = "/gpfs/cfel/cxi/labs/MLL-Sigray/scan-logs/"
scanssdir = "/gpfs/cfel/cxi/labs/MLL-Sigray/Scans/"
mask="/gpfs/cfel/cxi/labs/MLL-Sigray/mask/lambda_mask1.h5"
maskpath="/data"


#Load mask
maskfile=h5.File(mask,"r")
mask=maskfile[maskpath][()].astype(np.int)

#Load log and read tings
lognam="{0}{1}.log".format(logdir,scan_name)
logfile=open(,"r")
for line in log:
    if line.beginswith("# Points count"):
        N_points=int(line.split(":")[1])
    if line.beginswith():


def bld_basis(whole_data):
    basis_vectors=np.zeros((whole_data.shape[0],2,3))
    pixel_vector=np.array((cf.pxsize_x,cf.pxsize_y,cf.pxsize_z))
    for i in range(0,whole_data.shape[0],1):
        basis_vectors[i,0,:]=np.multiply(cf.unitvector_f,pixel_vector)
        basis_vectors[i,1,:]=np.multiply(cf.unitvector_s,pixel_vector)
    return(basis_vectors)


#Loading the data
ldir=scanssdir+scan_name
allfiles=os.listdir(ldir)
files=[]
for file in allfiles:
    if file.endswith(".nxs") and file.startswith(scan_name):
        files.append(file)
files.sort()
print("Detected {} files.".format(len(files)))
if N_points!=len(files):
    print("WARNING! Number of expected files does not match number of detected files.")
for i in range(0,len(files)):
    if i==0:
        example_f=h5.File(ldir+"/"+files[i],"r")
        if orientation=="v":
            example_data=example_f["/entry/instrument/detector/data"][0,:,:]
        else:
            example_data = example_f["/entry/instrument/detector/data"][0, :, :]
        data_full=np.zeros((len(files),example_data.shape[0],example_data.shape[1]))
    f_now=h5.File(ldir+"/"+files[i],"r")
    if orientation=="h":
        data_now=np.sum(f_now["/entry/instrument/detector/data"][0,cf.roi[0]:cf.roi[1],:],axis=0)
        data_insert = np.zeros_like(example_data)
        for i1 in range(0,data_insert.shape[0],1):
            data_insert[i1,:]=data_now
    else:
        data_now = np.sum(f_now["/entry/instrument/detector/data"][0, :, cf.roi[2]:cf.roi[3]],axis=1)
        data_insert = np.zeros_like(example_data)
        for i1 in range(0, data_insert.shape[1], 1):
            data_insert[:, i1] = data_now
    data_full[i,:,:]=data_insert
#now writing the cxifile
cxifile=h5.File(i_sdir+"/"+scan_prefix+".cxi","w")
#now creating the subbranches:
#entry_1
entry_1 = cxifile.create_group("entry_1")
    #data_1
print("Saving data...")
data_1 = entry_1.create_group("data_1")
        #data
data=data_1.create_dataset("data",data=data_full)
    #instrument_1
print("Saving metadata...")
instrument_1=entry_1.create_group("instrument_1")
        #detector_1
detector_1=instrument_1.create_group("detector_1")
                #basis vectors
basis_vectors_cxi=detector_1.create_dataset("basis_vectors",data=bld_basis(whole_data=data_full))
detector_dist=float(cf.detector_distance)
detector_distance_cxi=detector_1.create_dataset("distance",data=detector_dist)
                #mask
print("Making mask...")
mask=np.ones((data_full.shape[1],data_full.shape[2]))
energy=cf.energy*1000*sc.e
print("Energy is %s J"%energy)
print("Making wavelength")
wavelength=(sc.h*sc.c)/energy
print("Wavelength is %s m" %wavelength)
############################################
#Now import the positions
translation_file=np.loadtxt(fpos,skiprows=24)
translation_vec=np.zeros((translation_file.shape[0],3))
if orientation=="v":
	translation_vec[:,0]=translation_file[:,1]
else:
	translation_vec[:,1]=translation_file[:,1]

N_images=data_full.shape[0]
############################################
#powder
print("Makint powder array...")
if ff_tck==True:
    f_powder=h5.File(ff_fnam,"r")
    powderarr=f_powder["/entry/instrument/detector/data"][0,:,:]
else:
    powderarr=np.ones((data_full.shape[1],data_full.shape[2]))
    powderarr=np.sum(data_full,axis=0)
print("Making good frames array...")
good_frames_arr=np.arange(0,N_images,1)
##############################################
#looking for the files and the numbers
##############################################
			#x_pixel_size
x_pixel_cxi=detector_1.create_dataset("x_pixel_size",data=cf.pxsize_x)
y_pixel_cxi=detector_1.create_dataset("y_pixel_size",data=cf.pxsize_y)
            #source_1
source_1=instrument_1.create_group("source_1")
                #wavelength
wavelength_cxi=source_1.create_dataset("wavelength",data=wavelength)
                #energy
energy_cxi=source_1.create_dataset("energy",data=energy)
    #sample_3
sample_2=entry_1.create_group("sample_2")
sample_3=entry_1.create_group("sample_3")
        #geometry
geometry=sample_3.create_group("geometry")
            #translation
translation=geometry.create_dataset("translation",data=translation_vec*10**-9)
#process_2
process_2=cxifile.create_group("process_2")
#powder
powder=process_2.create_dataset("powder",data=powderarr)
#process_3
process_3=cxifile.create_group("process_3")
#good frames
good_frames=process_3.create_dataset("good_frames",data=good_frames_arr)
#mask
mask=mask.astype("bool")
mask_cxi_3=process_3.create_dataset("mask",data=mask)
print("Done.")