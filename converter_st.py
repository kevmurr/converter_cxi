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
year="2019"


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
savedir=""
#savedir= "/gpfs/cfel/cxi/labs/MLL-Sigray/Processed/{0}/".format(year)
mask="/gpfs/cfel/cxi/labs/MLL-Sigray/mask/lambda_mask1.h5"
maskpath="/data"

#Old permanent config things
pxsize_x=55*10**-6
pxsize_y=55*10**-6
pxsize_z=0
unit_scanmotor=10**-3
unitvector_f=np.array((-1,0,0))
unitvector_s=np.array((0,-1,0))




#Load mask
maskfile=h5.File(mask,"r")
mask=maskfile[maskpath][()].astype(np.int)

#Load log and read tings
lognam="{0}{1}.log".format(logdir,scan_name)
logfile=open(lognam,"r")
scanmotor=None
for line in logfile:
    if line.beginswith("# Points count"):
        N_points=int(line.split(":")[1])
    if line.beginswith("# Device:") and line.endswith("Scanner")==False and line.endswith("Lambda")==False:
        scanmotor=str(line.split(":")[1][1:])
if scanmotor==None:
    print("WARNING: No proper scanmotor detected!")
else:
    print("Scanmotor is {}".format(scanmotor))

if scanmotor=="SAMX":
    orientation="h"
elif scanmotor=="SAMY":
    orientation="v"
else:
    print("Scanmotor is not SAMX or SAMY. Don't know orientation. Assuming h.")
    orientation="h"




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
useframes=np.ones((N_points)).astype(np.bool)
startframe=0
for i in range(0,N_points):
    try:
        if i==startframe:
            example_f=h5.File(ldir+"/"+files[i],"r")
            example_data=example_f["/entry/instrument/detector/data"][0,:,:]
            if np.sum(example_data)==0:
                startframe+=1
                print("First frame was only zeros")
                useframes[i]=False
            else:
                db_coord=np.unravel_index(np.argmax(np.multiply(mask,example_data)),example_data.shape)
                print("Direct beam pixel coordinate is {0}.".format(db_coord))
                roi=(db_coord[0]-int(roi_size/2),db_coord[0]+int(roi_size/2),db_coord[1]-int(roi_size/2),db_coord[1]+int(roi_size/2))
                print("Startframe is {0}. Creating full data array.".format(i))
                data_full=np.zeros((N_points,example_data.shape[0],example_data.shape[1]))
        else:
            f_now=h5.File(ldir+"/"+files[i],"r")
            if orientation=="h":
                data_now=np.sum(f_now["/entry/instrument/detector/data"][0,roi[0]:roi[1],:],axis=0)
                data_now[(db_coord[1]-int(dbeam_cut_range/2)):(db_coord[1]+int(dbeam_cut_range/2))]=np.zeros_like(data_now[(db_coord[1]-int(dbeam_cut_range/2)):(db_coord[1]+int(dbeam_cut_range/2))])
                data_insert = np.zeros_like(example_data)
                for i1 in range(0,data_insert.shape[0],1):
                    data_insert[i1,:]=data_now
            else:
                data_now = np.sum(f_now["/entry/instrument/detector/data"][0, :, roi[2]:roi[3]],axis=1)
                data_now[(db_coord[0] - int(dbeam_cut_range / 2)):(db_coord[0] + int(dbeam_cut_range / 2))] = np.zeros_like(data_now[(db_coord[0] - int(dbeam_cut_range / 2)):(db_coord[0] + int(dbeam_cut_range / 2))])
                data_insert = np.zeros_like(example_data)
                for i1 in range(0, data_insert.shape[1], 1):
                    data_insert[:, i1] = data_now
            data_full[i,:,:]=data_insert
    except (KeyError,OSError):
        if i==startframe:
            startframe+=1
        print("Didnt find file {0} or proper data in file.".format(ldir+"/"+files[i]))
        useframes[i]=False

#writing the cxifile
if os.path.isdir(savedir+scan_name)==False:
    os.mkdir(savedir+scan_name)
cxifile=h5.File("{0}{1}/{1}.cxi".format(savedir,scan_name),"w")
entry_1 = cxifile.create_group("entry_1")
print("Saving data...")
data_1 = entry_1.create_group("data_1")
data=data_1.create_dataset("data",data=data_full[useframes])
instrument_1=entry_1.create_group("instrument_1")
detector_1=instrument_1.create_group("detector_1")

def bld_basis(whole_data,useframes=useframes):
    basis_vectors=np.zeros((whole_data.shape[0],2,3))
    pixel_vector=np.array((pxsize_x,pxsize_y,pxsize_z))
    for i in range(0,whole_data.shape[0],1):
        basis_vectors[i,0,:]=np.multiply(unitvector_f,pixel_vector)
        basis_vectors[i,1,:]=np.multiply(unitvector_s,pixel_vector)
    return(basis_vectors[useframes,:,:])

basis_vectors_cxi=detector_1.create_dataset("basis_vectors",data=bld_basis(whole_data=data_full,useframes=useframes))
detector_distance_cxi=detector_1.create_dataset("distance",data=float(detector_dist))
print("Energy is %s J"%energy)
print("Making wavelength")
wavelength=(sc.h*sc.c)/energy
print("Wavelength is %s m" %wavelength)

unit_dict={
" m":1,
"mm":10**-3,
"µm":10**-6,
"nm":10**-9,
"pm":10**-12,
}

#Now import the positions
translation_file=logfile
i=0
positions=np.zeros((N_points)).astype("float")
for line in translation_file and line.beginswith("#")==False:
    current_entry=line.split(";")[2]
    if i==0:
        unit=current_entry[-2:]
    unit_scaling=float(unit_dict[unit])
    current_pos=float(current_entry[:-2])*unit_scaling
    positions[i]=current_pos
    i+=1

translation_vec=np.zeros((N_points,3))
if orientation=="v":
	translation_vec[:,0]=positions
else:
	translation_vec[:,1]=positions
translation_vec=translation_vec[useframes,:]


print("Makint powder array...")
powderarr=np.ones((data_full.shape[1],data_full.shape[2]))
powderarr=np.sum(data_full[useframes],axis=0)
print("Making good frames array...")
good_frames_arr=np.arange(0,N_points,1)
x_pixel_cxi=detector_1.create_dataset("x_pixel_size",data=pxsize_x)
y_pixel_cxi=detector_1.create_dataset("y_pixel_size",data=pxsize_y)
source_1=instrument_1.create_group("source_1")
wavelength_cxi=source_1.create_dataset("wavelength",data=wavelength)
energy_cxi=source_1.create_dataset("energy",data=energy)
sample_2=entry_1.create_group("sample_2")
sample_3=entry_1.create_group("sample_3")
geometry=sample_3.create_group("geometry")
translation=geometry.create_dataset("translation",data=translation_vec)
make_whitefield=cxifile.create_group("make_whitefield")
whitefield=make_whitefield.create_dataset("whitefield",data=powderarr)
frame_selector=cxifile.create_group("frame_selector")
good_frames=frame_selector.create_dataset("good_frames",data=good_frames_arr[useframes])
mask=mask.astype("bool")
mask_maker=cxifile.create_group("mask_maker")
mask_cxi_3=mask_maker.create_dataset("mask",data=mask)
print("Done.")
cxifile.close()