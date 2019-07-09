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
import ini_maker
#Additional parameters
dbeam_cut_range=10
roi_size=20
year="2019"
manual_db=None

#input arguments
scanno=int(input("Please enter the scan number (e.g. 10):"))
targetnam=str(input("Please enter the target (Mo, Cu or Rh):"))
detector_dist=float(input("Please enter the detector distance in meters(e.g. 0.3):"))
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
#logdir="Log/"
scanssdir = "/gpfs/cfel/cxi/labs/MLL-Sigray/scan-frames/"
#scanssdir="Data/"
#savedir="output/"
savedir= "/gpfs/cfel/cxi/labs/MLL-Sigray/Processed/{0}/".format(year)
mask="/gpfs/cfel/cxi/labs/MLL-Sigray/mask/lambda_mask1.h5"
#mask="lambda_mask1.h5"

maskpath="/data"

#Old permanent config things
pxsize_x=55*10**-6
pxsize_y=55*10**-6
pxsize_z=0
unit_scanmotor=10**-3
unitvector_f=np.array((-1,0,0))
unitvector_s=np.array((0,-1,0))
db_coord=manual_db

unit_dict={
" m\n":1,
"mm\n":10**-3,
"µm\n":10**-6,
"nm\n":10**-9,
"pm\n":10**-12,
}


#Load mask
maskfile=h5.File(mask,"r")
mask=maskfile[maskpath][()].astype(np.int)

#Load log and read tings
lognam="{0}{1}.log".format(logdir,scan_name)
logfile=open(lognam,"r")
scanmotor=None
i=0
for line in logfile:
    if line.startswith("# Points count"):
        N_points=int(line.split(":")[1])
        positions = np.zeros((N_points)).astype("float")
    if line.startswith("# Device:") and line.endswith("Scanner")==False and line.endswith("Lambda\n")==False:
        scanmotor=str(line.split(":")[1][1:])
    if line.startswith("#")==False:
        current_entry=line.split(";")[2]
        if i==0:
            unit=current_entry[-3:]
            print(unit)
        unit_scaling=float(unit_dict[unit])
        current_pos=float(current_entry[:-3])*unit_scaling
        positions[i]=current_pos
        i+=1
if scanmotor==None:
    print("WARNING: No proper scanmotor detected!")
else:
    print("Scanmotor is {}".format(scanmotor))

if scanmotor=="X-SAM\n":
    orientation="h"
elif scanmotor=="Y-SAM\n":
    orientation="v"
else:
    print("Scanmotor is not X-SAM or Y-SAM. Don't know orientation. Assuming h.")
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
frames_arr=np.arange(0,N_points,1)
for i in range(0,N_points):
    try:
        if i==startframe:
            file_now="/{0}_".format(scan_name)+"{:05.0f}_Lambda.nxs".format(i)
            example_f=h5.File(ldir+"/"+file_now,"r")
            example_data=example_f["/entry/instrument/detector/data"][0,:,:]
            print(np.sum(example_data))
            if np.sum(example_data)<10:
                startframe+=1
                print("First frame was only zeros")
                useframes[i]=False
            else:
                if db_coord==None:
                    db_coord=np.unravel_index(np.argmax(np.multiply(mask,example_data)),example_data.shape)
                print("Direct beam pixel coordinate is {0}.".format(db_coord))
                roi=(db_coord[0]-int(roi_size/2),db_coord[0]+int(roi_size/2),db_coord[1]-int(roi_size/2),db_coord[1]+int(roi_size/2))
                print("Startframe is {0}. Creating full data array.".format(i))
                data_full=np.zeros((N_points,example_data.shape[0],example_data.shape[1]))
                raw_frames=np.zeros_like(data_full)
        else:
            file_now = "/{0}_".format(scan_name) + "{:05.0f}_Lambda.nxs".format(i)
            f_now=h5.File(ldir+"/"+file_now,"r")
            if orientation=="h":
                data_now=np.sum(f_now["/entry/instrument/detector/data"][0,roi[0]:roi[1],:],axis=0)
                frame_now=f_now["/entry/instrument/detector/data"][0,:,:]
                raw_frames[i,:,:]=frame_now
                data_now[(db_coord[1]-int(dbeam_cut_range/2)):(db_coord[1]+int(dbeam_cut_range/2))]=np.zeros_like(data_now[(db_coord[1]-int(dbeam_cut_range/2)):(db_coord[1]+int(dbeam_cut_range/2))])
                data_insert = np.zeros_like(example_data)
                for i1 in range(0,data_insert.shape[0],1):
                    data_insert[i1,:]=data_now
            else:
                data_now = np.sum(f_now["/entry/instrument/detector/data"][0, :, roi[2]:roi[3]],axis=1)
                frame_now=f_now["/entry/instrument/detector/data"][0,:,:]
                raw_frames[i,:,:]=frame_now
                data_now[(db_coord[0] - int(dbeam_cut_range / 2)):(db_coord[0] + int(dbeam_cut_range / 2))] = np.zeros_like(data_now[(db_coord[0] - int(dbeam_cut_range / 2)):(db_coord[0] + int(dbeam_cut_range / 2))])
                data_insert = np.zeros_like(example_data)
                for i1 in range(0, data_insert.shape[1], 1):
                    data_insert[:, i1] = data_now
            data_full[i,:,:]=data_insert
    except (KeyError,OSError):
        if i==startframe:
            startframe+=1
        print("Didnt find file {} or proper data in file.".format(file_now))
        useframes[i]=False

#writing the cxifile
if os.path.isdir(savedir+scan_name)==False:
    os.mkdir(savedir+scan_name)
cxifile=h5.File("{0}{1}/{1}.cxi".format(savedir,scan_name),"w")
entry_1 = cxifile.create_group("entry_1")
print("Saving data...")
data_1 = entry_1.create_group("data_1")
data=data_1.create_dataset("data",data=data_full[useframes])
print("Saving raw frames...")
data2=data_1.create_dataset("raw_frames",data=raw_frames)

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
good_frames_arr=np.arange(0,translation_vec.shape[0],1)
x_pixel_cxi=detector_1.create_dataset("x_pixel_size",data=pxsize_x)
y_pixel_cxi=detector_1.create_dataset("y_pixel_size",data=pxsize_y)
source_1=instrument_1.create_group("source_1")
wavelength_cxi=source_1.create_dataset("wavelength",data=wavelength)
db_coordinate=source_1.create_dataset("direct_beam_coordinate",data=db_coord)
energy_cxi=source_1.create_dataset("energy",data=energy)
sample_2=entry_1.create_group("sample_2")
sample_3=entry_1.create_group("sample_3")
geometry=sample_3.create_group("geometry")
translation=geometry.create_dataset("translation",data=translation_vec)
make_whitefield=cxifile.create_group("make_whitefield")
whitefield=make_whitefield.create_dataset("whitefield",data=powderarr)
frame_selector=cxifile.create_group("frame_selector")
good_frames=frame_selector.create_dataset("good_frames",data=good_frames_arr)
mask=mask.astype("bool")
mask_maker=cxifile.create_group("mask_maker")
mask_cxi_3=mask_maker.create_dataset("mask",data=mask)
print("Done.")
cxifile.close()
print("Making ini files")
ini_maker.mk_make_whitefield_ini(path="{0}{1}".format(savedir,scan_name))
ini_maker.mk_make_speckle_gui_ini(path="{0}{1}".format(savedir,scan_name))
ini_maker.mk_stitch_ini(path="{0}{1}".format(savedir,scan_name),roi=roi)
ini_maker.mk_update_pixel_map(path="{0}{1}".format(savedir,scan_name),roi=roi)
ini_maker.mk_zernike_ini(path="{0}{1}".format(savedir,scan_name),roi=roi)
print("Done.")