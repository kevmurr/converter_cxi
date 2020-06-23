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
from scipy.ndimage import median_filter
#Additional parameters
extrude=True #Set if data should be extruded according to old ST software
roi_size = 60
year = "2020"
manual_db = None # (378, 825) #None 370
manual_orientation = "h" #Is used in case if the scan motor differs from SAMX or SAMY

#input arguments
scanno = int(input("Please enter the scan number (e.g. 10) or e.g. -1 for the last scan:"))
targetnam = str(input("Please enter the target (Mo, Cu or Rh):"))
detector_dist = float(input("Please enter the detector distance in meters(e.g. 0.3):"))
print("All arguments are set.")

#Process input arguments
target_dict={
    "Mo": 17.48,
    "Cu": 8.,
    "Rh": 20.2}
energy = target_dict[targetnam] * 1e3 * sc.e


#Paths

logdir = "/gpfs/cfel/cxi/labs/MLL-Sigray/scan-logs/"
scanssdir = "/gpfs/cfel/cxi/labs/MLL-Sigray/scan-frames/"
savedir= "/gpfs/cfel/cxi/labs/MLL-Sigray/Processed/{0}/".format(year)
mask="/gpfs/cfel/cxi/labs/MLL-Sigray/mask/lambda_mask1.h5"

maskpath="/data"

#Old permanent config things
pxsize_x = 55*10**-6
pxsize_y = 55*10**-6
pxsize_z = 0
unit_scanmotor = 10**-3
unitvector_f = np.array([-1, 0, 0])
unitvector_s = np.array([0, -1, 0])
db_coord = manual_db

def bld_basis(whole_data,useframes):
    basis_vectors=np.zeros((whole_data.shape[0],2,3))
    pixel_vector=np.array((pxsize_x,pxsize_y,pxsize_z))
    for i in range(0,whole_data.shape[0],1):
        basis_vectors[i,0,:]=np.multiply(unitvector_f,pixel_vector)
        basis_vectors[i,1,:]=np.multiply(unitvector_s,pixel_vector)
    return(basis_vectors[useframes,:,:])

unit_dict = {
    "m":1.,
    "mm":10**-3,
    "µm":10**-6,
    "nm":10**-9,
    "pm":10**-12}

#Load mask
maskfile = h5.File(mask,"r")
mask = maskfile[maskpath][()].astype(np.int)

#Load log and read tings
if scanno<0:
    list_lognams=sorted(os.listdir(logdir))
    list_lognams_scans=[i for i in list_lognams if i.startswith('Scan_')]
    scanno_use=int(list_lognams_scans[scanno].split("_")[1].split(".")[0])
    print("Detected scan number:",scanno_use)
    scanno=scanno_use

scan_name = "Scan_{}".format(scanno)
lognam = "{0}{1}.log".format(logdir,scan_name)
logfile = open(lognam,"r")
scanmotor = None
i = 0
if scanno>1002:
    for line in logfile:
        if line.startswith("# Points count"):
            N_points=int(line.split(":")[1])
            positions = np.zeros((N_points)).astype("float")
        if line.startswith("# Device:") and not line.endswith("Scanner") and not line.endswith("Lambda\n"):
            scanmotor=str(line.split(":")[1][1:])
        if line.startswith("# Start point:"):
            unit=line.split(" ")[-1]
            try:
                unit_scaling = float(unit_dict[unit])
            except KeyError:
                if unit.startswith("n"):
                    unit_scaling = 1E-9
                else:
                    unit_scaling = 1.0
        if not line.startswith("#"):
            current_entry = line.split(";")[2]
            current_pos = float(current_entry.split(" ")[0])*unit_scaling
            positions[i] = current_pos
            i += 1
else:
    for line in logfile:
        if line.startswith("# Points count"):
            N_points=int(line.split(":")[1])
            positions = np.zeros((N_points)).astype("float")
        if line.startswith("# Device:") and not line.endswith("Scanner") and not line.endswith("Lambda\n"):
            scanmotor=str(line.split(":")[1][1:])
        if not line.startswith("#"):
            current_entry = line.split(";")[2]
            if i == 0:
                #unit = current_entry[-3:].strip()
                unit = current_entry.split(" ")[1]
                print("The unit is",unit)
            try:
                unit_scaling = float(unit_dict[unit])
            except KeyError:
                if unit.startswith("n"):
                    unit_scaling = 1E-9
                else:
                    unit_scaling = 1.0
            current_pos = float(current_entry.split(" ")[0])*unit_scaling
            positions[i] = current_pos
            i += 1
print("The scaling is ", unit_scaling)
if scanmotor==None:
    print("WARNING: No proper scanmotor detected!")
    orientation=manual_orientation
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
frames_arr=np.arange(0,N_points,1)
startframe=0
for i in range(0,N_points):
    try:
        if i==startframe:
            file_now="/{0}_".format(scan_name)+"{:05.0f}_Lambda.nxs".format(i)
            example_f=h5.File(ldir+"/"+file_now,"r")
            example_data=example_f["/entry/instrument/detector/data"][0,:,:]
            print(np.sum(example_data))
            if db_coord==None:
                db_coord=np.unravel_index(np.argmax(np.multiply(mask,median_filter(example_data,size=10))),example_data.shape)
            print("Direct beam pixel coordinate is {0}.".format(db_coord))
            roi=(db_coord[0]-int(roi_size/2),db_coord[0]+int(roi_size/2),db_coord[1]-int(roi_size/2),db_coord[1]+int(roi_size/2))
            print("Startframe is {0}. Creating full data array.".format(i))
            if orientation=="h":
                if extrude==True:
                    data_full = np.zeros((N_points, example_data.shape[0], example_data.shape[1]))
                else:
                    data_full=np.zeros((N_points,1,example_data.shape[1]))
            else:
                if extrude==True:
                    data_full = np.zeros((N_points, example_data.shape[0], example_data.shape[1]))
                else:
                    data_full = np.zeros((N_points, example_data.shape[0], 1))
            raw_frames=np.zeros((data_full.shape[0],example_data.shape[0],example_data.shape[1]))
        file_now = "/{0}_".format(scan_name) + "{:05.0f}_Lambda.nxs".format(i)
        f_now=h5.File(ldir+"/"+file_now,"r")
        if orientation=="h":
            data_now=np.sum(f_now["/entry/instrument/detector/data"][0,roi[0]:roi[1],:],axis=0)
            frame_now=f_now["/entry/instrument/detector/data"][0,:,:]
            raw_frames[i,:,:]=frame_now
            if extrude==True:
                data_insert=np.zeros_like(example_data)
                for irow in range(0,example_data.shape[0]):
                    data_insert[irow,:]=data_now
                data_full[i,:,:]=data_insert
            else:
                data_insert = data_now
                data_full[i, 0, :] = data_insert
        else:
            data_now = np.sum(f_now["/entry/instrument/detector/data"][0, :, roi[2]:roi[3]],axis=1)
            frame_now=f_now["/entry/instrument/detector/data"][0,:,:]
            raw_frames[i,:,:]=frame_now
            if extrude==True:
                data_insert = np.zeros_like(example_data)
                for icol in range(0,example_data.shape[1]):
                    data_insert[:,icol]=data_now
                data_full[i, :, :] = data_insert
            else:
                data_insert = data_now
                data_full[i, :, 0] = data_insert
    except (KeyError,OSError):
        if i==startframe:
            startframe+=1
        print("Didnt find file {} or proper data in file.".format(file_now))
        useframes[i]=False

#Make Ptychogram
if orientation=="h":
    ptychogram=np.zeros((data_full.shape[0],data_full.shape[2]))
    midy=int(data_full.shape[1]/2)
    for i in range(0,ptychogram.shape[0],1):
        ptychogram[i,:]=data_full[i,midy,:]
else:
    ptychogram = np.zeros((data_full.shape[0], data_full.shape[1]))
    midx=int(data_full.shape[2]/2)
    for i in range(0,ptychogram.shape[0],1):
        ptychogram[i,:]=data_full[i,:,midx]

#writing the cxifile
if os.path.isdir(savedir+scan_name)==False:
    os.mkdir(savedir+scan_name)
cxifile=h5.File("{0}{1}/{1}.cxi".format(savedir,scan_name),"w")
entry_1 = cxifile.create_group("entry_1")
print("Saving data...")
data_1 = entry_1.create_group("data_1")
data=data_1.create_dataset("data",data=data_full[useframes])
print("Data shape",data_full[useframes].shape)
print("Saving raw frames...")
data2=data_1.create_dataset("raw_frames",data=raw_frames[useframes])
print("Raw frames shape",data2.shape)
ptychogram_data=data_1.create_dataset("ptychogram",data=ptychogram)

instrument_1=entry_1.create_group("instrument_1")
detector_1=instrument_1.create_group("detector_1")



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
#powderarr=np.ones((data_full.shape[1],data_full.shape[2]))
powderarr=np.sum(data_full[useframes],axis=0)
print("whitefield shape", powderarr.shape)
print("Making good frames array...")
good_frames_arr=np.arange(0,translation_vec.shape[0],1)
x_pixel_cxi=detector_1.create_dataset("x_pixel_size",data=pxsize_x)
y_pixel_cxi=detector_1.create_dataset("y_pixel_size",data=pxsize_y)
source_1=instrument_1.create_group("source_1")
wavelength_cxi=source_1.create_dataset("wavelength",data=wavelength)
db_coordinate=source_1.create_dataset("direct_beam_coordinate",data=db_coord)
energy_cxi=source_1.create_dataset("energy",data=energy)
if extrude==True:
    sample_1=entry_1.create_group("sample_3")
else:
    sample_1 = entry_1.create_group("sample_1")
geometry=sample_1.create_group("geometry")
translation=geometry.create_dataset("translation",data=translation_vec)
make_whitefield=cxifile.create_group("make_whitefield")
whitefield=make_whitefield.create_dataset("whitefield",data=powderarr)
frame_selector=cxifile.create_group("frame_selector")
good_frames=frame_selector.create_dataset("good_frames",data=good_frames_arr)
if extrude==True:
    mask_save = np.ones((data_full.shape[1], data_full.shape[2])).astype("bool")
else:
    if orientation=="h":
        mask_save=np.ones((1,data_full.shape[2])).astype("bool")
    else:
        mask_save = np.ones((data_full.shape[1],1)).astype("bool")
mask_maker=cxifile.create_group("mask_maker")
mask_cxi_3=mask_maker.create_dataset("mask",data=mask_save)
print("Done.")
cxifile.close()
print("Making ini files")
ini_maker.mk_make_whitefield_ini(path="{0}{1}".format(savedir,scan_name))
ini_maker.mk_speckle_gui_ini(path="{0}{1}".format(savedir,scan_name))
ini_maker.mk_stitch_ini(path="{0}{1}".format(savedir,scan_name),roi=roi)
ini_maker.mk_update_pixel_map_ini(path="{0}{1}".format(savedir,scan_name),roi=roi)
ini_maker.mk_zernike_ini(path="{0}{1}".format(savedir,scan_name),roi=roi)
print("Done.")

