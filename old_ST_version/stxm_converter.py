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
import math

#Additional parameters
roi_size = 60
year = "2020"
manual_db = (431, 825) #None
manual_orientation = "h" #Is used in case if the scan motor differs from SAMX or SAMY
stxm_roi = (100,200,1190,1479)

#input arguments
scanno = int(input("Please enter the scan number (e.g. 10):"))
targetnam = str(input("Please enter the target (Mo, Cu or Rh):"))
detector_dist = float(input("Please enter the detector distance in meters(e.g. 0.3):"))
print("All arguments are set.")

#Process input arguments
target_dict={
    "Mo": 17.48,
    "Cu": 8.,
    "Rh": 20.2}
energy = target_dict[targetnam] * 1e3 * sc.e
scan_name = "Scan_{}".format(scanno)

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
lognam = "{0}{1}.log".format(logdir,scan_name)
logfile = open(lognam,"r")
scanmotor = None


#first figure out if scan is 1D or 2D.
for line in logfile:
    if line.startswith("# Datetime"):
        if len(line.split(";"))==3:
            scandim="1d"
        elif len(line.split(";"))==4:
            scandim="2d"
print("Scan dimension is",scandim)

if scandim=="1d":
    i = 0
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
elif scandim=="2d":
    i = 0
    i_ss= None
    N_ss= None
    i_fs= None
    N_fs= None
    logfile.close()
    logfile = open(lognam, "r")
    for line in logfile:
        if line.startswith("# Device:") and not line.endswith("Lambda\n"): # and not line.endswith("Scanner")
            print(line)
            if i_ss==None:
                i_ss = i
                print("i_ss",i_ss)
            else:
                i_fs=i
                print("i_fs", i_fs)
        if i_ss != None and i == i_ss + 6:
            #ss_scanmotor = line.split(":")[1]
            N_ss = int(line.split(":")[1])
        if i_fs != None and i == i_fs + 6:
            N_fs = int(line.split(":")[1])
            #fs_scanmotor = line.split(":")[1]
        i +=1
    #print("Detected slow scanmotor:", ss_scanmotor)
    #print("Detected fast scanmotor",fs_scanmotor)
    N_logged_points=i
    N_points=int(N_ss*N_fs)
    print(N_points)
    positions_2d = np.zeros((N_points, 2)).astype("float")
    i=0
    for line in logfile:
        if not line.startswith("#"):
            current_entry_ss = line.split(";")[2]
            current_entry_fs = line.split(";")[3]
            if i == 0:
                unit_ss = current_entry_ss.split(" ")[1]
                unit_fs = current_entry_fs.split(" ")[1]
                print("The unit of the slow scan motor is", unit_ss)
                print("The unit of the fast scan motor is", unit_fs)
            try:
                unit_scaling_ss = float(unit_dict[unit_ss])
            except KeyError:
                if unit_ss.startswith("n"):
                    unit_scaling_ss = 1E-9
                elif unit_ss.startswith("u"):
                    unit_scaling_ss = 1E-6
                else:
                    unit_scaling_ss = 1.0
            try:
                unit_scaling_fs = float(unit_dict[unit_fs])
            except KeyError:
                if unit_fs.startswith("n"):
                    unit_scaling_fs = 1E-9
                elif unit_fs.startswith("u"):
                    unit_scaling_fs = 1E-6
                else:
                    unit_scaling_fs = 1.0
            current_pos = (float(current_entry_ss.split(" ")[0]) * unit_scaling_ss, float(current_entry_fs.split(" ")[0]) * unit_scaling_fs)
            positions_2d[i,:] = current_pos
            i += 1

if scandim=="1d":
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
elif scandim=="2d":
    orientation= manual_orientation

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

if scandim=="1d":
    stxm_arr=np.zeros((N_points))
elif scandim=="2d":
    stxm_arr=np.zeros((N_ss,N_fs))

for i in range(0,N_points):
    if scandim=="1d":
        i_ssn=i
    else:
        i_ssn=int(math.floor(i/N_ss))
        i_fsn=int(i-N_ss*i_ssn)
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
                raw_frames[i,:,:]=example_data
                data_insert = np.zeros_like(example_data)
                if orientation=="h":
                    data_now = np.sum(example_data[roi[0]:roi[1], :], axis=0)
                    for i1 in range(0, data_insert.shape[0], 1):
                        data_insert[i1, :] = data_now
                else:
                    data_now = np.sum(example_data[:, roi[2]:roi[3]],axis=1)
                    for i1 in range(0, data_insert.shape[1], 1):
                        data_insert[:, i1] = data_now
                data_full[i, :, :] = data_insert
                if scandim=="1d":
                    stxm_arr[i_ssn]=np.sum(data_insert[stxm_roi[0]:stxm_roi[1],stxm_roi[2]:stxm_roi[3]])
                else:
                    try:
                        stxm_arr[i_ssn,i_fsn]=np.sum(data_insert[stxm_roi[0]:stxm_roi[1],stxm_roi[2]:stxm_roi[3]])
                    except IndexError:
                        pass
        else:
            file_now = "/{0}_".format(scan_name) + "{:05.0f}_Lambda.nxs".format(i)
            f_now=h5.File(ldir+"/"+file_now,"r")
            if orientation=="h":
                data_now=np.sum(f_now["/entry/instrument/detector/data"][0,roi[0]:roi[1],:],axis=0)
                frame_now=f_now["/entry/instrument/detector/data"][0,:,:]
                raw_frames[i,:,:]=frame_now
                data_insert = np.zeros_like(example_data)
                for i1 in range(0,data_insert.shape[0],1):
                    data_insert[i1,:]=data_now
            else:
                data_now = np.sum(f_now["/entry/instrument/detector/data"][0, :, roi[2]:roi[3]],axis=1)
                frame_now=f_now["/entry/instrument/detector/data"][0,:,:]
                raw_frames[i,:,:]=frame_now
                data_insert = np.zeros_like(example_data)
                for i1 in range(0, data_insert.shape[1], 1):
                    data_insert[:, i1] = data_now
            data_full[i,:,:]=data_insert
            if scandim == "1d":
                stxm_arr[i_ssn] = np.sum(data_insert[stxm_roi[0]:stxm_roi[1],stxm_roi[2]:stxm_roi[3]])
            else:
                try:
                    stxm_arr[i_ssn, i_fsn] = np.sum(data_insert[stxm_roi[0]:stxm_roi[1],stxm_roi[2]:stxm_roi[3]])
                except IndexError:
                    pass
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
print("Saving raw frames...")
data2=data_1.create_dataset("raw_frames",data=raw_frames)
ptychogram_data=data_1.create_dataset("ptychogram",data=ptychogram)
stxm_data=data_1.create_dataset("stxm",data=stxm_arr)
instrument_1=entry_1.create_group("instrument_1")
detector_1=instrument_1.create_group("detector_1")



basis_vectors_cxi=detector_1.create_dataset("basis_vectors",data=bld_basis(whole_data=data_full,useframes=useframes))
detector_distance_cxi=detector_1.create_dataset("distance",data=float(detector_dist))
print("Energy is %s J"%energy)
print("Making wavelength")
wavelength=(sc.h*sc.c)/energy
print("Wavelength is %s m" %wavelength)



translation_vec=np.zeros((N_points,3))
if scandim=="1d":
    if orientation=="v":
        translation_vec[:,0]=positions
    else:
        translation_vec[:,1]=positions
else:
    translation_vec[:,:2]= positions_2d
translation_vec=translation_vec[useframes,:]

if scandim=="2d":
    data_1.create_dataset("N_ss", data=N_ss)
    data_1.create_dataset("N_fs",data=N_fs)
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
x = np.ones((5,5))
print("Original array:")
print(x)
print("1 on the border and 0 inside in the array")
x[1:-1,1:-1] = 0
print(x)
ini_maker.mk_make_whitefield_ini(path="{0}{1}".format(savedir,scan_name))
ini_maker.mk_speckle_gui_ini(path="{0}{1}".format(savedir,scan_name))
ini_maker.mk_stitch_ini(path="{0}{1}".format(savedir,scan_name),roi=roi)
ini_maker.mk_update_pixel_map_ini(path="{0}{1}".format(savedir,scan_name),roi=roi)
ini_maker.mk_zernike_ini(path="{0}{1}".format(savedir,scan_name),roi=roi)
print(N_points)
print("Done.")

