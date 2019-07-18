import configparser
def mk_make_whitefield_ini(path):
    make_whitefield = configparser.ConfigParser()
    make_whitefield["make_whitefield"] = {
        "frames" : "/entry_1/data_1/data ;str, location of diffraction data",
        "sigma_t" : "10 ;None, or int, number of adjascent frames to average (None disables)",
    }
    make_whitefield["make_whitefield-advanced"] = {
        "h5_group" : "make_whitefield ;str, name of h5 group to write to",
    }
    with open('{0}/make_whitefield.ini'.format(path), 'w') as make_whitefield_file:
        make_whitefield.write(make_whitefield_file)

def mk_speckle_gui_ini(path):
    speckle_gui = configparser.ConfigParser()
    speckle_gui["speckle-gui"] = {
        "data_paths" : "['/entry_1/data_1/data']",
        "mask_paths" : "['/mask_maker/mask', 'entry_1/instrument_1/detector_1/mask']",
        "translation_paths" : "['/pos_refine/translation', '/entry_1/sample_3/geometry/translation']",
        "whitefield_paths" : "['/make_whitefield/whitefield', '/process_3/whitefield']",
        "good_frames_paths" : "['/frame_selector/good_frames', '/process_3/good_frames']"
    }
    with open('{0}/speckle-gui.ini'.format(path), 'w') as speckle_gui_file:
        speckle_gui.write(speckle_gui_file)

def mk_stitch_ini(path,roi):
    stitch = configparser.ConfigParser()
    stitch["stitch"] = {
        "roi" : "{0} ;tuple of length 4, (ss_start, ss_end, fs_start, fs_end)".format(roi),
        "whitefield" : "/make_whitefield/whitefield ;str, location of whitefield",
        "good_frames" : "/frame_selector/good_frames ;str, location of good frames list",
        "defocus" : "0.005  ;float, distance between focus and sample in meters",
        "reg" : "50 ;float, regularisation parameter in pixels (biases diffraction near the centre of the pupil)",
        "pixel_shifts" : "None ;str, location of pixel shift arrays (due to phase gradient effects)",
        "sub_pixel" : "True ;True or False, if True then sub pixel precision is enabled (slow)",
    }
    stitch["stitch-advanced"] = {
        "h5_group" : "stitch",
        "frames" : "/entry_1/data_1/data",
        "y_pixel_size": "/entry_1/instrument_1/detector_1/y_pixel_size",
        "x_pixel_size": "/entry_1/instrument_1/detector_1/x_pixel_size",
        "distance" : "/entry_1/instrument_1/detector_1/distance",
        "mask" : "/mask_maker/mask",
        "energy" : "/entry_1/instrument_1/source_1/energy",
        "translation" : "/entry_1/sample_3/geometry/translation",
        "basis_vectors" : "/entry_1/instrument_1/detector_1/basis_vectors",
    }
    with open('{0}/stitch.ini'.format(path), 'w') as stitch_file:
        stitch.write(stitch_file)

def mk_update_pixel_map_ini(path,roi):
    update_pixel_map = configparser.ConfigParser()
    update_pixel_map["update_pixel_map"] = {
        "roi" : "{0}".format(roi),
        "whitefield" : "/make_whitefield/whitefield",
        "good_frames" : "/frame_selector/good_frames",
        "defocus" : "0.005",
        "max_step" : "4.0",
        "pixel_shifts" : "None",
        "sub_pixel" : "True",
        "atlas" : "/stitch/O",
    }
    update_pixel_map["update_pixel_map-advanced"] = {
        "h5_group" : "update_pixel_map",
        "frames" : "/entry_1/data_1/data",
        "y_pixel_size": "/entry_1/instrument_1/detector_1/y_pixel_size",
        "x_pixel_size": "/entry_1/instrument_1/detector_1/x_pixel_size",
        "distance" : "/entry_1/instrument_1/detector_1/distance",
        "mask" : "/mask_maker/mask",
        "energy" : "/entry_1/instrument_1/source_1/energy",
        "translation" : "/entry_1/sample_3/geometry/translation",
        "basis_vectors" : "/entry_1/instrument_1/detector_1/basis_vectors",
    }
    with open('{0}/update_pixel_map.ini'.format(path), 'w') as upm_file:
        update_pixel_map.write(upm_file)

def mk_zernike_ini(path,roi):
    zernike = configparser.ConfigParser()
    zernike["zernike"] = {
        "roi" : "{0}".format(roi),
        "whitefield" : "/make_whitefield/whitefield",
        "pixel_shifts": "/update_pixel_map/pixel_shifts",
        "orders" : "36",
        "defocus" : "0.005",
        "remove_tilt" : "True",
        "remove_piston" : "True",
        "remove_defocus" : "True",
        "remove_astigmatism" : "True",
    }
    zernike["zernike-advanced"] = {
        "h5_group" : "zernike",
        "y_pixel_size": "/entry_1/instrument_1/detector_1/y_pixel_size",
        "x_pixel_size": "/entry_1/instrument_1/detector_1/x_pixel_size",
        "distance" : "/entry_1/instrument_1/detector_1/distance",
        "mask" : "/mask_maker/mask",
        "energy" : "/entry_1/instrument_1/source_1/energy",
        "basis_vectors" : "/entry_1/instrument_1/detector_1/basis_vectors",
    }
    with open('{0}/zernike.ini'.format(path), 'w') as zern_file:
        zernike.write(zern_file)

