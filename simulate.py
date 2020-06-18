## Simulates the images


###### IMPORTS #########

import os, sys

import numpy as np

import subprocess
import random
import time

import pickle

import sh

from astropy.io import fits, ascii
from astropy.table import Table, Column, MaskedColumn, hstack, vstack
import astropy.wcs as wcs
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.nddata import Cutout2D

import matplotlib.pyplot as plt
import matplotlib as mpl

from scipy import stats
from scipy import ndimage

######## FUNCTIONS ############

## This function creates config file for SkyMaker from an input parameter list.
def create_skymaker_config(filename,params):

    ## load config text
    config_txt = '''# Default configuration file for SkyMaker 3.10.5
# EB 2020-05-18
#
 
#--------------------------------- Image -------------------------------------
 
IMAGE_NAME      *IMAGE_NAME*        # Name of the output frame
IMAGE_SIZE      *IMAGE_SIZE*            # Width,[height] of the output frame
IMAGE_TYPE      *IMAGE_TYPE*             # PUPIL_REAL,PUPIL_IMAGINARY,PUPIL_MODULUS,
                                # PUPIL_PHASE,PUPIL_MTF,PSF_MTF,PSF_FULLRES,
                                # PSF_FINALRES,SKY_NONOISE,SKY,GRID
                                # or GRID_NONOISE
GRID_SIZE       64              # Distance between objects in GRID mode
IMAGE_HEADER    *IMAGE_HEADER*        # File name or INTERNAL
LISTCOORD_TYPE  *LISTCOORD_TYPE*           # Coordinates in input lists: PIXEL or WORLD
 
#-------------------------------- Detector -----------------------------------
 
GAIN            *GAIN*             # gain (e-/ADU)
WELL_CAPACITY   0               # full well capacity in e- (0 = infinite)
SATUR_LEVEL     *SATUR_LEVEL*           # saturation level (ADU)
READOUT_NOISE   *READOUT_NOISE*             # read-out noise (e-)
EXPOSURE_TIME   *EXPOSURE_TIME*           # total exposure time (s)
MAG_ZEROPOINT   *MAG_ZEROPOINT*            # magnitude zero-point ("ADU per second")

#-------------------------------- Sampling -----------------------------------
 
PIXEL_SIZE      *PIXEL_SIZE*           # pixel size in arcsec.
MICROSCAN_NSTEP 1               # number of microscanning steps (1=no mscan)
 
#---------------------------------- PSF --------------------------------------
 
PSF_TYPE        *PSF_TYPE*        # INTERNAL or FILE
PSF_NAME        *PSF_NAME*        # Name of the FITS image containing the PSF
PSFCENTER_TYPE  *PSFCENTER_TYPE*       # UPPERHALF, LOWERHALF, HALF, CENTROID,
                                # CENTROID_COMMON or PEAK
SEEING_TYPE     NONE   # (NONE, LONG_EXPOSURE or SHORT_EXPOSURE)
SEEING_FWHM     0.1             # FWHM of seeing in arcsec (incl. motion)
AUREOLE_RADIUS  0             # Range covered by aureole (pix) 0=no aureole
AUREOLE_SB      16.0            # SB (mag/arcsec2) at 1' from a 0-mag star
PSF_OVERSAMP    *PSF_OVERSAMP*               # Oversampling factor / final resolution
PSF_MAPSIZE     *PSF_MAPSIZE*            # PSF mask size (pixels): must be a power of 2
TRACKERROR_TYPE *TRACKERROR_TYPE*            # Tracking error model: NONE, DRIFT or JITTER
TRACKERROR_MAJ  *TRACKERROR_MAJ*             # Tracking RMS error (major axis) (in arcsec)
TRACKERROR_MIN  *TRACKERROR_MIN*             # Tracking RMS error (minor axis) (in arcsec)
TRACKERROR_ANG  *TRACKERROR_ANG*             # Tracking angle (in deg, CC/horizontal)
 
#----------------------------- Pupil features --------------------------------
 
M1_DIAMETER     2.5             # Diameter of the primary mirror (in meters)
M2_DIAMETER     0.5             # Obstruction diam. from the 2nd mirror in m.
ARM_COUNT       0               # Number of spider arms (0 = none)
ARM_THICKNESS   0.0            # Thickness of the spider arms (in mm)
ARM_POSANGLE    0.0             # Position angle of the spider pattern / AXIS1
DEFOC_D80       0.0             # Defocusing d80% diameter (arcsec)
DEFOC_CENTER    0.5,0.5         # Relative center of PSF focus variations
SPHER_D80       0.0             # Spherical d80% diameter (arcsec)
SPHER_CENTER    0.5,0.5         # Center of PSF spherical aber. variations
COMAX_D80       0.0             # Coma along X d80% diameter (arcsec)
COMAY_D80       0.0             # Coma along Y d80% diameter (arcsec)
COMA_CENTER     0.5,0.5         # Center of PSF coma variations
AST00_D80       0.0             # 0 deg. astigmatism d80% diameter (arcsec)
AST45_D80       0.0             # 45 deg. astigmatism d80% diameter (arcsec)
AST_CENTER      0.5,0.5         # Center of PSF astigmatism variations
TRI00_D80       0.0             # 0 deg. triangular d80% diameter (arcsec)
TRI30_D80       0.0             # 30 deg. triangular d80% diameter (arcsec)
TRI_CENTER      0.5,0.5         # Center of PSF triangular aber. variations
QUA00_D80       0.0             # 0 deg. quadratic d80% diameter (arcsec)
QUA22_D80       0.0             # 22.5 deg. quadratic d80% diameter (arcsec)
QUA_CENTER      0.5,0.5         # Center of PSF quad. aber. variations
 
#--------------------------------- Signal ------------------------------------
 
WAVELENGTH      *WAVELENGTH*             # average wavelength analysed (microns)
BACK_MAG        *BACK_MAG*            # background surface brightness (mag/arcsec2)
 
#------------------------------ Stellar field --------------------------------
 
STARCOUNT_ZP    *STARCOUNT_ZP*             # nb of stars /deg2 brighter than MAG_LIMITS
STARCOUNT_SLOPE *STARCOUNT_SLOPE*             # slope of differential star counts (dexp/mag)
MAG_LIMITS      *MAG_LIMITS*       # stellar magnitude range allowed
 
#------------------------------ Random Seeds ---------------------------------
 
SEED_MOTION     *SEED_MOTION*               # rand. seed for PSF turbulent motion (0=time)
SEED_STARPOS    *SEED_STARPOS*               # random seed for star positions (0=time)
 
#----------------------------- Miscellaneous ---------------------------------
 
VERBOSE_TYPE    *VERBOSE_TYPE*          # QUIET, NORMAL or FULL
NTHREADS        *NTHREADS*               # Number of simultaneous threads for
                                # the SMP version of SkyMaker
'''
    ## replace the keywords with the param set
    for key,val in params.items():
        config_txt = config_txt.replace("*%s*" % key , str(val))

    ## Save to file
    with open(filename , "w") as f:
        f.write(config_txt)
    
    return(True)



# save a header object to a file
def save_hdr_to_file(hdr_object,filename):
    with open(filename , "w") as f:
        f.write(str(hdr_object))
    return(True)

# reads a header from a file
def read_hdr_from_file(filename):
    with open(filename, "r") as f:
        txt = f.read()
    hdr = fits.Header.fromstring(txt)
    return( hdr )


# save a wcs object to a file
def save_wcs_to_file(wcs_object,filename):
    wcs_text = wcs_object.to_header_string()
    with open(filename , "w") as f:
        f.write(wcs_object)
    return(True)

# reads a wcs from a file
def read_wcs_from_file(filename):
    with open(filename, "r") as f:
        txt = f.read()
    return( wcs.WCS(txt) )

# Returns uniform random values
# Range must be [min,max] and size is number of random numbers
def get_random_uniform(range,size):
    return( np.random.uniform(low=range[0],high=range[1],size=size) )

## Creates half-normal magnitude distribution with max at "mag_max" and a tail with sigma of "mag_sig" towards brighter magnitudes
# This allows to simulate a simple magnitude distribution
def get_mag_distribution(mag_max,mag_sig,size):
    
    b = np.repeat(0,1)
    while(len(b) < size): # need to have this while loop to make sure that we have more galaxies than the size. It is easier to remove galaxies than to add some.
        a = np.random.normal(loc=0,scale=mag_sig, size=2*size)
        a = a + mag_max
        b = a[a < mag_max]
    
    b = b[:size] # since we know that there are more galaxies than requested, we can simply cut some random ones at the end.
    
    return(b)


## Convert magnitude to fluxes
def convert_mag_to_flux(mag,zp):
    return 10**(-0.4*(mag - zp))


## Create Mask Header for HSC to make fake mask
def create_hsc_mask_header(base_header):
    
    # additional keywords to add
    mask_header = {"EXTTYPE":"MASK",
                   "EXTNAME":"MASK",
                   "MP_BAD":0,
                   "MP_BRIGHT_OBJECT":9,
                   "MP_CLIPPED":13,
                   "MP_CR":3,
                   "MP_CROSSTALK":10,
                   "MP_DETECTED":5,
                   "MP_DETECTED_NEGATIVE":6,
                   "MP_EDGE":4,
                   "MP_INTRP":2,
                   "MP_NOT_DEBLENDED":11,
                   "MP_NO_DATA":8,
                   "MP_SAT":1,
                   "MP_SUSPECT":7,
                   "MP_UNMASKEDNAN":12
              }
    
    for key in mask_header:
        base_header[key] = mask_header[key]
        
    return base_header

##################################


def simulate(world_input,image_inputs):

    ## Create output directory --------------
    output_directory = os.path.join(world_input["output_directory"] , world_input["base_name"])
    if not os.path.exists(output_directory):
        print("Creating output directory: %s" % output_directory)
        sh.mkdir(output_directory)
    else:
        print("Output directory exists")

    ## Define some file names (easier to do it here for book keeping)
    FILES = dict()
    FILES["source_list"] = os.path.join(output_directory , "sources.csv")

    for image_id,image_input in enumerate(image_inputs):
        FILES["skymaker_list_%g" % image_id] = os.path.join(output_directory , "%s_input.list" % image_input["image_name"])
        FILES["header_file_%g" % image_id] = os.path.join(output_directory,"%s_header.txt" % image_input["image_name"])
        FILES["image_output_%g" % image_id] = os.path.join(output_directory,"%s.fits" % image_input["image_name"])
        FILES["skymaker_config_%g" % image_id] = os.path.join(output_directory,"%s.config" % image_input["image_name"])

    if not os.path.exists(FILES["source_list"]):
        CREATENEWCATALOG = True
    elif os.path.exists(FILES["source_list"]) & (world_input["overwrite_source_catalog"] == True):
        CREATENEWCATALOG = True
    else:
        CREATENEWCATALOG = False

    ### THE FOLLOWING DOES NOT DEPEND ON THE IMAGE INPUTS ########

    ## check if external catalog available. If so, load, else create new sources -------------
    GALDATA = dict()
    if CREATENEWCATALOG: # create new source catalog
        print("Creating new source catalog")

        ## a) Create galaxy positions
        GALDATA["nbr_gals"] = int( world_input["source_density"] * world_input["image_size_arcmin"]**2 )

        ras = np.random.uniform(low=world_input["field_center_ra"] - (world_input["image_size_arcmin"]/2./60. - 3/3600),
                                high=world_input["field_center_ra"] + (world_input["image_size_arcmin"]/2./60. - 3/3600),
                                size=int(GALDATA["nbr_gals"]))

        decs = np.random.uniform(low= world_input["field_center_dec"] - (world_input["image_size_arcmin"]/2./60. - 3/3600),
                                high=world_input["field_center_dec"] + (world_input["image_size_arcmin"]/2./60. - 3/3600),
                                size=int(GALDATA["nbr_gals"]))

        ids = np.arange(1,GALDATA["nbr_gals"]+1)

        ## b) magnitudes for sources
        mags = get_mag_distribution(mag_max=world_input["mag"][0],mag_sig=world_input["mag"][1],size=GALDATA["nbr_gals"]) # this generates a half-normal distribution with tail to brigther magnitudes.

        ## c) create structural parameters of sources for SkyMaker
        BTRs = get_random_uniform(world_input["BTR"],GALDATA["nbr_gals"])
        R_disks = get_random_uniform(world_input["R_disk"],GALDATA["nbr_gals"])
        AB_disks = get_random_uniform(world_input["AB_disk"],GALDATA["nbr_gals"])
        PA_disks = get_random_uniform(world_input["PA_disk"],GALDATA["nbr_gals"])
        R_bulges = R_disks * get_random_uniform(world_input["R_bulge_rel"],GALDATA["nbr_gals"])
        AB_bulges = get_random_uniform(world_input["AB_bulge"],GALDATA["nbr_gals"])

        if world_input["PA_bulge_rel"][0] > 0: world_input["PA_bulge_rel"][0] = world_input["PA_bulge_rel"]*(-1) # just check...
        PA_bulges = PA_disks + get_random_uniform(world_input["PA_bulge_rel"],GALDATA["nbr_gals"])
        
        ## d) Turn some sources into stars
        categories = np.repeat(200,GALDATA["nbr_gals"]) # 200 = galaxy
        npointsource = int(np.floor(world_input["fraction_stars"] * GALDATA["nbr_gals"]))
        sel_pointsource = np.random.choice(np.arange(0,GALDATA["nbr_gals"]),size=npointsource,replace=False).astype("int")
        categories[sel_pointsource] = 100 # 100 = star
        R_disks[sel_pointsource] = -1
        AB_disks[sel_pointsource] = -1
        PA_disks[sel_pointsource] = -1
        R_bulges[sel_pointsource] = -1
        AB_bulges[sel_pointsource] = -1
        PA_bulges[sel_pointsource] = -1


        ## e) put all in catalog ---------
        # ID RA DEC category mtot BTR R_bulge AB_bulge PA_bulge R_disk AB_disk PA_disk
        galtab = Table( [ids,ras,decs,categories,mags,BTRs,R_bulges,AB_bulges,PA_bulges,R_disks,AB_disks,PA_disks],
                        names=["id","ra","dec","category","magtot","BTR","R_bulge","AB_bulge","PA_bulge","R_disk","AB_disk","PA_disk"],
                        dtype=[np.int,np.float,np.float,np.int,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float])

       

        ## f) save catalog ---------
        galtab.write(FILES["source_list"] , format="csv" , overwrite=True)



    else:
        print("Loading existing source catalog; %s" % FILES["source_list"])

        # load catalog
        galtab = ascii.read(FILES["source_list"] , format="csv")
        GALDATA["nbr_gals"] = len(galtab)

    print("Total number of sources: %g" % len(galtab))
    print("Number of galaxies to be simulated: %g" % len(np.where(galtab["category"] == 200)[0]) )
    print("Number of point sources: %g" % len(np.where(galtab["category"] == 100)[0]))


    ######## THE FOLLOWING *DOES* DEPEND ON IMAGE INPUT #############
    ## iterate over the image inputs here

    for image_id,image_input in enumerate(image_inputs):


        ## Once catalog is created (or loaded) make SkyMaker catalog file  --------------

        ## a) Need to get the X and Y positions.
        ## For this, I want to create a placeholder header.

        # create new base WCS header
        hdu_wcs = fits.PrimaryHDU()

        hdu_wcs.header["EXTTYPE"] = "IMAGE"
        hdu_wcs.header["EXTNAME"] = "IMAGE"

        hdu_wcs.header["NAXIS"] = 2

        hdu_wcs.header["CD1_1"] = (-1)*image_input["pixscale"]/3600
        hdu_wcs.header["CD2_2"] = image_input["pixscale"]/3600
        hdu_wcs.header["CD1_2"] = 0
        hdu_wcs.header["CD2_1"] = 0

        hdu_wcs.header["CRVAL1"] = world_input["field_center_ra"] + world_input["image_size_arcmin"]/2/60
        hdu_wcs.header["CRVAL2"] = world_input["field_center_dec"] - world_input["image_size_arcmin"]/2/60

        hdu_wcs.header["CRPIX1"] = 0
        hdu_wcs.header["CRPIX2"] = 0

        hdu_wcs.header["CTYPE1"] = "RA---TAN"
        hdu_wcs.header["CTYPE2"] = "DEC--TAN"

        hdu_wcs.header["CDELT1"] = (-1)*image_input["pixscale"]/3600# -1.0
        hdu_wcs.header["CDELT2"] = image_input["pixscale"]/3600#1.0

        hdu_wcs.header["NAXIS1"] = int(world_input["image_size_arcmin"]*60/image_input["pixscale"])
        hdu_wcs.header["NAXIS2"] = int(world_input["image_size_arcmin"]*60/image_input["pixscale"])

        hdu_wcs.header["ZP"] = image_input["zp"]
        hdu_wcs.header["PIXSCALE"] = image_input["pixscale"]
        hdu_wcs.header["NOISE"] = image_input["noise_per_pixel"]


        # Just add some zeros for now to the data
        #hdu_wcs.data = np.zeros((hdu_wcs.header["NAXIS1"],hdu_wcs.header["NAXIS2"]))

        # get the WCS
        wcs_sim = wcs.WCS(hdu_wcs.header)

        # save the header information to a file
        save_hdr_to_file(hdr_object=hdu_wcs.header,filename=FILES["header_file_%g" % image_id] )

        # now convert the RA and DEC to X and Y
        tmp = [ wcs_sim.all_world2pix([ [galtab["ra"][ii],galtab["dec"][ii]] ] , 0) for ii in range(len(galtab["dec"])) ]
        Xs = [tmp[ii][0][0] for ii in range(len(galtab["dec"]))]
        Ys = [tmp[ii][0][1] for ii in range(len(galtab["dec"]))]

        ## b) Save the SkyMaker catalog list
        # ID X Y mtot BTR R_bulge AB_bulge PA_bulge R_disk AB_disk PA_disk
        skymakerlist = Table( [galtab["category"],Xs,Ys,galtab["magtot"],galtab["BTR"],galtab["R_bulge"],galtab["AB_bulge"],
                                galtab["PA_bulge"],galtab["R_disk"],galtab["AB_disk"],galtab["PA_disk"]],
                        names=["id","X","Y","magtot","BTR","R_bulge","AB_bulge","PA_bulge","R_disk","AB_disk","PA_disk"],
                        dtype=[np.int,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float])

        skymakerlist.write(FILES["skymaker_list_%g" % image_id] ,
                            format="ascii.commented_header" ,
                            overwrite=True,
                            formats={"X":"%4.4f",
                                    "Y":"%4.4f",
                                    "magtot":"%4.2f",
                                    "BTR":"%4.2f",
                                    "R_bulge":"%4.2f",
                                    "AB_bulge":"%4.2f",
                                    "PA_bulge":"%4.2f",
                                    "R_disk":"%4.2f",
                                    "AB_disk":"%4.2f",
                                    "PA_disk":"%4.2f"})

        ## Create the SkyMaker configuration file --------------------
        params = {"IMAGE_NAME":FILES["image_output_%g" % image_id],
                "IMAGE_SIZE":int(world_input["image_size_arcmin"]*60/image_input["pixscale"]),
                "IMAGE_TYPE":"SKY",
                "IMAGE_HEADER":"INTERNAL",
                "LISTCOORD_TYPE":"PIXEL",
                "GAIN":1000.0,
                "SATUR_LEVEL":100000,
                "READOUT_NOISE": image_input["noise_per_pixel"]*1000.0, # noise x gain
                "EXPOSURE_TIME": 1.0, # Always 1
                "MAG_ZEROPOINT": image_input["zp"],
                "PIXEL_SIZE": image_input["pixscale"],
                "PSF_TYPE":"FILE",
                "PSF_NAME":image_input["psf_file_name"],
                "PSFCENTER_TYPE":"CENTROID",
                "PSF_OVERSAMP":1,
                "PSF_MAPSIZE":1024,
                "TRACKERROR_TYPE":"NONE",
                "TRACKERROR_MAJ":0.0,
                "TRACKERROR_MIN":0.0,
                "TRACKERROR_ANG":0.0,
                "WAVELENGTH":0.8, # in microns
                "BACK_MAG":35.0, # background surface brightness (mag/arcsec2) (set to something small)
                "STARCOUNT_ZP":0,
                "STARCOUNT_SLOPE":0.2,
                "MAG_LIMITS":"17.0,26.0",
                "SEED_MOTION":1,
                "SEED_STARPOS":1,
                "VERBOSE_TYPE":"NORMAL",
                "NTHREADS":0
                }
        create_skymaker_config(filename=FILES["skymaker_config_%g" % image_id],params=params)

        ## Finally run SkyMaker -------------------------
        cmd = "sky %s -c %s" % (FILES["skymaker_list_%g" % image_id] , FILES["skymaker_config_%g" % image_id])
        print("Running SKYMAKER  . . . " , end="")
        subprocess.run(cmd , shell=True)
        print("Done!")


        ## Add WCS to header (and Mask and Variance if requested) and save -------------------

        # First read in the file
        with fits.open(FILES["image_output_%g" % image_id]) as hdul:
            hdu_sim = hdul[0].copy()

        # Prepare IMAGE header
        keys_to_copy = ["EXTTYPE","EXTNAME",
                        "CD1_1","CD2_2","CD1_2","CD2_1",
                        "CRVAL1","CRVAL2",
                        "CRPIX1","CRPIX2",
                        "CTYPE1","CTYPE2",
                        "CDELT1","CDELT2",
                        "ZP","PIXSCALE","NOISE"
                        ]

        for key in keys_to_copy:
            hdu_sim.header[key] = hdu_wcs.header[key]

        if image_input["fake_HSC"] == True:
            print("Adding Mask and Variance extension to FITS")

            hdu_sim.header["EXTTYPE"] = "PRIMARY"
            hdu_sim.header["EXTNAME"] = "PRIMARY"

            hdu_image = fits.ImageHDU()
            hdu_image.data = hdu_sim.data.copy()
            hdu_image.header = hdu_sim.header.copy()
            hdu_image.header["EXTTYPE"] = "IMAGE"
            hdu_image.header["EXTNAME"] = "IMAGE"

            hdu_mask = fits.ImageHDU()
            hdu_mask.data = np.zeros((hdu_sim.header["NAXIS1"],hdu_sim.header["NAXIS2"])).astype("int")
            hdu_mask.header = create_hsc_mask_header(base_header = hdu_sim.header.copy())
            hdu_mask.header["EXTTYPE"] = "MASK"
            hdu_mask.header["EXTNAME"] = "MASK"

            hdu_var = fits.ImageHDU()
            hdu_var.data = np.ones_like(hdu_sim.data) * (image_input["noise_per_pixel"]**2) # variance (NOT inverse variance)
            hdu_var.header = hdu_sim.header.copy()
            hdu_var.header["EXTTYPE"] = "VARIANCE"
            hdu_var.header["EXTNAME"] = "VARIANCE"

            hdul_final = fits.HDUList([hdu_sim,hdu_image,hdu_mask,hdu_var])
        
        else:

            hdul_final = fits.HDUList([hdu_sim])
        
        # verify and save again
        hdul_final.verify("silentfix")
        hdul_final.writeto(FILES["image_output_%g" % image_id] , overwrite=True)

    ## End iteration over image inputs

    # print list of files
    #print("List of Files:")
    #print(FILES)

