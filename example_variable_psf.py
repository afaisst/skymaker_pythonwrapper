## This script serves as example to add simulated sources to an existing image with an attached PSF.
# The code looks up the PSF for an existing image and uses it to add simulated sources. 
# The goal is to add the simulated sources to a high-resolution and low-resolution image.
# The images are from HSC (low-resolution) and ACS (high-resolution)
# All the paths need to be changed.
########

######## IMPORT ###########
import glob

import_file_list = ["simulate.py","simulate_to_existing.py"]
for file in import_file_list:
    exec(compile(open(file, "rb").read(), file, 'exec'))
#############################


## Path definitions
lr_image_path = "/stage/irsa-jointproc-data03/ACS_COSMOS/grizli-stacks/HSC_corresponding_calexp_oct152019/" #calexp-HSC-I-9813-5_7.fits/"
hr_image_path = "/stage/irsa-jointproc-data03/ACS_COSMOS/from_irsa/downloads/" #calexp-HSC-I-9813-5_7-2128_acs_I_mosaic_30mas_sci.fits"
lr_psf_path = "/stage/irsa-jointproc-data02/PSFex/allcalexp_psfex_pipeline/output_variableMag/psfs/" #calexp-HSC-I-9813-5_7-9813_2128_m21.fits
hr_psf_path = "/stage/irsa-jointproc-data03/ACS_COSMOS/from_irsa/sextractor/psf/" #HSC-I-9813-5_4-2812_psf.fits"

## List of images to add simulated sources
lr_image_list = ["calexp-HSC-I-9812-0_1.fits","calexp-HSC-I-9812-0_2.fits"]
hr_image_list = ["calexp-HSC-I-9812-0_1-4546_acs_I_mosaic_30mas_sci.fits","calexp-HSC-I-9812-0_2-4402_acs_I_mosaic_30mas_sci.fits"]

## Translation from image to PSF (* is placeholder)
lr_psf_name_template = "calexp-HSC-I-%s-%s-%s_*_m21.fits" # tract , patch (format x_y), tract
hr_psf_name_template = "HSC-I-9813-5_4-2812_psf.fits" # just one PSF



## For each image, create a dictionary for low and high resolution image and create simulated image.

# just check if lists are of same length
if len(lr_image_list) != len(hr_image_list):
    print("Image list have not same length. Check.")
    quit()



for hr_image , lr_image in zip(hr_image_list , lr_image_list):

    ## Get name of simulation. This changes for each patch because different location on sky.
    simulation_name = "sim_%s" % lr_image.split(".fits")[0]

    ## World properties
    world_input = {"base_name":simulation_name , # base simulation name (directory with this name will be created)
                    "output_directory":"./sim_output", # Directory in which Simulations are saved (in sub-folder named [base_name])
                    "overwrite_source_catalog":False, # if TRUE, overwrite source catalog and create new one
                    "source_density":100, # sources per arcmin2
                    "mag":[26.5,3], # in AB [max,half-normal std]
                    "BTR":[0.2,0.9], # bulge-to-total ratio range
                    "R_disk":[0.1,1], # disk length-scale range in arcsec
                    "AB_disk":[0.1,1], # disk aspect ratio range
                    "PA_disk":[0,180], # disk PA range (deg)
                    "R_bulge_rel":[0.2,0.8], # bulge length scale range (relative to R_disk)
                    "AB_bulge":[0.1,1], # bulge aspect ratio range
                    "PA_bulge_rel":[-20,20], # relative to PA_disk [-deg,+deg] in degrees
                    "fraction_stars":0.3, # fraction of sources that are stars (= point source)
                    }


    # for high-resolution image
    hr_psf_file_name = os.path.join( hr_psf_path , hr_psf_name_template)
    image_input_hr = {"image_name": os.path.join(hr_image_path , hr_image),
                    "zp":25.94734,
                    "psf_file_name": hr_psf_file_name,
                    "extensions":[0],
                    "delete_noiseless_image":True
                        }

    #print(image_input_hr)

    # for low-resolution image
    tract = lr_image.split("-")[3]
    patch = lr_image.split("-")[4].split(".fits")[0]
    lr_psf_file_name = glob.glob( os.path.join( lr_psf_path , lr_psf_name_template % (str(tract) , patch , str(tract)) ) )[0]
    
    image_input_lr = {"image_name": os.path.join(lr_image_path , lr_image),
                    "zp":27.0,
                    "psf_file_name": lr_psf_file_name,
                    "extensions":["IMAGE"],
                    "delete_noiseless_image":False
                        }

    print(image_input_lr)


    # simulate image
    #simulate_to_existing(world_input=world_input,
    #    image_inputs=[image_input_hr,image_input_lr])
    simulate_to_existing(world_input=world_input,
        image_inputs=[image_input_lr])