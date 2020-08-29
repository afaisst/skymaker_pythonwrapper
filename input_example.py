## Create a simulated image using SkyMaker


######## IMPORT ###########
import_file_list = ["simulate.py","simulate_to_existing.py"]
for file in import_file_list:
    exec(compile(open(file, "rb").read(), file, 'exec'))
#############################


##### 1. CREATING A NEW SIMULATED IMAGE FROM SCRATCH ############

## World properties
world_input = {"base_name":"sim1", # base simulation name (directory with this name will be created)
                "output_directory":"./sim_output/", # Directory in which Simulations are saved (in sub-folder named [base_name])
                "overwrite_source_catalog":False, # if TRUE, overwrite source catalog and create new one
                "source_density":100, # sources per arcmin2
                "image_size_arcmin":1, # image size in arcmin
                "field_center_ra":150.0, # field center RA
                "field_center_dec":2.0, # field center DEC
                "radec_distribution_type":"grid", # Distribution of galaxies: "random" or "grid"
                "mag_distribution_type":"halfnormal", # Type of the magnitude distribution: "uniform" or "halfnormal"
                "mag":[26.5,3], # in AB [max,half-normal std] for type "halfnormal", [min,max] for type "uniform" 
                "BTR":[0.2,0.9], # bulge-to-total ratio range
                "R_disk":[0.1,1], # disk length-scale range in arcsec
                "AB_disk":[0.1,1], # disk aspect ratio range
                "PA_disk":[0,180], # disk PA range (deg)
                "R_bulge_rel":[0.2,0.8], # bulge length scale range (relative to R_disk)
                "AB_bulge":[0.1,1], # bulge aspect ratio range
                "PA_bulge_rel":[-20,20], # relative to PA_disk [-deg,+deg] in degrees
                "fraction_stars":0.3, # fraction of sources that are stars (= point source)
                }

## Image input (Can create multiple images)
image_input1 = {"image_name":"hr", # name of the image
                "noise_per_pixel":0.0029295, # noise per pixel
                "zp":25.94734, # zero point AB
                "pixscale":0.03, # pixel scale (arcsec/px)
                "fake_HSC":False, # Fake HSC mask and other FITS extensions
                "psf_file_name":"./example_images/acs_psf.fits", # PSF file name
                "astro_offset":[[0,0],[0,0]] # [ [ra_mean,ra_std],[dec_mean,dec_std]] astrometry in MAS offset applied (coords = catalog_coords + offset)
                }
image_input2 = {"image_name":"lr", # name of the image
                "noise_per_pixel":0.02637168, # noise per pixel
                "zp":27, # zero point AB
                "pixscale":0.15, # pixel scale (arcsec/px)
                "fake_HSC":True, # Fake HSC mask and other FITS extensions
                "psf_file_name":"./example_images/hsc_psf.fits", # PSF file name
                "astro_offset":[[0,0],[0,0]] # [ [ra_mean,ra_std],[dec_mean,dec_std]] astrometry offset in MAS applied (coords = catalog_coords + offset)
                }


## Create simulated image
simulate(world_input=world_input,
        image_inputs=[image_input1,image_input2])



##### 2. ADD SIMULATED STARS/GALAXIES TO AN EXISTING IMAGE ############

## Word properties
world_input = {"base_name":"sim2", # base simulation name (directory with this name will be created)
                "output_directory":"./sim_output", # Directory in which Simulations are saved (in sub-folder named [base_name])
                "overwrite_source_catalog":False, # if TRUE, overwrite source catalog and create new one
                "source_density":100, # sources per arcmin2
                "radec_distribution_type":"grid", # Distribution of galaxies: "random" or "grid"
                "mag_distribution_type":"halfnormal", # Type of the magnitude distribution: "uniform" or "halfnormal"
                "mag":[26.5,3], # in AB [max,half-normal std] for type "halfnormal", [min,max] for type "uniform" 
                "BTR":[0.2,0.9], # bulge-to-total ratio range
                "R_disk":[0.1,1], # disk length-scale range in arcsec
                "AB_disk":[0.1,1], # disk aspect ratio range
                "PA_disk":[0,180], # disk PA range (deg)
                "R_bulge_rel":[0.2,0.8], # bulge length scale range (relative to R_disk)
                "AB_bulge":[0.1,1], # bulge aspect ratio range
                "PA_bulge_rel":[-20,20], # relative to PA_disk [-deg,+deg] in degrees
                "fraction_stars":0.3, # fraction of sources that are stars (= point source)
                }

## Image input (Can create multiple images)
image_input1 = {"image_name":"./example_images/acs_image.fits", # path to the image onto which simulated sources should be added
                "zp":25.94734, # zero point AB
                "psf_file_name":"./example_images/acs_psf.fits", # PSF file name
                "extensions":[0], # list of extensions to add simulated galaxies. For example [0,"IMAGE"] adds to Primary extension and extension called "IMAGE"
                "delete_noiseless_image":True  # set to True to delete the noiseless simulated image
                }
image_input2 = {"image_name":"./example_images/hsc_image.fits", # path to the image onto which simulated sources should be added
                "zp":27, # zero point AB
                "psf_file_name":"./example_images/hsc_psf.fits", # PSF file name
                "extensions":[0,"IMAGE"], # list of extensions to add simulated galaxies. For example [0,"IMAGE"] adds to Primary extension and extension called "IMAGE"
                "delete_noiseless_image":True  # set to True to delete the noiseless simulated image
                }

## Create simulated image
simulate_to_existing(world_input=world_input,
        image_inputs=[image_input1,image_input2])