## Create a simulated image using SkyMaker


######## IMPORT ###########
import_file_list = ["simulate.py"]
for file in import_file_list:
    exec(compile(open(file, "rb").read(), file, 'exec'))
#############################


## Source properties
source_input = {"base_name":"1sqarcmin", # base simulation name (directory with this name will be created)
                "overwrite_source_catalog":False, # if TRUE, overwrite source catalog and create new one
                "source_density":100, # sources per arcmin2
                "image_size_arcmin":1, # image size in arcmin
                "field_center_ra":150.0, # field center RA
                "field_center_dec":2.0, # field center DEC
                "mag":[26.5,3], # in AB [min,half-normal std]
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
                "psf_file_name":"../img_data/calexp-HSC-I-9813-5_4-9813_2812.fits" # PSF file name
                }
image_input2 = {"image_name":"lr", # name of the image
                "noise_per_pixel":0.02637168, # noise per pixel
                "zp":27, # zero point AB
                "pixscale":0.15, # pixel scale (arcsec/px)
                "fake_HSC":True, # Fake HSC mask and other FITS extensions
                "psf_file_name":"../img_data/calexp-HSC-I-9813-5_4-9813_2812.fits" # PSF file name
                }


## Create simulated image
simulate(source_input=source_input,
        image_inputs=[image_input1,image_input2])

