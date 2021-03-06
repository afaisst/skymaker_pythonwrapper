# Python Wrapper for SkyMaker

This is a Python wrapper for SkyMaker (https://www.astromatic.net/software/skymaker) to simulate astronomical images.
The wrapper allows to create a new simulated image or add simulated galaxies/stars to an existing image.

## 1. Installation

### a) SkyMaker Installation 
You need to have SkyMaker installed on your computer. The latest version can be downloaded from the astromatic website (https://www.astromatic.net/software/skymaker). To install SkyMaker, simply go to the directory and run
```
./configure
make
make install
```

If you have troubles installing SkyMaker, it's likely because of the FFTW dependency that is needed. You can download FFTW from http://www.fftw.org/download.html. Install it with
```
./configure
make
make install
```

After FFTW is installed, try to install SkyMaker again using the above steps. If you get an error like "configure: error: FFTW double precision library files not found in /sw/lib! Exiting.", then you need to do some additional steps.

First install FFTW with threads in double precision mode:
```
./configure --enable-threads
make
make install
```

Second, install FFTW with threads in single precision mode:
```
./configure --enable-single --enable-threads
make
install
```

Then install skymaker again using 
```
./configure
make
make all
```

This should work.

### b) Install the Python Wrapper

Simply clone this repository to your computer:
```
git clone https://github.com/afaisst/skymaker_pythonwrapper.git
```


## 2. Usage


### a) Creating a new simulated image

An example is included in the script "input_example.py". Simply run
```
python input_example.py
```
(Note that I am using Python 3.7, I haven't tested the code with other versions of Python.)



The simulation is created using the function (main function in "simulate.py")
```
simulate(world_input,image_input)
```


Note that "image_input" can be a list, for example
```
simulate(world_input,image_input=[image_input1,image_input2,image_input3, ...])
```

The "world_input" contains parameters to create the sources (stars and galaxies) while the "image_input" contains parameters defining the image itself (e.g., pixel size, zero point, etc). In the above example, the same galaxies/stars will used to simulate images 1, 2, 3, etc.

The "world_input" and "image_input" are both in dictionary format as explained below.

#### World Dictionary
Contains all the information about the sources. The parameters are:

- base_name: base simulation name (directory with this name will be created inside the "output_directory") [string]
- output_directory: Directory in which Simulations are saved (in sub-folder named according to "base_name") [string]
- overwrite_source_catalog: if TRUE, overwrite source catalog and create new one [True/False]
- source_density: sources per arcmin2 [float]
- image_size_arcmin: image size in arcmin [float]
- field_center_ra: field center RA (for creating WCS) [float]
- field_center_dec: field center DEC (for creating WCS) [float]
- radec_distribution_type:"random", # Distribution of galaxies on simulated sky: "random" or "grid" (evenly spaced)
- mag_distribution_type: Type of the magnitude distribution: "uniform" or "halfnormal" [string]
- mag: Magnitude range for objects in AB magnitudes. [max=faintest , half-normal std] for "halfnormal" magnitude distribution and [min=brightest , max=faintest] for "uniform" magnitude distribution [float , float]
- BTR: bulge-to-total ratio range [min=0 , max=1]
- R_disk: disk length-scale range in arcsec [min , max]
- AB_disk: disk aspect ratio range [min=0 , max=1]
- PA_disk: disk position angle (PA) range in degrees [min=0 , max=360]
- R_bulge_rel: bulge length scale range (relative to R_disk)  [min=0 , max=1]
- AB_bulge: bulge aspect ratio range [min=0 , max=1]
- PA_bulge_rel: Bulge PA relative to PA_disk in degrees [-deg,+deg]
- fraction_stars: fraction of sources that are stars (= point source) [min=0,max=1]

See the file "input_example.py" for example values.

#### Images Dictionary
Several image dictionaries can be created and fed to the simulator as a list.
Each dictionary must contain the following parameters:

- image_name: name of the image (for example "hr" for high resolution) [string]
- noise_per_pixel: noise per pixel in image units [float]
- zp: zero point AB [float]
- pixscale: pixel scale (arcsec/px) [float]
- fake_HSC: If TRUE, a fake HSC FITS extensions are created. This includes a MASK extension (here just zeros) and a VARIANCE extension (here constant given by the square of "noise_per_pixel"). [True/False]
- psf_file_name: File name of the PSF (in FITS format) to be used [string]
- astro_offset: Astrometry offset in MAS (coords = source_catalog_coords + offset) [ [ra_mean,ra_std],[dec_mean,dec_std]]

See the file "input_example.py" for example values.

Note on astro_offset: This allows the user to add astrometric offsets to the galaxies (e.g., RA_offset = RA + offset). This is simply done by changing the RA/DEC by the offset (requested in milli-arcseconds, mas). The source catalog (sources.csv) contains the original RA/DEC coordinates. In addition, a catalog "sources_[image_name].csv" is created including the same columns as "sources.csv" but in addition the X and Y coordinates (image coordinates), the offsets in RA/DEC as well as the final coordinates (with offset applied).

### b) Adding simulated galaxies/stars to an existing image

Instead of creating a new image, the wrapper also allows to add simulated galaxies and stars to an existing image.
The work flow looks something like this:
- First, the image is loaded and the pixel scale and dimension and RA/DEC center is determined (by determining the min/max RA/DEC). For now, it is assumed that the image is of rectangular shape.
- Galaxies and stars are created according to the "world" dictionary. The number is determined by the galaxy/star density and the area of sky covered by the image.
- A simulated image is created with zero noise (this image can be stored by setting "delete_noiseless_image" to FALSE (see below).
- The noiseless image is added to the input image (not that the simulated galaxies/stars will overlap with the real sources on the input image as there is no "avoidance" criteria implemented).

An example is included in the script "input_example.py". Simply run
```
python input_example.py
```
(Note that I am using Python 3.7, I haven't tested the code with other versions of Python.)


The simulated galaxies/stars are added to an existing image using the function
```
simulate_to_existing(world_input,image_input)
```

Again, "image_input" can be a list, for example
```
simulate_to_existing(world_input,image_input=[image_input1,image_input2,image_input3, ...])
```

The "world_input" contains parameters to create the sources (stars and galaxies) while the "image_input" contains parameters of the images.
Note that in the above example, the same galaxies will used to simulate images 1, 2, 3, etc. Therefore, the input images must cover the same region on the sky! If this is not the case, separate functions have to be run, for example:
```
simulate_to_existing(world_input1,image_input=[image_input1])
simulate_to_existing(world_input2,image_input=[image_input2])
...
```

Note that you need to change the "base_name" in the "world" dictionary else the same source catalog is used and the galaxies/stars RA/DEC may not be covered by the image! 


Since images are provided directly, the "world_input" and "image_input" look slightly different:

#### World Dictionary
Contains all the information about the sources. The parameters are:

- base_name: base simulation name (directory with this name will be created inside the "output_directory") [string]
- output_directory: Directory in which Simulations are saved (in sub-folder named according to "base_name") [string]
- overwrite_source_catalog: if TRUE, overwrite source catalog and create new one [True/False]
- source_density: sources per arcmin2 [float]
- radec_distribution_type:"random", # Distribution of galaxies on simulated sky: "random" or "grid" (evenly spaced)
- mag_distribution_type: Type of the magnitude distribution: "uniform" or "halfnormal" [string]
- mag: Magnitude range for objects in AB magnitudes. [max=faintest , half-normal std] for "halfnormal" magnitude distribution and [min=brightest , max=faintest] for "uniform" magnitude distribution [float , float]
- BTR: bulge-to-total ratio range [min=0 , max=1]
- R_disk: disk length-scale range in arcsec [min , max]
- AB_disk: disk aspect ratio range [min=0 , max=1]
- PA_disk: disk position angle (PA) range in degrees [min=0 , max=360]
- R_bulge_rel: bulge length scale range (relative to R_disk)  [min=0 , max=1]
- AB_bulge: bulge aspect ratio range [min=0 , max=1]
- PA_bulge_rel: Bulge PA relative to PA_disk in degrees [-deg,+deg]
- fraction_stars: fraction of sources that are stars (= point source) [min=0,max=1]

See the file "input_example.py" for example values.

#### Images Dictionary
Several image dictionaries can be created and fed to the simulator as a list.
Each dictionary must contain the following parameters:

- image_name: path to the image to which simulated galaxies/stars should be added [string]
- zp: zero point AB of the input image (it is important that this is correct!) [float]
- psf_file_name: File name of the PSF (in FITS format) to be used [string]
- extensions: A list of extensions to which the simulated galaxies/stars are added (for example [0,"IMAGE"] to add them to the primary (0) and "IMAGE" extension. [list]
- delete_noiseless_image: if FALSE, then the noiseless simulated image is kept, if TRUE, it is deleted.
- astro_offset: Astrometry offset in MAS (coords = source_catalog_coords + offset) [ [ra_mean,ra_std],[dec_mean,dec_std]]

See the file "input_example.py" for example values.


### c) Bonus example: How to add galaxies/stars to images that have different PSFs

The script "example_variable_psf.py" contains a small use case example on how to add simulated galaxies/stars to images with different PSFs.
In principle, this is a very simple "for-loop", iterating over the images and feeding the correct PSFs.

In this use case, the setup is the following. We have a list of high-resolution (from HST/ACS) and low-resolution (from Hyper Suprime-Cam) images. The PSF for each low-resolution image is different, however, the PSF of the high-resolution images is the same.
The low-resolution PSFs and images are linked via their file names. So the hardest part is to code the translation from image to PSF file name. For this, we create a PSF template name for the low-resolution image (which is in this case the Hyper Suprime-Cam image for a given tract and patch):

```
lr_psf_name_template = "calexp-HSC-I-%s-%s-%s_*_m21.fits" # tract , patch (format x_y), tract
```

Also, a wild card ( * ) is added and we use the "glob" package to search for the file we need.

We then simply loop over the list of images to which we want to add the galaxies/stars. Note that we need a new "world" dictionary for each of the images as they cover different area on the sky. Also, note that we want to add the simulated images to the primary FITS extension for the high-resolution image and the "IMAGE" extension for the low-resolution image.
Here is the pseudo-code (see "example_variable_psf.py" for the actual code example):

```

hr_image_path = ...
lr_image_path = ...
hr_psf_path = ...
lr_psf_path = ...

lr_psf_name_template = "calexp-HSC-I-%s-%s-%s_*_m21.fits" # tract , patch (format x_y), tract
hr_psf_name = ...

hr_image_list = [...]
lr_image_list = [...]

for hr_image , lr_image in zip(hr_image_list , lr_image_list):

    ## Get name of simulation. This changes for each patch because different location on sky.
    simulation_name = "sim_%s" % lr_image.split(".fits")[0]

    ## World properties
    world_input = {"base_name":simulation_name , # base simulation name (directory with this name will be created)
                    ...
                    }

    # for high-resolution image
    image_input_hr = {"image_name": os.path.join(hr_image_path , hr_image),
                    "zp":25.94734,
                    "psf_file_name": os.path.join(hr_psf_path , hr_psf_name),
                    "extensions":[0],
                    "delete_noiseless_image":True,
                    "astro_offset":[[0,0],[0,0]]
                        }

    # for low-resolution image
    tract = lr_image.split("-")[3]
    patch = lr_image.split("-")[4].split(".fits")[0]
    lr_psf_file_name = glob.glob( os.path.join( lr_psf_path , lr_psf_name_template % (str(tract) , patch , str(tract)) ) )[0]
    image_input_lr = {"image_name": os.path.join(lr_image_path , lr_image),
                    "zp":27.0,
                    "psf_file_name": lr_psf_file_name,
                    "extensions":["IMAGE"],
                    "delete_noiseless_image":True,
                    "astro_offset":[[0,0],[0,0]]
                        }

    # simulate image
    simulate_to_existing(world_input=world_input,
        image_inputs=[image_input_hr,image_input_lr])
    
```


