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

Simply clone this repository to your computer.


## 2. Usage


### a) Creating a new simulated image

To simulate images, modify and run "input_example.py".
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
- mag: Magnitude range for objects in AB magnitudes [max=faintest , half-normal std]
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

See the file "input_example.py" for example values.


### b) Adding simulated galaxies/stars to an existing image

Instead of creating a new image, the wrapper also allows to add simulated galaxies and stars to an existing image.
The work flow looks something like this:
- First, the image is loaded and the pixel scale and dimension and RA/DEC center is determined (by determining the min/max RA/DEC). For now, it is assumed that the image is of rectangular shape.
- Galaxies and stars are created according to the "world" dictionary. The number is determined by the galaxy/star density and the area of sky covered by the image.
- A simulated image is created with zero noise (this image can be stored by setting "delete_noiseless_image" to FALSE (see below).
- The noiseless image is added to the input image (not that the simulated galaxies/stars will overlap with the real sources on the input image as there is no "avoidance" criteria implemented).

An example input is also included in "input_example.py". Modify and run that script
```
python input_example.py
```
(Note that I am using Python 3.7, I haven't tested the code with other versions of Python.)

The add simulated galaxies/stars to an existing image, run
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
- mag: Magnitude range for objects in AB magnitudes [max=faintest , half-normal std]
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

See the file "input_example.py" for example values.



