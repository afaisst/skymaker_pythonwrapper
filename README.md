# Python Wrapper for SkyMaker

This is a Python wrapper for SkyMaker (https://www.astromatic.net/software/skymaker) to simulate astronomical images.

## 1. Installation

### a) SkyMaker Installation 
You need to have SkyMaker installed on your computer. The latest version can be downloaded from the astromatic website (https://www.astromatic.net/software/skymaker). To install SkyMaker, simply go to the directory and run
> ./configure
> make
> make install

If you have troubles installing SkyMaker, it's likely because of the FFTW dependency that is needed. You can download FFTW from http://www.fftw.org/download.html. Install it with
> ./configure
> make
> make install

After FFTW is installed, try to install SkyMaker again using the above steps. If you get an error like "configure: error: FFTW double precision library files not found in /sw/lib! Exiting.", then you need to do some additional steps.

First install FFTW with threads in double precision mode:
> ./configure --enable-threads
> make
> make install

Second, install FFTW with threads in single precision mode:
> ./configure --enable-single --enable-threads
> make
> install

Then install skymaker again using 
> ./configure
> make
> make all

This should work.
