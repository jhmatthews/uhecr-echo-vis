# uhecr-echo-vis
Visualisation, data and scripts for the Taylor et al. UHECR CoG echoes paper

## Directory Structure

* Movies/ -- all movies of skymaps and frames that create the movies
* Scripts/ -- all scripts for creating figures
* Data/ -- saved data. Does not include the raw data output from the file (email James if you want access to this)
* Figures/ -- figures used in the paper

## Description of Figures 

## Description of movies 



## Making all figures

Making figures from saved data:
* Run ```python Scripts/Make_Figures.py```

Making figures from raw data:
* To remake all figures from the raw data, the raw output directories must be obtained and then placed in the Data/ directory
* The scripts will look for subfolders ModelA, ModelB, ModelC, or appropriate symbolic links
* Run ```python Scripts/Make_Figures.py --raw``` to remake the figures from raw data
* Run ```python Scripts/Healpy_Figures.py``` to remake the skymaps from raw data

## Remaking all movies 
Making movies from png frames:
* make sure you have ffmpeg installed
* Run ```make_movie.sh/make_movie.sh```

Making movies from raw data:
* To remake all figures from the raw data, the raw output directories must be obtained and then placed in the Data/ directory
* The scripts will look for subfolders ModelA, ModelB, ModelC, or appropriate symbolic links
* Run ```python Scripts/Make_Figures.py --raw --movies --no-figures``` to remake the movie frames
* Run ```python Scripts/Healpy_Figures.py --movies --no-figures``` to remake the skymap movie frames

## Contact

James Matthews, james.matthews [at] physics [dot] ac [dot] uk

