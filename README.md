# uhecr-echo-vis
Visualisation, data and scripts for the Taylor et al. UHECR CoG echoes paper

## Directory Structure

* Movies/ -- all movies of skymaps and frames that create the movies
* Scripts/ -- all scripts for creating figures
* Data/ -- saved data. Does not include the raw data output from the file (email James if you want access to this)
* Figures/ -- figures used in the paper

## Description of Figures 
Figures from the paper are named according to figure number for figures 1 to 5. 

The skymaps are instead named 
```map_XX_YY_s20.png```
where XX is the model name, and YY is the time stamp in units of 0.49 Myr (the time resolution of the outputs). Composition skymaps follow the same convention but have the prefix light or heavy. 

## Description of movies 
* skymap movies are simply named ```ModelA.mp4```, etc. 
* composition movies have a prefix ```heavy``` or ```light``` 
* summary movies contain positional and local flux animations as well as skymaps and have the prefix ```summary```

## Making all figures

Making figures from saved data:
* Run ```python Scripts/Make_Figures.py```

Making figures from raw data:
* To remake all figures from the raw data, the raw output directories must be obtained and then placed in the Data/ directory
* The scripts will look for subfolders ModelA, ModelB, ModelC, or appropriate symbolic links
* Run ```python Scripts/Make_Figures.py --raw``` to remake the figures from raw data
* Run ```python Scripts/Healpy_Figures.py``` to remake the skymaps from raw data

## Remaking all movies 
Making movies from png frames (recommended):
* make sure you have ffmpeg installed
* Run ```cd Movies/; ./make_movie.sh```

Making movies from raw data:
* To remake all figures from the raw data, the raw output directories must be obtained and then placed in the Data/ directory
* The scripts will look for subfolders ModelA, ModelB, ModelC, or appropriate symbolic links
* Run ```python Scripts/fig3_alltimestamps.py ``` and ```Scripts/fig4_alltimestamps.py``` to remake the movie frames for figs 3 and 4
* Run ```python Scripts/Healpy_Figures.py --movies --no-figures``` to remake the skymap movie frames
* Run ```python Scripts/combine_frames.py``` to remake the summary movie frames
* Run ```cd Movies/; ./make_movie.sh```

## Contact

James Matthews, james.matthews [at] physics [dot] ac [dot] uk

