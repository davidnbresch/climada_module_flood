### Welcome to the *climada additional module*  **flood** 

The flood module encompasses the hazard flood and landslides, that are triggered through heavy precipitation events. For now, we are still working on the flood module. But the good news is, that the landslide part is ready to use.

Use the all-in-one function **climada_ls_hazard_set**

<br>
**Example**

We create a landslide hazard set for [Las Cañas](https://www.google.ch/maps/place/Ilopango,+El+Salvador) neighborhood in the city San Salvador, Central America. The edge coordinates of Las Cañas neighborhood are [-89.145 -89.1 13.692 13.727].
Type **[hazard, centroids, fig] = climada_ls_hazard_set([-89.145 -89.1 13.692 13.727]);** into the command line to create the hazard set.

![Landslide areas for different return periods in Las Cañas](docs/Landslide_hazard_las_canas.png)


<br>
**Get started**

* Get started reading the [landslides manual](/docs/climada_module_landslides.pdf) and the [climada manual](../../../climada/blob/master/docs/climada_manual.pdf).
* Download the climada modules, i.e. [climada core](https://github.com/davidnbresch/climada),  [climada advanced](https://github.com/davidnbresch/climada_advanced) and [flood](https://github.com/davidnbresch/climada_module_flood) into the required folder structure (see figure below).
![climada folder structure](https://raw.githubusercontent.com/wiki/davidnbresch/climada/images/climada_folder_structure.png)

<br>
<br>
<br>
**Maps for Las Cañas, a neighborhood in San Salvador, Central America**

![Digital elevation in Las Cañas](https://raw.githubusercontent.com/davidnbresch/climada_module_flood/master/docs/DEM_las_canas.png)
![Slope in Las Cañas](https://raw.githubusercontent.com/davidnbresch/climada_module_flood/master/docs/Slope_las_canas.png)
![Topographical wetness index in Las Cañas](https://raw.githubusercontent.com/davidnbresch/climada_module_flood/master/docs/TWI_las_canas.png)
![Aspect in Las Cañas](https://raw.githubusercontent.com/davidnbresch/climada_module_flood/master/docs/Aspect_las_canas.png)






