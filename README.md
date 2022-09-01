# MLDL_metalographic_nicrostructure_classification

The project aims to test the possibilities of teaching a network to recognize metallographic microstructures generated using random cellular automata grains grow model. Two types of structures were generated using own library for simulation of grains growth. Then the data set was divided into training and testing sets. The network was learned and then the model was tested for few example microstructures.

Example InputData:
<br />
	Desired Microstructures:
  <br />
       ![alt text](https://github.com/kpawliko/MLDL_metalographic_nicrostructure_classification/blob/main/c2.bmp?raw=true)
       ![alt text](https://github.com/kpawliko/MLDL_metalographic_nicrostructure_classification/blob/main/c3.bmp?raw=true)
  <br />   
	Other Microstructures:
  <br />
       ![alt text](https://github.com/kpawliko/MLDL_metalographic_nicrostructure_classification/blob/main/c4.bmp?raw=true)
       ![alt text](https://github.com/kpawliko/MLDL_metalographic_nicrostructure_classification/blob/main/d2.bmp?raw=true)
<br />
Data preprocessing:
Black spots are irrelevant and result from the ways in which the Random Cellular Automata models operate. The colors represent different crystallographic orientations, but in this analysis we are only interested in the shape and size of the grains. For this reason, so that the files do not have unnecessary information, they have been simplified by filling the empty spaces and creating a binary image with marked grain boundaries.
<br />
Preprocessing example:
   <br />
       ![alt text](https://github.com/kpawliko/MLDL_metalographic_nicrostructure_classification/blob/main/c1.bmp?raw=true)
       ![alt text](https://github.com/kpawliko/MLDL_metalographic_nicrostructure_classification/blob/main/gray.bmp?raw=true)
       ![alt text](https://github.com/kpawliko/MLDL_metalographic_nicrostructure_classification/blob/main/dilated.bmp?raw=true)
       ![alt text](https://github.com/kpawliko/MLDL_metalographic_nicrostructure_classification/blob/main/resu.bmp?raw=true)
 
