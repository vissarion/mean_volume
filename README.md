MeanVol is a software to compute the mean volume of all k-subsets 
of N points in D-dimensional space.

MeanVol is using CGAL and Eigen libraries. All written in C++.

To compile and use MeanVol, you need first to compile the CGAL library, or
download the precompiled library (this software has been tested with CGAL
4.7). You can follow these steps:

------------------------
1. Compile CGAL sources
------------------------

Following the CGAL installation manual, execute:

$ cmake .
$ make


--------------------------
2. Compile MeanVol sources
--------------------------

In folder mean_volume/examples execute:

$ cmake -DCGAL_DIR=_YOUR_CGAL_PATH_ .
$ make

where _YOUR_CGAL_PATH_ is the path where CGAL library was compiled. 

--------------
3. Use MeanVol
--------------

Then, you can test the code with some random iputs:

./test_mean_vol D N k b

where 

  D: dimension of points
  N: the number of points
  k: th size of the subsets
  b: if(!0) test the brute-force algorithm 


