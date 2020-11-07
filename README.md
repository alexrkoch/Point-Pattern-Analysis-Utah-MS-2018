# Overview
I wrote the MATLAB scripts in this repository to analyze my thesis data
 during my M.S. in Geology at the University of Utah. As such, they are
  extremely tailored to the work I was doing, and to my particular data set
  . They were also my first step into coding, and so are not optimized or
   very well organized.

Much of these scripts were heavily modified from the work of Wassim Benhallam
 who was a previous student in my research group. His work and original
  scripts can be found in his thesis:
  
  _Benhallam, W., 2015, "Spatial analysis of channel-belt
stacking patterns: Metrics to discriminate between local and regional
controls on deposition in the fluvial John Henry Member of the 
Straight Cliffs Formation, southern Utah"_

My thesis work involved spatial point pattern analysis of channel belt features 
from the Cretaceous John Henry member of the Straight Cliffs Formation
, outcropping in souther Utah. 


## Amalgamation.m
Measure amalgamation length and compare it to bottom lengths of 
channel belts over moving window range. To use this, you need two 
sets of polylines; the amalgamation contacts, and the channel belt 
bases. Shapefiles of amalgamation contacts must be meter or cm 
interpolated (see VerticeInterpolator.m)

## CellColorMap.m
This creates cells that can be colored by a map value, but don't need
to occur in a regular grid. The function in colorCell.m is called in this
script.

## colorCell.m
This script is the same as from Cell.m from Benhallam (2015), except
with an added field 'z' which is used to store a property for a cell.
This function is called on in the CellColorMaps.m script.

## KfnEdgeCorrectionAnalyzer.m
to analyze percentage of centroid points that require edge-correction 
in the K-function, by K-radius scale and by moving window scale.
Basically, this will pretend to run the K-function, but won't perform
the actual calculation. Instead it will just count the centroids that
would need correction.

## KfunctionLatWindVE.m
This script performs Ripley's K-function over a lateral moving window 
but with no vertical moving window. It can introduce varying  
anisotropy to the analysis by changing vertical exaggeration of the 
dataset. 

## MovingWindowCorrelation.m
Use this code to statistically compare the results of the Quadrat
analysis to any geologic statistic calculated over the same set of
moving windows (both vertically and laterally). It calculates a
Pearson's R^2 correlation coefficient between two properties for each
moving window.

## NTGMovingWindows.m
This script utilizes shapefiles of channel belt polygons, and
calculates the net-to-gross (NTG), which is net polygon to gross data
area (representing net sand to gross rock area). It calculates NTG
over lateral and vertical moving windows. It also can calculate the 
1D NTG of a vertical line in the center of any given window. This is 
to simulate the NTG that a single well would encounter, and compare
it to the 2D NTG in the surrounding subsurface. 

## QuadratLateralWindow.m
Performs a vertical moving window quadrat point pattern analysis, 
over any set of lateral moving windows.

## Vertice Interpolator.m
Interpolate vertices of a set of polygons so there is one for every 
meter or decimeter of X values. This is necessary in order to use
the width, thickness, NTG, and amalgamation scripts of this thesis.

## WidthThickWindowLimits.m
Calculates the maximum width and thickness of polygons within 
lateral and vertical moving windows.
