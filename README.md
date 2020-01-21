# Weighted running mean /spline APWP calculation scripts:
 Matlab code for creating a combined weighted running mean and spline apparent polar wander path

 This repository contains MATLAB scripts to generate apparent polar wander (APW) paths of Gondwana (WeightAPWPgond.m), Laurentia (WeightAPWPlau.m), Baltica (WeightAPWPbalt.m) and Laurussia (WeightAPWPlauru.m).

 Individual paleopoles are listed in an excel file (Paleopoles_Pangea_Formation.xlsx). For instance, Gondwanan poles are shown in Tab "gondwana".

 Results in Table 2 in the manuscript can be produced by running the above scripts which yield output=[1Age 2lonRM 3latRM 4e95a 5e95b 6omega 7Kx 8Ky 9N 10lonSP 11latSP].

 lonRM latRM: location of a running mean pole;
 e95a e95b omega: size and orientation of error ellispse of a running mean pole;
 Kx Ky: precision parameters along the error ellispse axes;
 N: number of input paleopoles used for calculating the running mean pole;
 lonSP latSP: location of a spline pole.

# GPlates reconstruction files of Pangea assembly:
All continents are reconstructed in the paleomagnetic frame (004).

Files included in Recon_Pangea_Formation.zip:
1. Reconstruction file: Pangea_Formation1.rot
2. Coastline file: Pangea_Formation_Coastlines1.shp (.dbf / .shx)
3. Input paleopole files for Gondwana (Plateid=701), Laurentia (Plateid=101), Baltica (Plateid=302) and Laurussia (Plateid=302): 
    ex: poleCGond.shp (.dbf / .shx): Gondwanan paleopole centers
          poleEGond.shp (.dbf / .shx): Gondwanan paleopole error ellipses
4. APWP files for Gondwana (Plateid=701), Laurentia (Plateid=101), Baltica (Plateid=302) and Laurussia (Plateid=302): 
    ex: rmSplpGond.shp (.dbf / .shx): Gondwanan running mean poles
          rmSpllGond.shp (.dbf / .shx): Gondwanan spline poles / curve
      
Three-step to show the mid-to-late Paleozoic continental reconstructions
1. Run GPlates, select 'File' > 'Open Project' > 'Pangea_Formation_Coastlines.prj' (all files will be loaded);
2. Select 'Reconstruction' > 'Specify Anchored Plate ID' > type '4' (paleomagnetic frame);
3. Select 'Reconstruction' > 'Reconstruction to Time' > type an integer between 300 and 460 Ma.

#  Please cite our paper:

 Wu, L., Murphy, J.B., Quesada, C., Li, Z-X., Waldron, J.W.F., Williams, S., Pisarevsky, S., Collins, W.J. The amalgamation of Pangea: Paleomagnetic and geological observations revisited. 
