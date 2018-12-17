# The Matlab OWC toolbox

This is a package of MATLAB routines for calibrating profiling float conductivity sensor drift. A description of the algorithms can be found in "An improved calibration method for the drift of the conductivity sensor on autonomous CTD profiling floats by θ-S climatology", by W.B. Owens and A.P.S. Wong, in Deep-Sea Research Part I: Oceanographic Research Papers, 56(3), 450-457, 2009.
Lately,  modifications suggested in “Improvement of bias detection in Argo float conductivity sensors and its application in the North Atlantic” , by C. Cabanes, V. Thierry and C. Lagadec, in Deep-Sea Research Part I, 114, 128-136, 2016  have been taken into account. 

# How to install the toolbox ?

Either clone the latest version of the git repository:

git clone https://github.com/ArgoDMQC/matlabow.git

or download and unzip the zip file (Clone or download button)

You can also access the different releases here: 
https://github.com/ArgoDMQC/matlabow/releases

# How to run the analysis?
Here is a summary of what should be done to run the analysis, please read the ./doc/README.doc file for more details


1. All files are to be used in MATLAB. The full package was tested with MATLAB R2014a. In addition, you will need:
a). The MATLAB Optimization Toolbox;
b). The ITS-90 version of the CSIRO SEAWATER library. The version 3\_3.0 of this library can be found in ./lib/seawater\_330\_its90. Please update if necessary.
c). The M_MAP toolbox. The version  1.4.c of this library can be found in ./lib/m\_map1.4. Please update if necessary.

2. Add the necessary path to your matlab path: addpath('./lib/seawater\_330\_its90';'./lib/m\_map1.4';'./matlab\_codes/')

3. Put your reference data in ./data/climatology/historical\_ctd, /historical\_bot, /argo\_profiles.

REFERENCE DATA can be obtain at ftp.ifremer.fr
cd /coriolis/data/DMQC-ARGO/   (if you need a login/pswd ask codac@ifremer.fr)

Then, create/update your ./data/constants/wmo\_boxes.mat file (more details in ./doc/README.doc, p3)

4. After you have decided where you want to install the package on your computer, edit ow\_config.txt at the following lines so the correct pathways are specified:

* HISTORICAL\_DIRECTORY =
* FLOAT\_SOURCE\_DIRECTORY =
* FLOAT\_MAPPED\_DIRECTORY =
* FLOAT\_CALIB\_DIRECTORY =
* FLOAT\_PLOTS\_DIRECTORY =
* CONFIG\_DIRECTORY =

5. The last section of ow\_config.txt below the heading "Objective Mapping Parameters" is where you set the various parameters (more details in .doc.README.doc, p4-6)

6.  If this is the first time you are using this system, then the 4 directories /data/float\_source, /float\_mapped, /float\_calib, and /float\_plots should be empty. Decide how you want to organise your floats, e.g. under different project names or different investigator names. Then make identical subdirectories under each of these 4 directories. For example:

/data/float\_source/project\_xx
/data/float\_mapped/project\_xx
/data/float\_calib/project\_xx
/data/float\_plots/project\_xx

7.  Create the float source file (./data/float\_source/project\_xx/$flt\_name$.mat) from the original netcdf files (more details in ./doc/README.doc,p6)

8. Open MATLAB in the top directory. List all the float files in a cell array "float\_names", with the corresponding subdirectories in another cell array "float_dirs". For example,
float_dirs = { 'project\_xx/'; 'project\_xx/'; 'jones/'; 'jones/' };
float_names = { 'float0001'; 'float0002'; 'myfloat\_a'; 'myfloat\_b' }.

Tips: If the files are not saved under a subdirectory and are only saved under ./float\_source/, specify float\_dirs = { ''; ''; ''; '' }, etc.

Run ow\_calibration.m. 


