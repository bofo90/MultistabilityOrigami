**Exploring multistability in prismatic metamaterials through local actuation**

by   
Agustin Iniguez-Rabago (1)  
Yun Li (1)  
Johannes T.B. Overvelde (1)\*
  

Affiliation:  
(1) AMOLF, Science Park 104, 1098 XG Amsterdam, The Netherlands

Corresponding author:  
\*overvelde@amolf.nl

Published in Nature Communications 10, 5577 (2019) https://doi.org/10.1038/s41467-019-13319-7

*If you use this code, please cite the article.*

---

## Purpose

This code scans the energy landscape of prismatic structures in order to find energy minima. 
It creates a prismatic structure based on uniform polyhedra by extruding faces. The structure is made out of faces that can stretch and hinges that can fold. 
The energy of these structures is based on linear and rotational springs at the faces and hinges, respectively. The minimization of such energy is done by using 
the SQP algorithm. The code also generates a 3D plot of these structures. For more detail on the structures, the energy, the algorithm and the results please 
refer to the paper. 

---

## Requirements

Matlab version 2018a with Optimization Toolbox.

---

## Files and directories

This repository has two directories and the file MAIN.m

1. **MAIN.m**: This is the main file of the code. For more information on how to use it see the **Usage** section and for an example see the section **Example**.
2. **Modules**: This directory has all the supporting code for the MAIN.m file. If there is something wrong that you want us to check, please contact the corresponding author. 
3. **hingeList_reduced**: This directory contains all hinge selections of different polyhedra. This can be used to scan through the energy landscapes of these structures in order to find energy minima.

---

## Usage

The **MAIN.m** file has an _option_ variable at the beginning. This variable controls the flow of the program. These are the main options: 

1. **inputType**: This option can be _individual_ or _material_:  
  * _individual_: It creates a single polyhedron that is extruded.  
  * _material_: It creates a polyhedron that is tessellated in space using periodic boundary conditions.
2. **template**: Here the type of polyhedron is defined. These are the options:  
  * _individual_: tetrahedron, triangular prism, cube, octahedron, truncated tetrahedron, pentagonal prism, hexagonal prism, heptagonal prism, octagonal prism, cuboctahedron, nonagonal prism, decagonal prism,
  dodecahedron, dodecagonal prism, truncated cube, truncated octahedron, rhombicuboctahedron, and truncated cuboctahedron.  
  * _material_: triangular prism_mat1, triangular prism_mat2, cubes_mat, hexagonal prism_mat, truncated tetrahedron_mat, octagonal prism_mat, cuboctahedron_mat, dodecagonal prism_mat, truncated cube_mat, truncated octahedron_mat,
  rhombicuboctahedron_mat1, rhombicuboctahedron_mat2, truncated cuboctahedron_mat1, truncated cuboctahedron_mat2, truncated cuboctahedron_mat3.  
3. **analysis**: There are four possible types of analysis:  
  * _info_: It plots the polyhedron and the extrusion. It also plots the numbers of the hinges on them in order for you to see which hinges you want to fold.
  * _result_: It folds the specified hinges in the option **hingeSet** by applying a momentum. Then it releases this momentum and let the structure relax to its local minimum. The code saves the data in a variable called result in the specified folder (option **saveFile**).  
  * _plot_: It reads from the specified file (option **saveFile**) the variable result and plots the deformation in a 3D graph.
  * _savedata_: It reads all results in the specified folder (option **saveFile**), obtains the energy, angles, position distribution and more, and saves this data in .csv files for further analysis.  
4. **readHingeFile**: This option can be either _on_ or _off_. If it is _on_, the code reads the file from the directory **hingeList_reduced** to obtain all possible hinge selections and it tests them. If it is _off_, the option **hingeSet** determines the hinges to fold.
5. **Kappa**: This variable is related to the stiffness of the springs. For more information see the original paper. 
6. **saveFile**: This is the sub-directory name where the files are saved. The simulation's results, plots and data are saved under a directory called **Result**, in a subdirectory with the name of the polyhedra, in this specicified sub-folder.  
7. **hingeSet**: This are the hinge numbers that the simulation folds. You can check which ones they are by selecting the option **analysis** as _info_.

There are more options in the file **initOpt.m**, however they are not needed to be changed. Only if you are curious enough, you can start changing them to see what happens. 

---

## Example

These are the necessary steps to follow a simple simulation and see the results:

1. These are the options to run the energy minimization:   

    ``opt=initOpt('inputType','individual', 'template','truncated tetrahedron', ...``  
    ``            'analysis', 'result', 'readHingeFile', 'off', 'Kappa',0.0001);``  
    `` ``   
    ``opt.saveFile = strcat('/',date,'_Example');``  
    ``opt.hingeSet = [1 2]; %This vector can be changed for the hinge numbers that you want to fold``  
    
2. These are the options to plot the results:   

    ``opt=initOpt('inputType','individual', 'template','truncated tetrahedron', ...``  
    ``            'analysis', 'plot', 'readHingeFile', 'off', 'Kappa',0.0001);``  
    `` ``   
    ``opt.saveFile = strcat('/',date,'_Example');``  
    ``opt.hingeSet = [1 2]; %This vector can be changed for the hinge numbers that you want to fold`` 