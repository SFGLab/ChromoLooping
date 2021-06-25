# ChromoLooping

Code and data in this repository comes from publication _Super-resolution visualization of chromatin loop folding in human lymphoblastoid cells using interferometric photoactivated localization microscopy 

You can use UCSF Chimera or any other viewer to visualize models and points in .pdb format.
 ## Getting Started
 
 Requirements for Python packages can be found in **requirements.txt**
 
 Requirements:
 * Python 3.6 or later
 * Jupyter Notebook
 * UCSF Chimera or other visualization software
 
 
 ## Jupyter Notebooks
 All examples of usage can be found in Jupyter notebooks. In notebook **image_statistics.ipynb** you can find 
 all the operations on ASCII files with image peaks exported from PeakSelector:
 * statistics
 * plotting
 * exporting image peaks to .pdb format
 
 In **image_driven_modeling.ipynb** we put whole modeling procedure:
 * Traveling Salesman Problem solver
 * Spline interpolation
 * Calculating distance maps from 3D models
 * Distance maps visualization 
 
 ## Command line usage
