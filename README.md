Please contact me at singh891@umn.edu / rsthakur922@gmail.com before you use this code. 
# hbond_elliptical

`hbond_elliptical` is a Python library for conducting H-bond analysis for molecular dynamics and Monte Carlo trajectories.

## Installation

`freud` [Link](https://freud.readthedocs.io/en/latest/) and `mdtraj` [Link](https://www.mdtraj.org/1.9.8.dev0/index.html) are required packages for this library. Those packages can be installed using the following commands.

```bash
conda install -c conda-forge freud
conda install -c conda-forge mdtraj
```

## Usage
`hbond_elliptical` expects the trajectories to be the in the `mdtraj` trajectory format, thus requires a `topology` file to load the traj (.xyz, .pdb, etc.) file. 

```python
from hbond_elliptical import calualateHBMap

#Making selections for Oxygens and Hydrogens
sel_oxygen_head = 'name OW' ;
sel_hydrogen = 'name HW1 or name HW2' ;
sel_oxygen_tail = 'name OW';
list_names_hydrogen = ["HW1", "HW2"] ;
list_names_oxygen_head = ["OW"] ;
list_names_oxygen_tail = ["OW"]

# Setting number of bins in r and theta directions
nbins_r = 200 ; nbins_a = 200 ; r_cutoff = 0.4 ; 

# Skipping frames if needed
skip_every_x_frames = 1
bonded_pdb_provided = True

# If your topology file does not contain the bonding information, set 
# bonded_pdb_provided = False and the code will generate bonds for you.
# See hbond_elliptical.create_bond_dict for more information.

rdf_output, inter_output, map_output,hbond,hbond_time = calualateHBMap(traj, r_cutoff, nbins_r, nbins_a, skip_every_x_frames, sel_oxygen_head, sel_oxygen_tail, sel_hydrogen, list_names_hydrogen, list_names_oxygen_head, list_names_oxygen_tail, bonded_pdb_provided)
```
After you run this analysis, you can generate a 2d-histgram of Hbond-freq vs r-theta (see the [jupyter notebook](https://github.com/ramanishsingh/hbond_elliptical_analysis/blob/master/compare_hbond_mdtraj/compare_md_hbondelliptical.ipynb) in the `compare_hbond_mdtraj` folder) and decide if you need to change the widths and center of the ellipse. If so, then you need to modify the values in `hbond_elliptical.ellipticalFun` and redo the analysis.

## Contributing
If you have any questions, or would like to contribute, email me at singh891@umn.edu
## License
[MIT](https://choosealicense.com/licenses/mit/)
