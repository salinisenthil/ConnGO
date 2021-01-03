# ConnGO
Workflow for CONNectivity preserving Geometry Optimization
==========================================================

For creating run files, do
> python ConnGO_master.py  filename.smi


For launching run files, do
> python ConnGO_master.py launch


For analysing results, do
> python  ConnGO_master.py  analyse


Files required in the same directory as filename.smi:
=====================================================
ConnGO_master.py,
ConnGO_xyz2sdf.py,
ConnGO_check_tautomers.py,
control.inp,
job.inp

Applications Required:
======================
obabel 2.3.1 (or above),
Gaussian 16,
Numpy

Preprint:
=========
Troubleshooting Unstable Molecules in Chemical Space\
Salini Senthil, Sabyasachi Chakraborty, Raghunathan Ramakrishnan\
https://arxiv.org/abs/2010.02635

@article{senthil2020troubleshooting,\
title={Troubleshooting Unstable Molecules in Chemical Space},\
author={Senthil, Salini and Chakraborty, Sabyasachi and Ramakrishnan, Raghunathan},\
journal={arXiv preprint arXiv:2010.02635},\
year={2020}\
}

Contact:
========
salini.1405@gmail.com, 
ramakrishnan@tifrh.res.in
