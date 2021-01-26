# ConnGO
Workflow for CONNectivity preserving Geometry Optimization

For creating run files, do
> python ConnGO_master.py  filename.smi


For launching run files, do
> python ConnGO_master.py launch


For analysing results, do
> python  ConnGO_master.py  analyse


Files required in the same directory as filename.smi:
ConnGO_master.py,
ConnGO_xyz2sdf.py,
ConnGO_check_tautomers.py,
control.inp,
job.inp

Applications Required:
obabel 2.3.2 (or above), 
Gaussian 16,
Numpy

For comments, please contact
salini.1405@gmail.com, 
ramakrishnan@tifrh.res.in
