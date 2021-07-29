# ConnGO

<a>
<img src="https://github.com/raghurama123/ConnGO/blob/master/conngo.gif"  height="350">
</a>

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
=======================
obabel 2.3.2(or above),
Gaussian 16,
Numpy

Citation:
=========
_Troubleshooting Unstable Molecules in Chemical Space_   
Salini Senthil, Sabyasachi Chakraborty, Raghunathan Ramakrishnan    
Chemical Science 12 (2021) 5566.    
DOI: https://doi.org/10.1039/D0SC05591C  

```
@article{senthil2021troubleshooting,
  title={Troubleshooting unstable molecules in chemical space},
  author={Senthil, Salini and Chakraborty, Sabyasachi and Ramakrishnan, Raghunathan},
  journal={Chemical Science},
  volume={12},
  number={15},
  pages={5566--5573},
  year={2021},
  publisher={Royal Society of Chemistry}
}
```


Data:
=====
curatedQM9 dataset can be downloaded from [https://github.com/moldis-group/curatedQM9](https://github.com/moldis-group/curatedQM9)
  
Contact:
========
salini.1405@gmail.com,  
ramakrishnan@tifrh.res.in


</div>
