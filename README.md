# 3D FDTD EM Simulator
3D Finite-Difference-Time-Domain electromagnetic simulator in python. Cython is used to speed up critical code and MayaVi is used for plotting. This is based on the theory outlined by this set of notes - https://eecs.wsu.edu/~schneidj/ufdtd/ (pages 244-248).

To use, first compile using:
```
python compile.py build_ext --inplace
``` 
Then run using
```
python FTDT.py
```
Example radiation from an infinitesimal dipole source:

![image](https://user-images.githubusercontent.com/15094591/224673047-cfafb9f4-5857-4961-8a87-f215a050f227.png)
