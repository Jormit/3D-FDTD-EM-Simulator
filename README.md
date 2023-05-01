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
## Example
A simulation of fields around 2 parallel conductors. Blue element is a sinusoidal voltage source and brown elements are conductors.

![image](https://user-images.githubusercontent.com/15094591/229079181-e978ea6f-7fac-49f4-bbae-f8232f0f7ce0.png)

Magnitude of E-field in the Z direction

![image](https://user-images.githubusercontent.com/15094591/229078004-dd4a748e-760e-468b-a18b-f0635a206349.png)
