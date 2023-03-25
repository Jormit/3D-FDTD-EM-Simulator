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
### Wave propagation between 2 parallel conductors:
#### Problem Setup
![image](https://user-images.githubusercontent.com/15094591/227667799-9311bc6d-fb7a-45e5-abdc-068fc99223c8.png)
#### Ez Field Magnitude
![image](https://user-images.githubusercontent.com/15094591/227667727-3e63bebb-82cc-4769-83a1-96f6de9db597.png)
### E Vector Field:
![image](https://user-images.githubusercontent.com/15094591/227667782-0787677e-476f-4385-abb4-46337218375c.png)
