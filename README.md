# PVmodel
A Python module that (1) generates the data (three-dimensional axes) of a position-velocity model for a sphere with conical cavities along either the major or minor axis and (2) plots the resulting axes. 

Here are the assumptions we employed:

1. We assume a rotating and infalling motion such that the specific angular momentum is conserved. 
2. We employ small angle approximation and assume that the cavity angle is very large, so that the envelope is almost a spatially thin disk. Thus, the centrifugal barrier is almost constant as a function of the sphere's height. To keep the error within 10%, the cavity angle must be >~ 130 deg.
3. We assume the surface density is proportional to R^p, where p is a free parameter. Default is to use p = -1.5.
4. We take into account the limited angular resolution of the telescope, so the intensity values on both axes of the plane of sky are averaged. 
5. Inside the centrifugal barrier, we assume rigid rotation with angular speed at the radius of the centrifugal barrier in order to avoid dividing by zero. This assumption is adjustable with one of the parameters (see docstring).

Additional notes:
1. The inclination angle is defined such that 0 deg is edge-on, and 90 deg is face-on.
2. The "max_vel" parameter refers to the maximum velocity observed on the PV diagram on the major axis. This is to derive the specific angular momentum and stellar mass parameters.
