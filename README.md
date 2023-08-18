# PVmodel
A Python module that (1) generates the data (three-dimensional axes) of a position-velocity model for a sphere with conical cavities along either the major or minor axis and (2) plots the resulting axes. 

Here are the assumptions we employed:
1. We assume a rotating and infalling motion such that the specific angular momentum is conserved. 
2. We employ small angle approximation and assume that the cavity angle is very large, so that the envelope is almost a spatially thin disk. Thus, the centrifugal barrier is almost constant as a function of the sphere's height. To keep the error within 10%, the cavity angle must be >~ 130 deg.
3. We assume the surface density is proportional to R^p, where p is a free parameter. Default is to use p = -1.5.
4. We take into account the limited angular resolution of the telescope, so the intensity values on both axes of the plane of sky are averaged. 
5. Inside the centrifugal barrier, we assume rigid rotation (constant angular speed) with an angular speed at the radius of the centrifugal barrier in order to avoid dividing by zero. This assumption is adjustable with one of the parameters (see docstring).
6. The thickness of the inner disk is equal to the height of the sphere with cavities at the centrifugal barrier. 

Features:
1. The data can be exported by setting save_txt=True and specifying a filename. If the filename already exists, the function does not override the existing file, and an Exception will be raised. You will have to manually delete the file. This is to prevent carelessness as the data generally take a long time to compute.
2. You can retrieve the v-, x-, and s-axes of the PV model from a generated text file using the 'read_txt' function.
3. The module is capable of simulating PV diagrams on both the major and minor axes within the plane of sky. 

Limitation:
1. Takes a long time to compute if linspace_num is high or int_step is low. However, low linspace_num results may differ from high linspace_num results (high linspace_num is more accurate). If you would like to simply see an overall picture, you can set the resol parameter to None, which assumes infinite (non-existent) angular resolution. This improves the performance time significantly, although the results may not be as accurate.

Additional notes:
1. The inclination angle is defined such that 0 deg is edge-on, and 90 deg is face-on.
2. The "max_vel" parameter refers to the maximum velocity observed on the PV diagram "on the major axis". This is to derive the specific angular momentum and stellar mass parameters.
3. This is a modification of the model in Oya et al. (2014, ApJ, 795, 152), which considered the case of an spatially-thin disk.

Derivation:
We consider a spherical symmetric envelope with conical cavities and an embedded circumstellar disk at an inclination angle $i$ relative to the edge-on configuration. We assume that the envelope and inner disk are rotating and infalling. Cartesian coordinates within the envelope-disk system are adopted with the x-axis along the line of sky, the y-axis along the line of sight, and the z-axis orthogonal to both the x- and y-axes. The origin is defined as the center of the system. 

We define the deprojected coordinates of $(X,\ Y,\ Z)$ as $(X',\ Y',\ Z')$. The relationship between the coordinates is given by:
\begin{equation}
X' = X,\ \ Y' = Y \cos i - Z \sin i,\ \ Z' = Y \sin i + Z \cos i
\end{equation}
The distances from (X', Y', Z') to the axis of rotation $R_{ax}$ and to the origin $R_{cen}$ can then be given by:
\begin{equation}
R_{ax} = \sqrt{(X')^{2} + (Y')^{2}},\ \ R_{cen} = \sqrt{(X')^{2}+(Y')^{2}+(Z')^{2}}
\end{equation}

If the angle of cavity $\theta_{cav}$ is very large, the envelope-disk system is almost a spatially-thin disk, specific angular momentum $l$ is roughly conserved, and the centrifugal barrier $R_{0}$ as a function of $Z'$ is almost constant. Hence, the magnitudes of the infall velocity $v_{inf}$ and rotational velocity $v_{rot}$ can be written as:
\begin{equation}
|\vec{v_{inf}}| = \sqrt{\frac{2GM}{R_{cen}} - (\frac{l}{R_{ax}})^{2}}
\end{equation}
and 
\begin{equation}
|\vec{v_{rot}}| = \frac{l}{R_{ax}}
\end{equation}
where $G$ and $M$ stand for the universal gravitational constant and the mass of the protostar (plus the inner disk), respectively. The centrifugal barrier (radius at which all kinetic energy is converted to rotational energy) is thus:
\begin{equation}
R_{0} = \frac{l^2}{2GM}
\end{equation}

Within any circular cross-section area at an arbitrary $Z'$, we define the azimuthal angle $\phi$ and polar angle $\theta$ such that:
\begin{equation}
\sin \phi = \frac{Y'}{R_{ax}}, \ \ 
\cos \phi = \frac{X'}{R_{ax}} 
\end{equation}
and 
\begin{equation}
\sin \theta = \frac{Z'}{R_{cen}}, \ \
\cos \theta = \frac{R_{ax}}{R_{cen}}
\end{equation}

To obtain the radial velocity, we first consider the component provided by the rotational velocity. The rotational velocity can be given by:
\begin{equation}
\vec{v_{rot}} = |\vec{v_{rot}}| \ (\sin \phi \ \hat{\mathbf{x'}} + \cos \phi \ \hat{\mathbf{y'}})
\end{equation}
where $\hat{\mathbf{x'}} = \langle 1,\ 0,\ 0 \rangle$ and $\hat{\mathbf{y'}} = \langle 0,\ \cos i,\ -\sin i \rangle$ are unit vectors along the $x'$- and $y'$-axes, respectively. 

By projecting the rotational velocity onto the y-axis, the component along the line of sight is derived to be:
\begin{equation}
v^{los}_{rot} = |\vec{v_{rot}}| \ \cos \phi \cos i
\end{equation}

Next, we consider the infall velocity, which is given by: 
\begin{equation}
\vec{v_{inf}} = |\vec{v_{inf}}|(\cos \theta \cos \phi \ \hat{\mathbf{x'}} 
+ \cos \theta \sin \phi \ \hat{\mathbf{y'}} - \sin \theta \ \hat{\mathbf{z'}}) 
\end{equation}
where $\hat{\mathbf{z'}} = \langle 0,\ \sin i,\ -\cos i \rangle$. 
Projecting the infall velocity onto the line of sight gives:
\begin{equation}
v^{los}_{inf} = |\vec{v_{inf}}| (\cos \theta \sin \phi \cos i - \sin \theta \sin i)
\end{equation}
The observed intensity $S_v$ is assumed to be proportional to the column density. We suppose that the disk surface density is dependent on the distance to the origin raised to a power law index $p = -1.5$ outside the centrifugal barrier (i.e., $\Sigma \propto R_{cen}^{-1.5}$). Inside the centrifugal barrier, we cap the disk surface density to $R_{0}^{-1.5}$ to avoid dividing by zero. This applies if the following conditions due to the geometric boundaries are satisfied; otherwise, we assume zero surface density:
1. $R_{cen}(X, Y, Z) \leq R_{env}$ and $\theta_{ref}=|\sin^{-1}(|\frac{Z'}{R_{cen}}|)| \leq 90^{\circ} - \frac{\theta_{cav}}{2} $   if $R_{ax} \leq R_0$
2.  $|Z'| \leq R_0 \cot{\frac{\theta_{cav}}{2}}$  if $R_{ax} > R_0$

We calculate the intensity of each pixel on the PV diagram using the following relation:
\begin{equation}
S_{v}(v, x) \propto \int \Sigma(x, y, z) \Delta y
\end{equation}
where $v$ is a function of $(x, y, z)$. 

Note: This is a modification of the model in Oya et al. (2014, ApJ, 795, 152), which considered the case of an optically-thin disk. 

Please feel free to contact me at ethanchou04@gmail.com if you have any questions about the code.
