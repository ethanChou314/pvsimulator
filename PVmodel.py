# We consider a spherical symmetric envelope with conical cavities and an embedded circumstellar disk at an inclination angle $i$ relative to the edge-on configuration. We assume that the envelope and inner disk are rotating and infalling. Cartesian coordinates within the envelope-disk system are adopted with the x-axis along the line of sky, the y-axis along the line of sight, and the z-axis orthogonal to both the x- and y-axes. The origin is defined as the center of the system. 
#
# We define the deprojected coordinates of $(X,\ Y,\ Z)$ as $(X',\ Y',\ Z')$. The relationship between the coordinates is given by:
# \begin{equation}
# X' = X,\ \ Y' = Y \cos i - Z \sin i,\ \ Z' = Y \sin i + Z \cos i
# \end{equation}
# The distances from (X', Y', Z') to the axis of rotation $R_{ax}$ and to the origin $R_{cen}$ can then be given by:
# \begin{equation}
# R_{ax} = \sqrt{(X')^{2} + (Y')^{2}},\ \ R_{cen} = \sqrt{(X')^{2}+(Y')^{2}+(Z')^{2}}
# \end{equation}
#
# If the angle of cavity $\theta_{cav}$ is very large, the envelope-disk system is almost a spatially-thin disk, specific angular momentum $l$ is roughly conserved, and the centrifugal barrier $R_{0}$ as a function of $Z'$ is almost constant. Hence, the magnitudes of the infall velocity $v_{inf}$ and rotational velocity $v_{rot}$ can be written as:
# \begin{equation}
# |\vec{v_{inf}}| = \sqrt{\frac{2GM}{R_{cen}} - (\frac{l}{R_{ax}})^{2}}
# \end{equation}
# and 
# \begin{equation}
# |\vec{v_{rot}}| = \frac{l}{R_{ax}}
# \end{equation}
# where $G$ and $M$ stand for the universal gravitational constant and the mass of the protostar (plus the inner disk), respectively. The centrifugal barrier (radius at which all kinetic energy is converted to rotational energy) is thus:
# \begin{equation}
# R_{0} = \frac{l^2}{2GM}
# \end{equation}
#
# Within any circular cross-section area at an arbitrary $Z'$, we define the azimuthal angle $\phi$ and polar angle $\theta$ such that:
# \begin{equation}
# \sin \phi = \frac{Y'}{R_{ax}}, \ \ 
# \cos \phi = \frac{X'}{R_{ax}} 
# \end{equation}
# and 
# \begin{equation}
# \sin \theta = \frac{Z'}{R_{cen}}, \ \
# \cos \theta = \frac{R_{ax}}{R_{cen}}
# \end{equation}
#
# To obtain the radial velocity, we first consider the component provided by the rotational velocity. The rotational velocity can be given by:
# \begin{equation}
# \vec{v_{rot}} = |\vec{v_{rot}}| \ (\sin \phi \ \hat{\mathbf{x'}} + \cos \phi \ \hat{\mathbf{y'}})
# \end{equation}
# where $\hat{\mathbf{x'}} = \langle 1,\ 0,\ 0 \rangle$ and $\hat{\mathbf{y'}} = \langle 0,\ \cos i,\ -\sin i \rangle$ are unit vectors along the $x'$- and $y'$-axes, respectively. 
#
# By projecting the rotational velocity onto the y-axis, the component along the line of sight is derived to be:
# \begin{equation}
# v^{los}_{rot} = |\vec{v_{rot}}| \ \cos \phi \cos i
# \end{equation}
#
# Next, we consider the infall velocity, which is given by: 
# \begin{equation}
# \vec{v_{inf}} = |\vec{v_{inf}}|(\cos \theta \cos \phi \ \hat{\mathbf{x'}} 
# + \cos \theta \sin \phi \ \hat{\mathbf{y'}} - \sin \theta \ \hat{\mathbf{z'}}) 
# \end{equation}
# where $\hat{\mathbf{z'}} = \langle 0,\ \sin i,\ -\cos i \rangle$. 
# Projecting the infall velocity onto the line of sight gives:
# \begin{equation}
# v^{los}_{inf} = |\vec{v_{inf}}| (\cos \theta \sin \phi \cos i - \sin \theta \sin i)
# \end{equation}
# The observed intensity $S_v$ is assumed to be proportional to the column density. We suppose that the disk surface density is dependent on the distance to the origin raised to a power law index $p = -1.5$ outside the centrifugal barrier (i.e., $\Sigma \propto R_{cen}^{-1.5}$). Inside the centrifugal barrier, we cap the disk surface density to $R_{0}^{-1.5}$ to avoid dividing by zero. This applies if the following conditions due to the geometric boundaries are satisfied; otherwise, we assume zero surface density:
# 1. $R_{cen}(X, Y, Z) \leq R_{env}$ and $\theta_{ref}=|\sin^{-1}(|\frac{Z'}{R_{cen}}|)| \leq 90^{\circ} - \frac{\theta_{cav}}{2} $   if $R_{ax} \leq R_0$
# 2.  $|Z'| \leq R_0 \cot{\frac{\theta_{cav}}{2}}$  if $R_{ax} > R_0$
#
# We calculate the intensity of each pixel on the PV diagram using the following relation:
# \begin{equation}
# S_{v}(v, x) \propto \int \Sigma(x, y, z) \Delta y
# \end{equation}
# where $v$ is a function of $(x, y, z)$. 
#
# Note: This is a modification of the model in Oya et al. (2014, ApJ, 795, 152), which considered the case of an optically-thin disk. 
#
# Please feel free to contact me at ethanchou04@gmail.com if you have any questions about the code.

def pv_sphere_model(xlim=(-2100, 2100, 350),    # line of sky
                    ylim=(-2100, 2100, 5),      # line of sight
                    vlim=(-0.795, 0.805, 0.1),
                    vresol=0.21,
                    resol=None,
                    linspace_num=5,
                    max_vel=1.5,
                    cent_barrier=20,
                    inc_angle=0,
                    cav_angle=180,
                    envelope_radius=2000,
                    int_step=10,
                    rigid_rot=True,
                    axis_unit="AU",
                    distance=403,
                    power_law=-1.5,
                    axis="major",
                    save_txt=False,
                    filename=None):
    """
    To-do: generate the PV diagram data (with intensity) using the spherical 
            symmetric assumption, with cones as upper & lower cavity and a disk within centrifugal barrier.
    Parameters:
        xlim (tuple[float]): the limits of the line of sky position axis [AU] (x_start, x_end, x_step). 
                             We recommend setting x_step as the linear/angular resolution.
                             This is always the position along the line of sky, no matter whether 
                             axis == "major" or "minor" is selected.
        ylim (tuple[float]): the limits of the line of sight axis [AU] (y_start, y_end, y_step).
                             y_step can be any small value.
        vlim (tuple[float]): the limits of the vaxis [AU] (v_start, v_end, v_step). 
                             v_step is the velocity width of the channel map. 
        vresol (float): the velocity resolution of the observation [km/s]
        resol (float): the resolution of the observation in the units of 'axis_unit' parameter. 
                       Default (None) is to consider infinite resolution.
        linspace_num (int): the number by which the other plane-of-sky axis is divided within -resol/2 to resol/2.
                            Decrease this number to get a less accurate result but faster program time.
        max_vel (float): the maximum observed (rotational) velocity [km/s] in the PV diagram.
        cent_barrier (float): the centrifugal barrier [AU].
        inc_angle (float): the inclination angle in deg. 0 is edge-on; 90 is face-on.
        cav_angle (float): the angle of cavity in deg. 180 is spatially thin disk. 0 is a sphere without any cavity.
                            We recommend setting the cavity angle to be greater than 150 deg to 
                            minimize small angle approximation error (<10%).
        envelope_radius (float): the radius of the outer envelope [AU].
        int_step (float): the integration step of xsmooth, a continuous positional axis. 
        rigid_rot (bool): True to use rigid rotation assumption inside centrifugal barrier (to avoid dividing by zero). 
        axis_unit (str): the unit of xlim and ylim. Linear or angular units (e.g. AU, arcsec) are allowed.
        distance (float): the distance of the stellar object from observation in parsec. 
                          This parameter will only be relevant if 'axis_unit' parameter is angular.
        power_law (float): the power law assumption for surface density. Default is R_center ** -1.5. 
        axis (str): the axis ('major' or 'minor') on which the PV model is created
        save_txt (bool): True to save output data as a text file. '.txt' can either be included or excluded. 
        filename (str): the name of the textfile if save_txt is set to True. 
                        An exception will be raised if save_txt is False and filename is not None.
        
    Returns:
        xaxis_ret / zaxis_ret (ndarray): the offset axis (in specified unit) of the PV diagram.
        vaxis (ndarray): the velocities along the line of sight (km/s).
        saxis (ndarray): the relative intensity axis, normalized from 0~1 (no unit).
    """
    
    # modules
    import numpy as np
    from astropy import constants as const, units as u
    import time
    import warnings
    import os
    
    start = time.time()  # start counting program time

    warnings.filterwarnings('ignore')  # ignore warnings (dividing by zero, etc.)
    
    # error checking
    if save_txt and (filename is None or filename == ""):
        raise Exception("Filename paremter cannot be None if save_text is set to True.")
    elif os.path.isfile(filename) or os.path.exists(filename) or os.path.isfile(filename+".txt") or os.path.exists(filename+".txt"):
        raise Exception(f"Filename '{filename}' already exists.")
    else:
        params = {"xlim": xlim,   
                  "ylim": ylim, 
                  "vlim": vlim,
                  "vresol": vresol,
                  "resol": resol,
                  "linspace_num": linspace_num,
                  "max_vel": max_vel,
                  "cent_barrier": cent_barrier,
                  "inc_angle": inc_angle,
                  "cav_angle": cav_angle,
                  "envelope_radius": envelope_radius,
                  "int_step": int_step,
                  "rigid_rot": rigid_rot,
                  "axis_unit": axis_unit,
                  "distance": distance,
                  "power_law": power_law,
                  "axis": axis}
        
    if linspace_num % 2 == 0:
        linspace_num += 1
    if linspace_num >= 11 and resol is not None:
        print("Warning: This program might take very long. You can decrease the 'linspace_num' parameter to shorten the program time, although the results will not be as accurate.")
    
    # physical constants
    G = const.G.cgs

    # assign units
    max_vel *= u.km / u.s
    cent_barrier *= u.AU
    envelope_radius *= u.AU
    inc_angle *= u.deg
    cav_angle *= u.deg
    
    # check error due to small angle approximation
    small_angle_approx_error = (1 - np.cos(90*u.deg-cav_angle/2))/np.cos(90*u.deg-cav_angle/2)
    if small_angle_approx_error > 0.1:
        print(f"WARNING: small angle approx error ~{np.round(100*small_angle_approx_error, 1)}% is greater than 10%")
    else:
        print(f"Note: small angle approx error ~{np.round(100*small_angle_approx_error, 1)}%")
    
    # functions to convert angular units to linear units and vice versa.
    def ang2lin_res(ang_resol, distance=distance*u.pc):
        return (distance*np.tan(ang_resol)).to(u.AU)

    def lin2ang_res(lin_resol, distance=distance*u.pc):
        return np.arctan(lin_resol/distance).to(u.arcsec)    

    # other derived values
    spec_ang_mom = (max_vel*cent_barrier / np.cos(inc_angle)).to(u.AU*u.km/u.s)
    stellar_mass = (spec_ang_mom**2 / 2 / G / cent_barrier).to(u.Msun)
    disk_thickness = (2*cent_barrier/np.tan(cav_angle/2)).to(u.AU)
    max_infall_velocity = (G*stellar_mass/spec_ang_mom).to(u.km/u.s)
    max_rot_velocity = 2 * max_infall_velocity
    
    if save_txt:
        derived_params = {"Specific angular momentum": spec_ang_mom,
                          "Stellar mass": stellar_mass,
                          "Disk thickness": disk_thickness,
                          "Max infall velocity": max_infall_velocity,
                          "Max rot velocity": max_rot_velocity
                          }
    
    # print all important values
    if axis != "Major" and axis != "major" and axis != "Minor" and axis != "minor":
        raise Exception(f"Axis parameter '{axis}' not recognized. Available options are 'Major' and 'Minor'.")
    print()
    print(f"############## Important Information ##############")
    print(f"AXIS: Major") if axis == "Major" or axis == "major" else print(f"AXIS: Minor")
    print()
    print(f"RESOLUTION: Infinite") if resol is None else print(f"RESOLUTION: {resol} {axis_unit} (linspace_num={linspace_num})")
    print()
    print(f"PARAMETERS:")
    print(f"Max observed velocity: {np.round(max_vel, 3)}")
    print(f"Centrifugal barrier: {np.round(cent_barrier, 3)} = {np.round(lin2ang_res(cent_barrier), 3)}")
    print(f"Envelope radius: {np.round(envelope_radius, 3)} = {np.round(lin2ang_res(envelope_radius), 3)}")
    print(f"Inclination angle: {np.round(inc_angle, 3)}")
    print(f"Cavity angle: {np.round(cav_angle, 3)}")
    print()
    print(f"DERIVED VALUES:")
    print(f"Stellar mass: {np.round(stellar_mass, 3)}")
    print(f"Specific angular momentum: {np.round(spec_ang_mom, 3)}")
    print(f"Inner disk thickness: {np.round(disk_thickness, 3)} = {np.round(lin2ang_res(disk_thickness), 3)}")
    print(f"Deprojected maximum infall: {np.round(max_infall_velocity, 3)}")
    print(f"Deprojected maximum rotation: {np.round(max_rot_velocity, 3)}")
    print(f"###################################################")
    print()

    # other functions
    def get_rotated_coordinates(x, y, z):
        x_prime = x  # checked
        y_prime = y*np.cos(inc_angle) - z*np.sin(inc_angle)  
        z_prime = y*np.sin(inc_angle) + z*np.cos(inc_angle)  
        return (x_prime, y_prime, z_prime) 

    def get_raxis(x, y, z):
        x_prime, y_prime, z_prime = get_rotated_coordinates(x, y, z)  
        return np.hypot(x_prime, y_prime)  

    def get_rcenter(x, y, z):
        x_prime, y_prime, z_prime = get_rotated_coordinates(x, y, z)  
        return np.sqrt(x_prime**2 + y_prime**2 + z_prime**2)

    def surf_dens(x, y, z):
        x *= u.AU
        y *= u.AU
        z *= u.AU
        rcenter = get_rcenter(x, y, z)
        raxis = get_raxis(x, y, z)
        return np.where(raxis>cent_barrier, np.power(rcenter, power_law), np.power(cent_barrier, power_law)).value

    def rot_vel(x, y, z):
        x_prime, y_prime, z_prime = get_rotated_coordinates(x, y, z)  
        r_axis = get_raxis(x, y, z)
        cos_phi = x_prime/r_axis
        return (spec_ang_mom / r_axis * cos_phi * np.cos(inc_angle)).to(u.km/u.s).value

    def rigid_rot_vel(x, y, z):
        x_prime, y_prime, z_prime = get_rotated_coordinates(x, y, z)
        raxis = get_raxis(x, y, z)
        cos_phi = x_prime/raxis
        return ((raxis*spec_ang_mom/np.square(cent_barrier)) * cos_phi * np.cos(inc_angle)).to(u.km/u.s).value

    def inf_vel(x, y, z):
        x_prime, y_prime, z_prime = get_rotated_coordinates(x, y, z) # x_prime: float; y_prime, z_prime: array
        r_axis = get_raxis(x, y, z)  
        r_center = get_rcenter(x, y, z)  
        sin_theta = z_prime / r_center  
        cos_theta = r_axis / r_center  
        sin_phi = y_prime / r_axis 
        v_magnitude = np.sqrt(2*G*stellar_mass/r_center - np.square(spec_ang_mom/r_axis))  # array
        return (v_magnitude.to(u.km/u.s) * (cos_theta*sin_phi*np.cos(inc_angle) - sin_theta*np.sin(inc_angle))).value

    def los_vel(x, y, z):
        x *= u.AU
        y *= u.AU
        z *= u.AU

        xprime, yprime, zprime = get_rotated_coordinates(x, y, z)  
        raxis = get_raxis(x, y, z) 
        rcenter = get_rcenter(x, y, z)  

        # check for condition inside the centrifugal barrier
        ang_cond = (np.abs(np.arcsin(np.abs(zprime/rcenter)).to(u.deg)) <= 90*u.deg - cav_angle/2)
        vel = inf_vel(x, y, z) + rot_vel(x, y, z)
        vel = np.where(rcenter<=envelope_radius, vel, np.nan)
        if rigid_rot:
            vel = np.where(raxis>=cent_barrier, vel, rigid_rot_vel(x, y, z))
        vel = np.where(raxis<cent_barrier, np.where(np.abs(zprime)<disk_thickness/2, vel, np.nan), 
                       np.where(ang_cond, vel, np.nan))
        return vel

    def int_surf_dens(xvalue, yarr, zvalue, v1_val, v2_val):
        mask = (v1_val <= los_vel(xvalue, yarr, zvalue)) & (los_vel(xvalue, yarr, zvalue) < v2_val) & (los_vel(xvalue, yarr, zvalue) != np.nan)
        surface_densities = np.where(mask, surf_dens(xvalue, yarr, zvalue), 0)
        return np.nansum(surface_densities)
    
    # y-axis:
    yaxis = np.arange(ylim[0], ylim[1]+ylim[2], ylim[2]) 
    
    # x and z axes:
    if axis == "major" or axis == "Major":
        x1 = np.arange(xlim[0]-xlim[2], xlim[1], xlim[2])
        x2 = x1 + xlim[2]*2
        xaxis = (x1 + x2)/2
        ret_xaxis = xaxis[:]  # keep original unit to be returned
        xsmooth = np.arange(xaxis[0], xaxis[-1]+int_step, int_step)
        zaxis = np.linspace(-resol/2, resol/2, linspace_num) if resol is not None else np.array([0])
        if np.round(xaxis[-1], 5) != xlim[1]:  # check for OBOB
            print("WARNING: off-by-one bug detected. We recommend setting x-axis upper & lower limits as multiples of step size.")
        # convert positional axes to AU
        if axis_unit != "AU":
            try:
                x1 = eval(f"(x1*u.{axis_unit}).to(u.AU).value")
                x2 = eval(f"(x2*u.{axis_unit}).to(u.AU).value")
                xaxis = eval(f"(xaxis*u.{axis_unit}).to(u.AU).value")
                xsmooth = eval(f"(xsmooth*u.{axis_unit}).to(u.AU).value")
                yaxis = eval(f"(yaxis*u.{axis_unit}).to(u.AU).value")
                zaxis = eval(f"(zaxis*u.{axis_unit}).to(u.AU).value")
            except:
                try:
                    x1 = eval(f"ang2lin_res(x1*u.{axis_unit}).value")
                    x2 = eval(f"ang2lin_res(x2*u.{axis_unit}).value")
                    xaxis = eval(f"ang2lin_res(xaxis*u.{axis_unit}).value")
                    xsmooth = eval(f"ang2lin_res(xsmooth*u.{axis_unit}).value")
                    yaxis = eval(f"ang2lin_res(yaxis*u.{axis_unit}).value")
                    zaxis = eval(f"ang2lin_res(zaxis*u.{axis_unit}).value")
                except:
                    raise Exception(f"Axis unit '{axis_unit}' is not recognized. Try angular or linear unit and check for whitespaces.")
    elif axis == "minor" or axis == "Minor":
        z1 = np.arange(xlim[0]-xlim[2], xlim[1], xlim[2])
        z2 = z1 + xlim[2]*2
        zaxis = (z1 + z2)/2
        ret_xaxis = zaxis[:]  # keep original unit to be returned
        zsmooth = np.arange(zaxis[0], zaxis[-1]+int_step, int_step)
        xaxis = np.linspace(-resol/2, resol/2, linspace_num) if resol is not None else np.array([0])
        if np.round(zaxis[-1], 5) != xlim[1]:  # check for OBOB
            print("WARNING: off-by-one bug detected. We recommend setting x-axis upper & lower limits as multiples of step size.")
        # convert positional axes to AU
        if axis_unit != "AU":
            try:
                z1 = eval(f"(z1*u.{axis_unit}).to(u.AU).value")
                z2 = eval(f"(z2*u.{axis_unit}).to(u.AU).value")
                zaxis = eval(f"(zaxis*u.{axis_unit}).to(u.AU).value")
                zsmooth = eval(f"(zsmooth*u.{axis_unit}).to(u.AU).value")
                yaxis = eval(f"(yaxis*u.{axis_unit}).to(u.AU).value")
                xaxis = eval(f"(xaxis*u.{axis_unit}).to(u.AU).value")
            except:
                try:
                    z1 = eval(f"ang2lin_res(z1*u.{axis_unit}).value")
                    z2 = eval(f"ang2lin_res(z2*u.{axis_unit}).value")
                    zaxis = eval(f"ang2lin_res(zaxis*u.{axis_unit}).value")
                    zsmooth = eval(f"ang2lin_res(zsmooth*u.{axis_unit}).value")
                    yaxis = eval(f"ang2lin_res(yaxis*u.{axis_unit}).value")
                    xaxis = eval(f"ang2lin_res(xaxis*u.{axis_unit}).value")
                except:
                    raise Exception(f"Axis unit '{axis_unit}' is not recognized. Try angular or linear unit and check for whitespaces.")

    # v-axis
    v1 = np.arange(vlim[0]-vresol/2, vlim[1]-vresol/2+vlim[2], vlim[2])
    v2 = v1 + vresol
    vaxis = (v1 + v2)/2
    if np.round(vaxis[-1], 5) != vlim[1]:
        print("WARNING: off-by-one bug detected. We recommend setting v-axis upper & lower limits as multiples of step size. Ignore this message if not applicable.")

    # s-axis        
    saxes = []
    print_every = 10  # percent
    dot_print_every = print_every * 0.25  # Print a dot every 25% of print_every (3 dots in between)
    next_print = progress = next_dot_print = 0
    total_progress = zaxis.size * vaxis.size * xaxis.size
    print("Progress: ", end="")
    
    if axis == "Major" or axis == "major":
        for z in zaxis:
            svals = np.zeros((vaxis.size, xaxis.size))  # Initialize saxis
            for i in range(vaxis.size):
                for j in range(xaxis.size):
                    x_int = xsmooth[(x1[j] <= xsmooth) & (xsmooth < x2[j])]
                    svals[i][j] = np.mean([int_surf_dens(x, yaxis, z, v1[i], v2[i]) for x in x_int])
                    # Progress bar:
                    if (100 * progress) // total_progress >= next_print:
                        print(f"{next_print}%", end="")
                        next_print += print_every
                        next_dot_print += dot_print_every 
                    elif (100 * progress) // total_progress >= next_dot_print:
                        print(".", end="")
                        next_dot_print += dot_print_every
                    progress += 1
            saxes.append(svals)
    elif axis == "Minor" or axis == "minor":
        for x in xaxis:
            svals = np.zeros((vaxis.size, zaxis.size))  # Initialize saxis
            for i in range(vaxis.size):
                for j in range(zaxis.size):
                    z_int = zsmooth[(z1[j] <= zsmooth) & (zsmooth < z2[j])]
                    svals[i][j] = np.mean([int_surf_dens(x, yaxis, z, v1[i], v2[i]) for z in z_int])
                    # Progress bar:
                    current_progress_percent = (100 * progress) // total_progress
                    if current_progress_percent >= next_print:
                        print(f"{next_print}%", end="")
                        next_print += print_every
                        next_dot_print += dot_print_every 
                    elif current_progress_percent >= next_dot_print:
                        print(".", end="")
                        next_dot_print += dot_print_every
                    progress += 1
            saxes.append(svals)
        
    saxis = np.mean(saxes, axis=0)
    saxis = (saxis - np.nanmin(saxis)) / (np.nanmax(saxis) - np.nanmin(saxis))  # normalize
    print("100%")

    # save results to a textfile
    if save_txt:
        filename += ".txt" if not filename.endswith(".txt") else ""
        print()
        print(f"Saving results to '{filename}'...")
        lines = ["params = {"] + [f"'{key}': '{value}'," for key, value in params.items()] + ["}"] + [" "]
        lines += ["derived_params = {"] + [f"'{key}': '{value},'" for key, value in derived_params.items()] + ["}"] + [" "]
        lines += [f"saxis = np.array({str(saxis.tolist())})"] + [" "]
        lines += [f"vaxis = np.array({str(vaxis.tolist())})"] + [" "]
        lines += [f"xaxis = np.array({str(ret_xaxis.tolist())})"]
        with open(filename, "w") as f:
            f.write("\n".join(lines))
    
    # end timer and print program time
    print()
    print(f"Program time: {time.strftime('%M:%S', time.gmtime(np.round(time.time()-start)))}")

    return vaxis, ret_xaxis, saxis


def plt_pv_data(vaxis, 
                xaxis, 
                saxis, 
                title="PV Model",
                cbaron=True, 
                cmap="inferno",
                plt_type="color",
                datatable=True,
                use_offset_as_x=False,
                figsize=(6, 4), 
                figlim=None,
                ref_lines=True,
                interpolate=False,
                clevels=None,
                offset_unit="AU",
                srange=None,
                transpose=False,
                flip=False,
                line_color="w",
                ):
    """
    To-do: plot the PV diagram given vaxis, xaxis, and saxis.
    Parameters:
        vaxis (1d array): the velocity axis of the PV diagram [km/s].
        xaxis (1d array): the offset axis of the PV diagram in the specified unit.
        saxis (2d array): the (relative) intensity axis.
        title (str): the title of the plot. 
        cbaron (bool): if True, plot the colorbar on the side.
        cmap (str): the color map of the plot.
        plt_type (str): the type of plot. Available options are "color", "contour," and "contourf".
        datatable (bool): if True, the function will plot the data in a datatable.
        use_offset_as_x (bool): if True, the offset axis label will be plotted on the x axis of the PV diagram.
        figsize (array-like with shape (2,)): the size of the figure. Default is to use (6, 4).
        figlim (array-like with shape (4,)): the limits of the figure given in (xmin, xmax, ymin, ymax). 
                                             Default is to use min & max of given axes.
        ref_lines (bool): True to plot reference lines at x = 0 and y = 0.
        interpolate (bool): True to interpolate the data. Only applicable for plt_type == "contour" or "contourf".
        clevels (1d array-like): list of contour levels to be plotted if plt_type == "contour" or "contourf".
        offset_unit (str): unit of offset axis.
        srange (array-like): the range of the intensity axis. Default to use min & max of given saxis.
        transpose (bool): if True, the saxis will be transposed before plotting.
        flip (bool): if True, the saxis will be flipped across x=0 before plotting.
        line_color (str): the color of the reference lines, if ref_lines == True.
    """
    # modules
    import numpy as np
    import pandas as pd
    from matplotlib import pyplot as plt
    from matplotlib import rcParams
    
    saxis = saxis.T if transpose else saxis
    
    # Plot results
    letterratio = 1.294
    nrows = ncols = 1
    width = height = 250
    fig_width_pt  = width*ncols
    fig_height_pt = height*nrows
    inches_per_pt = 1.0/72.27  # Convert pt to inch
    fig_width     = fig_width_pt*inches_per_pt  # width in inches
    fig_height    = fig_height_pt*inches_per_pt # height in inches
    fig_size      = [fig_width, fig_height]
    params = {'axes.labelsize': 12,
              'axes.titlesize': 12,
              'font.size': 12,
              'legend.fontsize': 12,
              'xtick.labelsize': 12,
              'ytick.labelsize': 12,
              'xtick.top': True,   # draw ticks on the top side
              'xtick.major.top': True,
              'figure.figsize': fig_size,
              'figure.dpi': 600,
              'font.family': 'Times New Roman',
              "mathtext.fontset": 'stix', #"Times New Roman"
              'mathtext.tt': 'Times New Roman',
              'axes.linewidth': 2.5,
              'xtick.major.width': 1.0,
              'ytick.major.width': 1.0,
              'xtick.minor.width': 1.0,
              'ytick.minor.width': 1.0,
               'xtick.major.size': 6,
               'ytick.major.size': 6,
                'xtick.minor.size': 4.,
                'ytick.minor.size': 4.,
    }
    rcParams.update(params)
    
    if not use_offset_as_x:
        saxis = saxis.T
    if flip:
        saxis = np.flipud(saxis)
        
    if interpolate and plt_type != "color":
        if use_offset_as_x:
            x, v = np.meshgrid(xaxis[:], vaxis[:])
            points = np.vstack((x.flatten(), v.flatten())).T
            values = saxis[:].flatten()
            grid_x, grid_v = np.mgrid[xaxis[:].min():xaxis[:].max():100j, vaxis[:].min():vaxis[:].max():100j]
            grid_s = griddata(points, values, (grid_x, grid_v), method='cubic')
            xaxis, vaxis, saxis = grid_x[:], grid_v[:], grid_s[:]
        else:
            v, x = np.meshgrid(vaxis[:], xaxis[:])
            points = np.vstack((v.flatten(), x.flatten())).T
            values = saxis[:].flatten()
            grid_v, grid_x = np.mgrid[vaxis[:].min():vaxis[:].max():100j, xaxis[:].min():xaxis[:].max():100j]
            grid_s = griddata(points, values, (grid_v, grid_x), method='cubic')
            vaxis, xaxis, saxis = grid_v[:], grid_x[:], grid_s[:]

    plt.figure(figsize=(6, 4))

    plt.tick_params(which='both',direction='in',bottom=True, top=True, left=True, right=True,
            colors='k',labelrotation=0, labelcolor='k')
    plt.title(title)

    if use_offset_as_x:
        if plt_type == "contourf":
            plt.contourf(xaxis, vaxis, saxis, cmap=cmap, levels=clevels)
        elif plt_type == "contour":
            saxis = saxis.T if interpolate else saxis
            if clevels is None:
                clevels = np.arange(0.2, 1., 0.2)*np.nanmax(saxis)
                print(f"Contour levels: {list(clevels)}")
            plt.contour(saxis, colors="k", origin='lower', extent=[xaxis[0], xaxis[-1], vaxis[0], vaxis[-1]],
                        levels=clevels, linewidths=1.25, alpha=None)
        elif plt_type == "color" or plt_type == "colour":
            plt.imshow(np.flipud(saxis), origin='lower', extent=[xaxis[0], xaxis[-1], vaxis[0], vaxis[-1]], cmap=cmap, aspect='auto')
        plt.xlabel(f'Offset ({offset_unit})') 
        plt.ylabel('Velocity (km/s)')
    else:  
        if plt_type == "contourf":
            plt.contourf(vaxis, xaxis, saxis, cmap=cmap, levels=clevels)
        elif plt_type == "contour":
            saxis = saxis.T if interpolate else saxis
            if clevels is None: 
                clevels = np.arange(0.2, 1., 0.2)*np.nanmax(saxis)
                print(f"Contour levels: {list(clevels)}")
            plt.contour(saxis, colors="k", origin='lower', extent=[vaxis[0], vaxis[-1], xaxis[0], xaxis[-1]],
                        levels=clevels, linewidths=1.25, alpha=None)
        elif plt_type == "color" or plt_type == "colour":
            plt.imshow(np.flipud(saxis), origin='lower', extent=[vaxis[0], vaxis[-1], xaxis[0], xaxis[-1]], cmap=cmap, aspect='auto')
        plt.xlabel('Velocity (km/s)')
        plt.ylabel(f'Offset ({offset_unit})')

    if cbaron: 
        plt.colorbar(label='Relative Intensity')

    if ref_lines:
        plt.axhline(y = 0, color = line_color, linestyle = '-', linewidth=1)
        plt.axvline(x = 0, color = line_color, linestyle = '-', linewidth=1)

    if figlim is not None:
        if figlim[0] > figlim[1]: 
            figlim[0], figlim[1] = figlim[1], figlim[0]
        if figlim[2] > figlim[3]:
            figlim[2], figlim[3] = figlim[3], figlim[2]
        plt.xlim([figlim[0], figlim[1]])
        plt.ylim([figlim[2], figlim[3]])
    else:
        if use_offset_as_x:
            plt.xlim([xaxis[0], xaxis[-1]])
            plt.ylim([vaxis[0], vaxis[-1]])
        else:
            plt.ylim([xaxis[0], xaxis[-1]])
            plt.xlim([vaxis[0], vaxis[-1]])

    if srange is not None: 
        if srange[0] > srange[1]:
            srange[0], srange[1] = srange[1], srange[0]
        plt.clim(srange[0], srange[1])
    
    if datatable:
        try:
            from IPython.display import display
        except:
            print("Your environment does not support 'IPython.display'. Datatable cannot be shown.")
        finally:
            pd.set_option('display.precision', 3)

            def highlight_cells(val):
                color = 'yellow' if val > 0 else ''
                return 'background-color: {}'.format(color)

            sdt = pd.DataFrame(np.round(saxis, 3))
            if use_offset_as_x:
                sdt.columns = ["%.1f" % x for x in xaxis]
                sdt.index = ["%3f" % v for v in vaxis]
            else:
                sdt.columns = ["%.3f" % v for v in vaxis]
                sdt.index = ["%.1f" % x for x in xaxis]
            sdt = sdt.style.format('{:.3f}', na_rep="-").applymap(highlight_cells)
            display(sdt)
                
    plt.show()
