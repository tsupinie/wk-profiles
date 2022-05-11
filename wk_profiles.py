
import numpy as np

from scipy.optimize import minimize
from scipy.linalg import norm

from metpy.calc import parcel_profile, lcl, cape_cin
from metpy.units import units


def wk_sounding(th_trop=343., t_trop=213., z_trop=12000., th_sfc=300., qv_bl=0.014):
    """ 
    Weisman-Klemp profiles
    
    The defaults for the free parameters are specified the same as in the Weisman and 
    Klemp (1982) paper (with the exception of qv_bl, which isn't specified in the paper, 
    but is in a reasonable range).
    
    Parameters:
    th_trop: Potential temperature at the tropopause (K)
    t_trop: Temperature at the tropopause (K)
    z_trop: Height of the tropopause (m)
    th_sfc: Potential temperature at the surface (K)
    qv_bl: Water vapor mixing ratio in the boundary layer (kg/kg)
    
    Returns:
    A function to which you can pass an array of heights in meters. This function returns
    profiles of temperature (K), dewpoint (K), and pressure (Pa).
    
    Example:
    >>> zs = np.arange(0, 20000, 500)
    >>> prof = wk_sounding()
    >>> temp, dewp, pres = prof(zs)
    """
    c_p = 1004.
    g = 9.806
    kappa = 2. / 7.
    p_0 = 100000.
    
    p_trop = p_0 * (t_trop / th_trop) ** (1 / kappa)
    
    def hypsometric_theta(p1, z1, z2, th1, th2):
        """ 
        Obtained by re-deriving the hypsometric equation, but substituting 
        potential temperature for sensible temperature 
        """
        th_bar = (th1 + th2) / 2
        return (p1 ** kappa - (g * p_0 ** kappa) / (c_p * th_bar) * (z2 - z1)) ** (1 / kappa)
        
    def create_sounding(z):
        # Weisman-Klemp profiles are defined in terms of potential temperature and relative 
        # humidity, which is mildly inconvenient
        theta = np.where(z < z_trop, 
                         th_sfc + (th_trop - th_sfc) * (z / z_trop) ** (5. / 4.),
                         th_trop * np.exp(g / (c_p * t_trop) * (z - z_trop)))
        rh = np.where(z < z_trop, 
                      1 - 3. / 4. * (z / z_trop) ** (5. / 4.),
                      0.25)
        
        pres = np.empty(theta.shape)
        k_trop = np.argmin(np.abs(z_trop - z))
        
        # Integrate pressure upwards and downwards from the tropopause using the theta-modified
        # version of the hypsometric equation
        pres[k_trop] = hypsometric_theta(p_trop, z_trop, z[k_trop], th_trop, theta[k_trop])
        
        for kz in range(k_trop, 0, -1):
            pres[kz - 1] = hypsometric_theta(pres[kz], z[kz], z[kz - 1], theta[kz], theta[kz - 1])
            
        for kz in range(k_trop, len(z) - 1):
            pres[kz + 1] = hypsometric_theta(pres[kz], z[kz], z[kz + 1], theta[kz], theta[kz + 1])
        
        # Now we can get temperature
        temp = theta * (pres / p_0) ** kappa
        
        # Get qv from relative humidity and temperature using Bolton's formula
        vapr_sat = 611.2 * np.exp((17.67 * (temp - 273.15)) / (temp - 273.15 + 243.5))
        qv = rh * 0.622 * vapr_sat / (pres - vapr_sat)
        
        # Enforce the well-mixed qv in the boundary layer
        qv = np.minimum(qv_bl, qv)
        
        # Convert qv to dewpoint
        vapr = (qv * pres) / (0.622 + qv)
        bolton_fac = np.log(vapr / 611.2)
        dwpt = (243.5 * bolton_fac) / (17.67 - bolton_fac) + 273.15
        
        return temp, dwpt, pres
        
    return create_sounding


def compute_opt_params(x):
    hght = np.arange(0, 15000, 50)
    
    th_trop, t_trop, th_sfc, qv_bl, z_trop = x
    wk_prof = wk_sounding(th_trop=th_trop, t_trop=t_trop, z_trop=z_trop, th_sfc=th_sfc, qv_bl=qv_bl)
    temp, dewp, pres = wk_prof(hght)

    temp = temp * units.K
    dewp = dewp * units.K
    pres = pres * units.Pa

    parcel_temp = parcel_profile(pres, temp[0], dewp[0])
    cape, cinh = cape_cin(pres, temp, dewp, parcel_temp)
    lclpres, lcltemp = lcl(pres[0], temp[0], dewp[0])

    lclhght = np.interp(lclpres.m, pres[::-1].m, hght[::-1])
    return temp[0].m, pres[0].m, cape.m, lclhght
    

def find_wk_params(tsfc, psfc, sbcape, sblcl, ztrop): 
    """
    Find parameters for Weisman-Klemp profiles based on some more intuitive values.

    Parameters:
    tsfc:   Requested surface temperature (K)
    psfc:   Requested surface pressure (Pa)
    sbcape: Requested convective available potential energy from the surface parcel (J/kg)
    sblcl:  Requested lifted condensation level from the surface parcel (m)
    ztrop:  Requested height for the tropopause (m)

    Returns:
    A dictionary containing the values for the parameters for the Weisman-Klemp profiles

    Examples
    >>> params = find_wk_params(293., 97000., 1000., 500., 12000.)
    >>> wk_prof = wk_sounding(**params)
    >>> hght = np.arange(0, 20000, 500)
    >>> temp, dewp, pres = wk_prof(hght)
    """   
    def do_optimize(x, ztrop):
        x_tsfc, x_psfc, x_sbcape, x_sblcl = compute_opt_params(np.append(x, ztrop))
        return norm([tsfc - x_tsfc, psfc - x_psfc, sbcape - x_sbcape, sblcl - x_sblcl], ord=2)
    
    initial_guess = np.array([343., 213., 300., 0.014])
    ret = minimize(do_optimize, initial_guess, args=(ztrop,), 
                   method='Nelder-Mead', options={'adaptive': True})

    ret = minimize(do_optimize, ret.x, args=(ztrop,), 
                   method='Nelder-Mead', options={'adaptive': True})
    
    return {param: v for param, v in zip(['th_trop', 't_trop', 'th_sfc', 'qv_bl', 'z_trop'], np.append(ret.x, ztrop))}

if __name__ == "__main__":
    params = find_wk_params(293., 97000., 1000., 500., 12000.)
    print(params)