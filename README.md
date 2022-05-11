## wk-profiles
Create Weisman-Klemp (WK82) profiles for initializing idealized simulations of severe convective storms.

### Dependencies
* Metpy (https://unidata.github.io/MetPy/latest/index.html)
* Numpy/Scipy

### Usage
There are two functions: `wk_sounding` and `find_wk_parameters`. `wk_sounding` is the main function for generating profiles. It takes the parameters defined in the WK82 paper:
| Parameter | Description                                 | Default value |
| --------- | --------------------------------------------| ------------- |
| th_trop   | Potential temperature at the tropopause (K) | 343           |
| t_trop    | Temperature at the tropopause (K)           | 213           |
| z_trop    | Height of the tropopause (m)                | 12000         |
| th_sfc    | Potential temperature at the surface (K)    | 300           |
| qv_bl     | Mixing ratio in the boundary layer (kg/kg)  | 0.014         |

It returns another function that you can call to generate a profile. For example:
```python
prof = wk_sounding() # Create WK82 profile object with the default parameters
hght = np.arange(0, 20000, 500) # Generate a height array in meters
temp, dewp, pres = prof(hght) # Get profiles of temperature (K), dewpoint (K), and pressure (Pa)
```

These parameters are a bit unintuitive, and it can be difficult to find two profiles that hold a quantity like CAPE constant. So for that, you can use `find_wk_parameters` to find parameters for a profile with specific values for various quantities. The arguments to `find_wk_parameters` are as follows:
| Argument | Description                                      |
| -------- | ------------------------------------------------ |
| tsfc     | Requested surface temperature (K)                |
| psfc     | Requested surface pressure (Pa)                  |
| sbcape   | Requested CAPE from the surface parcel (J/kg)    |
| sblcl    | Requested LCL height from the surface parcel (m) |
| ztrop    | Requested height for the tropopause (m)          |

`find_wk_parameters` returns a dictionary of parameters you can pass straight to `wk_sounding`. For example:

```python
params = find_wk_params(293., 97000., 1000., 500., 12000.) # Find the WK82 profile parameters corresponding to
                                                           # a surface temperature of 293 K, a surface pressure
                                                           # of 970 hPa, an SBCAPE of 1000 J/kg, an SBLCL height 
                                                           # of 500 m, and a tropopause height of 12 km.
prof = wk_sounding(**params) # Create the WK82 profile object with those parameters
hght = np.arange(0, 20000, 500) # Generate a height array in meters
temp, dewp, pres = prof(hght) # Get profiles of temperature (K), dewpoint (K), and pressure (Pa)
```
