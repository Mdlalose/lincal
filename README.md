# lincal
lincal.py is a scipt which that got a function that compute best fits gains & true sky using linear redundant calibration alogorith.

lincal.py is test scrpt, to run it , you to specify args as follows : xdim ydim noise_frac_gains noise_frac_sky.
get_chi2_func.py contains chi squared, gradient & curvature functions.
Ant_esp.apy is a script that generate beam errors from gaussain fake data.
get_anta.py is a script that simulate antenna positions given xdim,ydim, latitude and longitude location of an array.
get_vis_sim.py this script simulate visibility given point source catalogue information, antenna positions and beam errors.
All oberve_8x8 is the data simulated from an 8x8 array.
