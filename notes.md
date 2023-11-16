# Notes on code

- There is a Calib struct variable that eventually get's saved as .mat file

- Calib file get's updated in the code, which introduces a lot of redundancy when assigning variables

- Code is very bloated, hard to follow large codes when there is not clearly defined functions that state output and inputs. 
- Some modularization is needed to make the code understandable.

## calibration.ini

- The edge finding section is a mess, needs to be modularized in function

- Edge calculation done on I_parallel_nobias_field, which is assigned I0

- The images are 696 x 520 pixels.         
N0 = dimension(1); 
L = dimension(2);

- One possible reason for inconsistent results is that ROI for calibration is hand-selected, which is open to interpretation by the user...


## Calculate_Efield.m

- Seems like a very important function that does performs distortion correction and Efield calculation

- How is Efield calculated in general? Efield profile is proportional to the intensity of each pixel
The integrated Efield is in units of V.
Efield is usually plot in 1-D, so the pixels must be summed or averaged across the line profile

- There is a difference between Intensity of image and E-field profile, the naming of variables have it confused sometimes

- Try and understand Func_distortion_correction.m, it does the main edge correction, but not sure if it's correct

- Calculate_Efield ...

- Why not apply distortion_correction on all the images rather than the Efield ???
