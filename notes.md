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




