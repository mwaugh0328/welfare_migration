Here is the some run down of this code. I need to fill this out more.

To compute this
```
load('calibration_allmoments_new.mat')
compute_outcomes(new_cal, 0)
```
Then it should compute everything and then spit out the moments. In ``compute_outcomes.m`` you can see each step (i) value function iteration (ii) simulation (iii) implementation of experiment and (iv) collect results.  

Outstanding issue is that I'm using `.mex` files within it. Its not clear that the non-mex files match up. I have to verify and/or search my computer on this. 
