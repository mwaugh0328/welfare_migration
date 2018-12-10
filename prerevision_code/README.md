Here is the some run down of this code. I need to fill this out more.

To compute this
```
load('calibration_allmoments_new.mat')
compute_outcomes(new_cal, 0)
```
Then it should compute everything and then spit out the moments. In ``compute_outcomes.m`` you can see each step (i) value function iteration (ii) simulation (iii) implementation of experiment and (iv) collect results.  The ``.mex`` and non-mex files now match up.

The calibration routine is
```
calibrate_wrap_tight
```
And then within it you can see how it works. The key to get this thing to fit was using the ``ga`` solver which is essentially a search of the entire parameter space in a smart way.

**Near-term goal** Use trick employed by Morten and Lyon and Waugh to make the problem smooth, then use a derivative based solver to calibrate the model.
  - This should just be about creating a smooth lottery in the simulation routine between the discrete choices. So rather than a hard go or not, some probability is assigned. This is some logit function over the value functions for each choice, for the appropriate states. 
