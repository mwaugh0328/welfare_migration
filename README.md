### [The Welfare Effects of Encouraging Rural-Urban Migration](http://www.waugheconomics.com/uploads/2/2/5/6/22563786/LMW.pdf)

---
This repository contains code to reproduce aspects of the paper ["The Welfare Effects of Encouraging Rural-Urban Migration"](http://www.waugheconomics.com/uploads/2/2/5/6/22563786/LMW.pdf). It also includes replication files (empirical and quantitative results) for the paper ["Underinvestment in a Profitable
Technology: The Case of Seasonal Migration in Bangladesh"](https://onlinelibrary.wiley.com/doi/abs/10.3982/ECTA10489)

**Complete explanations of the repository are currently under construction.**

##### Basic Calls
The most basic call to start inside the ``\revision_code\calibration`` folder. From there to generate outcomes from the model and the partial equilibrium welfare numbers is here:

```
>> load('calibration_final.mat')

>> analyze_outcomes_prefshock(exp(new_val), 1)
```
Then it should compute everything and then spit out the moments. In ``analyze_outcomes_prefshock.m`` you can see each step (i) value function iteration (ii) simulation (iii) implementation of experiment and (iv) collect results. The results should mimic (or come very close) to those in Table 2, 6, and 8 in the January 2020 version of the paper.

The calibration routine is implemented by
```
>> calibrate_wrap_tight
```
And then within it you can see how it works. It calls ``compute_outcomes_prefshock.m`` which is similar to the ``analyze...`` file above but is optimized for the calibration routine.  The key to get this thing to fit was using the ``ga`` solver which is essentially a search of the entire parameter space in a smart way. Alternative approaches with different minimizers are in the ``graveyard`` folder.
