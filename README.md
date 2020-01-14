### The Welfare Effects of Encouraging Rural-Urban Migration

---
This repository contains code to reproduce aspects of the paper ["The Welfare Effects of Encouraging Rural-Urban Migration"]((http://www.waugheconomics.com/uploads/2/2/5/6/22563786/LMW.pdf)).

Complete explanations of the repository are currently under construction, but the most basic call to generate outcomes from the model and the partial equilibrium welfare numbers is here:

```
>> load('calibration_final.mat')

>> analyze_outcomes_prefshock(exp(new_val), 1)
```
