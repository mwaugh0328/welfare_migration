**GE + Taxation Setup**

This is setup so I can talk about results and code.

First, what did I do. In the economy, the government levies a labor income tax. I set it up like this:
```
aftertax = tax.rate.*(laborincome).^(1-tax.prog);
```
This all shows up in [``\revision2_code\calibration\labor_income_tax.m``](https://github.com/mwaugh0328/welfare_migration/blob/master/revision2_code/calibration/labor_income_tax.m). Here ``tax.rate`` is the tax rate and I put in the HSV progressivity parameter ``.^(1-tax.prog)`` but per the discussion below, just set this thing to one.

I explored several things. One was a urban, permanent residents only. This is less distortionary and hence gives the highest welfare gains. The baseline is that it hits everybody. It is probably less realistic, but this is the most distortionary since it will affect both the permanent location decision and the migration decision.

A second one was progressivity. This was less interesting than I thought. I think the key issue here is (i) now the tax system starts to play a public insurance role and confused what's up with migration and (ii) while labor supply in a sense is elastic across space, it's very inelastic within a location. What this last point means is that the gov can really exploit people without much behavioral cost.

The second issue is the government budget constraint. The key here is that we have two periods within a year and one question is...does the government run a balanced budget every periods or over the year. I opted for over the year and set it up like this.
```
params.R.*(accounting.all.monga.tax - accounting.all.monga.fiscalcost) + (accounting.all.notmonga.tax - accounting.all.notmonga.fiscalcost);
```
So the government gets some revenues during the Monga and pays out some moving costs but saves it (hence the R showing up there), then the next period gets some more revenue and pays out more moving costs. So at the end of the fiscal year, this is what it has. We want this to be zero. I played with other versions of this and it did not matter a whole lot. What is not possible is period by period zero. The issue is that the since migration is varying by season, then the government would have to have a seasonal varying tax rate todo this period by period.

---
**Results GE + Taxation**

The main driver file here is [``revision2_code/ge_taxation/tax_eq.m``](https://github.com/mwaugh0328/welfare_migration/blob/master/revision2_code/ge_taxation/tax_eq.m). As described in the file, I set things up a little bit differently. The idea here was to be able to **fix policy functions** and compute welfare and then let actions change.

1. First, let's look at the initial allocation that corresponds with the calibrated economy.

```
Aggregate Statistics
Average Rural Population
    0.5983

Migrants, Control Group, Mushfiqs Sample
    0.3729

Fraction of Rural with Access to Migration Subsity
    0.5012

Experince, Control Group, Mushfiqs Sample
    0.2334

Consumption, Mushfiqs Sample
    0.6284

Control Group, Welfare by Income Quintile: Welfare, Migration Rate, Experience, Consumption
         0   49.3615   22.9460    0.3952
         0   42.1718   22.7238    0.4734
         0   33.7160   21.7077    0.6445
         0   31.8674   21.5502    0.6933
         0   29.3306   27.7570    0.9356

Social Welfare: All, Rural, Urban
     0     0     0
```  

2. Next, let's put in the free moving costs, **but do not let people change their behavior**

```
Average Rural Population
    0.5983

Migrants, Control Group, Mushfiqs Sample
    0.3729

Fraction of Rural with Access to Migration Subsity
    0.5012

Experince, Control Group, Mushfiqs Sample
    0.2334

Consumption, Mushfiqs Sample
    0.6576

Control Group, Welfare by Income Quintile: Welfare, Migration Rate, Experience, Consumption
    2.5395   49.3615   22.9460    0.4340
    2.1170   42.1718   22.7238    0.5067
    1.7746   33.7160   21.7077    0.6707
    1.7091   31.8674   21.5502    0.7181
    1.2541   29.3306   27.7570    0.9583

Social Welfare: All, Rural, Urban
    0.9296    1.4827    0.1057

```
So note how the moving quantities did not change, but welfare went up. So this basically was like relaxing the budget constraint for all the guys and it gets pushed into consumption (note too, assets are fixed here). Obviously everyone gains here, it's just free money. This is kind of the problem with what we had in the paper now, it has to go up when you give out free money, its only a question of how much. With the tax below, it's not obvious because the tax will hurt people too. So it's a tension between the tax and migration subsidy. And if the allocation were efficient, this kind of policy has to hurt. So if it helps, then this is also evidence speaking to the economy not being efficient has misallocation.

3. Now put the tax in and find the tax to the government budget constraint is satisfied, **but do not let people change their behavior**

```
Aggregate Statistics
Average Rural Population
    0.5983

Migrants, Control Group, Mushfiqs Sample
    0.3729

Fraction of Rural with Access to Migration Subsity
    0.5012

Experince, Control Group, Mushfiqs Sample
    0.2334

Consumption, Mushfiqs Sample
    0.6535

Control Group, Welfare by Income Quintile: Welfare, Migration Rate, Experience, Consumption
    2.1191   49.3615   22.9460    0.4319
    1.6983   42.1718   22.7238    0.5041
    1.3561   33.7160   21.7077    0.6668
    1.2905   31.8674   21.5502    0.7138
    0.8354   29.3306   27.7570    0.9511

Social Welfare: All, Rural, Urban
    0.5197    1.0749   -0.3073
```
The tax rate here is ``0.41`` percent. So about 40 cents on 100 dollars.

Again, not how moving quantities did not change. But here welfare is going up, but not as much as before. This is about the tax and it shows up with about everybody's welfare going down by about ``0.40`` percent. But from a social perspective, this improves things, there is still surplus even after taxing things.

4. Now let behavior change. And then clear labor markets and the governments budget constraint.

```
Aggregate Statistics
Average Rural Population
    0.6566

Migrants, Control Group, Mushfiqs Sample
    0.6996

Fraction of Rural with Access to Migration Subsity
    0.7215

Experince, Control Group, Mushfiqs Sample
    0.4178

Consumption, Mushfiqs Sample
    0.7129

Control Group, Welfare by Income Quintile: Welfare, Migration Rate, Experience, Consumption
    2.8954   89.7291   42.1217    0.4538
    2.3294   82.8340   42.3158    0.5261
    1.7096   64.2755   41.0246    0.7374
    1.6620   63.3134   40.8558    0.7802
    1.1382   49.6498   42.5859    1.0669

Social Welfare: All, Rural, Urban
    0.6254    1.5510   -1.1439
```
The tax rate here is ``1.18`` percent. So about 1.20 on 100 dollars.

This is where it gets interesting. Urban guys are hurt by this, again basically by the same amount as the tax. But here, social welfare is actually larger when people can change their actions (0.63 vs 0.52).

What is going on here. I think the key issue is the much of the rural population substitutes into migration in two ways (1) migration rates go up a bunch (2) experience goes up a bunch too (from about 0.23 to 0.41), so they a migrant now receives higher utility from the migration experience than before and (3) notice how the fraction of rural with access to the subsidy increased (50 percent in baseline to 72 percent in GE). So guys are **decumulating** assets to get under the asset threshold and take advantage of the subsidy.

This is not obvious that it would have had to happen, you could have argued that since more people with migrate -> higher fiscal burden -> lower welfare (relative to fixed allocation). Here we are finding the opposite. What is hard to figure out is the role of better insurance (the free migration option) or that the migration experience improved. In a sense, these are fundamentally linked...the more you do it, the better of an option that it becomes.

One interesting point here too is the political economy of this. The urban guys would be against this. Not because of congestion externalities, but that the fiscal burden (while modest) hit them directly. They get no benefit from this policy, just higher taxes. But from a social welfare perspective it's welfare enhancing. 
