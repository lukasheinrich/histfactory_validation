# histfactory_validation

The HistFactory p.d.f. template is per-se independent of its implementation in ROOT and sometimes, it's useful to be able to run statistical analysis outside
of ROOT, RooFit, RooStats framework.

This repo has some example code for multi-bin histogram-based analysis based on the asymptotic formulas of arxiv:1007.1727 and  

So far it only implements a simple model of

* one signal historgam with a NormFactor
* one background histogram with a ShapeSys (uncorrelated bin-by-bin uncertainties)

and validates it against output computed from HistFactory workspaces


<img src="img/manual_1bin_55_50_7.png" alt="manual" style="width: 200px;"/>
<img src="img/hfh_1bin_55_50_7.png" alt="manual" style="width: 200px;"/>
