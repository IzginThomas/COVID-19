# Numerical Experiments

This directory contains code to reproduce the numerical experiments described
in the manuscript.

This code is developed with Matlab version R2025b. To reproduce the
results, start Matlab in this directory and execute the following commands in
the Matlab command window to create the figures and tables shown in the paper.

The reference solutions are produced using 

```bash
ODESolverTestSuite();
```
where the respective test case is set in line 66, e.g.

```bash
test ='npzd';
```
To create the tables in the Appendix use

```bash
WriteTables();
```
which will generate latex files in the parent folder. Use
```bash
read_standard_costs();
```
for Table 1. Table 2 is made by hand.

Finally, to produce the Work-Precision diagrams use
```bash
WorkPrecison();
```
There, choose the method in line 24, e.g. with
```bash
solverstr = 'ROS2';
```
In addition, the chosen parameters of the controller can be adjusted in lines 48-50. The Work-Precision diagrams for the comparison of the solvers with the built-in Matlab schemes will be also done, if
```bash
Matlab_comp_bool = 1;
```
in line 4. In this case, choose the respective test case in line 152. E.g., for the Porous Medium Equation, write
```bash
str = {'PME'};
```

Finally, if you wish to run the Basian Optimization again, use

```bash
Bayes();
```



