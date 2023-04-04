### +smi_stat/@ChangeDetection

ChangeDetection contains methods for change detection analysis.
This class contains several methods used to detect change points 
from a sequence of intensity data, e.g. intensity trace from 
single-particle tracking.

See reference:
Daniel L. Ensign and Vijay S. Pande, Bayesian (2010), 
Detection of Intensity Changes in Single Molecule and 
Molecular Dynamics Trajectories, J. Phys. Chem. B 2010
https://doi.org/10.1021/jp906786b

---

```
properties:
   % Input variables

   Data;   % Input data given to constructor or set with setData
   LogBayesThreshold; % Input threshold for deciding on signficance level of
                      % Bayes factors
   Nobservations;     % Length of data

   % Output variables

   NchangePoints;   % The number of change points detected
   ChangePoints;    % The discrete times at which a change point was detected
                    % (size=NchangePoints)
   Intensity;       % The mean intensities for each subsequence
                    % (size=NchangePoints+1)
   IntensityModel;  % Intensity of model sequence
   LogBayesFactors; % The Bayes factors for each accepted change points
                    % (size=NchangePoints)
   NrejectedChangePoints;   % The number of rejected change points
   RejectedChangePoints;    % The rejected change points
                            % (size=NrejectedChangePoints)
   RejectedLogBayesFactors; % The Bayes factors for each rejected change points
                            % (size=NrejectedChangePoints)
```
---

Function **estimateChangePoints**(Data,LogBayesFactorThreshold)
is the main recursive method for identifiying change points
since it and other helper methods work on parts of the data set
recursively we make it a static method that takes in only the
portion of the data that it is to operate on

For a detailed example, see
[Example_ChangeDetection.m](../../examples/Example_ChangeDetection.m).
