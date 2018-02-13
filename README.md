# PenalizedBTPB
Penalized Bivariate Tensor Product B-Spline Code
This code is for the paper An Investigation of Actuarial Fair Crop Insurance Rates Using Partial Derivatives of Penalized Bivariate Tensor Product B-splines (BTPB).

There are five R scripts in this project and two simulated data frames to show how the scripts and function work.  

The scripts are: 

FunctionsBTPS.R

GraphSetUp.R

SimulationSetUP.R

Workflow_Script.R

SmallWorkingExample.R

The data frames are:

one

two

########FunctionsBTPS.R###########

The FunctionsBTPS.R r script contains five functions.  These functions are used for simulating data and and fitting data using the penalized bivariate tensor product B-spline.  There is also a function for estimating the variance of the partial derivatives of the penalized BTPB, the main focas of the paper.

The functions are: fitBTPS, yieldDataSim, coupute_bandwidth, predictDensity, and VarianceEstimator.

########fitBTPS

fitBTPS is a function used for fitting a penalized bivariate tensor product B-spline. In addition to fitting the B-spline, which uses the function gam in the mgcv package, it returns estimates, first partial derivative wrt the x value estimates, and second partial derivative wrt the x value estimate for a given set of points.  Knots, degree of the function and degree of the penalty can all be adjusted.  the knots represent the number of interior knots. The residuals of the penalized BTPB are also fitted in this function in order to save a step when getting variance estimates in the function VarianceEstimator.  

########yieldDataSim

This function simulates data according to the simulation study specified in chapter 4.  This function simulates coverage rate (x_i), APH (average historical/actual production history) yield (z_i), current yield (w_i), premium values (mu(x,z)), and observed premium values (y_i). There are three main components, the mean structure (mu(z_i)), the error structure (epsilon_i), and the variance structure (sigma) used to estimate the current yield.  The mean structure is either Linear or Quad, Linear: mu(z_i)=-25+1.3z_i, Quad: mu(z_i)=1.2(z_i-50)+(z_i-150)^2/200. The error structure is Normal or Beta. Normal:10*N(0,1), Beta: 50*(Beta(5,3)-5/8). The variance structure is NonConst or Const. NonConst: abs(z_i)^.2 Const: 10 if error structure in Normal or 50 if error structure is Beta.


########compute_bandwidth

The compute_bandwidth function is a helper function for predictDensity. It is used for calculating bandwidths for bivariate kernal density estimation.   

########predictDensity

The predictDensity function is a helper function for VarianceEstimator.  It is used for calculating the bivariate kernel density for f_hat.  

########VarianceEstimator

This function calculates the variance for any partial derivative up to the fourth using the penalized BTPB.  This function uses the object returned by the fitBTPS function.  The equation used is (22) from the paper.

########kernelDensity

########fitKernel
