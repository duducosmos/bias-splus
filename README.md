# bias-splus

Name: bias_splus  
Version 0.2  
Description: Provide statistical analysis for s-plus bias images over periods of time  
Author: Walter Santos  
Created on: Set 12th 2016  
Last Updated: Nov 26th 2016  
Latest Changes:  
- Try/exception clauses to ignore bad bias files
- Added a new input parameter: master bias
- Results include biases division by this master, including thefinal plot
- When no dates are given, only new downaloaded (updated in the server
    are considered for the analysis)
- If a bias for a given night is off by a certain percentage compared to the master,
    a warning is printed out and the point gets redin the plot, instead of the default blue  
    
Instructions: Run the code with a -h/--help option to list all the arguments,
necessary and optional, on the command line.  
Requirements: Numpy, Astropy, Matplotlib

EXAMPLES:

`python bias_splus_v02.py --initdate "2016/09/05" --finaldate "2016/09/11" -o "bias_splus_05_10" -m "splus_master_bias.fits"`  

or for new mode, ony takes into account new bias/dates:  

`python bias_splus_v02.py -o "bias_splus_26_sept" -m "splus_master_bias.fits"`  


NOTE:
For human input, unless stated otherwise, dates should be given as a string by yyyy/mm/dd, 
for example, '2016/09/09'.  
For function input, are given as python date objects
