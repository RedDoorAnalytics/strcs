
Name: Hannah Bower
Date: 20181121
Address: Department for Medical Epidemiology and Biostatistics (MEB)
	 Karolinska Institutet
	 Box 281
	 171 77 Stockholm
	hannah.bower@ki.se

------------------------------------------------------------------------------------
This document provides explanations regarding the folders included in the folder strcs version 1.85.


Project: strcs: A Stata command for fitting flexible parametric survival models on the log
		hazard scale
Authors: Hannah Bower(1), Michael J. Crowther(1,2), Paul C. Lambert(1,2)
	(1) Department of Medical Epidemiology and Biostatistics, Karolinska Institutet,
		Stockholm, Sweden.
	(2) Department of Health Sciences, University of Leicester, Leicester, UK.


------------------------------------------------------------------------------------
DESCRIPTION:
This folder contains the strcs command where the two-step integration process has been removed,
i.e., Gaussian quadrature is used to numerically integrate over the entire time scale. Score and
Hessian equations have been calculated to increase computational speed. 

Version 1.85 adds a fix for delayed entry. When delayed entry is used, MLMethod
gf0 is used (problem lies in that the delayed entry could be after the first knots
and that this isn't incorporated in the analytical score and Hessian.)
