options ls=150 formdlim="#";


/**********
***********
**** #2 ***
***********
**********/

data receptor;
  infile 'c:\My Documents\stat611\receptor2.txt';
  input Ethylene Propane nButane x2Mp x3Mp Benzene Cyclohexane x2Mh x224Tmp Acetylene Ethane;
run;


title '#3a (Receptor Data): Principal Factor Method with Varimax Rotation';
proc factor method=principal 
            priors=smc   /* COMMENT: smc tells Proc Factor to use the squared multiple 
                            correlation with all the other vars as the prior 
			                communality estimate. */
            mineigen=1
            rotate=promax   
            reorder 
            plot;
  var Ethylene Propane nButane x2Mp x3Mp Benzene Cyclohexane x2Mh x224Tmp Acetylene Ethane;
run;


title '#3b (Receptor Data): Maximum Likelihood Method (k=1 factor)';
proc factor method=ml
            heywood   /* sets to 1 any communality greater than var{x_i}=1  */
			n=1;   /* use 1 factor */
  var Ethylene Propane nButane x2Mp x3Mp Benzene Cyclohexane x2Mh x224Tmp Acetylene Ethane;
run;


title '#3b (Receptor Data): Maximum Likelihood Method (k=2 factor)';
proc factor method=ml
            heywood   /* sets to 1 any communality greater than var{x_i}=1  */
			n=2;   /* use 2 factors */
  var Ethylene Propane nButane x2Mp x3Mp Benzene Cyclohexane x2Mh x224Tmp Acetylene Ethane;
run;


title '#3b (Receptor Data): Maximum Likelihood Method (k=3 factor)';
proc factor method=ml
            heywood   /* sets to 1 any communality greater than var{x_i}=1  */
			n=3;   /* use 3 factors */
  var Ethylene Propane nButane x2Mp x3Mp Benzene Cyclohexane x2Mh x224Tmp Acetylene Ethane;
run;


title '#3b (Receptor Data): Maximum Likelihood Method (k=2 factor) with Promax Rotation';
proc factor method=ml
            heywood 
			n=2   /* use 2 factors */
            rotate = promax
            plot;
  var Ethylene Propane nButane x2Mp x3Mp Benzene Cyclohexane x2Mh x224Tmp Acetylene Ethane;
run;





/**********
***********
*** 14.5 ***
***********
**********/

data seishu;
  infile 'c:\Documents and Settings\William Christensen\My Documents\stat666\MULTIVAR\SEISHU.DAT';
  input taste odor ph acid1 acid2 sake drsugar totsugar alcohol fn;
run;

title '#4 Seishu: PROC CALIS--Evaluating a specific hypothesis about the structure';
proc calis method=ml cov maxiter=5000 maxfunc=5000;
  lineqs
     taste    =        f1                                     + e01,
	 odor     = lam021 f1 + lam022 f2 + lam023 f3 + lam024 f4 + e02,
	 ph       =                    f2                         + e03,
     acid1    = lam041 f1 + lam042 f2 + lam043 f3 + lam044 f4 + e04,
     acid2    = lam051 f1 + lam052 f2 + lam053 f3 + lam054 f4 + e05,
     sake     = lam061 f1 + lam062 f2 + lam063 f3 + lam064 f4 + e06,
     drsugar  = lam071 f1 + lam072 f2 + lam073 f3 + lam074 f4 + e07,
     totsugar =                                f3             + e08,
     alcohol  =                                            f4 + e09,
     fn       = lam101 f1 + lam102 f2 + lam103 f3 + lam104 f4 + e10;
  std
     e01-e10= psi01-psi10,
	 f1-f4= phi1-phi4;
  cov
     f1-f4= phi11-phi16;
  bounds
     0 <= phi1-phi4,
	 0 <= psi01-psi10;
run;
/* NOTE:  Chi^2 = 14.8718, df = 12, p-value = 0.2485, CFI = .9824, RMSEA = .0908, AIC = -9.1282 */



title '#4 Seishu: Set unimportant lam051 to 0 due to non-significant t-value';
proc calis method=ml cov maxiter=5000 maxfunc=5000;
  lineqs
     taste    =        f1                                     + e01,
	 odor     = lam021 f1 + lam022 f2 + lam023 f3 + lam024 f4 + e02,
	 ph       =                    f2                         + e03,
     acid1    = lam041 f1 + lam042 f2 + lam043 f3 + lam044 f4 + e04,
     acid2    =             lam052 f2 + lam053 f3 + lam054 f4 + e05,
     sake     = lam061 f1 + lam062 f2 + lam063 f3 + lam064 f4 + e06,
     drsugar  = lam071 f1 + lam072 f2 + lam073 f3 + lam074 f4 + e07,
     totsugar =                                f3             + e08,
     alcohol  =                                            f4 + e09,
     fn       = lam101 f1 + lam102 f2 + lam103 f3 + lam104 f4 + e10;
  std
     e01-e10= psi01-psi10,
	 f1-f4= phi1-phi4;
  cov
     f1-f4= phi11-phi16;
  bounds
     0 <= phi1-phi4,
	 0 <= psi01-psi10;
run;
/* NOTE:  Chi^2 = 15.1252, df = 13, p-value = 0.2996, CFI = .9870, RMSEA = .0751, AIC = -10.8748 */



title 'BAD IDEA...set an important lambda to 0: lam073';
proc calis method=ml cov maxiter=5000 maxfunc=5000;
  lineqs
     taste    =        f1                                     + e01,
	 odor     = lam021 f1 + lam022 f2 + lam023 f3 + lam024 f4 + e02,
	 ph       =                    f2                         + e03,
     acid1    = lam041 f1 + lam042 f2 + lam043 f3 + lam044 f4 + e04,
     acid2    = lam051 f1 + lam052 f2 + lam053 f3 + lam054 f4 + e05,
     sake     = lam061 f1 + lam062 f2 + lam063 f3 + lam064 f4 + e06,
     drsugar  = lam071 f1 + lam072 f2             + lam074 f4 + e07,
     totsugar =                                f3             + e08,
     alcohol  =                                            f4 + e09,
     fn       = lam101 f1 + lam102 f2 + lam103 f3 + lam104 f4 + e10;
  std
     e01-e10= psi01-psi10,
	 f1-f4= phi1-phi4;
  cov
     f1-f4= phi11-phi16;
  bounds
     0 <= phi1-phi4,
	 0 <= psi01-psi10;
run;
/* NOTE:  Chi^2 = 34.9273, df = 13, p-value = 0.0009, CFI = .8660, RMSEA = .2412, AIC = 8.9273 */


title '#4 Seishu: Set lam022, 023, 041, 051, 053, 101, and 103 (all that had |t|<1) 
  to 0 due to non-significant t-value';
proc calis method=ml cov maxiter=5000 maxfunc=5000;
  lineqs
     taste    =        f1                                     + e01,
	 odor     = lam021 f1                         + lam024 f4 + e02,
	 ph       =                    f2                         + e03,
     acid1    =             lam042 f2 + lam043 f3 + lam044 f4 + e04,
     acid2    =             lam052 f2             + lam054 f4 + e05,
     sake     = lam061 f1 + lam062 f2 + lam063 f3 + lam064 f4 + e06,
     drsugar  = lam071 f1 + lam072 f2 + lam073 f3 + lam074 f4 + e07,
     totsugar =                                f3             + e08,
     alcohol  =                                            f4 + e09,
     fn       =             lam102 f2             + lam104 f4 + e10;
  std
     e01-e10= psi01-psi10,
	 f1-f4= phi1-phi4;
  cov
     f1-f4= phi11-phi16;
  bounds
     0 <= phi1-phi4,
	 0 <= psi01-psi10;
run;
/* NOTE:  Chi^2 = 17.1967, df = 19, p-value = 0.5765, CFI = 1.0000, RMSEA = .0000, AIC = -20.8033 */



title '#4 Seishu: Additionally, set lam024, 042, 043, 061, 064, 071, 072, and 074 (all that had |t|<2) 
  to 0 due to non-significant t-value';
proc calis method=ml cov maxiter=5000 maxfunc=5000;
  lineqs
     taste    =        f1                                     + e01,
	 odor     = lam021 f1                                     + e02,
	 ph       =                    f2                         + e03,
     acid1    =                                   + lam044 f4 + e04,
     acid2    =             lam052 f2             + lam054 f4 + e05,
     sake     =             lam062 f2 + lam063 f3             + e06,
     drsugar  =                         lam073 f3             + e07,
     totsugar =                                f3             + e08,
     alcohol  =                                            f4 + e09,
     fn       =             lam102 f2             + lam104 f4 + e10;
  std
     e01-e10= psi01-psi10,
	 f1-f4= phi1-phi4;
  cov
     f1-f4= phi11-phi16;
  bounds
     0 <= phi1-phi4,
	 0 <= psi01-psi10;
run;
/* NOTE:  Chi^2 = 36.6472, df = 27, p-value = 0.1018, CFI = .9410, RMSEA = .1110, AIC = -17.3528 */



