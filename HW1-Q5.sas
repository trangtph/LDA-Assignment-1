libname m '/home/u62794741/LDA/';


proc contents data = m.hemodialysis;
run;

/* ---------- Data manipulation ---------------- */
/* Create new DOSE2 variable with 1 month lag */ 

proc sort data=m.hemodialysis;
by ID month;
run;

data m.hemodialysis2;
   set m.hemodialysis;
   by ID month;

   /* Retain the previous month's DOSE for each ID */
   retain prev_DOSE;
   
   /* Initialize DOSE2 */
   if first.ID then DOSE2 = 0;
   else DOSE2 = prev_DOSE;
   
   /* Update prev_DOSE with the current DOSE for the next observation */
   prev_DOSE = DOSE;

run;

/* reparameterize month */
/* Dummy variables for sex */
data m.hemodialysis2;
   set m.hemodialysis2;
   month_class = month;
   month0 = month - 1;
   month0_class = month0;
   SEX_male = (SEX = 1);
run;




/* ---------------- First fit a most elaborated OLS model  -------- */

proc glm data=m.hemodialysis2;
   model hb = DOSE2 AGE SEX_male iron month 
   	DOSE2*AGE DOSE2*SEX_male DOSE2*iron DOSE2*month iron*month age*month SEX_male*month / solution;
   output out=m.residuals_data1 r=residual;
   title "Elaborated OLS model for Hb";
run;
quit;

proc glm data=m.hemodialysis2;
	class SEX iron;
   model hb = DOSE2 AGE SEX iron month 
   	DOSE2*AGE DOSE2*SEX DOSE2*iron DOSE2*month iron*month age*month SEX*month / noint solution;
   title "Elaborated OLS model for Hb - SEX as factor";
run;
quit;

proc glm data=m.hemodialysis2;
   model hb = DOSE2 AGE SEX_male iron month0 
   	DOSE2*AGE DOSE2*SEX_male DOSE2*iron DOSE2*month0 iron*month0 age*month0 SEX_male*month0 / solution;
      output out=m.residuals_data r=residual;
   title "Elaborated OLS model for Hb with month - 1";
run;
quit;


/* Plot the residuals */
/* Select 50 random patients */

proc surveyselect data=m.residuals_data out=sampled_residual
   method=srs n=50 seed=1003; /*simple random sampling */
   id ID;
run;

proc sql;
   create table subset_50_res as
   select a.*
   from m.residuals_data as a
   inner join sampled_residual as b
   on a.ID = b.ID;
quit;


/* Plot OLS profiles of 50 random patients */
proc sort data=subset_50_res;
by ID month;
run;

goptions reset=all ftext=swiss device=psepsf gsfname=fig1
gsfmode=replace rotate=landscape i=join;
proc gplot data=subset_50_res;
plot residual*month=id / haxis=axis1 vaxis=axis2 nolegend;
axis1 label=(h=2 'Time in month') value=(h=1.5)
order=(0 to 6 by 1) minor=none;
axis2 label=(h=2 A=90 'Residuals') value=(h=1.5)
order=(-5 to 5 by 1) minor=none;
title h=3 'OLS residual profiles for 50 random subjects';
run;quit;

/* Plot the smoothed average of squared OLS residuals */

data residuals_data_squared;
   set m.residuals_data;
   squared_residual = residual**2;
run;

proc loess data=residuals_data_squared;
ods output scoreresults=out;
model squared_residual=month;
score data=residuals_data_squared;
run;

proc sort data=out;
by month;
run;

data out;
	set out;
run;

goptions reset=all ftext=swiss device=psepsf
gsfname=fig1 gsfmode=replace rotate=landscape;
proc gplot data=out;
plot squared_residual*month=1 p_squared_residual*month=2
/ overlay haxis=axis1 vaxis=axis2;
symbol1 c=red v=dot h=0.2 mode=include;
symbol2 c=black i=join w=2 mode=include;
axis1 label=(h=2 'Time in month')
value=(h=1.5) order=(0 to 6 by 1) minor=none;
axis2 label=(h=2 A=90 'Squared residual') value=(h=1.5)
order=(0 to 8 by 1) minor=none;
title h=3 'Smoothed variance function';
run;quit;

/* ------------------------------------------------------------------------------------------*/
/* -------- Fit elaborated mixed model with random intercept and random slope ---------------*/


proc mixed data=m.hemodialysis2 method=reml;
   	class ID SEX iron month_class;
   	model hb = DOSE2 AGE SEX iron month 
   	DOSE2*AGE DOSE2*SEX DOSE2*iron DOSE2*month iron*month age*month SEX*month / solution;
   	random intercept month / type=un subject=ID g gcorr v vcorr;
	repeated month_class / type=simple subject=ID r rcorr;
run;
quit;

/* Dummy for sex */
proc mixed data=m.hemodialysis2 method=reml;
   	class ID month_class;
   	model hb = DOSE2 AGE SEX_male iron month 
   	DOSE2*AGE DOSE2*SEX_male DOSE2*iron DOSE2*month iron*month age*month SEX_male*month / solution;
   	random intercept month / type=un subject=ID g gcorr v vcorr;
	repeated month_class / type=simple subject=ID r rcorr;
run;
quit;

/* Month 0-5 */
proc mixed data=m.hemodialysis2 method=reml;
   	class ID SEX iron month0_class;
   	model hb = DOSE2 AGE SEX iron month0 
   	DOSE2*AGE DOSE2*SEX DOSE2*iron DOSE2*month0 iron*month0 age*month0 SEX*month0 / solution;
   	random intercept month0 / type=un subject=ID g gcorr v vcorr;
	repeated month0_class / type=simple subject=ID r rcorr;
run;
quit;

proc mixed data=m.hemodialysis2 method=reml;
   	class ID month0_class;
   	model hb = DOSE2 AGE SEX_male iron month0 
   	DOSE2*AGE DOSE2*SEX_male DOSE2*iron DOSE2*month0 iron*month0 AGE*month0 SEX_male*month0 / solution;
   	random intercept month0 / type=un subject=ID g gcorr v vcorr;
	repeated month0_class / type=simple subject=ID r rcorr;
	title "Elaborated mixed model with month0 & dummy SEX";
run;
quit;

/* ------------------------------------------------------------------------------------------*/
/* -------- Model with serial correlation ---------------*/

/* Exponential */

proc mixed data=m.hemodialysis2 method=reml;
   	class ID month0_class;
   	model hb = DOSE2 AGE SEX_male iron month0 
   	DOSE2*AGE DOSE2*SEX_male DOSE2*iron DOSE2*month0 iron*month0 AGE*month0 SEX_male*month0 / solution;
   	random intercept month0 / type=un subject=ID g gcorr v vcorr;
	repeated month0_class / type=sp(gau)(month0) local subject=ID r rcorr;
	title "Elaborated mixed model with Exponential serial correlation";
run;
quit;

/* Gaussian */
proc mixed data=m.hemodialysis2 method=reml;
   	class ID month0_class;
   	model hb = DOSE2 AGE SEX_male iron month0 
   	DOSE2*AGE DOSE2*SEX_male DOSE2*iron DOSE2*month0 iron*month0 AGE*month0 SEX_male*month0 / solution;
   	random intercept month0 / type=un subject=ID g gcorr v vcorr;
	repeated month0_class / type=sp(gau)(month0) local subject=ID r rcorr;
	title "Elaborated mixed model with Gaussian serial correlation";
	ods output FitStatistics=FitFull;
run;
quit;

/* Model with serial correlation fits significantly better than without. But cannot differentiate between
Gaussian and Exponential -> Go with Gaussian */

/* ------------------------------------------------------------------------------------------*/
/* -------- Reduce random effect ---------------*/

/* Only intercept */

proc mixed data=m.hemodialysis2 method=reml;
   	class ID month0_class;
   	model hb = DOSE2 AGE SEX_male iron month0 
   	DOSE2*AGE DOSE2*SEX_male DOSE2*iron DOSE2*month0 iron*month0 AGE*month0 SEX_male*month0 / solution;
   	random intercept / type=un subject=ID g gcorr v vcorr;
	repeated month0_class / type=sp(gau)(month0) local subject=ID r rcorr;
	title "Random Intercept-only mixed model";
	ods output FitStatistics=FitReduced;
run;
quit;

/* LR test random slope + intercept vs random intercept */
data LRT_Results;
   merge FitFull(where=(Descr="-2 Res Log Likelihood") rename=(Value=LL_full))
         FitReduced(where=(Descr="-2 Res Log Likelihood") rename=(Value=LL_reduced));
   LRT_stat = LL_reduced - LL_full;  /* Test statistic */
   df = '1:2';  /* Degrees of freedom: difference in the number of parameters */
   p_value = 0.5*(1 - probchi(LRT_stat, 1) + 1 - probchi(LRT_stat, 2)); /* p-value for chi-square test */
run;

proc print data=LRT_Results;
   var LRT_stat df p_value;
   title "LR test: full model vs intercept only";
run;

/* Only slope */
proc mixed data=m.hemodialysis2 method=reml;
   	class ID month0_class;
   	model hb = DOSE2 AGE SEX_male iron month0 
   	DOSE2*AGE DOSE2*SEX_male DOSE2*iron DOSE2*month0 iron*month0 AGE*month0 SEX_male*month0 / solution;
   	random month0 / type=un subject=ID g gcorr v vcorr;
	repeated month0_class / type=sp(gau)(month0) local subject=ID r rcorr;
	title "Random slope-only mixed model";
	ods output FitStatistics=FitReduced2;
run;
quit;

/* LR test random slope + intercept vs random slope */
data LRT_Results2;
   merge FitFull(where=(Descr="-2 Res Log Likelihood") rename=(Value=LL_full))
         FitReduced2(where=(Descr="-2 Res Log Likelihood") rename=(Value=LL_reduced));
   LRT_stat = LL_reduced - LL_full;  /* Test statistic */
   df = '1:2';  /* Degrees of freedom: difference in the number of parameters */
   p_value = 0.5*(1 - probchi(LRT_stat, 1) + 1 - probchi(LRT_stat, 2)); /* p-value for chi-square test */
run;

proc print data=LRT_Results2;
   var LRT_stat df p_value;
   title "LR test: full model vs slope only";
run;

/* Model with random intercept only was not significantly different from full model 
-> Keep random intercept only */

/* ------------------------------------------------------------------------------------------*/
/* -------- Reduce fixed effect ---------------*/

proc mixed data=m.hemodialysis2 method=reml;
   	class ID month0_class;
   	model hb = DOSE2 AGE SEX_male iron month0 
   	DOSE2*AGE DOSE2*SEX_male DOSE2*iron DOSE2*month0 iron*month0 AGE*month0 SEX_male*month0 / solution;
   	random intercept / type=un subject=ID g gcorr v vcorr;
	repeated month0_class / type=sp(gau)(month0) local subject=ID r rcorr;
	title "Random intercept model";
	contrast 'Reduce interaction of DOSE' DOSE2*iron 1,
										DOSE2*AGE 1,
										DOSE2*SEX_male 1, /chisq;
run;
quit;

proc mixed data=m.hemodialysis2 method=reml;
   	class ID month0_class;
   	model hb = DOSE2 AGE SEX_male iron month0 
   	DOSE2*AGE DOSE2*SEX_male DOSE2*iron DOSE2*month0 iron*month0 AGE*month0 SEX_male*month0 / solution;
   	random intercept / type=un subject=ID g gcorr v vcorr;
	repeated month0_class / type=sp(gau)(month0) local subject=ID r rcorr;
	title "Random intercept model";
	contrast 'Reduce intera Dose & month' DOSE2*iron 1,
										DOSE2*AGE 1,
										DOSE2*SEX_male 1,
										AGE*month0 1,
										SEX_male*month0 1,/chisq;
run;
quit;

proc mixed data=m.hemodialysis2 method=reml;
   	class ID month0_class;
   	model hb = DOSE2 AGE SEX_male iron month0 
   	DOSE2*AGE DOSE2*SEX_male DOSE2*iron DOSE2*month0 iron*month0 AGE*month0 SEX_male*month0 / solution;
   	random intercept / type=un subject=ID g gcorr v vcorr;
	repeated month0_class / type=sp(gau)(month0) local subject=ID r rcorr;
	title "Random intercept model";
	contrast 'Reduce intera Dose & month' DOSE2*iron 1,
										DOSE2*AGE 1,
										DOSE2*SEX_male 1,
										AGE*month0 1,
										SEX_male*month0 1,
										iron*month0 1,/chisq;
run;
quit;