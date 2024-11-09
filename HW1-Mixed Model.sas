libname m '/home/u62794741/LDA/';


proc contents data = m.hemodialysis;
run;

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

/* ---------------- First fit a most elaborated OLS model  -------- */

proc glm data=m.hemodialysis2;
   class SEX iron;
   model hb = DOSE2 AGE SEX iron month DOSE2*AGE DOSE2*SEX DOSE2*iron DOSE2*month iron*month / solution;
   output out=m.residuals_data r=residual;
   title "Elaborated OLS model for Hb";
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
   set residuals_data;
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

/* -------- Fit elaborated mixed model with random effect and random slope */

data m.hemodialysis2;
   set m.hemodialysis2;
   month_class = month;
run;

proc mixed data=m.hemodialysis2 method=reml;
   	class ID SEX iron month_class;
   	model hb = DOSE2 AGE SEX iron month DOSE2*AGE DOSE2*SEX DOSE2*iron DOSE2*month iron*month /noint solution;
   	random intercept month / type=un subject=ID g gcorr v vcorr;
	repeated month_class / type=simple subject=ID r rcorr;
	ods output solutionr= m.lmm1_out;
run;
quit;