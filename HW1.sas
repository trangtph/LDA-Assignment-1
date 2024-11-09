libname m '/home/u62794741/LDA/';

proc contents data = m.hemodialysis;
run;

/* --------- Data exploration ---------------------------------- */

/* Description of the sample */


proc sql; /* Age and sex summary */
   create table summary_hemo as
   select count(distinct ID) as num_patients,
   		mean(AGE) as age_avg, 
   		std(AGE) as age_sd, 
   		mean(sex=1) as male_proportion format=percent8.2
   from (select distinct ID, AGE, SEX
         from m.hemodialysis);
quit;

proc print data=summary_hemo;
   title "Patient characteristics summary";
run;

/* Level of missing data */

proc sql;
   create table id_level_missing as
   select distinct ID, SEX, AGE
         from m.hemodialysis
quit;

proc means data=id_level_missing nmiss n noprint;
   var SEX AGE;
   output out= id_level_missing2 (drop=_type_ _freq_)
      nmiss=miss_SEX miss_AGE
      n = total_SEX total_AGE;
run;

data missing_sex_age_summary;
   set id_level_missing2;
   pct_miss_sex = miss_SEX / (total_SEX + miss_SEX) ;
   pct_miss_age = miss_AGE / (total_AGE + miss_AGE);
   format pct_miss_: percent8.2;
run;

proc print data=missing_sex_age_summary;
   title "Overall Missing Data  for age and sex";
run;

proc means data=m.hemodialysis nmiss n noprint;
   var DOSE hb iron;
   output out=overall_missing (drop=_type_ _freq_)
      nmiss=miss_DOSE miss_hb miss_iron
      n=total_DOSE total_hb total_iron;
run;

data missing_summary;
   set overall_missing;
   pct_miss_DOSE = miss_DOSE / (total_DOSE + miss_DOSE) ;
   pct_miss_hb = miss_hb / (total_hb + miss_hb);
   pct_miss_iron = miss_iron / (total_iron +miss_iron);
   format pct_miss_: percent8.2;
run;


proc print data=missing_summary;
   title "Overall Missing Data Summary for DOSE, hb, and iron";
run;

/* average number of measures per patient */
proc sql;
	create table avg_measure_patient as
   select avg(measures_per_patient) as avg_measures_per_patient
   from (select ID, count(*) as measures_per_patient
         from m.hemodialysis
         group by ID);
quit;

proc print data=avg_measure_patient;
   title "Average number of measures per patient";

/* Nr of patients, avg EPO doses, Hb and Iron deficiency status per month*/
proc sql;
   create table patients_per_month as
   select month, count(distinct ID) as number_patients
   from m.hemodialysis
   group by month;
quit;

proc means data=m.hemodialysis noprint;
   class month;
   var hb DOSE;
   output out=monthly_summary (drop=_type_ _freq_)
      mean(hb)=avg_hb std(hb)=std_hb
      mean(DOSE)=avg_DOSE std(DOSE)=std_DOSE;
run;

proc sql;
   create table iron_proportion as
   select month,
          mean(IRON=0) as prop_iron_1 format=percent8.2
   from m.hemodialysis
   group by month;
quit;

data final_monthly_summary;
   merge patients_per_month monthly_summary iron_proportion;
   by month;
run;

proc print data=final_monthly_summary;
   title "Monthly Summary";
run;

/* ----------- Individual profile ------------ */

/* Select 50 random patients */

proc surveyselect data=m.hemodialysis out=sampled_patients
   method=srs n=50 seed=1003; /*simple random sampling */
   id ID;
run;

proc sql;
   create table subset_50_hemo as
   select a.*
   from m.hemodialysis as a
   inner join sampled_patients as b
   on a.ID = b.ID;
quit;

/*proc print data=subset_50_hemo;
   title "Data from Random Sample of 50 Patients";
run; */

/* Plot Individual profiles of 50 random patients */
proc sort data=subset_50_hemo;
by ID month;
run;

goptions reset=all ftext=swiss device=psepsf gsfname=fig1
gsfmode=replace rotate=landscape i=join;
proc gplot data=subset_50_hemo;
plot hb*month=id / haxis=axis1 vaxis=axis2 nolegend;
axis1 label=(h=2 'Time in month') value=(h=1.5)
order=(0 to 6 by 1) minor=none;
axis2 label=(h=2 A=90 'Hemoglobin level (g/dl)') value=(h=1.5)
order=(4 to 18 by 2) minor=none;
title h=3 'Individual profiles for 50 random subjects';
run;quit;


/* ----------Mean structure ----------*/

goptions reset=all ftext=swiss device=psepsf gsfname=fig2 gsfmode=replace
rotate=landscape;


proc gplot data=m.hemodialysis;
plot hb*month / haxis=axis1 vaxis=axis2;
symbol c=brown i=std1mjt w=2 mode=include;
axis1 label=(h=2 'Time (months)') value=(h=1.5) order=(1 to 6 by 1) minor=none;
axis2 label=(h=2 A=90 'Hemoglobin level (g/dl)') value=(h=1.5) order=(10 to 12 by 0.5)
minor=none;
title h=3 'Average evolution of Hb level, with standard errors of means';
run;quit;

/* Mean structure by sex */
goptions reset=all ftext=swiss device=psepsf gsfname=fig3 gsfmode=replace
   rotate=landscape;

/* Define different symbol settings for males and females */
symbol1 c=cadetblue i=std1mjt w=2 v=none mode=include; /* male */
symbol2 c=brown i=std1mjt w=2 v=none mode=include; /* female */

legend1 label=none 
         value=('Male' 'Female')
         position=(top right) mode=share across=1;

proc gplot data=m.hemodialysis;
   plot hb*month=SEX / haxis=axis1 vaxis=axis2 legend=legend1;
   axis1 label=(h=2 'Time (months)') value=(h=1.5) order=(1 to 6 by 1) minor=none;
   axis2 label=(h=2 A=90 'Hemoglobin level (g/dl)') value=(h=1.5) order=(9 to 15 by 0.5)
         minor=none;
   title h=3 'Average evolution of Hb level by sex, with standard errors of means';
run;
quit;


/* -------- Variance structure --------- */

/* Calculate monthly mean Hb */
proc means data=m.hemodialysis noprint;
   class month;
   var hb;
   output out=mean_hb mean=mean_hb;
run;

/* Calculate variance by month */

proc sort data=m.hemodialysis;
   by month;
run;

proc sort data=mean_hb;
   by month;
run;

data variance_calc;
   merge m.hemodialysis (in=a) mean_hb (in=b);
   by month;
   if a and b; /*Keeps only records where month is present in both datasets */
   /* Calculate squared deviation for each observation */
   sq_residual = (hb - mean_hb)**2;
run;


/* Plot var(t) and its standard error by month */


goptions reset=all ftext=swiss device=psepsf gsfname=fig4 gsfmode=replace
rotate=landscape;

proc gplot data=variance_calc;
plot sq_residual*month / haxis=axis1 vaxis=axis2;
symbol c=brown i=std1mjt w=2 mode=include;
axis1 label=(h=2 'Time (months)') value=(h=1.5) order=(1 to 6 by 1) minor=none;
axis2 label=(h=2 A=90 'Squared residuals') value=(h=1.5) order=(0 to 4 by 0.5)
minor=none;
title h=3 'Variance structure of Hb values over time';
run;quit;

/* ---- Correlation structure --- */

proc sort data=m.hemodialysis;
   by ID month;
run;

proc transpose data=m.hemodialysis out=transposed_hemo prefix=hb_;
   by ID;
   id month;
   var hb;
run;

data transposed_hemo;
   set transposed_hemo;
   label hb_1 = "Month 1"
         hb_2 = "Month 2"
         hb_3 = "Month 3"
         hb_4 = "Month 4"
         hb_5 = "Month 5"
         hb_6 = "Month 6";
run;

proc sgscatter data=transposed_hemo;
   matrix hb_1 hb_2 hb_3 hb_4 hb_5 hb_6 / diagonal=(histogram);
   title "Scatter Plot Matrix of monthly Hb values";
run;

