/*using data from SAS LIFEREG documentation*/

title 'Motorette Failures With Operating Temperature as a Covariate';
data motors;
input time censor temp @@;
if _N_=1 then
do;
temp=130;
time=.;
control=1;
z=1000/(273.2+temp);
output;
temp=150;
time=.;
control=1;
z=1000/(273.2+temp);
output;
end;
if temp>150;
control=0;
z=1000/(273.2+temp);
output;
datalines;
8064 0 150 8064 0 150 8064 0 150 8064 0 150 8064 0 150
8064 0 150 8064 0 150 8064 0 150 8064 0 150 8064 0 150
1764 1 170 2772 1 170 3444 1 170 3542 1 170 3780 1 170
4860 1 170 5196 1 170 5448 0 170 5448 0 170 5448 0 170
408 1 190 408 1 190 1344 1 190 1344 1 190 1440 1 190
1680 0 190 1680 0 190 1680 0 190 1680 0 190 1680 0 190
408 1 220 408 1 220 504 1 220 504 1 220 504 1 220
528 0 220 528 0 220 528 0 220 528 0 220 528 0 220
;

run;

/* RUN REGRESSION
*/

proc lifereg data=motors;
model time*censor(0)=z /dist=gamma;
output out=out_gamma predicted=pred_ xbeta=xbeta_;
ods output ParameterEstimates=parm_gamma;
run;

/*analytical calculation of GENGAMMA quantile that uses PROC LIFEREG parametrization is broken down into three steps:

1) obtain coefficients from PROC LIFEREG (stable parametrization https://go.documentation.sas.com/doc/en/pgmsascdc/9.4_3.4/statug/statug_lifereg_details04.htm)
2) map obtained coefficients to original GENGAMMA parametrization
3) map original gengamma parametrization to standard Gamma parameters based on https://blogs.sas.com/content/iml/2021/03/15/generalized-gamma-distribution.html

NOTE: calculations in step 2, 3 must be performed by hand, result are below
*/

/* extract estimates based on SAS documentation stable parametrization */

data _null_;
set work.parm_gamma;
where Parameter eq "Scale";
call symput('sigma', Estimate);
run;
%put &=sigma;

/* Q = delta in stable parametrization */
data _null_;
set work.parm_gamma;
where Parameter eq "Shape";
call symput('Q', Estimate);
run;
%put &=Q;

/*map stable parametrization to parameters of original GENGAMMA distribution
and to PROC LIFEREG quantile using relation to GAMMA distribution (step 2, 3).
Make use of Xbeta from out dataset
*/

data res_gamma;
set work.out_gamma;
k = 1/(&Q.*&Q.);
b = &Q./&sigma.;
log_a = xbeta_ - log(k)/b;
log_up = log(quantile("GAMMA", 0.5, k)); /* use percentile = 0.5 as by default in lifereg*/
pred_reproduce = exp(log_a + (1/b)*log_up); /*compare to pred_*/
run;

proc print data=res_gamma;
