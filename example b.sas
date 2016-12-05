* Result for clustered survival times with informative censoring in Huang and Wolfe (2002);
* Event=1: failure, Event=2: informative censoring, Event=0: indepdent censoring;

data one;
input id stoptime event trt;
aa=1;
cards;
1 5.0 1 1	/* failure at time 5.0 */
1 6.7 1 2	/* informative censoring at time 6.7 */
1 10.0 2 0	/* independent censoring at time 10.0 */
2 3.0 1 0	/* independent censoring at time 3.0 */
2 5.4 1 1	/* failure at time 5.4 */
2 7.9 0 0	/* independent censoring at time 7.9 */
3 4.5 1 1
3 6.0 0 1
......
;
run;
* event=1: death;
* Get all the death times;
data two;
set one;
if event=1;
run;
* Calculate the quantiles of the failure time;
proc univariate data=two noprint;
var stoptime; 
output out=quant_1 pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=qa; 
run;
data quant_1;
set quant_1;
aa=1;
run;
* Get all the informative censoring event times;
* event=2: informative censoring;
data three;
set one;
if event=2;
run;
* Calculate the quantiles of the informative censoring event time;
proc univariate data=three noprint;
var stoptime; 
output out=quant_2 pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=qb; 
run;
data quant_2;
set quant_2;
aa=1;
run;

* Merge data with the quantiles;
data four;
merge one quant_1 quant_2;
by aa;
run;

* Calculate the duration in each quantile interval, together with the indicator of event in each interval;
data five;
set four;
array quant_1 {11} qa0 qa10 qa20 qa30 qa40 qa50 qa60 qa70 qa80 qa90 qa100;
array quant_2 {11} qb0 qb10 qb20 qb30 qb40 qb50 qb60 qb70 qb80 qb90 qb100;

array dur_a {10} dur_a1-dur_a10;
array dur_b {10} dur_b1-dur_b10;

array event_a {10} event_a1-event_a10;
array event_b {10} event_b1-event_b10;

do i=1 to 10;
	dur_a{i}=0;
	dur_b{i}=0;
end;

do i=1 to 10;
	event_a{i}=0;
	event_b{i}=0;
end;

* For failure;
if event=1 then do;
	do i=2 to 11;
		if stoptime<=quant_a{i} then do;
			event_a{i-1}=1;		/* Failure indicator */
			i=11;
		end;
	end;
end;

* For informative censoring;

if event=2 then do;
	do i=2 to 11;
		if stoptime<=quant_b{i} then do;
			event_b{i-1}=1;		/* Informative censoring indicator */
			i=11;
		end;
	end;
end;

do i=2 to 11;
	if stoptime<=quant_a{i} then do;
		dur_a{i-1}=max(0, stoptime-quant_a{i-1});	/* Duration in each quantile interval for failure */
		i=11;
	end;
	else dur_a{i-1}=quant_a{i}-quant_a{i-1};
end;

do i=2 to 11;
	if stoptime<=quant_b{i} then do;
		dur_b{i-1}=max(0, stoptime-quant_b{i-1});	/* Duration in each quantile interval for informative censoring*/
		i=11;
	end;
	else dur_b{i-1}=quant_b{i}-quant_b{i-1};
end;


run;

* This is the code for Normal frailty with variance theta;
proc nlmixed data=five qpoints=5;

parms r01=1 r02=1 r03=1 r04=1 r05=1 r06=1 r07=1 r08=1 r09=1 r10=1 beta1=.1 beta2=-1.4 theta=1
	  h01=1 h02=1 h03=1 h04=1 h05=1 h06=1 h07=1 h08=1 h09=1 h10=1 alpha1=.2 alpha2=1.2 gamma=-1;
bounds r01 r02 r03 r04 r05 r06 r07 r08 r09 r10 h01 h02 h03 h04 h05 h06 h07 h08 h09 h10 theta >=0;

/* baseline hazard for failure */

base_haz_a=h01 * event_a1 + h02 * event_a2 + h03 * event_a3 + h04 * event_a4 + h05 * event_a5 + 
h06 * event_a6 + h07 * event_a7 + h08* event_a8 +h09 * event_a9 + h10 * event_a10;

/* cumulative baseline hazard for failure */

cum_base_haz_a=h01 * dur_a1 + h02 * dur_a2 + h03 * dur_a3 + h04 * dur_a4 + h05 * dur_a5 + h06 * dur_a6 + 
h07 * dur_a7 + h08* dur_a8 +h09 * dur_a9 + h10 * dur_a10;

/* baseline hazard for informative censoring */

base_haz_b=r01 * event_b1 + r02 * event_b2 + r03 * event_b3 + r04 * event_b4 + r05 * event_b5 + 
r06 * event_b6 + r07 * event_b7 + r08* event_b8 +r09 * event_b9 + r10 * event_b10;

/* cumulative baseline hazard for informative censoring */

cum_base_haz_b=r01 * dur_b1 + r02 * dur_b2 + r03 * dur_b3 + r04 * dur_b4 + r05 * dur_b5 + r06 * dur_b6 + 
r07 * dur_b7 + r08* dur_b8 +r09 * dur_b9 + r10 * dur_b10;

mu1= beta1 * X1 + beta2 * X2  + a;	/* for failure event */

mu2= alpha1 * X1 + alpha2 * X2 +  gamma * a;	/* for informative censoring event */


loglik1=-exp(mu1) * cum_base_haz_a;

loglik2=-exp(mu2) * cum_base_haz_b;

if event=1 then loglik= log(base_haz_a) + mu1 + loglik1 + loglik2; 	/*log likelihood for faiure */
if event=2 then loglik= log(base_haz_d) + mu2 + loglik1 + loglik2;	/*log likelihood for informative censoring */
if event=0 then loglik=loglik1 + loglik2;							/*log likelihood for indepedent censoring */
model stoptime ~ general(loglik);
random a ~ normal(0, theta) subject=id;
ods output ParameterEstimates=est1_all FitStatistics=fit1_all; 
run;
