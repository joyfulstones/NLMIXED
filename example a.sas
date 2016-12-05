* Result for recurrent event with a terminal event in Liu, Wolfe, and Huang (2004);
* Event=1: recurrent event, Event=2: death, Event=0: independent censoring;
data one;
input id stoptime event trt;
aa=1;
cards;
1 5.0 1 1	/* recurrent event at time 5.0 */
1 6.7 1 1	/* recurrent event at time 6.7 */
1 10.0 2 1	/* death event at time 10.0 */
2 3.0 1 0
2 5.4 1 0
2 7.9 0 0	/* independent censoring event at time 7.9 */
3 4.5 1 1
3 6.0 0 1
......
;
run;
* dataset two includes only recurrent events;
data two;
set one;
if event=1;	/* Keep all the recurrent event */
run;
* Get the quantiles for the recurrent events, store in dataset quant_r;
proc univariate data=two noprint;
var stoptime; 
output out=quant_r pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=qr; 
run;
data quant_r;
set quant_r;
aa=1;
run;
* Get the quantiles for the death/censoring event;
data three;
set one;
if event ne 1;
run;
* Get the quantiles for the death event, store in dataset quant_d;
proc univariate data=three noprint;
var stoptime; 
output out=quant_d pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=qd; 
where event=2;
run;
data quant_d;
set quant_d;
aa=1;
run;
* Merge data with the quantiles;
data four;
merge one quant_r quant_d;
by aa;
run;

* Calculate the duration in each quantile interval, together with the indicator of event in each interval;
data five;
set four;
array quant_r {11} qr0 qr10 qr20 qr30 qr40 qr50 qr60 qr70 qr80 qr90 qr100;
array quant_d {11} qd0 qd10 qd20 qd30 qd40 qd50 qd60 qd70 qd80 qd90 qd100;

array dur_r {10} dur_r1-dur_r10;
array dur_d {10} dur_d1-dur_d10;

array event_r {10} event_r1-event_r10;
array event_d {10} event_d1-event_d10;

do i=1 to 10;
	dur_r{i}=0;
	dur_d{i}=0;
end;

do i=1 to 10;
	event_r{i}=0;
	event_d{i}=0;
end;

* For recurrent event;
if event=1 then do;
	do i=2 to 11;
		if stoptime<=quant_r{i} then do;
			event_r{i-1}=1;		/* indicator of recurrent event in each interval */
			i=11;
		end;
	end;
end;

else do; /* If death or censored observation */
	do i=2 to 11;	/* duration in each recurrent event quantiles */
		if stoptime<=quant_r{i} then do;
			dur_r{i-1}=stoptime-quant_r{i-1};
			i=11;
		end;
		else dur_r{i-1}=quant_r{i}-quant_r{i-1};
	end;

	do i=2 to 11; /* duration in each death event quantiles */
		if stoptime<=quant_d{i} then do;
			event_d{i-1}=(event=2);		/* indicator of death event in each interval */
			dur_d{i-1}=stoptime-quant_d{i-1};
			i=11;
		end;
		else dur_d{i-1}=quant_d{i}-quant_d{i-1};
	end;
end;

run;
* This is the code for Gamma frailty with variance theta;
proc nlmixed data=five qpoints=30 noad;

parms r01=1 r02=1 r03=1 r04=1 r05=1 r06=1 r07=1 r08=1 r09=1 r10=1 /* Starting values for recurrent baseline hazard */
	  h01=1 h02=1 h03=1 h04=1 h05=1 h06=1 h07=1 h08=1 h09=1 h10=1 /* Starting values for death baseline hazard */ 
	  beta1=0 theta=.2 alpha1=-.2 gamma=0;
bounds r01 r02 r03 r04 r05 r06 r07 r08 r09 r10 h01 h02 h03 h04 h05 h06 h07 h08 h09 h10 theta >=0;

* base hazard and cumulative baseline hazard for recurrent events;

base_haz_r=r01 * event_r1 + r02 * event_r2 + r03 * event_r3 + r04 * event_r4 + r05 * event_r5 + 
r06 * event_r6 + r07 * event_r7 + r08* event_r8 +r09 * event_r9 + r10 * event_r10;
cum_base_haz_r=r01 * dur_r1 + r02 * dur_r2 + r03 * dur_r3 + r04 * dur_r4 + r05 * dur_r5 + r06 * dur_r6 + 
r07 * dur_r7 + r08* dur_r8 +r09 * dur_r9 + r10 * dur_r10;

* base hazard and cumulative baseline hazard for death events;

base_haz_d=h01 * event_d1 + h02 * event_d2 + h03 * event_d3 + h04 * event_d4 + h05 * event_d5 + 
h06 * event_d6 + h07 * event_d7 + h08* event_d8 +h09 * event_d9 + h10 * event_d10;
cum_base_haz_d=h01 * dur_d1 + h02 * dur_d2 + h03 * dur_d3 + h04 * dur_d4 + h05 * dur_d5 + h06 * dur_d6 + 
h07 * dur_d7 + h08* dur_d8 +h09 * dur_d9 + h10 * dur_d10;

p=cdf('NORMAL', a);
if p > .999999 then p=.999999;
g2=quantile('GAMMA', p, 1/theta); 
g=g2 * theta;		/* g is gamma distributed with mean 1 and variance theta */

mu1= beta1 * trt +  log(g);			/* for recurrent event */

mu2= alpha1 * trt  + gamma * log(g);	/* for death event */

loglik1=-exp(mu1) * cum_base_haz_r;

loglik2=-exp(mu2) * cum_base_haz_d;

if event=1 then loglik= log(base_haz_r) + mu1 ; 					/*log likelihood for recurrent event */
if event=2 then loglik=loglik1 + log(base_haz_d) + mu2 + loglik2;	/*log likelihood for death */
if event=0 then loglik=loglik1 + loglik2;							/*log likelihood for censoring */

model stoptime ~ general(loglik);
random a ~ normal(0, 1) subject=id;
ods output ParameterEstimates=est1 FitStatistics=fit1; 
run;

* This is the code for Normal frailty with variance theta;

proc nlmixed data=five qpoints=5;

parms r01=1 r02=1 r03=1 r04=1 r05=1 r06=1 r07=1 r08=1 r09=1 r10=1 /* Starting values for recurrent baseline hazard */
	  h01=1 h02=1 h03=1 h04=1 h05=1 h06=1 h07=1 h08=1 h09=1 h10=1 /* Starting values for death baseline hazard */ 
	  beta1=0 theta=.2 alpha1=-.2 gamma=0;
bounds r01 r02 r03 r04 r05 r06 r07 r08 r09 r10 h01 h02 h03 h04 h05 h06 h07 h08 h09 h10 theta >=0;

* base hazard and cumulative baseline hazard for recurrent events;

base_haz_r=r01 * event_r1 + r02 * event_r2 + r03 * event_r3 + r04 * event_r4 + r05 * event_r5 + 
r06 * event_r6 + r07 * event_r7 + r08* event_r8 +r09 * event_r9 + r10 * event_r10;
cum_base_haz_r=r01 * dur_r1 + r02 * dur_r2 + r03 * dur_r3 + r04 * dur_r4 + r05 * dur_r5 + r06 * dur_r6 + 
r07 * dur_r7 + r08* dur_r8 +r09 * dur_r9 + r10 * dur_r10;

* base hazard and cumulative baseline hazard for death events;

base_haz_d=h01 * event_d1 + h02 * event_d2 + h03 * event_d3 + h04 * event_d4 + h05 * event_d5 + 
h06 * event_d6 + h07 * event_d7 + h08* event_d8 +h09 * event_d9 + h10 * event_d10;
cum_base_haz_d=h01 * dur_d1 + h02 * dur_d2 + h03 * dur_d3 + h04 * dur_d4 + h05 * dur_d5 + h06 * dur_d6 + 
h07 * dur_d7 + h08* dur_d8 +h09 * dur_d9 + h10 * dur_d10;

mu1= beta1 * trt +  a;			/* for recurrent event */

mu2= alpha1 * trt  + gamma * a;	/* for death event */

loglik1=-exp(mu1) * cum_base_haz_r;

loglik2=-exp(mu2) * cum_base_haz_d;

if event=1 then loglik= log(base_haz_r) + mu1 ; 					/*log likelihood for recurrent event */
if event=2 then loglik=loglik1 + log(base_haz_d) + mu2 + loglik2;	/*log likelihood for death */
if event=0 then loglik=loglik1 + loglik2;							/*log likelihood for censoring */

model stoptime ~ general(loglik);
random a ~ normal(0, theta) subject=id;
ods output ParameterEstimates=est1 FitStatistics=fit1; 
run;
