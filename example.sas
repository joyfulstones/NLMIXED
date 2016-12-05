* Section 4.1;
* Input the data;
data frailty;
input center pid followup status X1 X2;
aa=1;
lines;
1 1 12.0 1 1 0.5
1 2 8.5 0 0 0.8
1 3 7.8 1 1 0.2
...
;
run;
* Get the event times only;
data frailty2;
set frailty;
if status=1;
run;
* Calculate the quantiles of the event times;
proc univariate data=frailty2 noprint;
var followup; 
output out=quantile pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=q; 
run;
data quantile;
set quantile;
aa=1;
run;
* Merge data with the quantiles;
data frailty3;
merge frailty quantile;
by aa;
run;

* Calculate the duration in each quantile interval, together with the indicator of event in each interval;
data frailty4;
set frailty3;
array quantile {11} q0 q10 q20 q30 q40 q50 q60 q70 q80 q90 q100;
array duration {10} dur1-dur10;
array event {10} event1-event10;

do i=1 to 10;
	duration{i}=0;
end;

do i=1 to 10;
	event{i}=0;
end;

do i=2 to 11;
	if followup<=quantile{i} then do;
		duration{i-1}=followup-quantile{i-1};
		event{i-1}=status;
		i=11;
	end;
	else duration{i-1}=quantile{i}-quantile{i-1};
end;

run;
* Code for Normal frailty with variance theta;
proc nlmixed data=frailty4 qpoints=5;
parms r01=1 r02=1 r03=1 r04=1 r05=1 r06=1 r07=1 r08=1 r09=1 r10=1 beta1=1 beta2=-1  theta=1;
bounds r01 r02 r03 r04 r05 r06 r07 r08 r09 r10 theta >=0;

base_haz=r01 * event1 + r02 * event2 + r03 * event3 + r04 * event4 + r05 * event5 + 
r06 * event6 + r07 * event7 + r08* event8 +r09 * event9 + r10 * event10;
cum_base_haz=r01 * dur1 + r02 * dur2 + r03 * dur3 + r04 * dur4 + r05 * dur5 + r06 * dur6 + 
r07 * dur7 + r08* dur8 +r09 * dur9 + r10 * dur10;

mu=  beta1 * x1 + beta2 * x2 + a;

loglik2=-exp(mu) * cum_base_haz;

if status=1 then loglik= log(base_haz) + mu + loglik2; 	/*log likelihood for failure */
if status=0 then loglik=loglik2;						/*log likelihood for censoring */
model followup ~ general(loglik);
random a ~ normal(0, theta) subject=center;
ods output ParameterEstimates=est1 FitStatistics=fit1; 
run;
* Code for Gamma frailty with variance theta;
proc nlmixed data=frailty4 qpoints=10 noad;

parms r01=1 r02=1 r03=1 r04=1 r05=1 r06=1 r07=1 r08=1 r09=1 r10=1 beta1=1 beta2=-1 theta=1;
bounds r01 r02 r03 r04 r05 r06 r07 r08 r09 r10 theta >=0;

base_haz=r01 * event1 + r02 * event2 + r03 * event3 + r04 * event4 + r05 * event5 + 
r06 * event6 + r07 * event7 + r08* event8 +r09 * event9 + r10 * event10;
cum_base_haz=r01 * dur1 + r02 * dur2 + r03 * dur3 + r04 * dur4 + r05 * dur5 + r06 * dur6 + 
r07 * dur7 + r08* dur8 +r09 * dur9 + r10 * dur10;

/* Code from Kerrie et al. 2006 to generate gamma random number */

p=cdf('NORMAL', a);
if p > .999999 then p=.999999;
g2=quantile('GAMMA', p, 1/theta); 
g=g2 * theta;		/* g is gamma distributed with mean 1 and variance theta */

mu=  beta1 * x1 + beta2 * x2 + log(g);

loglik2=-exp(mu) * cum_base_haz;

if status=1 then loglik= log(base_haz) + mu + loglik2; 	/*log likelihood for failure */
if status=0 then loglik=loglik2;						/*log likelihood for censoring */
model followup ~ general(loglik);
random a ~ normal(0, 1) subject=center;
ods output ParameterEstimates=est2 FitStatistics=fit2; 
run;
