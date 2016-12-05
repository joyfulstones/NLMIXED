* Event_type=1: local recurrent event, Event_type=2: distant recurrent event;
* Event_type=3: death, Event_type=0: independent censoring;
data one;
input id stoptime trt event_type;
aa=1;
cards;
1 5.0 1 1	/* local recurrent event at time 5.0 */
1 6.7 1 1	/* local recurrent event at time 6.7 */
1 3.0 1 2	/* distant recurrent event at time 3.0 */
1 5.4 1 2	/* distant recurrent event at time 5.4 */
1 10.0 1 3	/* death event at time 10.0 */
2 7.9 0 0	/* independent censoring event at time 7.9 */
3 4.5 1 1
3 6.0 1 1
3 5.0 1 2
3 6.0 1 4
......
;
run;
* dataset two includes only local recurrent events;
data two;
set one;
if event_type=1;	/* Keep all local recurrent events */
run;
* Get the quantiles for the recurrent events, store in dataset quant_r;
proc univariate data=two noprint;
var stoptime; 
output out=quant_rl pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=qrl; 
run;
data quant_rl;
set quant_rl;
aa=1;
run;
* dataset two includes only distant recurrent events;
data three;
set one;
if event_type=2;	/* Keep all distant recurrent events */
run;
* Get the quantiles for the recurrent events, store in dataset quant_r;
proc univariate data=three noprint;
var stoptime; 
output out=quant_rd pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=qrd; 
run;
data quant_rd;
set quant_rd;
aa=1;
run;
* Get the quantiles for the death/censoring event;
data four;
set one;
if event_type=3 or event_type=4;	
run;
* Get the quantiles for the death event, store in dataset quant_d;
proc univariate data=four noprint;
var stoptime; 
output out=quant_d pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=qd; 
where event_type=3;
run;
data quant_d;
set quant_d;
aa=1;
run;
* Merge data with the quantiles;
data five;
merge one quant_rl quant_rd quant_d;
by aa;
run;

* Calculate the duration in each quantile interval, together with the indicator of event in each interval;
data six;
set five;
array quant_rl {11} qrl0 qr10 qr20 qr30 qr40 qr50 qr60 qr70 qr80 qr90 qr100;
array quant_rd {11} qrd0 qr10 qr20 qr30 qr40 qr50 qr60 qr70 qr80 qr90 qr100;
array quant_d {11} qd0 qr10 qr20 qr30 qr40 qr50 qr60 qr70 qr80 qr90 qr100;

array dur_rl {10} dur_rl1-dur_rl10;
array dur_rd {10} dur_rd1-dur_rd10;
array dur_d {10} dur_d1-dur_d10;

array event_rl {10} event_rl1-event_rl10;
array event_rd {10} event_rd1-event_rd10;
array event_d {10} event_d1-event_d10;

do i=1 to 10;
	dur_rl{i}=0;
	dur_rd{i}=0;
	dur_d{i}=0;
	event_rl{i}=0;
	event_rd{i}=0;
	event_d{i}=0;
end;

* For local recurrent event;
if event_type=1 then do;
	do i=2 to 11;
		if stoptime<=quant_rl{i} then do;
			event_rl{i-1}=1;		/* indicator of local recurrent event in each interval */
			i=11;
		end;
	end;
end;

* For distant recurrent event;
if event_type=2 then do;
	do i=2 to 11;
		if stoptime<=quant_rd{i} then do;
			event_rd{i-1}=1;		/* indicator of distant recurrent event in each interval */
			i=11;
		end;
	end;
end;

else do; /* If death or censored observation */
	do i=2 to 11;	/* duration in each local recurrent event quantiles */
		if stoptime<=quant_rl{i} then do;
			dur_rl{i-1}=stoptime-quant_rl{i-1};
			i=11;
		end;
		else dur_rl{i-1}=quant_rl{i}-quant_rl{i-1};
	end;

	do i=2 to 11;	/* duration in each distant recurrent event quantiles */
		if stoptime<=quant_rd{i} then do;
			dur_rd{i-1}=stoptime-quant_rd{i-1};
			i=11;
		end;
		else dur_rd{i-1}=quant_rd{i}-quant_rd{i-1};
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
* This is the code for Normal frailty;
proc nlmixed data=five qpoints=30 noad;

parms rl1=1 rl2=1 rl3=1 rl4=1 rl5=1 rl6=1 rl7=1 rl8=1 rl9=1 rl10=1 /* Starting values for local recurrent baseline hazard */
	  rd1=1 rd2=1 rd3=1 rd4=1 rd5=1 rd6=1 rd7=1 rd8=1 rd9=1 rd10=1 /* Starting values for distant recurrent baseline hazard */
	  h01=1 h02=1 h03=1 h04=1 h05=1 h06=1 h07=1 h08=1 h09=1 h10=1 /* Starting values for death baseline hazard */ 
	  beta1=0 beta2=0 sigma11=.2 sigma12=0 sigma22=.2 alpha1=-.2 gamma1=0 gamma2=0;
bounds rl1 rl2 rl3 rl4 rl5 rl6 rl7 rl8 rl9 rl10 rd1 rd2 rd3 rd4 rd5 rd6 rd7 rd8 rd9
h01 h02 h03 h04 h05 h06 h07 h08 h09 h10 theta11 theta22 >=0;

* base hazard and cumulative baseline hazard for local recurrent events;

base_haz_rl=rl1 * event_rl1 + rl2 * event_rl2 + rl3 * event_rl3 + rl4 * event_rl4 + rl5 * event_rl5 + 
rl6 * event_rl6 + rl7 * event_rl7 + rl8* event_rl8 +rl9 * event_rl9 + rl10 * event_rl10;
cum_haz_rl=rl1 * dur_rl1 + rl2 * dur_rl2 + rl3 * dur_rl3 + rl4 * dur_rl4 + rl5 * dur_rl5 + 
rl6 * dur_rl6 + rl7 * dur_rl7 + rl8* dur_rl8 +rl9 * dur_rl9 + rl10 * dur_rl10;

* base hazard and cumulative baseline hazard for distant recurrent events;

base_haz_rd=rd1 * event_rd1 + rd2 * event_rd2 + rd3 * event_rd3 + rd4 * event_rd4 + rd5 * event_rd5 + 
rd6 * event_rd6 + rd7 * event_rd7 + rd8* event_rd8 +rd9 * event_rd9 + rd10 * event_rd10;
cum_haz_rd=rd1 * dur_rd1 + rd2 * dur_rd2 + rd3 * dur_rd3 + rd4 * dur_rd4 + rd5 * dur_rd5 + 
rd6 * dur_rd6 + rd7 * dur_rd7 + rd8* dur_rd8 +rd9 * dur_rd9 + rd10 * dur_rd10;

* base hazard and cumulative baseline hazard for death events;

base_haz_d=h01 * event_d1 + h02 * event_d2 + h03 * event_d3 + h04 * event_d4 + h05 * event_d5 + 
h06 * event_d6 + h07 * event_d7 + h08* event_d8 +h09 * event_d9 + h10 * event_d10;
cum_base_haz_d=h01 * dur_d1 + h02 * dur_d2 + h03 * dur_d3 + h04 * dur_d4 + h05 * dur_d5 + h06 * dur_d6 + 
h07 * dur_d7 + h08* dur_d8 +h09 * dur_d9 + h10 * dur_d10;

mu1= beta1 * trt +  a;			/* for local recurrent event */

mu2= beta2 * trt +  b;			/* for distant recurrent event */

mu2= alpha1 * trt  + gamma1 * a + gamma2 *b;	/* for death event */

loglik1=-exp(mu1) * cum_base_haz_rl;

loglik2=-exp(mu2) * cum_base_haz_rd;

loglik2=-exp(mu2) * cum_base_haz_d;

if event_type=1 then loglik= log(base_haz_rl) + mu1 ; 					/*log likelihood for recurrent event */
if event_type=2 then loglik= log(base_haz_rd) + mu2 ; 					/*log likelihood for recurrent event */
if event_type=3 then loglik= loglik1 + loglik2 + log(base_haz_d) + mu + loglik3; 	/*log likelihood for recurrent event *//*log likelihood for death */
if event_type=4 then loglik= loglik1 + loglik2 +loglik3;							/*log likelihood for censoring */

model stoptime ~ general(loglik);
random a b~ normal([0,0], [theta11, theta21, theta22]) subject=id;
ods output ParameterEstimates=est1 FitStatistics=fit1; 
run;
