function mpc = case5laurent
%CASE5LAURENT

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	0	  0	  0	0	1	1	0	230	1	1.05	0.95;
	2	2	0	  0	  0	0	1	1	0	230	1	1.05	0.95;
	3	2	0	  0	  0	0	1	1	0	230	1	1.05	0.95;
	4	2	0	  0	  0	0	1	1	0	230	1	1.05	0.95;
	5	1	250	80	0	0	1	1	0	230	1	1.05	0.95;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	0	0	900	-900	1.0	100	1	999	0	0	0	0	0	0	0	0	0	0	0	0;
	2	70	25	900	-900	1.0	100	1	999	0	0	0	0	0	0	0	0	0	0	0	0;
	3	110	30	900	-900	1.0	100	1	999	0	0	0	0	0	0	0	0	0	0	0	0;
  4	70	25	900	-900	1.0	100	1	999	0	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	2	0.01	0.10	0	999	999	999	0	0	1	-360	360;
  2	3	0.01	0.18	0	999	999	999	0	0	1	-360	360;
  3	4	0.01	0.10	0	999	999	999	0	0	1	-360	360;
  4	5	0.01	0.10	0	999	999	999	0	0	1	-360	360;
  5	2	0.01	0.10	0	999	999	999	0	0	1	-360	360;
];

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	0	0	3	0.00533	11.669	213.1;
	2	0	0	3	0.00889	10.333	200;
	2	0	0	3	0.00741	10.833	240;
	2	0	0	3	0.00741	10.833	240;
];
