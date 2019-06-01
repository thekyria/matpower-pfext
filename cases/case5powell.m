function mpc = powell5
%POWELL9    Power flow data for 5 bus, 1 slack and 4 PQ buses case.
%   Please see CASEFORMAT for details on the case file format.
%
%   Based on data from Lynn Powell's book "Power system load flow
%   analysis", p. 24.

%   This file is part of MATPOWER-PFEXT
%   $Id: powell5.m,v 1.00 2013/01/30$
%   by Theodoros Kyriakidis, EPFL

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
  1 3 0  0  0 0 1 1 0 132 1 1.1 0.9;
  2 1 40 20 0 0 1 1 0 132 1 1.1 0.9;
  3 1 25 15 0 0 1 1 0 132 1 1.1 0.9;
  4 1 40 20 0 0 1 1 0 132 1 1.1 0.9;
  5 1 50 20 0 0 1 1 0 132 1 1.1 0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
  1	0	0	999	-999 1 100	1	999	0	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
  1 2 0.05 0.11 0.02 999 999 999 0 0 1 -360 360;
  1 3 0.05 0.11 0.02 999 999 999 0 0 1 -360 360;
  1 5 0.03 0.08 0.02 999 999 999 0 0 1 -360 360;
  2 3 0.04 0.09 0.02 999 999 999 0 0 1 -360 360;
  2 5 0.04 0.09 0.02 999 999 999 0 0 1 -360 360;
  3 4 0.06 0.13 0.03 999 999 999 0 0 1 -360 360;
  4 5 0.04 0.09 0.02 999 999 999 0 0 1 -360 360;
];

%%-----  OPF Data  -----%%
%% area data
%	area	refbus
mpc.areas = [
	1	1;
];

%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
  2	1500	0	3	0.11	5	150;
];

end
