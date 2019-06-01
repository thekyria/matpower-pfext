function [Bp, Bpp] = makeB_BX(baseMVA, bus, branch)
%MAKEB_BX   Builds the FDPF matrices, B prime and B double prime.
%   [BP, BPP] = MAKEB_BX(BASEMVA, BUS, BRANCH, ALG) returns the two
%   matrices B prime and B double prime used in the fast decoupled power
%   flow.
%
%   Example:
%       [Bp, Bpp] = makeB_BX(baseMVA, bus, branch);
%
%   This file is part of MATPOWER-PFEXT
%   by Theodoros Kyriakidis, EPFL

%% constants
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of lines

%% define named indices into bus, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%%-----  form Bp (B prime)  -----
temp_branch = branch;                       %% modify a copy of branch
temp_bus = bus;                             %% modify a copy of bus
temp_bus(:, BS) = zeros(nb, 1);             %% zero out shunts at buses
temp_branch(:, BR_B) = zeros(nl, 1);        %% zero out line charging shunts
temp_branch(:, TAP) = ones(nl, 1);          %% cancel out taps
Bp = -imag( makeYbus(baseMVA, temp_bus, temp_branch) );

%%-----  form Bpp (B double prime)  -----
temp_branch = branch;                       %% modify a copy of branch
temp_branch(:, SHIFT) = zeros(nl, 1);       %% zero out phase shifters
temp_branch(:, BR_R) = zeros(nl, 1);        %% zero out line resistance
Bpp = -imag( makeYbus(baseMVA, bus, temp_branch) );

end