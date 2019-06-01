function [dVm2_dE, dVm2_dF] = dVm2_dV(V, Vissparse)
%DVM2_DV   Computes partial derivatives of bus voltage magnitude squared
% w.r.t. voltage.
%   [DVM2_DE, DVM2_DF] = DVM2_DV(V, VISSPARSE) returns 2 matrices containng
%   partial derivatives of the bus voltage magnitude squared w.r.t real (E)
%   and imaginary (F) bus voltage parts (for all buses). If SPARSE is true
%   the return values will be sparse. The following explains the
%   expressions used to form the matrices:
%
%   size of vectors is nx1
%   V = Vm < Va (polar)
%   V = E + j * F (cartesian)
%   Vm = sqrt(diag(E) * E + diag(F) * F)
%
%   ---------- Voltage magnitude squared ----------
%   Partials of Vm^2 w.r.t. real (E) & imaginary (F) voltage part
%     dVm2/dE = 2*diag(E)
%     dVm2/dF = 2*diag(F)
%   These derivatives are useful for cartesian coordinates formulations of
%   the NR PF
%
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [dVm2_dE, dVm2_dF] = dVm2_dV(Ybus, V);

%   This file is part of MATPOWER-PFEXT
%   $Id: dVm2_dV.m,v 1.00 2013/02/01$
%   by Theodoros Kyriakidis, EPFL

n = length(V);

if Vissparse                   %% sparse version
  diagV = sparse(1:n, 1:n, V, n, n);
else                        %% dense version
  diagV = diag(V);
end

dVm2_dE   = 2*real(diagV);
dVm2_dF   = 2*imag(diagV);

end
