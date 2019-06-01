function [V, success, i, FNorm, eNorm, p, mu] = alg_dcpf(baseMVA, bus, gen, branch, opt, optAn, Vsol)

%% analysis options
if nargin < 6
  residual  = false;
  errorNorm = false;
else
  residual  = optAn{1};
  errorNorm = optAn{2};
%   rate % NOT SUPPORTED
%   conv_monitor % NOT SUPPORTED
end

%% initialize
% named indices
GS = 5;
% build B matrices and phase shift injections
[B, ~, Pbusinj, ~] = makeB_dc(bus, branch);
% compute complex bus power injections (generation - load)
% adjusted for phase shifters and real shunts
Pbus = real(makeSbus(baseMVA, bus, gen)) - Pbusinj - bus(:, GS) / baseMVA;
% bus initial voltage
V0 = getInitialVoltage(bus, gen);
% get bus index lists of each type of bus
[ref, pv, pq] = bustypes(bus, gen);

i = 0;
success = true; % DC loadflow is always "successful"
V = V0;         % bus voltage
Va = angle(V);  % bus voltage angle
Vm = abs(V);    % bus voltage magnitude

%% evaluate F
% build admittance matrices
[Ybus, ~, ~] = makeYbus(baseMVA, bus, branch);
% compute complex bus power injections (generation - load)
Sbus = makeSbus(baseMVA, bus, gen);
% evaluate F proper
Ibus = Ybus * V;
S = V .* conj(Ibus);
Smis = S - Sbus;
Pmis = real(Smis);
Qmis = imag(Smis);
F = [ Pmis(pv);
      Pmis(pq);
      Qmis(pq)  ];
normF  = norm(F, inf); % residual

%% analysis output arguments 
if residual
  FNorm(i+1) = normF;
end
if errorNorm
  eNorm(i+1) = norm(V-Vsol,inf);
end
p  = [];
mu = [];

%% calculate angles for non-reference buses
Va([pv; pq]) = B([pv; pq],[pv; pq]) \ (Pbus([pv; pq]) - B([pv; pq], ref)*Va(ref));

%% update output arguments
V = Vm .* exp(1j * Va);
i = 1;          % DCPF has no iterations

%% evaluate F
Ibus = Ybus * V;
S = V .* conj(Ibus);
Smis = S - Sbus;
Pmis = real(Smis);
Qmis = imag(Smis);
F = [ Pmis(pv);
      Pmis(pq);
      Qmis(pq)  ];
normF  = norm(F, inf); % residual

%% analysis output arguments 
if residual
  FNorm(i+1) = normF;
end
if errorNorm
  eNorm(i+1) = norm(V-Vsol,inf);
end

end
