function [V, success, i, FNorm, eNorm, p, mu] = alg_gauss(baseMVA, bus, gen, branch, opt, optAn, Vsol)

%% options
term_crit  = opt{1};
tolF       = opt{2};
tol        = opt{3};
max_it     = opt{4};
alpha      = opt{5}; % SOR acceleration factor

%% analysis options
if nargin < 6
  residual  = false;
  errorNorm = false;
  rate      = false;
  conv_monitor = false;
else
  residual  = optAn{1};
  errorNorm = optAn{2};
  rate      = optAn{3};
  conv_monitor = optAn{4};
end

%% initialize
% build admittance matrices
[Ybus, ~, ~] = makeYbus(baseMVA, bus, branch);
% compute complex bus power injections (generation - load)
Sbus = makeSbus(baseMVA, bus, gen);
% bus initial voltage
V0 = getInitialVoltage(bus, gen);
% get bus index lists of each type of bus
[ref, pv, pq] = bustypes(bus, gen);

success = 0;   % success flag
i = 0;         % iteration counter
V = V0;        % bus voltage
Vnew = V;      % new bus voltage
dV = Vnew - V; % voltage correction
Vm = abs(V);   % bus voltage magnitude

%% set up indexing for updating V
npv = length(pv);
npq = length(pq);
% the following used only for analysis
j1 = 1;         j2 = npv;           % j1:j2 - V angle of pv buses
j3 = j2 + 1;    j4 = j2 + npq;      % j3:j4 - V angle of pq buses
j5 = j4 + 1;    j6 = j4 + npq;      % j5:j6 - V mag of pq buses

%% evaluate F(x0)
Ibus = Ybus * V;
S    = V .* conj(Ibus);
Smis = S - Sbus;
Pmis = real(Smis);
Qmis = imag(Smis);
F    = [ Pmis(pv);
         Pmis(pq);
         Qmis(pq)  ];

%% check tolerance
normF = norm(F, inf); % residual
normF0 = normF;       % residual of the initial guess
normFp = normF;       % residual of the previous iteration
if normF < tolF
  success = 1;
end

%% analysis: initialize extra output arguments 
if residual
  FNorm(i+1) = normF;
end
if errorNorm
  eNorm(i+1) = norm(V-Vsol,inf);
end
if rate
  p(i+1) = rateOfConvergence(V, 0, Vsol);
end
if conv_monitor
  mu(i+1) = NaN;
end

%% do Gauss-Seidel iterations
while (~success && i < max_it)
  % update iteration counter
  i = i + 1;
  
  %% analysis: evaluate true Jacobian
  [dSbus_dVmT, dSbus_dVaT] = dSbus_dV(Ybus, V);
  J11T = real(dSbus_dVaT([pv; pq], [pv; pq]));
  J12T = real(dSbus_dVmT([pv; pq], pq));
  J21T = imag(dSbus_dVaT(pq, [pv; pq]));
  J22T = imag(dSbus_dVmT(pq, pq));
  JT = [ J11T J12T;
         J21T J22T; ];
  % factor true Jacobian
  [LJT, UJT, PJT] = lu(JT);
  
  %% update voltage
  % analysis: updates could be performed directly to the V vector, but the
  % intermediate Vnew is used for analysis puproses
  % at PQ buses
  for k = pq(1:npq)'
    Vnew(k) = Vnew(k) + alpha*(conj(Sbus(k) / Vnew(k)) - Ybus(k,:) * Vnew ) / Ybus(k,k);
  end
  
  % at PV buses
  for k = pv(1:npv)'
    Sbus(k) = real(Sbus(k)) + 1j * imag( Vnew(k) .* conj(Ybus(k,:) * Vnew));
    Vnew(k) = Vnew(k) + alpha*(conj(Sbus(k) / Vnew(k)) - Ybus(k,:) * Vnew ) / Ybus(k,k);
  end
  % Correct voltage magnitude back to setpoint values for pv buses
  Vnew(pv) = Vm(pv) .* (Vnew(pv)./abs(Vnew(pv)));
  
  % analysis: calculate intermediate updates
  dV = Vnew - V;
  Va    = angle(V);
  Vnewa = angle(Vnew);
  Vm    = abs(V);
  Vnewm = abs(Vnew);
  dx(j1:j2) = Vnewa(pv) - Va(pv);
  dx(j3:j4) = Vnewa(pq) - Va(pq);
  dx(j5:j6) = Vnewm(pq) - Vm(pq);
  % actually update voltage
  V = Vnew;
  
  %% evalute F(x)
  Ibus = Ybus * V;
  S    = V .* conj(Ibus);
  Smis = S - Sbus;
  Pmis = real(Smis);
  Qmis = imag(Smis);
  F    = [ Pmis(pv);
           Pmis(pq);
           Qmis(pq)  ];
  
  %% check for convergence
  normF  = norm(F, inf);
  normdV = norm(dV, inf);
  normV  = norm(V, inf);
  switch term_crit
    case 1 % absolute residual
      if normF < tol
        success = 1;
      end 
    case 2 % relative residual against 0
      if normF/normF0 < tol
        success = 1;
      end
    case 3 % relative residual against i-1
      if normF/normFp < tol
        success = 1;
      end
    case 4 % absolute correction
      if normdV < tol
        success = 1;
      end
    case 5 % relative correction (against i-1)
      if normdV/normV < tol
        success = 1;
      end
  end
  normFp = normF;
  
  %% analysis
  % residual error norm
  if residual
    FNorm(i+1) = normF;
  end
  % error norm
  if errorNorm
    eNorm(i+1) = norm(V-Vsol,inf);
  end
  % rate of convergence
  if rate
    p(i+1) = rateOfConvergence(V,i); 
  end
  % convergence monitor
  if conv_monitor
    dxtemp  = -(UJT \ (LJT \ (PJT * F)));
    mu(i+1) = norm(dxtemp) / norm(dx);
  end
end

end
