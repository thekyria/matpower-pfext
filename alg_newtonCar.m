function [V, success, i, FNorm, eNorm, p, mu] = alg_newtonCar(baseMVA, bus, gen, branch, opt, optAn, Vsol)

%% options
term_crit  = opt{1};
tolFS      = opt{2};
tolFVm     = opt{3};
tolS       = opt{4};
tolVm      = opt{5};
max_it     = opt{6};
retain_jac = opt{7};

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

success = 0;
i = 0;
jaci = 0; % shows how many more iterations the Jacobian will be retained
V = V0;
Vnew = V;
Vm = abs(V);
Vmsp = abs(V); % setpoint Vm

%% set up indexing for updating V
npv = length(pv);
npq = length(pq);
j1 = 1;         j2 = npv + npq;      %% j1:j2 - real voltage (E) of pv & pq buses
j3 = j2 + 1;    j4 = j2 + npv + npq; %% j3:j4 - imag voltage (F) of pv & pq buses

%% evaluate F(x0)
% power mismatch
Ibus = Ybus * V;
S = V .* conj(Ibus);
Smis = S - Sbus;
Pmis = real(Smis);
Qmis = imag(Smis);
F = [ Pmis(pv);
      Pmis(pq);
      Qmis(pq)  ];
% voltage magnitude mismatch
Vmmis = Vm.^2 - Vmsp.^2;
FVm = [ Vmmis(pv) ];

%% check tolerance
normF = norm(F, inf);
normF0 = normF;
normFp = normF;
normFVm = norm(FVm, inf);
normFVm0 = normFVm;
normFVmp = normFVm;
if normF < tolFS && normFVm < tolFVm
  success = 1;
end

%% analysis: initialize extra output arguments 
if residual
  FNorm(i+1,1) = normF;
  FNorm(i+1,2) = normFVm;
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

%% do cartesian Newton iterations
while (~success && i < max_it)

  %% evaluate Jacobian
  if jaci == 0
    [~, ~, dSbus_dE, dSbus_dF] = dSbus_dV(Ybus, V);  % power Jacobian
    [dVm2_dE, dVm2_dF] = dVm2_dV(V, issparse(Ybus)); % Vm^2 Jacobian
    J11 = real(dSbus_dE([pv; pq], [pv; pq])); % dP/dE
    J12 = real(dSbus_dF([pv; pq], [pv; pq])); % dP/dF
    J21 = imag(dSbus_dE(pq, [pv; pq]));       % dQ/dE
    J22 = imag(dSbus_dF(pq, [pv; pq]));       % dQ/dF
    JE  = dVm2_dE(pv, [pv; pq]);              % dVm^2/dE
    JF  = dVm2_dF(pv, [pv; pq]);              % dVm^2/dF
    J = [ J11, J12;
          J21, J22;
          JE,  JF];
    % factor Jacobian
    [LJ, UJ, PJ] = lu(J);
    % reset jaci
    jaci = retain_jac;    
  end
  
  %% analysis: evaluate true Jacobian
  [~, ~, dSbus_dET, dSbus_dFT] = dSbus_dV(Ybus, V);  % power Jacobian
  [dVm2_dET, dVm2_dFT] = dVm2_dV(V, issparse(Ybus)); % Vm^2 Jacobian
  J11T = real(dSbus_dET([pv; pq], [pv; pq])); % dP/dE
  J12T = real(dSbus_dFT([pv; pq], [pv; pq])); % dP/dF
  J21T = imag(dSbus_dET(pq, [pv; pq]));       % dQ/dE
  J22T = imag(dSbus_dFT(pq, [pv; pq]));       % dQ/dF
  JET  = dVm2_dET(pv, [pv; pq]);              % dVm^2/dE
  JFT  = dVm2_dFT(pv, [pv; pq]);              % dVm^2/dF
  JT = [ J11T, J12T;
         J21T, J22T;
         JET,  JFT];
  % factor Jacobian
  [LJT, UJT, PJT] = lu(JT);

  %% compute update step
  dx = -(UJ \ (LJ \ (PJ * [F; FVm])));
  % update iteration counter
  i = i + 1;
  % decrease jaci by one (Jacobian has been used this iteration)
  jaci = jaci-1;

  %% update voltage
  Vnew([pv; pq]) = Vnew([pv; pq]) + dx(j1:j2);
  Vnew([pv; pq]) = Vnew([pv; pq]) + 1j* dx(j3:j4);
  dV = Vnew - V;
  V = Vnew;
  % update Vm
  Vm = abs(V);

  %% evalute F(x)
  % power mismatch
  Ibus = Ybus * V;
  S = V .* conj(Ibus);
  Smis = S - Sbus;
  Pmis = real(Smis);
  Qmis = imag(Smis);
  F = [ Pmis(pv);
        Pmis(pq);
        Qmis(pq)  ];
  % voltage magnitude mismatch
  Vmmis = Vm.^2 - Vmsp.^2;
  FVm = [ Vmmis(pv) ];

  %% check for convergence
  normF  = norm(F, inf);
  normFVm = norm(FVm, inf);
  normdV = norm(dV, inf);
  normV  = norm(V, inf);
  switch term_crit
    case 1 % absolute residual
      if normF < tolS && normFVm < tolVm
        success = 1;
      end 
    case 2 % relative residual against 0
      if normF/normF0 < tolS && normFVm/normFVm0 < tolVm
        success = 1;
      end
    case 3 % relative residual against i-1
      if normF/normFp < tolS && normFVm/normFVmp < tolVm
        success = 1;
      end
    case 4 % absolute correction
      if normdV < tolVm
        success = 1;
      end
    case 5 % relative correction (against i-1)
      if normdV/normV < tolVm
        success = 1;
      end
  end
  normFp = normF;
  normFVmp = normFVm;
  
  %% analysis
  % residual error norm
  if residual
    FNorm(1,i+1) = normF;
    FNorm(2,i+1) = normFVm;
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
    dxtemp  = -(UJT \ (LJT \ (PJT * [F; FVm])));
    mu(i+1) = norm(dxtemp) / norm(dx);
  end
end

end
