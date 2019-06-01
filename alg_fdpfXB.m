function [V, success, i, FNorm, eNorm, p, mu] = alg_fdpfXB(baseMVA, bus, gen, branch, opt, optAn, Vsol)

%% options
term_crit = opt{1};
tolFP     = opt{2};
tolFQ     = opt{3};
tolP      = opt{4}; % used also as tol for norm(V)
tolQ      = opt{5};
max_it    = opt{6};
vm1       = opt{7};

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
% build FDPF B-matrices (XB version)
[Bp, Bpp] = makeB_XB(baseMVA, bus, branch);
% compute complex bus power injections (generation - load)
Sbus = makeSbus(baseMVA, bus, gen);
% bus initial voltage
V0 = getInitialVoltage(bus, gen);
% get bus index lists of each type of bus
[ref, pv, pq] = bustypes(bus, gen);

success = 0;   % success flag
i = 0;         % iteration counter
doPiter = true;
doQiter = true;

V = V0;        % bus voltage
Vnew = V;      % new bus voltage
dV = Vnew - V; % voltage correction
Va = angle(V); % bus voltage angle
Vm = abs(V);   % bus voltage magnitude

%% set up indexing for updating V
npv = length(pv);
npq = length(pq);
j1 = 1;         j2 = npv;           % j1:j2 - V angle of pv buses
j3 = j2 + 1;    j4 = j2 + npq;      % j3:j4 - V angle of pq buses
j5 = j4 + 1;    j6 = j4 + npq;      % j5:j6 - V mag of pq buses

%% evaluate initial mismatch
Ibus = Ybus * V;
S = V .* conj(Ibus);
Smis = S - Sbus;
Pmis = real(Smis);
Qmis = imag(Smis);

%% check tolerance
normP = norm(Pmis([pv; pq]), inf);
normP0 = normP;        % active power residual of the initial guess
normPp = normP;        % active power residual of the previous iteration
normQ = norm(Qmis(pq), inf);
normQ0 = normQ;        % reactive power residual of the initial guess
normQp = normQ;        % reactive power residual of the previous iteration
if normP < tolFP && normQ < tolFQ
  success = 1;
end

%% analysis: initialize extra output arguments 
if residual
  FNorm(i+1,1) = normP;
  FNorm(i+1,2) = normQ;
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

%% factor (reduced) B matrices
[Lp,  Up,  Pp ] = lu(Bp([pv; pq], [pv; pq]));
[Lpp, Upp, Ppp] = lu(Bpp(pq, pq));

%% do P and Q iterations
while (~success && i < max_it)
  % update iteration counter
  i = i + 1;
  
  %% do P iteration
  if doPiter
    %% analysis: evaluate true Jacobian
    [dSbus_dVmT, dSbus_dVaT] = dSbus_dV(Ybus, V);
    J11T = real(dSbus_dVaT([pv; pq], [pv; pq]));
    J12T = real(dSbus_dVmT([pv; pq], pq));
    J21T = imag(dSbus_dVaT(pq, [pv; pq]));
    J22T = imag(dSbus_dVmT(pq, pq));
    JT = [ J11T J12T;
           J21T J22T; ];
    % factor Jacobian
    [LJT, UJT, PJT] = lu(JT);
    
    %% update Voltage
    % update voltage angle
    Pmis_dV = Pmis;
    if ~vm1
      Pmis_dV = Pmis_dV./Vm; % Pmis divided by Voltage magnitude
    end
    dVa = -( Up \  (Lp \ (Pp * Pmis_dV([pv; pq]))));

    %% update voltage
    Va([pv; pq]) = Va([pv; pq]) + dVa;
    dx(j1:j4) = dVa;
    dx(j5:j6) = zeros(npq,1);
    Vnew = Vm .* exp(1j * Va);
    dV = Vnew - V;
    V = Vnew;
    % update Vm and Va again in case we wrapped around with a negative Vm
    Vm = abs(V);
    Va = angle(V);

    %% evalute mismatch
    Ibus = Ybus * V;
    S = V .* conj(Ibus);
    Smis = S - Sbus;
    Pmis = real(Smis);
    Qmis = imag(Smis);
    % analysis: calculate F
    F = [ Pmis(pv);
          Pmis(pq);
          Qmis(pq) ];

    %% check tolerance
    normP = norm(Pmis([pv; pq]), inf);
    normQ = norm(Qmis(pq), inf);
    normdV = norm(dV, inf);
    normV  = norm(V, inf);
    switch term_crit
      case 1 % absolute residual
        if normP < tolP && normQ < tolQ
          success = 1;
          doQiter = false;
        end 
      case 2 % relative residual against 0
        if normP/normP0 < tolP && normQ/normQ0 < tolQ
          success = 1;
          doQiter = false;
        end
      case 3 % relative residual against i-1
        if normP/normPp < tolP && normQ/normQp < tolQ
          success = 1;
          doQiter = false;
        end
      case 4 % absolute correction
        if normdV < tolP
          success = 1;
          doQiter = false;
        end
      case 5 % relative correction (against i-1)
        if normdV/normV < tolP
          success = 1;
          doQiter = false;
        end
    end
    normPp = normP;
    normQp = normQ;
    
    %% analysis
    % residual error norm
    if residual
      FNorm(i+1,1) = normP;
      FNorm(i+1,2) = normQ;
    end
    % error norm
    if errorNorm
      eNorm(i+1,1) = norm(V-Vsol,inf);
    end
    % rate of convergence
    if rate
      p(i+1,1) = rateOfConvergence(V,i); 
    end
    % convergence monitor
    if conv_monitor
      dxtemp  = -(UJT \ (LJT \ (PJT * F)));
      mu(i+1,1) = norm(dxtemp) / norm(dx);
    end
  end
  
  %% do Q iteration
  if doQiter
    %% analysis: evaluate true Jacobian
    [dSbus_dVmT, dSbus_dVaT] = dSbus_dV(Ybus, V);
    J11T = real(dSbus_dVaT([pv; pq], [pv; pq]));
    J12T = real(dSbus_dVmT([pv; pq], pq));
    J21T = imag(dSbus_dVaT(pq, [pv; pq]));
    J22T = imag(dSbus_dVmT(pq, pq));
    JT = [ J11T J12T;
           J21T J22T; ];
    % factor Jacobian
    [LJT, UJT, PJT] = lu(JT);
    
    %% update Voltage 
    % update voltage magnitude
    Qmis_dV = Qmis;
    if ~vm1
      Qmis_dV = Qmis_dV./Vm; % Qmis divided by Voltage norm
    end
    dVm = -(Upp \ (Lpp \ (Ppp * Qmis_dV(pq))));

    % update voltage
    Vm(pq) = Vm(pq) + dVm;
    dx(j1:j4) = zeros(npv+npq,1);
    dx(j5:j6) = dVm;
    Vnew = Vm .* exp(1j * Va);
    dV = Vnew - V;
    V = Vnew;
    % update Vm and Va again in case we wrapped around with a negative Vm
    Vm = abs(V);
    Va = angle(V);

    %% evalute mismatch
    Ibus = Ybus * V;
    S = V .* conj(Ibus);
    Smis = S - Sbus;
    Pmis = real(Smis);
    Qmis = imag(Smis);
    % analysis: calculate F
    F = [ Pmis(pv);
          Pmis(pq);
          Qmis(pq) ];

    %% check tolerance
    normP = norm(Pmis([pv; pq]), inf);
    normQ = norm(Qmis(pq), inf);
    normdV = norm(dV, inf);
    normV  = norm(V, inf);
    switch term_crit
      case 1 % absolute residual
        if normP < tolP && normQ < tolQ
          success = 1;
          doQiter = false;
        end 
      case 2 % relative residual against 0
        if normP/normP0 < tolP && normQ/normQ0 < tolQ
          success = 1;
          doQiter = false;
        end
      case 3 % relative residual against i-1
        if normP/normPp < tolP && normQ/normQp < tolQ
          success = 1;
          doQiter = false;
        end
      case 4 % absolute correction
        if normdV < tolP
          success = 1;
          doQiter = false;
        end
      case 5 % relative correction (against i-1)
        if normdV/normV < tolP
          success = 1;
          doQiter = false;
        end
    end
    normPp = normP;
    normQp = normQ;
    
    %% analysis
    % residual error norm
    if residual
      FNorm(i+1,3) = normP;
      FNorm(i+1,4) = normQ;
    end
    % error norm
    if errorNorm
      eNorm(i+1,2) = norm(V-Vsol,inf);
    end
    % rate of convergence
    if rate
      p(i+1,2) = rateOfConvergence(V,i); 
    end
    % convergence monitor
    if conv_monitor
      dxtemp  = -(UJT \ (LJT \ (PJT * F)));
      mu(i+1,2) = norm(dxtemp) / norm(dx);
    end
  end
end

end
