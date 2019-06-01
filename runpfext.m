function [res] = runpfext(casedata, alg, opt, optAn)

%% named indices
% branch matrix 
PF = 14;
QF = 15;
PT = 16;
QT = 17;
% gen matrix 
PG = 2;
QG = 3;

%% initialize
% read data
mpc = loadcase(casedata);
% add zero columns to branch for flows if needed
if size(mpc.branch,2) < QT
  mpc.branch = [ mpc.branch zeros(size(mpc.branch, 1), QT-size(mpc.branch,2)) ];
end
% convert to internal indexing
mpc = ext2int(mpc);
[baseMVA, bus, gen, branch] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch);

%% prepare algorithm bank
alg_bank  = { ...
  @alg_newton;      % 1  {term_crit tolF tol max_iter retain_jac}
  @alg_gauss;       % 2  {term_crit tolF tol max_iter sor_acceleration}
  @alg_jacobi;      % 3  {term_crit tolF tol max_iter sor_acceleration}
  @alg_dNewton;     % 4  {term_crit tolF tol max_iter retain_jac}
  @alg_newtonCar;   % 5  {term_crit tolFS tolFVm tolS tolVm max_iter retain_jac}
  @alg_dcpf;        % 6  {}
  @alg_dcpfvm;      % 7  {formulation approximation}
  @alg_zgauss;      % 8  {term_crit tolF tol max_iter sor_acceleration}
  @alg_zjacobi;     % 9  {term_crit tolF tol max_iter sor_acceleration}
  @alg_fdpfBX;      % 10 {term_crit tolFP tolFQ tolP tolQ max_iter vm1}
  @alg_fdpfXB;      % 11 {term_crit tolFP tolFQ tolP tolQ max_iter vm1}
  @alg_fdpfICP;     % 12 {term_crit tolFP tolFQ tolP tolQ max_iter fillin}
  @alg_fdpfICD;     % 13 {term_crit tolFP tolFQ tolP tolQ max_iter fillin}
  @alg_halley;      % 14 {term_crit tolF tol max_iter}
  @alg_lanz;        % 15 {term_crit tolF tol max_iter beta1 beta2}
  @alg_icNewton;    % 16 {term_crit tolFP tolFQ tolP tolQ max_iter fillinH fillinL}
  @dummy; };

%% determine a solution for the problem
optSol = {1, 10^-12, 10^-12, 100, 1};
[Vsol, successSol, ~] = alg_newton(baseMVA, bus, gen, branch, optSol);
% TODO: if solution fails do something
if ~successSol
  error('Failed!');
end

%% execute
t0 = clock;
% run algorithm proper
[V, success, iterations, FNorm, eNorm, p, mu] = alg_bank{alg}(baseMVA, bus, gen, branch, opt, optAn, Vsol);
% update solution matrices
[bus, gen, branch] = updateSolution(baseMVA, bus, gen, branch, V);

%% results
% retrieve timing and success status
mpc.et = etime(clock, t0);
mpc.success = success;
mpc.iter = iterations;
mpc.FNorm = FNorm;
mpc.eNorm = eNorm;
mpc.p = p;
mpc.mu = mu;

% convert back to original bus numbering & print results
[mpc.bus, mpc.gen, mpc.branch] = deal(bus, gen, branch);
res = int2ext(mpc);

% zero out result fields of out-of-service gens & branches
if ~isempty(res.order.gen.status.off)
  res.gen(res.order.gen.status.off, [PG QG]) = 0;
end
if ~isempty(res.order.branch.status.off)
  res.branch(res.order.branch.status.off, [PF QF PT QT]) = 0;
end

end
