function [p] = rateOfConvergence(V, i, Vsol)

%% declare persistent variables (static in the c++ sense :P)
persistent Vsol_
persistent ekm1
persistent ek
persistent ekp1

%% Initialize
if nargin == 3 && i == 0
  Vsol_ = Vsol;
  ekp1 = V - Vsol; % error (k+1 iter, current)
  ek   = ekp1;     % error (k iter, previous)
  ekm1 = ekp1;     % error (k-1 iter, second to last)
end

%% update ek's
ekm1 = ek;
ek   = ekp1;
ekp1 = V - Vsol_;

%% calculate order of convergence
if i >= 2
  norm_ekm1 = norm(ekm1,inf);
  norm_ek   = norm(ek,inf);
  norm_ekp1 = norm(ekp1, inf);
  p = (log(norm_ekp1)-log(norm_ek))/(log(norm_ek)-log(norm_ekm1));
else
  p = NaN;
end

end
