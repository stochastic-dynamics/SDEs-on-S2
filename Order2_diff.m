function [Df2] = Order2_diff(driftdiffusion_term, x, y, z)
%ORDER2_DIFF Computes the second-order partial derivatives (Hessian matrix)
%
%   Df2 = Order2_diff(driftdiffusion_term, x, y, z)
%
%   Inputs:
%       driftdiffusion_term - symbolic scalar expression
%       x, y, z             - symbolic variables (e.g., syms x y z)
%
%   Output:
%       Df2 - 3x3 Hessian matrix:
%             Df2(i,j) = ∂²(driftdiffusion_term)/∂var_i∂var_j
%             where var_i, var_j ∈ {x, y, z}
%
%   This function computes the second-order partial derivatives of the
%   input expression, corresponding to Eq. (4.14) of the manuscript.

% Compute first-order partial derivatives
Df = [diff(driftdiffusion_term, x), ...
      diff(driftdiffusion_term, y), ...
      diff(driftdiffusion_term, z)];

% Compute second-order partial derivatives (Hessian)
Df2 = [diff(Df(1), x), diff(Df(1), y), diff(Df(1), z);
       diff(Df(2), x), diff(Df(2), y), diff(Df(2), z);
       diff(Df(3), x), diff(Df(3), y), diff(Df(3), z)];

end
