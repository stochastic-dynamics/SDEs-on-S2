function [Df1] = Order1_diff(driftdiffusion_term, x, y, z)
%ORDER1_DIFF Computes the first-order partial derivatives of a symbolic expression
%
%   Df1 = Order1_diff(driftdiffusion_term, x, y, z)
%
%   Inputs:
%       driftdiffusion_term - symbolic scalar expression
%       x, y, z             - symbolic variables (e.g., syms x y z)
%
%   Output:
%       Df1 - 1x3 row vector of first-order partial derivatives:
%             [∂/∂x, ∂/∂y, ∂/∂z] of the input expression
%
%   This function corresponds to the first-order derivative computation 
%   in Eq. (4.14) of the manuscript.

Df1 = [diff(driftdiffusion_term, x), ...
       diff(driftdiffusion_term, y), ...
       diff(driftdiffusion_term, z)];
end
