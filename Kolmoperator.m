function [L0a, L0b, L1a, L1b, L0L0a, L0L0b, L0L1a, L1L0a, L1L1a, L1L0b, L0L1b, L1L1b] = Kolmoperator(drift_term, diffusion_term, a, b, z)
%KOLMOPERATOR Computes nested Kolmogorov operators for a given system
%
%   [L0a, L0b, L1a, L1b, L0L0a, L0L0b, L0L1a, L1L0a, L1L1a, L1L0b, L0L1b, L1L1b] = ...
%       Kolmoperator(drift_term, diffusion_term, a, b, z)
%
%   Inputs:
%       drift_term     - symbolic expression or function for drift
%       diffusion_term - symbolic expression or function for diffusion
%       a              - scalar coefficient for drift
%       b              - scalar coefficient for diffusion
%       z              - 3x1 state vector [z1; z2; z3]
%
%   Outputs:
%       L0a, L0b       - L0 applied to drift_term and diffusion_term
%       L1a, L1b       - L1 applied to drift_term and diffusion_term
%       L0L0a, L0L0b   - L0 composed twice on drift_term and diffusion_term
%       L0L1a, L0L1b   - L0 applied to L1a and L1b
%       L1L0a, L1L0b   - L1 applied to L0a and L0b
%       L1L1a, L1L1b   - L1 applied to L1a and L1b
%
%   This function uses symbolic first- and second-order differential operators
%   (Order1_diff and Order2_diff) to compute nested applications of Kolmogorov
%   operators L0 and L1, used in stochastic analysis.

% Extract state variables
z1 = z(1);
z2 = z(2);
z3 = z(3);

% First-order Kolmogorov operators
L0a = Order1_diff(drift_term, z1, z2, z3) * a * z + ...
      0.5 * (b * z).' * Order2_diff(drift_term, z1, z2, z3) * (b * z) + ...
      0.5 * Order1_diff(drift_term, z1, z2, z3) * (b^2 * z);

L0b = Order1_diff(diffusion_term, z1, z2, z3) * a * z + ...
      0.5 * (b * z).' * Order2_diff(diffusion_term, z1, z2, z3) * (b * z) + ...
      0.5 * Order1_diff(diffusion_term, z1, z2, z3) * (b^2 * z);

L1a = Order1_diff(drift_term, z1, z2, z3) * b * z;
L1b = Order1_diff(diffusion_term, z1, z2, z3) * b * z;

% Second-order nested Kolmogorov operators
L0L0a = Order1_diff(L0a, z1, z2, z3) * a * z + ...
        0.5 * (b * z).' * Order2_diff(L0a, z1, z2, z3) * (b * z) + ...
        0.5 * Order1_diff(L0a, z1, z2, z3) * (b^2 * z);

L0L0b = Order1_diff(L0b, z1, z2, z3) * a * z + ...
        0.5 * (b * z).' * Order2_diff(L0b, z1, z2, z3) * (b * z) + ...
        0.5 * Order1_diff(L0b, z1, z2, z3) * (b^2 * z);

L0L1a = Order1_diff(L1a, z1, z2, z3) * a * z + ...
        0.5 * (b * z).' * Order2_diff(L1a, z1, z2, z3) * (b * z) + ...
        0.5 * Order1_diff(L1a, z1, z2, z3) * (b^2 * z);

L1L0a = Order1_diff(L0a, z1, z2, z3) * b * z;
L1L1a = Order1_diff(L1a, z1, z2, z3) * b * z;
L1L0b = Order1_diff(L0b, z1, z2, z3) * b * z;

L0L1b = Order1_diff(L1b, z1, z2, z3) * a * z + ...
        0.5 * (b * z).' * Order2_diff(L1b, z1, z2, z3) * (b * z) + ...
        0.5 * Order1_diff(L1b, z1, z2, z3) * (b^2 * z);

L1L1b = Order1_diff(L1b, z1, z2, z3) * b * z;

end
