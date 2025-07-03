function [Phi] = so3_wedge(phi)
%SO3_WEDGE Wedge operator: Converts a 3x1 vector to a 3x3 skew-symmetric matrix
%
%   Phi = so3_wedge(phi)
%
%   Inputs:
%       phi - 3x1 vector
%
%   Outputs:
%       Phi - 3x3 skew-symmetric matrix (element of so(3))
%
%   This function maps a vector in R^3 to its corresponding element
%   in the Lie algebra so(3) (the set of 3x3 skew-symmetric matrices).

% Input validation
if ~isvector(phi) || length(phi) ~= 3
    error('Input must be a 3x1 vector.');
end

Phi = [    0,   -phi(3),  phi(2);
        phi(3),     0,   -phi(1);
       -phi(2),  phi(1),     0 ];

end
