function [Rot] = so3_exp_new(phi)
%SO3_EXP_NEW Exponential map from so(3) to SO(3)
%
%   Rot = so3_exp_new(phi)
%
%   Inputs:
%       phi - 3x1 vector representing rotation in axis-angle form
%
%   Output:
%       Rot - 3x3 rotation matrix in SO(3)
%
%   This function implements the exponential map from so(3) to SO(3)
%   using Rodrigues' rotation formula. For small ||phi||, a first-order
%   approximation is used.

TOL = 1e-9;

angle = norm(phi);
if size(phi,1) > 1 || size(phi,2) > 1
    if angle < TOL
        % Near |phi|==0, use first order Taylor expansion
        Rot = eye(3) + phi;
    else
        axis = phi / angle;
%         Rot = eye(3) - (1-cos(angle))*(axis*axis') + ...
%             sin(angle) * axis;
        Rot = eye(3) + (1-cos(angle))*(axis^2) + ...
            sin(angle) * axis;
    end
else
    if angle < TOL
        % Near |phi|==0, use first order Taylor expansion
        Rot = eye(3) + so3_wedge(phi);
    else
        axis = phi / angle;
        Rot = cos(angle) * eye(3) + (1-cos(angle))*(axis*axis') + ...
            sin(angle) * so3_wedge(axis);
    end
end
