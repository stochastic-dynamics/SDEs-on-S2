function [v] = logmap_sphere(p, x)
%LOGMAP_SPHERE Logarithmic map on the unit 2-sphere (S^2)
%
%   v = logmap_sphere(p, x)
%
%   Inputs:
%       p - 3x1 vector on the unit sphere (base point), ||p|| = 1
%       x - 3x1 vector on the unit sphere, ||x|| = 1
%
%   Output:
%       v - 3x1 tangent vector at point p pointing toward x
%
%   This function computes the logarithmic map from x back to the
%   tangent space at point p on the unit sphere S^2.
%   The result is the minimal tangent vector v such that:
%       expmap_sphere(p, v) ≈ x

% Check if x is numerically the same as p
if norm(p - x) < 1e-10
    v = zeros(size(p));
    return;
end

% Compute the angle between p and x
theta = acos(dot(p, x));

% Compute the unit tangent direction from p to x
direction = (x - dot(p, x) * p) / sqrt(1 - (dot(p, x))^2);

% Scale by the angle to get the logarithmic map
v = theta * direction;

end
