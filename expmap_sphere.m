function [x] = expmap_sphere(p, v)
%EXPMAP_SPHERE Exponential map on the 2-sphere (S^2)
%
%   x = expmap_sphere(p, v)
%
%   Inputs:
%       p - 3x1 vector on the unit sphere (base point), ||p|| = 1
%       v - 3x1 tangent vector at point p, v' * p = 0
%
%   Output:
%       x - 3x1 vector on the unit sphere resulting from exponential map
%
%   This function computes the exponential map of a tangent vector v 
%   at a point p on the unit sphere S^2. It maps v from the tangent 
%   space at p onto the manifold.

% Check if the tangent vector is negligible
if norm(v) < 1e-10
    x = p;
    return;
end

% Normalize the tangent vector and apply exponential map
theta = norm(v);
x = cos(theta) * p + sin(theta) * v / theta;

end
