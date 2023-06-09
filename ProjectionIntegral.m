function [ projint ] = ProjectionIntegral(wf, nf, jacobian, g, f)
% SubdomainIntegral Numerical computation of a projection integral or inner
% product of a function f with the function g, <g,f>.
%--------------------------------------------------------------------------
% Input:
%   wf          - Quadrature weights.
%   nf          - Number of integration points.
%   jacobian    - Jacobian determinant (extents of the subdomain).
%   g           - Testing or projector function.
%   f           - Source function that needs to be tested.
%--------------------------------------------------------------------------
% Output:
%   projint     - The projection integral <g,f>.
%--------------------------------------------------------------------------
% Copyright: Klearchos A. Samaras, kleasama@gmail.com, June 2023
%--------------------------------------------------------------------------
  projint = 0;
  for t = 1 : nf
    projint = projint + wf(t) .* g(t,:) * transpose(f(t,:));
  end
  projint = jacobian .* projint;
end
