function [ subint ] = SubdomainIntegral(ws, ns, jacobian, f, kern, doextract, kern_extr, subint_sing)
%SubdomainIntegral Numerical computation of a subdomain integral using
% singularity extraction techniques.
%--------------------------------------------------------------------------
% Input:
%   ws          - Quadrature weights.
%   ns          - Number of integration points.
%   jacobian    - Jacobian determinant (extents of the subdomain).
%   f           - Integrand function as evaluated at the integration points.
%   kern        - Kernel function as evaluated at the field point and the
%                 collection of the quadrature points.
%   doextract   - Control flag to perform the singularity extraction.
%   kernel_extr - Extracted part of the kernel function as evaluated at the
%                 field point and the collection of the quadrature points.
%   subint_sing - The analytically evaluated integral of the extracted part
%                 of the kernel as evaluated at the feild point.
%--------------------------------------------------------------------------
% Output:
%   subint      - The integral of the integrand function f on the subdomain
%                 using the kernel function kern.
%--------------------------------------------------------------------------
% Copyright: Klearchos A. Samaras, kleasama@gmail.com, June 2023
%--------------------------------------------------------------------------
  subint = 0;
  if doextract ~= 1
    for s = 1 : ns
      subint = subint + ws(s) .* kern(s,:) .* f(s,:);
    end
    subint = jacobian .* subint;
  else
    for s = 1 : ns
      subint = subint + ws(s) .* (kern(s,:) - kern_extr(s,:)) .* f(s,:);
    end
    subint = jacobian .* subint;
    subint = subint + subint_sing;
  end
end
