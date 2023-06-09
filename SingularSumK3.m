function [ intsum ] = SingularSumK3(KER, T, j, nf)
%%SingularSumK3 Summation of the nl analytic integrals K3 derived from the
% singular terms extracted from the free space Helmholtz Green's function
% gradient kernel.
%--------------------------------------------------------------------------
% Input:
%   KER    - Kernel object.
%   T      - Triangle object.
%   j      - The index of the local node.
%   nf     - The nf-th integration point index.
%--------------------------------------------------------------------------
% Output:
%   intsum - The summation of the analytic integrals K3 of the extracted
%            terms from the free space Helmholtz Green's function gradient
%            kernel.
%--------------------------------------------------------------------------
% Copyright: Klearchos A. Samaras, kleasama@gmail.com, June 2023
%--------------------------------------------------------------------------
  intsum = 0;
  for l = 1 : KER.nl
    n = 2 * (l - 1);
    intsum = intsum + KER.dGFsingCoeff(n) ...
           .* T.getSingularTermK3(KER.rf(nf,:), n - 1);
  end
  intsum = intsum .* 0.5 .* T.getEdgeLength(j) ./ T.getArea;
end
