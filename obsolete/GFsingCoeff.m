function [ coeff ] = GFsingCoeff(k, n)
%%GFsingCoeff Coefficient of the extracted term of the free space Helmholtz
% Green's function intended for singularity extraction techniques.
%--------------------------------------------------------------------------
% Input:
%   k     - The wavenumber of operation.
%   n     - The exponent of the extrated term using MacLaurin series.
%           We increment n by one to be in sync with the singular term
%           power.
%--------------------------------------------------------------------------
% Output:
%   coeff - The coefficient of the singular term raised to the n-th power.
%--------------------------------------------------------------------------
% Copyright: Klearchos A. Samaras, kleasama@gmail.com, June 2023
%--------------------------------------------------------------------------
  oneover4pi = 0.25 ./ pi;
  mjk = -1i .* k;
  coeff = oneover4pi .* mjk.^(n) ./ factorial(n);
end

