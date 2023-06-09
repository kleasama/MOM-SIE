function [ dGext ] = dGFext(k, R, nl)
%%dGFext Extracted portion of the scalar part of the free space Helmholtz
% Green's function gradient intended for singularity extraction techniques.
%--------------------------------------------------------------------------
% Input:
%   k     - The wavenumber of operation.
%   R     - The distance between the field and source point.
%   nl    - The number of extrated terms using MacLaurin series.
%--------------------------------------------------------------------------
% Output:
%   dGext - The overall extracted part of the Green's function kernel.
%--------------------------------------------------------------------------
% Copyright: Klearchos A. Samaras, kleasama@gmail.com, June 2023
%--------------------------------------------------------------------------
  dGext = 0;
  for l = 1 : nl
    n = 2 * (l - 1);
    dGext = dGext + dGFsingCoeff(k, n) .* R.^(n - 3);
  end
end
