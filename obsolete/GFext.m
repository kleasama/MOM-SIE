function [ Gext ] = GFext(k, R, nl)
%%GFext Extracted portion of the free space Helmholtz Green's function
% intended for singularity extraction techniques.
%--------------------------------------------------------------------------
% Input:
%   k    - The wavenumber of operation.
%   R    - The distance between the field and source point.
%   nl   - The number of extrated terms using MacLaurin series.
%--------------------------------------------------------------------------
% Output:
%   Gext - The overall extracted part of the Green's function kernel.
%--------------------------------------------------------------------------
% Copyright: Klearchos A. Samaras, kleasama@gmail.com, June 2023
%--------------------------------------------------------------------------
  Gext = 0;
  for l = 1 : nl
    n = 2 * (l - 1);
    Gext = Gext + GFsingCoeff(k, n) .* R.^(n - 1);
  end
end

