function [ Aschur, Bschur ] = SC( A, B, m, opt )
% Schur Complement domain decomposition for the matrix equation A * X = B, 
% where
%
%          /          \            /  \            /  \
%         | A11    A12 |          | X1 |          | B1 |
%     A = |            | ,    X = |    | ,    B = |    |
%         | A21    A22 |          | X2 |          | B2 |
%          \          /            \  /            \  /
%
% Targeting the top-left block matrix or the unknowns X1 we obtain:
%   Aschur * X1  = Bschur, 
% where
%   Aschur = A11 - A12 * inv(A22) * A21
%   Bschur = B1  - A12 * inv(A22) * B2
%
% Targeting the bottom-right block matrix or the unknowns X2 we obtain:
%   Aschur * X2  = Bschur, 
% where
%   Aschur = A22 - A21 * inv(A11) * A12
%   Bschur = B2  - A21 * inv(A11) * B1
%--------------------------------------------------------------------------
% Usage:
%    [ Aschur, Bschur ] = SC( A, B, m, opt );
%--------------------------------------------------------------------------
% Input:
%    A      - The square matrix of the matrix equation.
%    B      - The RHS vector of the matrix equation.
%    m      - The truncation size for the matrix block.
%    opt    - The option for the Schur complement target.
%             opt = 0 --> top-left block matrix (X1).
%             opt = 1 --> bottom-right block matrix (X2).
%--------------------------------------------------------------------------
% Output:
%    Aschur - The Schur complement matrix.
%    Bschur - The equivalent RHS vector that corresponds to Aschur.
%--------------------------------------------------------------------------
% Copyright: Klearchos A. Samaras, kleasama@gmail.com, June 2023
%--------------------------------------------------------------------------
  
  % Schur Complement on the upper left block of Z-matrix.
  if opt == 0
    P = A(1:m,m+1:end) / A(m+1:end,m+1:end);
    Aschur = A(1:m,1:m) - P * A(m+1:end,1:m);
    Bschur = B(1:m)  - P * B(m+1:end);
    % Schur Complement on the lower right block of Z-matrix.
  elseif opt == 1
    P = A(m+1:end,1:m) / A(1:m,1:m);
    Aschur = A(m+1:end,m+1:end) - P * A(1:m,m+1:end);
    Bschur = B(m+1:end)  - P * B(1:m);
  end
end