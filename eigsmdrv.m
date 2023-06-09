function [ v ] = eigsmdrv( A, M, k )
  sigma = 0;
  n = size( A, 1 );
  maxit = 300;
  tol = 1e-10;
  p = min( max( 2 * k, 20 ), n );
  method = 'largestabs';
  displ = 0;
  randStr = randn( n, 1 );
  shiftAndInvert = 0;
  OP = inv( A - sigma * M ) * M;
  v0 = randn(n,1) + 1i * randn(n,1);
  [ U, v, isNotConverged, X ] = eigsm( OP, eye( n ), n, k, v0, maxit, tol, p, method, displ, randStr );
  v = 1./v + sigma;
  [ ~, i ] = sort( abs( v ) );
end
