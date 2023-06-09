function [ V, d, isNotConverged, VV ] = eigsm( OP, M, n, k, v0, maxit, tol, p, method, displ, randStr )
if displ
    disp( '--- StartKrylovSchur ---' );
end
% Normalize starting vector
v = v0;
v = v / sqrt( abs( v' * ( M * v ) ) );
if ~( abs( sqrt( v' * ( M * v ) ) - 1 ) < 1e-14 )
    % Normalization of starting vector failed
    v2 = v0 / norm( v0 );
    if ~( abs( norm( v2 ) - 1 ) < 1e-14 )
        error( 'InvalidStartingVector' )
    else
        error( 'InvalidStartingVectorBnorm' )
    end
end
v0 = v;
[ V, d, isNotConverged, stopAlgorithm, VV ] = KSnonherm( OP, M, n, k, v0, maxit, tol, p, method, displ, randStr );
if stopAlgorithm
    % Failure to build an orthogonal subspace
    error( 'NoOrthogonalSubspace' );
end
end
%-------------------------------------------------------------------------------
function [ U, d, isNotConverged, stopAlgorithm, V ] = KSnonherm( OP, M, n, k, v0, maxit, tol, p, method, displ, randStr )
% Applies Krylov Schur algorithm for non-Hermitian matrices (reference 1)

% Initialize variables
V = zeros( n, p );
v = v0;
stopAlgorithm = false;
k0 = k; % k0 is original k. (variable k will be adaptively increased)
H = [];
nconv = 0;
sizeV = 1;
% Begin Main Algorithm
for mm = 1 : maxit
    % Build Krylov subspace in V, H:
    for jj = sizeV : p
        V(:,jj) = v;
        r = OP * ( M * v );
        Vjj = V(:,1:jj);
        w = Vjj' * ( M * r );
        r = r - Vjj * w;
        % Reorthogonalize
        [ r, normRes, stopAlgorithm, w ] = robustReorthogonalize( V, r, M, jj, randStr, w );
        if stopAlgorithm
            U = [];
            d = [];
            isNotConverged = false( 0, 1 );
            return;
        end
        % Save data
        v = r;
        H = [ H, w; zeros(1,jj-1), normRes ];
    end
    % Should we expect conjugate pairs?
    isrealprob = isreal( H );
    % Returns 2x2 block form if H is real
    [ X, T ] = schur( H(1:end-1,:) );
    % Compute eigenvalues
    [ U, d ] = eig( T, 'vector' );
    U = X * U;
    % Implicitly calculate residuals
    res = abs( H(end,:) * U );
    % Sort eigenvalues and residuals
    ind = whichEigenvalues( d, method );
    d = d(ind);
    res = res(ind);
    % Number of converged eigenpairs:
    nconvold = nconv;
    isNotConverged = ~( res(1:k0)' < tol * max( eps^(2/3), abs( d(1:k0) ) ) );
    nconv = nnz( ~isNotConverged );
    if displ
        minrelres = min( res(isNotConverged)' ./ max( eps^(2/3), abs( d(isNotConverged) ) ) );
        if ~isempty( minrelres )
            disp( [ 'Iter ', num2str(mm), ' of ', num2str(maxit),...
                ', Ritz values converged ', num2str(nconv), ' of ', num2str(k0),...
                ', resid ', num2str(minrelres) ] );
        else
            disp( [ 'Converged at iter ', num2str(mm), ' of ', num2str(maxit),...
                ', Ritz values converged ', num2str(nconv), ' of ', num2str(k0),...
                ' against tol ', num2str(tol) ] );
        end
    end
    if nconv >= k0 || mm == maxit
        % Stop the algorithm now
        break;
    else
        % Adjust k to prevent stagnating (see reference 2)
        k = k0 + min(nconv, floor( ( p - k0 ) / 2 ) );
        if k == 1 && p > 3
            k = floor( p / 2 );
        end
        % Lola's heuristic
        if  k + 1 < p && nconvold > nconv
            k = k + 1;
        end
    end
    % Get original ordering of eigenvalues back
    d = ordeig( T );
    % Choose desired eigenvalues in d to create a Boolean select vector
    ind = whichEigenvalues( d, method );
    ind = ind(1:k);
    select = false( 1, p );
    select(ind) = true;
    % Make sure both parts of a conjugate pair are present
    if isrealprob
        for i = ind'
            if ( i < p && T(i+1,i) ~= 0 && ~select(i+1) )
                select(i+1) = true;
                k = k + 1;
            end
            if ( i > 1 && T(i, i-1) ~= 0 && ~select(i-1) )
                select(i-1) = true;
                k = k + 1;
            end
        end
    end
    % Reorder X and T based on select
    [ X, T ] = ordschur( X, T, select );
    % Store variables for next iteration
    Xk = X(:,1:k);
    H = [ T(1:k, 1:k); H(end, :) * Xk ];
    V(:,1:k) = V * Xk;
    sizeV = k + 1;
end
U = U(:,ind(1:k0));
d = d(1:k0);
U = V * U;
end
%-------------------------------------------------------------------------------
function [ r, normRes, stopAlgorithm, w]  = robustReorthogonalize( V, r, M, index, randStr, wIn )

normr0 = sqrt( abs( r' * ( M * r ) ) );
if nargin < 6
    wIn = zeros( index, 1 );
end
w = wIn;
stopAlgorithm = false;
% Reorthogonalize:
Vjj = V(:,1:index);
dw = Vjj' * ( M * r );
r = r - Vjj * dw;
w = w + dw;
normRes = sqrt( abs( r' * ( M * r ) ) );
numReorths = 1;
while ( normRes <= ( 1 / sqrt( 2 ) ) * normr0 && numReorths < 5 )
    dw = Vjj' * ( M * r );
    r = r - Vjj * dw;
    w = w + dw;
    normr0 = normRes;
    normRes = sqrt( abs( r' * ( M * r ) ) );
    numReorths = numReorths + 1;
end

if ( normRes <= ( 1 / sqrt( 2 ) ) * normr0 )
    % Cannot Reorthogonalize, invariant subspace found.
    normRes = 0;
    w = wIn;
    % Try a random restart
    stopAlgorithm = true;
    for restart = 1 : 3
        % Do a random restart: Will try at most three times
        r = randn( randStr, size( r, 1 ), 1 );
        % Orthogonalize r
        Vindex = V(:,1:index);
        r = r - Vindex * ( Vindex' * ( M * r ) );
        rMr = sqrt( abs( r'* ( M * r ) ) );
        r = r / rMr;
        % Re-orthogonalize if necessary
        stopAlgorithm = true;
        for reorth = 1 : 5
            % Check orthogonality
            Mr = M * r;
            VMr = Vindex' * Mr;
            rMr = sqrt( abs( r' * Mr ) );
            if ( abs( rMr - 1 ) <= 1e-10 && all( abs( VMr ) <= 1e-10 ) )
                stopAlgorithm = false;
                break;
            end
            % Re-orthogonalize
            r = r - Vindex * VMr;
            r = r / sqrt( abs( r' * ( M * r ) ) );
        end
        if ~stopAlgorithm
            break;
        end
    end
else
    r = r / normRes;
end
end
%-------------------------------------------------------------------------------
function [ ind ] = whichEigenvalues( d, method )

switch method
    case 'largestabs'
        [ ~, ind ] = sort( abs( d ), 'descend' );
    case 'largestreal'
        [ ~, ind ] = sort( real( d ), 'descend' );
    case 'smallestreal'
        [ ~, ind ] = sort( real( d ), 'ascend' );
    case 'largestimag'
        [ ~, ind ] = sort( imag( d ), 'descend' );
    case 'smallestimag'
        [ ~, ind ] = sort( imag( d ), 'ascend' );
    case 'bothendsreal'
        [ ~, ind ] = sort( real( d ), 'descend' );
        ind2 = [ ind, flip( ind ) ]';
        ind2 = ind2(:);
        ind = ind2(1:size(d,1));
    case 'bothendsimag'
        [ ~, ind ] = sort( imag( d ), 'descend' );
        ind2 = [ ind, flip( ind ) ]';
        ind2 = ind2(:);
        ind = ind2(1:size(d,1));
    case 'smallestabs'
        [ ~,ind ] = sort( abs( d ), 'ascend' );
    case 'smallestimagabs'
        [ ~, ind ] = sort( abs( imag( d ) ), 'ascend' );
end
end

