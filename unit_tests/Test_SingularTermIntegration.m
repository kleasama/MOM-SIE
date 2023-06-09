clear;clc;close all;addpath('..\');
% Node coordinate list.
% p = 0.1*rand(3);
p = [-0.30, 0.00, 0.00;...
      0.70, 0.00, 0.00;...
      0.00, 0.87, 0.00];
% Triangle-to-node connectivity.
t = [1, 2, 3, 1, 1];
% Triangle object.
T = Triangle(p(t(:,1),:), p(t(:,2),:), p(t(:,3),:));
% Field point.
rf = [0.3, 0.2, 0.2];
% Number of integration points.
Ns = 33;
% Singular term selection n = -1,1,3,...
n = -1;
% Source-quadrature points & weights.
[rs, ws] = GQ(T.getP(1), T.getP(2), T.getP(3), Ns);
% Source-to-field distance vector.
D = rf - rs;
% Source-to-field distance.
R = sqrt(sum(D .* D, 2));
% Gradient ∇Rⁿ
gradRn = n .* R.^(n - 2) .* D;
% Analytic evaluation of singular terms (use recursion to obtain the n-th term).
Im   = T.getSingularTermIm(rf, n);
K1   = T.getSingularTermK1(rf, n);
K2v1 = T.getSingularTermK2(rf, n, 1);
K2v2 = T.getSingularTermK2(rf, n, 2);
K2v3 = T.getSingularTermK2(rf, n, 3);
K3   = T.getSingularTermK3(rf, n);
K4v1 = T.getSingularTermK4(rf, n, 1);
K4v2 = T.getSingularTermK4(rf, n, 2);
K4v3 = T.getSingularTermK4(rf, n, 3);
% Numerical approximation of singular terms.
K1_approx   = T.getArea .* sum(ws .* R.^n, 1);
K2v1_approx = T.getArea .* sum(ws .* R.^n .* (rs - T.getP(1)), 1);
K2v2_approx = T.getArea .* sum(ws .* R.^n .* (rs - T.getP(2)), 1);
K2v3_approx = T.getArea .* sum(ws .* R.^n .* (rs - T.getP(3)), 1);
K3_approx   = T.getArea .* sum(ws .* gradRn, 1);
K4v1_approx = T.getArea .* sum(ws .* cross(-gradRn, rs - T.getP(1), 2), 1);
K4v2_approx = T.getArea .* sum(ws .* cross(-gradRn, rs - T.getP(2), 2), 1);
K4v3_approx = T.getArea .* sum(ws .* cross(-gradRn, rs - T.getP(3), 2), 1);
% Results.
format longeng;
disp('K1')
disp([K1_approx; K1])
disp('K2v1')
disp([K2v1_approx; K2v1])
disp('K2v2')
disp([K2v2_approx; K2v2])
disp('K2v3')
disp([K2v3_approx; K2v3])
disp('K3')
disp([K3_approx; K3])
disp('K4v1')
disp([K4v1_approx; K4v1])
disp('K4v2')
disp([K4v2_approx; K4v2])
disp('K4v3')
disp([K4v3_approx; K4v3])
format;
% Visualisation.
viewer(p, t); hold on;
plot3(rf(1), rf(2), rf(3), 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
plot3(rs(:,1), rs(:,2), rs(:,3),'.r');
hold off
