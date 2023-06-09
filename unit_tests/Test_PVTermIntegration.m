clear;clc;close all;addpath('..\');
% Principal value term Ii,mn = ∫∫fm*(nxfn)ds.
% Principal value term In,mn = ∫∫(nxfm)*(nxfn)ds = ∫∫fm*fnds.
p = randn(3);
t = [1, 2, 3, 1, 1];
T = Triangle(p(t(:,1),:), p(t(:,2),:), p(t(:,3),:));
ord = 7;
[qp, w] = GQ(T.getP(1), T.getP(2), T.getP(3), ord);
Ii = zeros(3); Iia = zeros(3); In = zeros(3); Inn = zeros(3); Ina = zeros(3);
for m=1:3
  for n=1:3
    Iia(m,n) = T.getSelfTermIi(m, n);
    Ina(m,n) = T.getSelfTermIn(m, n);
    for q = 1:ord
      Ii(m,n)   = Ii(m,n)   + T.getArea .* sum(w(q) .* T.getF(qp(q,:), m)   .* T.getNxF(qp(q,:), n), 2);
      In(m,n)   = In(m,n)   + T.getArea .* sum(w(q) .* T.getNxF(qp(q,:), m) .* T.getNxF(qp(q,:), n), 2);
      Inn(m,n)  = Inn(m,n)  + T.getArea .* sum(w(q) .* T.getF(qp(q,:), m)   .* T.getF(qp(q,:), n), 2);
    end
  end
end
Ii, Iia
In, Inn, Ina

