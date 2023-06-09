classdef Triangle
%TRIANGLE Planar triangle object.
% Complete set of geometric calculations for planar triangles.
%
%
%           p3
%            O
%           / \
%          /   \
%     e2  / A   \ e1
%        /   .   \
%       /     c   \
%   p1 O - - > - - O p2
%            e3
%               [T]
%
%--------------------------------------------------------------------------
% REFERENCES:
% [1] - S. Järvenpää, M. Taskinen, P. Ylä-Oijala, "Singularity Subtraction
%       Technique for High-Order Polynomial Vector Basis Functions on
%       Planar Triangles", IEEE Transactions on Antennas and Propagation,
%       Vol. 54, no. 1, January 2006.
%
% [2] - P. Ylä-Oijala, M. Taskinen, "Calculation of CFIE Impedance Matrix
%       Elements With RWG and nxRWG Functions", IEEE Transactions on
%       Antennas and Propagation, Vol. 51, No. 8, August 2003.
%
% [3] - W. C. Gibson, "The Method of Moments in Electromagnetics", Second
%       Edition, CRC Press, 2015.
%
% [4] - P. Arcioni, M. Bressan, L. Perregrini, "On the Evaluation of the
%       Double Surface Integrals Arising in the Application of the
%       Boundary Integral Method to 3-D Problems", IEEE Transactions on
%       Microwave Theory and Techniques, Vol. 45, No. 3, March 1997.
%
%--------------------------------------------------------------------------
% Copyright: Klearchos A. Samaras, kleasama@gmail.com, June 2023
%--------------------------------------------------------------------------

  properties
    % The coordinates of all points in columns 1 through 9.
    p
    % The number of triangles processed by the triangle object.
    nt
    % Container for the edge vectors.
    e
    % Container for the edge lengths.
    L
    % Container for the edge unit vectors.
    s
    % Container for the trianglular facet area.
    A
    % Container for the face normal unit vectors.
    n
    % Container for the edge normal unit vectors.
    m
    % Container for the triangle centroid.
    rc
  end

  methods
    function [ obj ] = Triangle(point1, point2, point3)
    %%Triangle Constructor of Triangle object.
      obj.p  = [point1, point2, point3];
      obj.nt = obj.calcTrianglesTotal;
      obj.e  = [obj.calcEdgeVector(1),...
                obj.calcEdgeVector(2),...
                obj.calcEdgeVector(3)];
      obj.L  = [obj.calcEdgeLength(1),...
                obj.calcEdgeLength(2),...
                obj.calcEdgeLength(3)];
      obj.s  = [obj.calcEdgeUnit(1),...
                obj.calcEdgeUnit(2),...
                obj.calcEdgeUnit(3)];
      obj.A  = obj.calcArea;
      obj.n  = obj.calcFaceNormal;
      obj.m  = [obj.calcEdgeNormal(1),...
                obj.calcEdgeNormal(2),...
                obj.calcEdgeNormal(3)];
      obj.rc = obj.calcCentroid;
    end
    
    function [ obj ] = processTriangle(obj)
    %%processTriangle Processes the triangle object.
      obj.nt = obj.calcTrianglesTotal;
      obj.e  = [obj.calcEdgeVector(1),...
                obj.calcEdgeVector(2),...
                obj.calcEdgeVector(3)];
      obj.L  = [obj.calcEdgeLength(1),...
                obj.calcEdgeLength(2),...
                obj.calcEdgeLength(3)];
      obj.s  = [obj.calcEdgeUnit(1),...
                obj.calcEdgeUnit(2),...
                obj.calcEdgeUnit(3)];
      obj.A  = obj.calcArea;
      obj.n  = obj.calcFaceNormal;
      obj.m  = [obj.calcEdgeNormal(1),...
                obj.calcEdgeNormal(2),...
                obj.calcEdgeNormal(3)];
      obj.rc = obj.calcCentroid;
    end
   
    function [ nt ] = calcTrianglesTotal(obj)
    %calcTrianglesTotal Calculates the number of triangles processed by the
    % triangle object.
      nt = length(obj.p(:,1));
    end
    
    function [ nt ] = getTrianglesTotal(obj)
    %%getTrianglesTotal Returns number of triangles processed by the
    % triangle object.
      nt = obj.nt;
    end

    function [ idx ] = getIdxN(obj, j)
    %%getIdxN Returns the extended local index to facilitate operations
    % that make use of the (j+1) and (j+2) permutations.
    % i.e. if j = 3, then (j+1) = 4 --> 1 and (j+2) = 5 --> 2.
      node_idx = [1, 2, 3, 1, 2];
      idx = node_idx(j);
    end

    function [ pj ] = getP(obj, j)
    %%getP Returns the coordinates of the j-th local node.
      pj = obj.p(:,(obj.getIdxN(j)-1)*3+1:(obj.getIdxN(j))*3);
    end

    function [ e ] = calcEdgeVector(obj, j)
    %%calcEdgeVector Calculates the j-th edge vector.
      e = obj.getP(j + 2) - obj.getP(j + 1);
    end
    
    function [ e ] = getEdgeVector(obj, j)
    %%getEdgeVector Returns the j-th edge vector.
      e = obj.e(:,3*(j-1)+1:3*j);
    end

    function [ l ] = calcEdgeLength(obj, j)
    %%calcEdgeLength Calculates the j-th edge length.
      l = sqrt(sum(obj.getEdgeVector(j).^2, 2));
    end
    
    function [ l ] = getEdgeLength(obj, j)
    %%getEdgeLength Returns the j-th edge length.
      l = obj.L(:,j);
    end

    function [ s ] = calcEdgeUnit(obj, j)
    %%calcEdgeUnit Calculates the j-th edge tangential unit vector.
      s = obj.getEdgeVector(j) ./ obj.getEdgeLength(j);
    end
    
    function [ s ] = getEdgeUnit(obj, j)
    %%getEdgeUnit Returns the j-th edge tangential unit vector.
      s = obj.s(:,3*(j-1)+1:3*j);
    end

    function [ mp ] = getEdgeMidPoint(obj, j)
    %%getEdgeMidPoint Returns the j-th edge mid point coordinates.
      mp = 0.5 * (obj.getP(j + 1) + obj.getP(j + 2));
    end

    function [ a ] = getAngle(obj, j)
    %%getAngle Returns the j-th angle of the planar triangle facet.
      a = obj.getEdgeLength(obj.getIdxN(j));
      b = obj.getEdgeLength(obj.getIdxN(j + 1));
      c = obj.getEdgeLength(obj.getIdxN(j + 2));
      a = acos((b .* b + c .* c - a .* a) ./ (2 * b .* c));
    end

    function [ pm ] = getPerimeter(obj)
    %%getPerimeter Returns the perimeter of the planar triangle facet.
      pm = obj.getEdgeLength(1) ...
         + obj.getEdgeLength(2) ...
         + obj.getEdgeLength(3);
    end

    function [ c ] = calcCentroid(obj)
    %%calcCentroid Returns the centroid coordinates of the planar
    % triangle facet.
      onethird = 1 / 3;
      c = onethird .* (obj.getP(1) + obj.getP(2) + obj.getP(3));
    end
    
    function [ c ] = getCentroid(obj)
    %%getCentroid Returns the centroid coordinates of the planar
    % triangle facet.
      c = obj.rc;
    end

    function [ a ] = calcArea(obj)
    %%calcArea Calculates the area of the planar triangle facet.
      a = 0.5 * sqrt(...
        sum(cross(obj.getEdgeVector(1), obj.getEdgeVector(2), 2).^2, 2));
    end
    
    function [ a ] = getArea(obj)
    % Returns the area of the planar triangle facet.
      a = obj.A;
    end

    function [ n ] = calcFaceNormal(obj)
    %%calcFaceNormal Calculates the unit normal vector of the planar
    % triangle facet.
      n = 0.5 * cross(obj.getEdgeVector(1), obj.getEdgeVector(2), 2) ...
        ./ obj.getArea;
    end

    function [ n ] = getFaceNormal(obj)
    %%getFaceNormal Returns the unit normal vector of the planar
    % triangle facet.
      n = obj.n;
    end

    function [ m ] = calcEdgeNormal(obj, j)
    %%calcEdgeNormal Calculates the j-th edge normal unit vector with
    % outward direction and on the plane defined by the planar triangle.
      m = cross(obj.getEdgeUnit(j), obj.getFaceNormal, 2);
    end
    
    function [ m ] = getEdgeNormal(obj, j)
    %%getEdgeNormal Returns the j-th edge normal unit vector with
    % outward direction and on the plane defined by the planar triangle.
      m = obj.m(:,3*(j-1)+1:3*j);
    end

    function [ h ] = getHeight(obj, j)
    %%getHeight Returns the j-th height.
      h = sum((obj.getP(j + 1) - obj.getP(j)) .* obj.getEdgeNormal(j), 2);
    end

    function [ md ] = getMedian(obj, j)
    %%getMedian Returns the j-th median.
      md = 0.5 * sqrt(...
           2 * obj.getEdgeLength(obj.getIdxN(j + 1)).^2 ...
         + 2 * obj.getEdgeLength(obj.getIdxN(j + 2)).^2 ...
             - obj.getEdgeLength(obj.getIdxN(j)).^2);
    end

    function [ u ] = getU(obj)
      u = obj.getEdgeUnit(3);
    end

    function [ v ] = getV(obj)
      v = cross(obj.getFaceNormal, obj.getU, 2);
    end

    function [ u3 ] = getU3(obj)
      u3 = -sum(obj.getEdgeVector(2) .* obj.getEdgeUnit(3), 2);
    end

    function [ v3 ] = getV3(obj)
      v3 = 2 * obj.getArea ./ obj.getEdgeLength(3);
    end

    function [ ic ] = getIncenter(obj)
    %%getIncenter Returns the incenter coordinates.
      u  = obj.getEdgeUnit(3);
      v  = obj.getV;
      l2 = obj.getEdgeLength(2);
      h3 = obj.getHeight(3);
      l3_over_per = obj.getEdgeLength(3) ./ obj.getPerimeter;
      x3 = sqrt(abs(l2.^2 - h3.^2));
      log_vec = (cos(obj.getAngle(1)) < 0);
      x3(log_vec) = -x3(log_vec);
      x_loc = (l2 + x3) .* l3_over_per;
      y_loc = h3 .* l3_over_per;
      ic = obj.getP(1) + u .* x_loc + v .* y_loc;
    end

    function [ ir ] = getInradius(obj)
    %%getInradius Returns the inradius.
      ir = 2 * obj.getArea ./ obj.getPerimeter;
    end

    function [ cc ] = getCircumcenter(obj)
    %%getCircumcenter Returns the circumcenter coordinates.
      u = obj.getEdgeUnit(3);
      v = obj.getV;
      l2 = obj.getEdgeLength(2);
      l3 = obj.getEdgeLength(3);
      h3 = obj.getHeight(3);
      sin2a = sin(2 * obj.getAngle(1));
      sin2b = sin(2 * obj.getAngle(2));
      sin2c = sin(2 * obj.getAngle(3));
      x3 = sqrt(abs(l2.^2 - h3.^2));
      log_vec = (cos(obj.getAngle(1)) < 0);
      x3(log_vec) = -x3(log_vec);
      w = sin2a + sin2b + sin2c;
      x_loc = (l3 .* sin2b + x3 .* sin2c) ./ w;
      y_loc = h3 .* sin2c ./ w;
      cc = obj.getP(1) + u .* x_loc + v .* y_loc;
    end

    function [ cr ] = getCircumradius(obj)
    %%getCircumradius Returns the circumradius.
      cr = 0.5 * obj.getEdgeLength(1) ./ sin(obj.getAngle(1));
    end
    
    function [ w0 ] = getW0(obj, r)
      w0 = sum((r - obj.getP(1)) .* obj.getFaceNormal, 2);
    end

    function [ t0j ] = getT0(obj, r, j)
      t0j = sum((obj.getP(obj.getIdxN(j + 1)) - r) ...
          .* obj.getEdgeNormal(j), 2);
    end

    function [ smj ] = getSm(obj, r, j)
      smj = sum((obj.getP(obj.getIdxN(j + 1)) - r) ...
          .* obj.getEdgeUnit(j), 2);
    end

    function [ spj ] = getSp(obj, r, j)
      spj = obj.getSm(r, j) + obj.getEdgeLength(j);
    end

    function [ rpj ] = getRp(obj, r, j)
      rpj = sqrt(sum((r - obj.getP(obj.getIdxN(j + 2))).^2, 2));
    end

    function [ rmj ] = getRm(obj, r, j)
      rmj = obj.getRp(r, obj.getIdxN(j + 2));
    end

    function [ r0j ] = getR0(obj, r, j)
      r0j = sqrt(obj.getT0(r, j).^2 + obj.getW0(r).^2);
    end

    function [ nmj ] = getNm(obj, r, s, j)
        nmj = (obj.getSp(r, j) - s) ./ obj.getEdgeLength(j);
    end

    function [ nmj ] = getNmNS(obj, r, j)
        nmj = obj.getN(r, obj.getIdxN(j + 1));
    end

    function [ npj ] = getNp(obj, r, s, j)
        npj = (s - obj.getSm(r, j)) ./ obj.getEdgeLength(j);
    end

    function [ npj ] = getNpNS(obj, r, j)
        npj = obj.getN(r, obj.getIdxN(j + 2));
    end

    function [ rsj ] = getRs(obj, r, s, j)
    %%getRs Returns the representation of R = sqrt(s^2 + R0j(r)^2) on the
    % j-th edge. Value s runs from smj to spj.
      rsj = sqrt(s.^2 + obj.getR0(r, j).^2);
    end

    function [ nj ] = getN(obj, r, j)
    %%getN Returns the value of the j-th node shape interpolator function
    % as evaluated at the point r.
        mj = obj.getEdgeNormal(j);
        nj = 1 - sum(mj .* (r - obj.getP(j)), 2) ./ obj.getHeight(j);
    end

    function [ gnj ] = getGradN(obj, r, j, s)
    %%getGradN Returns the surface gradient of the j-th node shape
    % interpolator function raised to s-th power the as evaluated
    % at the point r.
        gnj = -s ./ obj.getHeight(j) .* obj.getEdgeNormal(j) ...
            .* obj.getN(r, j).^(s - 1);
    end

    function [ rpj ] = getProjection(obj, r)
    %getProjection Returns the coordinates of the projection point on the
    % plane defined by the planar triangle.
      rpj = r - obj.getW0(r) .* obj.getFaceNormal;
    end

    function [ rhoj ] = getRho(obj, r, j)
    %%getRho Returns the local surface radial vector that starts from the
    % j-th vertex and ends to the projection point of r on the plane
    % defined the planar triangle.
        rhoj = (obj.getProjection(r) - obj.getP(j));
    end

    function [ rhoj ] = getRhoNS(obj, r, j)
    %%getRhoNS Returns the local surface radial vector that starts from
    % the j-th vertex and ends to the projection point of r on the plane
    % defined the planar triangle. In this routine the result is derived
    % using the nodal shape interpolator functions Nn+1 and Nn+2.
    % ρj(r) = Nj+1(r)ej+2 - Nj+2(r)ej+1.
        rhoj = obj.getN(r, obj.getIdxN(j + 1)) ...
             .* obj.getEdgeVector(obj.getIdxN(j + 2)) ...
             - obj.getN(r, obj.getIdxN(j + 2)) ...
             .* obj.getEdgeVector(obj.getIdxN(j + 1));
    end

    function [ fj ] = getF(obj, r, j)
    % Returns the positive signed half RWG-linear element that corresponds
    % to the j-th vertex of the triangle.
        fj = 0.5 * obj.getEdgeLength(j) ./ obj.getArea .* obj.getRho(r, j);
    end

    function [ fj ] = getFNS(obj, r, j)
    % Returns the positive signed half RWG-linear element that corresponds
    % to the j-th vertex of the triangle. In this routine the result is
    % derived using the nodal shape interpolator functions Nj+1 and Nj+2.
        fj = 0.5 * obj.getEdgeLength(j) ./ obj.getArea ...
           .* obj.getRhoNS(r, j);
    end

    function [ divfj ] = getDivF(obj, j)
    % Returns the positive signed surface divergence of the RWG-linear
    % element that corresponds to the j-th vertex of the triangle.
        divfj = obj.getEdgeLength(j) ./ obj.getArea;
    end

    function [ nxfj ] = getNxF(obj, r, j)
    % Returns the positive signed half nxRWG-linear element that
    % corresponds to the j-th vertex of the triangle.
        f = obj.getF(r, j);
        nxfj = [obj.n(:,2) .* f(:,3) - obj.n(:,3) .* f(:,2),...
                obj.n(:,3) .* f(:,1) - obj.n(:,1) .* f(:,3),...
                obj.n(:,1) .* f(:,2) - obj.n(:,2) .* f(:,1)];
    end

    function [ nxfj ] = getNxFNS(obj, r, j)
    % Returns the positive signed half nxRWG-linear element that
    % corresponds to the j-th vertex of the triangle.In this routine the result
    % is derived using the nodal shape interpolator functions Nj+1 and Nj+2.
        f = obj.getFNS(r, j);
        nxfj = [obj.n(:,2) .* f(:,3) - obj.n(:,3) .* f(:,2),...
                obj.n(:,3) .* f(:,1) - obj.n(:,1) .* f(:,3),...
                obj.n(:,1) .* f(:,2) - obj.n(:,2) .* f(:,1)];
    end

    function [ im1 ] = getSingularTermIm1(obj, r, j)
    % Initial value for contour singular integral recursive formulas
    % associated with the j-th local edge of the triangle.
      im1 = log((obj.getRp(r, j) + obj.getSp(r, j)) ...
          ./(obj.getRm(r, j) + obj.getSm(r, j)));
    end

    function [ bj ] = getB(obj, r, j)
    % Returns the b value used for the evaluation of the initial value of
    % the surface singular surface intergral recursive formulas associated
    % with the j-th node.
      t0j = obj.getT0(r, j);
      r0j = obj.getR0(r, j);
      w0 = obj.getW0(r);
      bj = atan2((t0j .* obj.getSp(r, j)), ...
          (r0j.^2 + abs(w0) .* obj.getRp(r, j))) ...
        - atan2((t0j .* obj.getSm(r, j)), ...
          (r0j.^2 + abs(w0) .* obj.getRm(r, j)));
    end

    function [ w0k1m3 ] = getSingularTermW0K1m3(obj, r)
    % Initial value for surface singular integral recursive formulas.
      w0 = obj.getW0(r);
      w0k1m3 = sign(w0) ...
             .* (obj.getB(r, 1) + obj.getB(r, 2) + obj.getB(r, 3));
      w0k1m3(w0==0) = 0;
    end

    function [ inj, term_count ] = getSingularTermI(obj, r, n, j)
    % Recursive formula for the contour integrals Im.
    % r - field point
    % n - power of distance R, n = -1,1,3,5...
    % j - edge index
      if (n == -1)
        % If n == -1 then we are returning the inital term Im1j.
        inj = obj.getSingularTermIm1(r, j);
        return;
      end
      % Recursion.
      for ni = 1 : 2 : n
        % Indexing for the singular term arrays.
        term_count = (ni + 1) / 2 + 1;
        inj = 1./(ni + 1) ...
            .* ( obj.getSp(r, j) .* obj.getRp(r, j).^ni ...
            - obj.getSm(r, j) .* obj.getRm(r, j).^ni ...
            + ni .* obj.getR0(r, j).^2 ...
            .* obj.getSingularTermI(r, ni - 2, j));
      end
    end

    function [ imn ] = getSingularTermIm(obj, r, n)
    % The contour integral along the boundary of the planar triangle facet.
    % Im = ∫ Rⁿ dl'
    % r - field point
    % n - power of distance R, n = -1,1,3,5...
      imn = obj.getEdgeNormal(1) .* obj.getSingularTermI(r, n, 1) ...
          + obj.getEdgeNormal(2) .* obj.getSingularTermI(r, n, 2) ...
          + obj.getEdgeNormal(3) .* obj.getSingularTermI(r, n, 3);
    end

    function [ k1n, term_count ] = getSingularTermK1(obj, r, n)
    % The surface integral K1 on the surface of the planar triangle.
    % K1 = ∫∫ Rⁿ ds'
    % r - field point
    % n - power of distance R, n = -1,1,3,5..
      w0 = obj.getW0(r);
      if n == -1
        w0k1m3 = obj.getSingularTermW0K1m3(r);
        k1n = obj.getT0(r, 1) .* obj.getSingularTermI(r, -1, 1) ...
            + obj.getT0(r, 2) .* obj.getSingularTermI(r, -1, 2) ...
            + obj.getT0(r, 3) .* obj.getSingularTermI(r, -1, 3);
        k1n(w0~=0,:) = k1n(w0~=0,:) - w0(w0~=0) .* w0k1m3(w0~=0);
        return;
      end
      for ni = 1 : 2 : n
        term_count = (ni + 1) / 2 + 1;
        k1n = 1./(ni + 2) ...
              .* ( obj.getT0(r, 1) .* obj.getSingularTermI(r, ni, 1) ...
                 + obj.getT0(r, 2) .* obj.getSingularTermI(r, ni, 2) ...
                 + obj.getT0(r, 3) .* obj.getSingularTermI(r, ni, 3) );
        k1nm2 = obj.getSingularTermK1(r, ni - 2);
        k1n(w0~=0,:) = k1n(w0~=0,:) ...
                     + 1./(ni + 2) .* (ni.*(w0(w0~=0)).^2.*k1nm2(w0~=0,:));
      end
    end

    function [ k2nj ] = getSingularTermK2(obj, r, n, j)
    % The surface integral K2j on the surface of the planar triangle
    % associated with the j-th node.
    % K2j = ∫∫ Rⁿ(r' - vj) ds'
    % r - field point
    % n - power of distance R, n = -1,1,3,5..
    % j - the local node index, j = 1,2,3
      k2nj = 1./(n + 2) .* obj.getSingularTermIm(r, n + 2) ...
           + obj.getRho(r, j) .* obj.getSingularTermK1(r, n);
    end

    function [ k3n ] = getSingularTermK3(obj, r, n)
    % The surface integral K3 on the surface of the planar triangle.
    % K3 = ∫∫ ∇Rⁿ ds'
    % r - field point
    % n - power of distance R, n = -1,1,3,5..
      if n == -1
        k3n = -obj.getSingularTermW0K1m3(r) .* obj.getFaceNormal ...
            - obj.getSingularTermIm(r, -1);
        return;
      end
      k3n = n .* obj.getW0(r) .* obj.getFaceNormal ...
          .* obj.getSingularTermK1(r, n - 2) ...
          - obj.getSingularTermIm(r, n);
    end

    function [ k4nj ] = getSingularTermK4(obj, r, n, j)
    % The surface integral K2j on the surface of the planar triangle
    % associated with the j-th node.
    % K4j = ∫∫ ∇Rⁿ×(r' - vj) ds'
    % r - field point
    % n - power of distance R, n = -1,1,3,5..
    % j - the local node index, j = 1,2,3
      k4nj = cross((r - obj.getP(j)), obj.getSingularTermK3(r, n), 2);
    end
    
    function [ i ] = getSelfTermI0(obj, a, b, c)
    % I0 = ∫∫ N1^a N2^b N3^c ds = 2A(a!b!c!)/(a + b + c + 2)!
      i = 2 * obj.getArea ...
        .* factorial(a) .* factorial(b) .* factorial(c) ...
        ./ factorial(a + b + c + 2);
    end

    function [ i ] = getSelfTermIi(obj, m, n)
    % Principal value term Ii,mn = ∫∫fm*(nxfn)ds, where m,n = 1,2,3.
    %
    %           /                     \
    %          |    0   -L1*L2  L1*L3  |
    % Ii = 1/6 |  L2*L1    0   -L2*L3  |
    %          | -L3*L1  L3*L2    0    |
    %           \                     /
      onesixth = 1/6;
      d = [0, -1,  1;
           1,  0, -1;
          -1,  1,  0];
      i = d(m,n) * onesixth * obj.getEdgeLength(m) .* obj.getEdgeLength(n);
    end

    function [ i ] = getSelfTermIn(obj, m, n)
    % Principal value term In,mn = ∫∫(nxfm)*(nxfn)ds = ∫∫fm*fnds,
    % where m,n = 1,2,3.
    % The analytical evaluation yields:
    % In,mn =
    % (LmLn/4A)( (1/12)(p1^2+p2^2+p3^2) + 3/4(c^2) - c*(pm + pn) + pm*pn )
      onetwelvth = 1 / 12;
      c  = obj.getCentroid;
      i = (0.25 ...
        * obj.getEdgeLength(m) .* obj.getEdgeLength(n) ./ obj.getArea) ...
        .* ( onetwelvth * (...
        sum(obj.getP(1).^2, 2) ...
        + sum(obj.getP(2).^2, 2) ...
        + sum(obj.getP(3).^2, 2)) ...
        + 0.75 * sum(c.*c, 2)...
        - sum(c .* (obj.getP(m) + obj.getP(n)), 2)...
        + sum(obj.getP(m) .* obj.getP(n), 2) );
    end

    function [ i ] = getSelfTermId(obj)
      % Id = ∫∫{∫∫(1/R)ds'}ds (Arcioni et al)
      sp = 0.5 .* obj.getPerimeter;
      i = -4/3 .* obj.getArea.^2 ...
        .* (1 ./ obj.getEdgeLength(1) ...
        .* log(1 - obj.getEdgeLength(1) ./ sp) ...
        + 1 ./ obj.getEdgeLength(2) ...
        .* log(1 - obj.getEdgeLength(2) ./ sp) ...
        + 1 ./ obj.getEdgeLength(3) ...
        .* log(1 - obj.getEdgeLength(3) ./ sp));
    end

    function [ i ] = getSelfTermIs(obj, m, n)
    % Is = ∫∫{∫∫(fm*fn/R)ds'}ds (Arcioni et al)
      sp = 0.5 * obj.getPerimeter;
      if m == n
        a = obj.getEdgeLength(m);
        b = obj.getEdgeLength(obj.getIdxN(m + 1));
        c = obj.getEdgeLength(obj.getIdxN(m + 2));
        i = obj.getArea.^2 ./ 30 ... 
          .* ((10 + 3.*(c.^2 - a.^2) ./ b.^2 - 3 .* (a.^2 - b.^2) ./ c.^2) .* a ...
          - (5 - 3 .* (a.^2 - b.^2) ./ c.^2 - 2 .* (b.^2 - c.^2) ./ a.^2) .* b ...
          - (5 + 3 .* (c.^2 - a.^2) ./ b.^2 + 2 .* (b.^2 - c.^2) ./ a.^2) .* c ...
          + (a.^2 - 3 .* b.^2 - 3 .* c.^2 - 8 .* obj.getArea.^2 ./ a.^2) .* 2 ./ a .* log(1 - a ./ sp) ...
          + (a.^2 - 2 .* b.^2 - 4 .* c.^2 + 6 .* obj.getArea.^2 ./ b.^2) .* 4 ./ b .* log(1 - b ./ sp) ...
          + (a.^2 - 4 .* b.^2 - 2 .* c.^2 + 6 .* obj.getArea.^2 ./ c.^2) .* 4 ./ c .* log(1 - c ./ sp));
      else
        kj = 6 - m - n;
        a = obj.getEdgeLength(obj.getIdxN(kj));
        b = obj.getEdgeLength(obj.getIdxN(kj + 1));
        c = obj.getEdgeLength(obj.getIdxN(kj + 2));
        i = obj.getArea.^2 ./ 60 ... 
          .* ((-10 + (c.^2 - a.^2) ./ b.^2 - (a.^2 - b.^2) ./ c.^2) .* a ...
          + (5 + (a.^2 - b.^2) ./ c.^2 - 6 .* (b.^2 - c.^2) ./ a.^2) .* b ...
          + (5 - (c.^2 - a.^2) ./ b.^2 + 6 .* (b.^2 - c.^2) ./ a.^2) .* c ...
          + (2 .* a.^2 - b.^2 - c.^2 + 4 .* obj.getArea.^2 ./ a.^2) .* 12 ./ a .* log(1 - a ./ sp) ...
          + (9 .* a.^2 - 3 .* b.^2 - c.^2 + 4 .* obj.getArea.^2 ./ b.^2) .* 2 ./ b .* log(1 - b ./ sp) ...
          + (9 .* a.^2 - b.^2 - 3 .* c.^2 + 4 .* obj.getArea.^2 ./ c.^2) .* 2 ./ c .* log(1 - c ./ sp));
      i = i .* 0.25 .* obj.getEdgeLength(m) .* obj.getEdgeLength(n) ./ obj.getArea.^2; 
      end
    end
  end
end
