classdef Edge
%EDGE Linear edge object.
% Complete set of geometric calculations for linear edges.
%
%    p1 o
%        \
%         \ L
%      E   \
%           \
%            o p2
%
%--------------------------------------------------------------------------
% Copyright: Klearchos A. Samaras, kleasama@gmail.com, June 2023
%--------------------------------------------------------------------------

  properties
    % The coordinates of points 1 and 2 in columns 1 - 6.
    p
    % Total number of edges processed by the Edge object.
    ne
    % Container for the edge vector.
    e
    % Container for the edge length.
    L
    % Container for the edge unit vector (from p1 to p2).
    s
    % Container for the midpoints of the edges.
    mp
  end
  methods
    function [ obj ] = Edge( point1, point2 )
    %%Edge Constructor of edge object.
      obj.p  = [point1, point2];
      obj.ne = obj.calcEdgesTotal;
      obj.e  = obj.calcEdgeVector;
      obj.L  = obj.calcEdgeLength;
      obj.s  = obj.calcEdgeUnit;
      obj.mp = obj.calcEdgeMidPoint;
    end
    
    function [ obj ] = processEdge( obj )
    %%processEdge Processes the edge object.
      obj.ne = obj.calcEdgesTotal;
      obj.e  = obj.calcEdgeVector;
      obj.L  = obj.calcEdgeLength;
      obj.s  = obj.calcEdgeUnit;
      obj.mp = obj.calcEdgeMidPoint;
    end
    
    function [ ne ] = calcEdgesTotal(obj)
    %calcEdgesTotal Calculates the number of edges processed by the
    % edge object.
      ne = length(obj.p(:,1));
    end
    
    function [ ne ] = getEdgesTotal(obj)
    %%getEdgesTotal Returns number of edges processed by the
    % edge object.
      ne = obj.ne;
    end
   
    function [ pj ] = getP(obj, j)
    %%getP Returns the coordinates of the j-th local node.
      pj = obj.p(:,(j-1)*3+1:(j*3));
    end
    
    function [ e ] = calcEdgeVector(obj)
    %%calcEdgeVector Calculates the edge vector.
      e = obj.getP(2) - obj.getP(1);
    end
    
    function [ e ] = getEdgeVector(obj)
    %%getEdgeVector Returns the edge vector.
      e = obj.e;
    end

    function [ l ] = calcEdgeLength(obj)
    %%calcEdgeLength Calculates the edge length.
      l = sqrt(sum(obj.getEdgeVector.^2, 2));
    end
    
    function [ l ] = getEdgeLength(obj)
    %%getEdgeLength Returns the edge length.
      l = obj.L;
    end

    function [ s ] = calcEdgeUnit(obj)
    %%calcEdgeUnit Calculates the edge tangential unit vector.
      s = obj.getEdgeVector ./ obj.getEdgeLength;
    end
    
    function [ s ] = getEdgeUnit(obj)
    %%getEdgeUnit Returns the edge tangential unit vector.
      s = obj.s;
    end

    function [ mp ] = calcEdgeMidPoint(obj)
    %%calcEdgeMidPoint Calculates the edge mid point coordinates.
      mp = 0.5 * (obj.getP(1) + obj.getP(2));
    end
    
    function [ mp ] = getEdgeMidPoint(obj)
    %%getEdgeMidPoint Returns the edge mid point coordinates.
      mp = obj.mp;
    end
    
    function [ sm ] = getSm(obj, r)
      sm = sum((obj.getP(1) - r) .* obj.getEdgeUnit, 2);
    end

    function [ sp ] = getSp(obj, r)
      sp = obj.getSm(r) + obj.getEdgeLength;
    end

    function [ rp ] = getRp(obj, r)
      rp = sqrt(sum((r - obj.getP(2)).^2, 2));
    end

    function [ rm ] = getRm(obj, r)
      rm = sqrt(sum((r - obj.getP(1)).^2, 2));
    end
    
    function [ r0 ] = getR0(obj, r)
      l0 = sum(obj.getEdgeUnit .* (r - obj.getP(1)),2);
      p0 = obj.getP(1) + obj.getEdgeUnit .* l0;
      r0 = sqrt(sum((r - p0).^2, 2));
    end
    
    function [ im1 ] = getSingularTermIm1(obj, r)
    % Initial value for contour singular integral recursive formulas
    % associated with the linear edge.
      im1 = log((obj.getSp(r) + obj.getRp(r)) ...
          ./(obj.getSm(r) + obj.getRm(r)));
    end
    
    function [ in, term_count ] = getSingularTermI(obj, r, n)
    % Recursive formula for the contour integrals Im.
    % r - field point
    % n - power of distance R, n = -1,1,3,5...
      if (n == -1)
        % If n == -1 then we are returning the inital term Im1.
        in = obj.getSingularTermIm1(r);
        return;
      end
      % Recursion.
      for ni = 1 : 2 : n
        % Indexing for the singular term arrays.
        term_count = (ni + 1) / 2 + 1;
        in = 1./(ni + 1) ...
            .* ( obj.getSp(r) .* obj.getRp(r).^ni ...
            - obj.getSm(r) .* obj.getRm(r).^ni ...
            + ni .* obj.getR0(r).^2 ...
            .* obj.getSingularTermI(r, ni - 2));
      end
    end
  end
end