function viewer(p,t)
  mtris = t(t(:,4)==1,1:3);
  dtris = t(t(:,4)==2,1:3);
  figure('Name','SurfaceMesh');hold on;
  if length(mtris)~=0
    trisurf(mtris, p(:,1), p(:,2), p(:,3),'FaceColor', [1.0,1.0,0.8]);
  end
  if length(dtris)~=0
    trisurf(dtris, p(:,1), p(:,2), p(:,3),'FaceColor', [0.0,1.0,1.0]);
  end
  axis equal; axis off; rotate3d on;hold off;
end
