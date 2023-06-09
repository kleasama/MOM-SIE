%% geometric precalculations over triangles
GP      = GPC(p(t(:,1),:),p(t(:,2),:),p(t(:,3),:));
[r ,w ] = GQ(p,t,nP);
[r_,w_] = GQ(p,t,nQ);
[R,D]   = DC(r,r_);
w_      = permute(w_,[2,1,4,3]);
F1      = permute(0.5*GP(:,2).*(r - p(t(:,1),:))./GP(:,1),[1,5,3,4,2]);
F2      = permute(0.5*GP(:,3).*(r - p(t(:,2),:))./GP(:,1),[1,5,3,4,2]);
F3      = permute(0.5*GP(:,4).*(r - p(t(:,3),:))./GP(:,1),[1,5,3,4,2]);
F1_     = permute(0.5*GP(:,2).*(r_- p(t(:,1),:))./GP(:,1),[5,1,4,3,2]);
F2_     = permute(0.5*GP(:,3).*(r_- p(t(:,2),:))./GP(:,1),[5,1,4,3,2]);
F3_     = permute(0.5*GP(:,4).*(r_- p(t(:,3),:))./GP(:,1),[5,1,4,3,2]);
DF1     = GP(:,2)./GP(:,1);
DF2     = GP(:,3)./GP(:,1);
DF3     = GP(:,4)./GP(:,1);
DxF1_   = cross(D,repmat(F1_,[length(t(:,1)),1,nP,1]));
DxF2_   = cross(D,repmat(F2_,[length(t(:,1)),1,nP,1]));
DxF3_   = cross(D,repmat(F3_,[length(t(:,1)),1,nP,1]));
clear r r_