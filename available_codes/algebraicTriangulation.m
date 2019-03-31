function [pts3D] = algebraicTriangulation(pts2D_1, pts2D_2, ProjMat_1, ProjMat_2)

ProjMat_1=[ProjMat_1(1,1), ProjMat_1(1,2), ProjMat_1(1,4); ProjMat_1(2,1), ProjMat_1(2,2), ProjMat_1(2,4); ProjMat_1(3,1), ProjMat_1(3,2), ProjMat_1(3,4)];
ProjMat_2=[ProjMat_2(1,1), ProjMat_2(1,2), ProjMat_2(1,4); ProjMat_2(2,1), ProjMat_2(2,2), ProjMat_2(2,4); ProjMat_2(3,1), ProjMat_2(3,2), ProjMat_2(3,4)];
pts3D=inv(ProjMat_1)*pts2D_1;
temp=inv(ProjMat_2)*pts2D_2;
pts3D=[pts3D temp];
%pts3D=transpose(pts3D);
end