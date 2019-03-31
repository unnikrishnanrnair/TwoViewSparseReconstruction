clear all
close all
clc

img1=imread('../images/img1.png');
img2=imread('../images/img2.png');

img1=rgb2gray(img1);
img2=rgb2gray(img2);

pts1=detectSURFFeatures(img1,'MetricThreshold',1000);
pts2=detectSURFFeatures(img2,'MetricThreshold',1000);

%disp(pts1.Location(:,1))

mu1x=mean(pts1.Location(:,1));
mu1y=mean(pts1.Location(:,2));
mu2x=mean(pts2.Location(:,1));
mu2y=mean(pts2.Location(:,2));
mux=(mu1x+mu2x)/2;
muy=(mu1y+mu2y)/2;

distances = hypot(pts1.Location(:,1),pts1.Location(:,2));
avg_dist1 = mean(distances);
distances = hypot(pts2.Location(:,1),pts2.Location(:,2));
avg_dist2 = mean(distances);
d=(avg_dist1+avg_dist2)/2;

scale1=sqrt(2)/avg_dist1;
scale2=sqrt(2)/avg_dist2;

T1=[scale1,0,-scale1*mu1x; 0,scale1,-scale1*mu1y; 0,0,1];
T2=[scale2,0,-scale2*mu2x; 0,scale2,-scale2*mu2y; 0,0,1];

X = ones(length(pts1.Location(:,1)),1);
hpts1=[pts1.Location X];
X = ones(length(pts2.Location(:,1)),1);
hpts2=[pts2.Location X];

spts1=T1*transpose(hpts1);
spts2=T2*transpose(hpts2);

[f1,vpts1] = extractFeatures(img1,pts1);
[f2,vpts2] = extractFeatures(img2,pts2);

indexPairs = matchFeatures(f1,f2);

mpts1 = vpts1(indexPairs(:,1));
mpts2 = vpts2(indexPairs(:,2));

figure; showMatchedFeatures(img1,img2,mpts1,mpts2);

X = ones(length(mpts1.Location(:,1)),1);
mhpts1=[mpts1.Location X];
X = ones(length(mpts2.Location(:,1)),1);
mhpts2=[mpts2.Location X];


tmhpts1=T1*transpose(mhpts1);
tmhpts2=T2*transpose(mhpts2);

F= estimateFundamentalMatrix(mpts1,mpts2,'Method','RANSAC','NumTrials',4000,'DistanceThreshold',1e-4);

[u d v] = svd(F);
new_d = diag([d(1,1) d(2,2) , 0]);
F = u * new_d * transpose(v);

K = [558.7087,0.0000,310.3210;0.0000,558.2827,240.2395;0.0000,0.0000,1.0000];

E=transpose(K)*F*K;

[u d v] = svd(E);
new_d = diag([(d(1,1)+d(2,2))/2, (d(1,1)+d(2,2))/2 , 0]);
E = u * new_d * transpose(v) ;

[R,t]= decomposeEssentialMatrix(E, transpose(mhpts1), transpose(mhpts2), K);

ProjMat_1 = K * [eye(3), [0; 0; 0]];
ProjMat_2 = K * [R,t];

[pts3D] = algebraicTriangulation(transpose(mhpts1), transpose(mhpts2), ProjMat_1, ProjMat_2);
%plotCameraFrustum(ProjMat_1,'Red',1);
%scatter3(pts3D(1,:),pts3D(2,:),pts3D(3,:))
%surf(pts3D)
%pcwrite(transpose(pts3D), 'OutputFileName.ply');

