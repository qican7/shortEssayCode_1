clear;close all;clc;
f=imread('基于贝叶斯估计NSCT域去噪增强.jpg');
f=double(f);
[m,n]=size(f);
leix=zeros(m,n-1);
for i=1:n-1 
    leix(:,i)=abs(f(:,i+1)-f(:,i));
end
% figure;
% imshow(mat2gray(leix));
w1=[1 1 1 1 1;1 1 1 1 1;1 1 1 1 1;1 1 1 1 1;1 1 1 1 1];
U1=imfilter(leix,w1);
avU1=mean(U1(:));
c=find(U1(:)>avU1);
