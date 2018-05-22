clear;close all;clc;
f=imread('FLIR0172.jpg');
f=rgb2gray(f);
f=imnoise(f,'gaussian',0.02,0.004);%对图像加入高斯噪声
figure,imshow(f);
ft=f';
[M,N]=size(f);
Data=zeros(M*N,1);
Data=ft(:);

tic;
%FCM聚类
% [U,P,Dist,Cluster_Res,Obj_Fcn,iter]=fuzzycm(Data,3,0,2,1.0e-4);
[U,P,Dist,Cluster_Res,Obj_Fcn,iter]=imhistfuzzycm(Data,3,0,2,1.0e-4);
% [U,P,Dist,Obj_Fcn,iter]=GGspconimhistfuzzycm(f,Data,2,0,2,1.0e-4);%//
% [U,P,Dist,Obj_Fcn,iter]=GGspconimhistfuzzycm_3(f,Data,3,0,2,1.0e-4);

%给它标签
UL=zeros(1,size(U,2));

for i=1:size(U,2)
    [va,UL(i)]=max(U(:,i));
end

fdb=double(f);
fnew=zeros(M,N);
for i=1:3 %//
    Idx=find(UL==i)-1;
    KIdx=ismember(fdb,Idx);
    fnew=fnew+i*KIdx;
end
figure,imshow(fnew,[]);

fU=zeros(3,M*N); %//
fUt=zeros(M,N);
for i=1:3 %//
    for j=1:256
        fUt(fdb==(j-1))=U(i,j);
    end
    fU(i,:)=fUt(:);
end

ImG=fnew;
[afaImS,afaSL]=afaMFCM(ImG,fU);
% [afaImS,afaSL]=ulsound_2cls_aMFCM(ImG,fU);

figure,imshow(afaImS,[]);
toc;




