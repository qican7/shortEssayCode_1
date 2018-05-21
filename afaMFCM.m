function [ImS,S]=afaMFCM(FImg,U)
[M,N]=size(FImg);
% X=zeros(M+2,N+2);
% X(2:M+1,2:N+1)=FImg;
X=padarray(FImg,[1,1],'symmetric');
% X=padarray(FImg,[2,2],'symmetric');
% X=zeros(M+4,N+4);X(3:M,3:N)=FImg;
d1=zeros(M+2,N+2);
% d1=zeros(M+4,N+4);
d1(X==1)=0;
d1(X~=1)=1;
d2=zeros(M+2,N+2);
% d2=zeros(M+4,N+4);
d2(X==2)=0;
d2(X~=2)=1;
d3=zeros(M+2,N+2);%//
% d3=zeros(M+4,N+4);
d3(X==3)=0;%//
d3(X~=3)=1;%//
w1=[1 1 1;1 0 1;1 1 1];
%  w1=[1 1 1 1 1;1 1 1 1 1;1 1 0 1 1;1 1 1 1 1;1 1 1 1 1];
U1=imfilter(d1,w1);
U2=imfilter(d2,w1);
U3=imfilter(d3,w1); %//
eU1=exp(-U1);
eU2=exp(-U2);
eU3=exp(-U3);%//
eU=eU1+eU2+eU3; %//
% eU=eU1+eU2; %//
P1=eU1./eU;
P2=eU2./eU;
P3=eU3./eU; %//
PP1=P1(2:M+1,2:N+1)';
PP2=P2(2:M+1,2:N+1)';
PP3=P3(2:M+1,2:N+1)'; %//
% PP1=P1(3:M+2,3:N+2)';
% PP2=P2(3:M+2,3:N+2)';
% PP3=P3(3:M+2,3:N+2)';
P=zeros(3,M*N);%//
P(1,:)=PP1(:);
P(2,:)=PP2(:);
P(3,:)=PP3(:); %//
S=zeros(3,M*N);%//
w2=ones(3)/9;
% w2=ones(5)/9;
mX=imfilter(X,w2);

mmX=mX(2:M+1,2:N+1);
%   mmX=mX(3:M+2,3:N+2);
temp=abs(FImg-mmX);

Af=zeros(M,N);
inAf=ones(M,N);
for i=2:M+1
    for j=2:N+1
% for i=3:M+2
%     for j=3:N+2
        Wd=X(i-1:i+1,j-1:j+1);
%           Wd=X(i-2:i+2,j-2:j+2);  
        fX(i-1,j-1)=(std(Wd(:),1))^2;
%           fX(i-2,j-2)=(std(Wd(:),1))^2;
        if(temp(i-1,j-1)< fX(i-1,j-1))
           Af(i-1,j-1)=fX(i-1,j-1);
        end
%          if(temp(i-2,j-2)< fX(i-2,j-2))
%              Af(i-2,j-2)=fX(i-2,j-2);
%          end
    end
end


% Af(temp<fX)=fX;
MMax=max(Af(:));
Aftt=Af/MMax;
% Aftt(Aftt==0)=1;
% Af(Af~=0)=1-Af;   
% mmX=mX(2:M+1,2:N+1);
% temp=abs(FImg-mmX);
% pp=0.85;

% Af=zeros(M,N);
% Af(temp<fX)=pp;
% Af(temp<fX)=;

% Aftt=rand(M,N);
Aft=Aftt';
Aff=Aft(:)';
% prem1=0.01;%噪声0.4,0.6模拟图像
prem1=0.01;
% prem2=100;
S(1,:)=prem1*Aff.*U(1,:)+(1-prem1*Aff).*P(1,:);
S(2,:)=prem1*Aff.*U(2,:)+(1-prem1*Aff).*P(2,:);
S(3,:)=prem1*Aff.*U(3,:)+(1-prem1*Aff).*P(3,:);%//

% S(1,:)=(1-Aff).*U(1,:)+Aff.*P(1,:);
% S(2,:)=(1-Aff).*U(2,:)+Aff.*P(2,:);
% S(3,:)=(1-Aff).*U(3,:)+Aff.*P(3,:);
% S(1,:)=0.3.*U(1,:)+0.7.*P(1,:);
% S(2,:)=0.3.*U(2,:)+0.7.*P(2,:);
% S(3,:)=0.3.*U(3,:)+0.7.*P(3,:);

%给它标签
SL=zeros(1,size(U,2));
for i=1:size(U,2)
    [va,SL(i)]=max(S(:,i));
end

%看看原图中标签分布
ImS=zeros(M,N);
for i=1:M
    for j=1:N
        ImS(i,j)=SL((i-1)*N+j);
    end
end

