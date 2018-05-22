% function [U,P,Dist,Cluster_Res,Obj_Fcn,iter]=GGspconimhistfuzzycm(f,Data,C,plotflag,M,epsm)
function [U,P,Dist,Obj_Fcn,iter]=GGspconimhistfuzzycm_3(f,Data,C,plotflag,M,epsm)
% 模糊 C 均值聚类 FCM: 从随机初始化划分矩阵开始迭代
% [U,P,Dist,Cluster_Res,Obj_Fcn,iter] = fuzzycm(Data,C,plotflag,M,epsm)
% 输入:
%     Data: N×S 型矩阵,聚类的原始数据,即一组有限的观测样本集,
%           Data 的每一行为一个观测样本的特征矢量,S 为特征矢量
%           的维数,N 为样本点的个数
%     C:    聚类数,1<C<N
%     plotflag: 聚类结果 2D/3D 绘图标记,0 表示不绘图,为缺省值
%     M:    加权指数,缺省值为 2
%     epsm: FCM 算法的迭代停止阈值,缺省值为 1.0e-6
% 输出:
%     U:    C×N 型矩阵,FCM 的划分矩阵
%     P:    C×S 型矩阵,FCM 的聚类中心,每一行对应一个聚类原型
%     Dist: C×N 型矩阵,FCM 各聚类中心到各样本点的距离,聚类中
%           心 i 到样本点 j 的距离为 Dist(i,j)
%     Cluster_Res: 聚类结果,共 C 行,每一行对应一类
%     Obj_Fcn: 目标函数值
%     iter: FCM 算法迭代次数
% See also: fuzzydist maxrowf fcmplot
if nargin<5
%     epsm=1.0e-6;
      epsm=1.0e-4;
end
if nargin<4
    M=2;
end
if nargin<3
    plotflag=0;
end
% f=255*mat2gray(f);
w=[1/(1+sqrt(2)),0.5,1/(1+sqrt(2)),0.5,0,0.5,1/(1+sqrt(2)),0.5,1/(1+sqrt(2))];
fdoub=double(f);

[f_m, f_n] = size(f);

start=min(fdoub(:));
finish=max(fdoub(:));
fdoub1=fdoub;
fdoub2=fdoub;
fdoub3=fdoub;%//

HData=imhist(uint8(Data));
NN=size(HData,1);
HData=HData';
HDatanew=HData(ones(C,1),:);

% m=1/(M-1);
iter=0;
Dist(C,NN)=0; U(C,NN)=0; P(C,1)=0;
 UTEMP(C,NN)=0;
% 随机初始化划分矩阵
U0 = rand(C,NN);
% U0=U0./(ones(C,1)*sum(U0));
% FCM 的迭代算法
KK=0:255;
KKN=KK(ones(C,1),:);

P = [0.3;0.5;0.7];

while true
    
     %初始目标中心
    row = 111;
    col = 99;
    diagonal_len = 59;
    Z1 = zeros(size(f));
    Z2 = zeros(size(f));
    Z3 = zeros(size(f));
    Z1_1 = zeros(1,256);
    Z2_1 = zeros(1,256);
    Z3_1 = zeros(1,256);
    for i = 1:f_m
        for j = 1:f_n
            Z1(i,j) = diagonal_len./sqrt((i - row).^2 + (j - col).^2);
            Z2(i,j) = (diagonal_len./sqrt((i - row).^2 + (j - col).^2)).^(0.5);
            Z3(i,j) = (diagonal_len./sqrt((i - row).^2 + (j - col).^2)).^(1./3);
        end
    end
    for i=0:255
       Idxx=find(fdoub==i);
        if(isempty(Idxx))
             Z1_1(i+1)=0;
             Z2_1(i+1)=0;
             Z3_1(i+1)=0;
        else
             Z1_1(i+1)=sum(Z1(Idxx));
             Z2_1(i+1)=sum(Z2(Idxx));
             Z3_1(i+1)=sum(Z3(Idxx));
        end
    end
    
    Z(1,:) = Z1_1(:);
    Z(2,:) = Z2_1(:);
    Z(3,:) = Z3_1(:);
    
    Z(find(isnan(Z)==1)) = 1;
    Z(Z==inf)=1;
    Z(1,:) = Z(1,:)./max(Z);
    Z(2,:) = Z(2,:)./max(Z);
    Z(3,:) = Z(3,:)./max(Z);
    
    
    % 迭代计数器
    iter=iter+1;
    % 计算或更新聚类中心 P
    Um=(U0.*Z).^M;
%     P=sum(Um.*HDatanew.*KKN,2)./sum(Um.*HDatanew,2);
    
   
    % 更新划分矩阵 U
    for i=1:C
        for j=1:NN
            %             Dist(i,j)=fuzzydist(P(i,:),Data(j,:));
            Dist(i,j)=fuzzydist(P(i,:),j-1);
        end
    end
  
    Uqufa=(ones(3,256)- Um).^2; %//
    U1=Uqufa(1,:);
    U2=Uqufa(2,:);
    U3=Uqufa(3,:);%//
    
    for i=start:finish
    Idx=find(fdoub==i);
    fdoub1(Idx)=U1(i+1);
    fdoub2(Idx)=U2(i+1);
    fdoub3(Idx)=U3(i+1);%//
    end
    
    diff_cen1=(fdoub-P(1,:)).^2;
    diff_cen2=(fdoub-P(2,:)).^2;
    diff_cen3=(fdoub-P(3,:)).^2;%//
     G1=fdoub1.*diff_cen1;
     G2=fdoub2.*diff_cen2; 
     G3=fdoub3.*diff_cen3;%//
     
    G1_sp=imfilter(G1,w);
    G2_sp=imfilter(G2,w);
    G3_sp=imfilter(G3,w);%//
    
%     bw1=edge(fdoub1);
%     [x1,y1]=find(bw1==1);
%     for ii=1:length(x1)
%       fdoub1_sp(x1(ii),y1(ii))=fdoub1(x1(ii),y1(ii));
%     end
%     
%      bw2=edge(fdoub2);
%     [x2,y2]=find(bw2==1);
%     for ii=1:length(x2)
%       fdoub2_sp(x2(ii),y2(ii))=fdoub2(x2(ii),y2(ii));
%     end
%     
%      bw3=edge(fdoub3);
%     [x3,y3]=find(bw3==1);
%     for ii=1:length(x3)
%       fdoub3_sp(x3(ii),y3(ii))=fdoub3(x3(ii),y3(ii));
%     end
%     U1new=U1;
%     U2new=U2;
%     U3new=U3;

    for i=0:255
       Idxx=find(fdoub==i);
        if(isempty(Idxx))
             sumG1(i+1)=0;
             sumG2(i+1)=0;
             sumG3(i+1)=0;%//
        else
             sumG1(i+1)=sum(G1_sp(Idxx));
             sumG2(i+1)=sum(G2_sp(Idxx));
             sumG3(i+1)=sum(G3_sp(Idxx));%//
        end
    end
%     fdoub1_spT=fdoub1_sp';
%     fdoub2_spT=fdoub2_sp';
%     fdoub3_spT=fdoub3_sp';
% % 
%     U1new=reshape(fdoub1_spT,1,41895);
%     U2new=reshape(fdoub2_spT,1,41895);
%     U3new=reshape(fdoub3_spT,1,41895);
    sumG(1,:)=sumG1(:);  
    sumG(2,:)=sumG2(:); 
    sumG(3,:)=sumG3(:);%//
%     sumG=sumG;
%    U=1./((Dist.^M+sumG).*(ones(C,1)*sum((Dist.^2+sumG).^(-1))));%//


    
    U=((Z.^M.*Dist.^2+sumG).*(ones(C,1)*sum((Z.^M.*Dist.^2+sumG).^(-1)))).^(-1./(M-1));
    
 
    % 目标函数值: 类内加权平方误差和
    if nargout>4 | plotflag
       Obj_Fcn(iter)=sum(sum(Um.*(Dist.^2).*HDatanew))+sum(sumG(:)) + 3*M;
    end
    fprintf('FCM:Iteration count = %d, obj. fcn = %f, expo=%f\n', iter, Obj_Fcn(iter),M); 
    
    P=sum(Um.*HDatanew.*KKN,2)./sum(Um.*HDatanew,2);
    M = sum(sum(log(-3./(log(U.*Z).*Dist.^2))))./sum(sum(log(U.*Z)));
    if isnan(M)
        M = 2;
    end
    % FCM 算法迭代停止条件
    if (norm(U-U0,Inf))<epsm||(iter==50)
        break
    end
    U0=U;
%     U=imfilter(U,w);
%     U=U./(ones(C,1)*sum(U));
%     U0=U;   
end
% 聚类结果
% if nargout > 3
%     res = maxrowf(U);
%     for c = 1:C
%         v = find(res==c);
%         Cluster_Res(c,1:length(v))=v;
%     end
% end


% 绘图
% if plotflag
%     fcmplot(Data,U,P,Obj_Fcn);
% end
