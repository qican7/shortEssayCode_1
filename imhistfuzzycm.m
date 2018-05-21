function [U,P,Dist,Cluster_Res,Obj_Fcn,iter]=imhistfuzzycm(Data,C,plotflag,M,epsm)
% ģ�� C ��ֵ���� FCM: �������ʼ�����־���ʼ����
% [U,P,Dist,Cluster_Res,Obj_Fcn,iter] = fuzzycm(Data,C,plotflag,M,epsm)
% ����:
%     Data: N��S �;���,�����ԭʼ����,��һ�����޵Ĺ۲�������,
%           Data ��ÿһ��Ϊһ���۲�����������ʸ��,S Ϊ����ʸ��
%           ��ά��,N Ϊ������ĸ���
%     C:    ������,1<C<N
%     plotflag: ������ 2D/3D ��ͼ���,0 ��ʾ����ͼ,Ϊȱʡֵ
%     M:    ��Ȩָ��,ȱʡֵΪ 2
%     epsm: FCM �㷨�ĵ���ֹͣ��ֵ,ȱʡֵΪ 1.0e-6
% ���:
%     U:    C��N �;���,FCM �Ļ��־���
%     P:    C��S �;���,FCM �ľ�������,ÿһ�ж�Ӧһ������ԭ��
%     Dist: C��N �;���,FCM ���������ĵ���������ľ���,������
%           �� i �������� j �ľ���Ϊ Dist(i,j)
%     Cluster_Res: ������,�� C ��,ÿһ�ж�Ӧһ��
%     Obj_Fcn: Ŀ�꺯��ֵ
%     iter: FCM �㷨��������
% See also: fuzzydist maxrowf fcmplot
if nargin<5
    epsm=1.0e-6;
end
if nargin<4
    M=2;
end
if nargin<3
    plotflag=0;
end

HData=imhist(uint8(Data));
NN=size(HData,1);
HData=HData';
HDatanew=HData(ones(C,1),:);

m=2/(M-1);iter=0;
Dist(C,NN)=0; U(C,NN)=0; P(C,1)=0;
w=ones(3,3)/9;
% �����ʼ�����־���
U0 = rand(C,NN);
U0=U0./(ones(C,1)*sum(U0));
% FCM �ĵ����㷨
KK=0:255;
KKN=KK(ones(C,1),:);
while true
    % ����������
    iter=iter+1;
    % �������¾������� P
    Um=U0.^M;
    P=sum(Um.*HDatanew.*KKN,2)./sum(Um.*HDatanew,2);
    % ���»��־��� U
    for i=1:C
        for j=1:NN
            %             Dist(i,j)=fuzzydist(P(i,:),Data(j,:));
            Dist(i,j)=fuzzydist(P(i,:),j-1);
        end
    end
    U=1./(Dist.^m.*(ones(C,1)*sum(Dist.^(-m))));
    % Ŀ�꺯��ֵ: ���ڼ�Ȩƽ������
    if nargout>4 | plotflag
        Obj_Fcn(iter)=sum(sum(Um.*(Dist.^2).*HDatanew));
    end
    % FCM �㷨����ֹͣ����
    if norm(U-U0,Inf)<epsm
        break
    end
    U0=U;
%     U=imfilter(U,w);
%     U=U./(ones(C,1)*sum(U));
%     U0=U;   
end
% ������
if nargout > 3
    res = maxrowf(U);
    for c = 1:C
        v = find(res==c);
        Cluster_Res(c,1:length(v))=v;
    end
end
% ��ͼ
if plotflag
    fcmplot(Data,U,P,Obj_Fcn);
end
