function mr=maxrowf(U,c)
% ����� U ÿ�е� c ��Ԫ��������,c ��ȱʡֵΪ 1
% ���ø�ʽ: mr = maxrowf(U,c)
% See also: addr
if nargin<2
    c=1;
end
N=size(U,2);mr(1,N)=0;
for j=1:N
    aj=addr(U(:,j),'descend');
    mr(j)=aj(c);
end
