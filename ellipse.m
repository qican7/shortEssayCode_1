function ellipse(a,b,center,style,c_3d)
% ����һ����Բ
% ����: ellipse(a,b,center,style,c_3d)
% ����:
%     a: ��Բ���᳤(ƽ���� x ��)
%     b: ��Բ���᳤(ƽ���� y ��)
%     center: ��Բ������ [x0,y0],ȱʡֵΪ [0,0]
%     style: ���Ƶ����ͺ���ɫ,ȱʡֵΪʵ����ɫ
%     c_3d:   ��Բ�������� 3D �ռ��е� z ������,��ȱʡ
if nargin<4
    style='b';
end
if nargin<3 | isempty(center)
    center=[0,0];
end
t=1:360;
x=a/2*cosd(t)+center(1);
y=b/2*sind(t)+center(2);
if nargin>4
    plot3(x,y,ones(1,360)*c_3d,style)
else
    plot(x,y,style)
end
