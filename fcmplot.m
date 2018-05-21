function fcmplot(Data,U,P,Obj_Fcn)
% FCM �����ͼ����
% See also: fuzzycm maxrowf ellipse
[C,S] = size(P); res = maxrowf(U);
str = 'po*x+d^v><.h'; 
% Ŀ�꺯����ͼ
figure(1),plot(Obj_Fcn)
title('Ŀ�꺯��ֵ�仯����','fontsize',8)
% 2D ��ͼ
if S==2 
    figure(2),plot(P(:,1),P(:,2),'rs'),hold on
    for i=1:C
        v=Data(find(res==i),:); 
        plot(v(:,1),v(:,2),str(rem(i,12)+1))      
        ellipse(max(v(:,1))-min(v(:,1)), ...
                max(v(:,2))-min(v(:,2)), ...
                [max(v(:,1))+min(v(:,1)), ...
                max(v(:,2))+min(v(:,2))]/2,'r:')    
    end
    grid on,title('2D ������ͼ','fontsize',8),hold off
end
% 3D ��ͼ
if S>2 
    figure(2),plot3(P(:,1),P(:,2),P(:,3),'rs'),hold on
    for i=1:C
        v=Data(find(res==i),:);
        plot3(v(:,1),v(:,2),v(:,3),str(rem(i,12)+1))      
        ellipse(max(v(:,1))-min(v(:,1)), ...
                max(v(:,2))-min(v(:,2)), ...
                [max(v(:,1))+min(v(:,1)), ...
                max(v(:,2))+min(v(:,2))]/2, ...
                'r:',(max(v(:,3))+min(v(:,3)))/2)   
    end
    grid on,title('3D ������ͼ','fontsize',8),hold off
end
