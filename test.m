f = imread('test_1.bmp');

[m,n] = size(f);
figure,imshow(f);
for i = 1:m
    for j = 1:n
        if f(i,j) > 200
            f(i,j) = 255;
        else
            f(i,j) = 0;
        end
    end
end

%se=strel('square',5');%���ͽṹԪ��
se=strel('disk',7');%Բ���ͽṹԪ��
f=imopen(f,se);%ֱ�ӿ�����
f=imopen(f,se);%ֱ�ӿ�����
% f=imopen(f,se);%ֱ�ӿ�����
% f=imopen(f,se);%ֱ�ӿ�����
% f=imclose(f,se);%ֱ�ӱ�����

% f=imopen(f,se);%ֱ�ӿ�����
% f=imopen(f,se);%ֱ�ӿ�����

figure,imshow(f);
imwrite(f,'test_1.bmp');