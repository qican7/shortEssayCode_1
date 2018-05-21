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

%se=strel('square',5');%方型结构元素
se=strel('disk',7');%圆盘型结构元素
f=imopen(f,se);%直接开运算
f=imopen(f,se);%直接开运算
% f=imopen(f,se);%直接开运算
% f=imopen(f,se);%直接开运算
% f=imclose(f,se);%直接闭运算

% f=imopen(f,se);%直接开运算
% f=imopen(f,se);%直接开运算

figure,imshow(f);
imwrite(f,'test_1.bmp');