%%
%选择目标中心
clear;close all;clc;
f = imread('noise_0359.bmp');

f = rgb2gray(f);

%高斯滤波
sigma = 1.8;      %sigma赋值  
N = 3;            %大小是（2N+1）×（2N+1）  
N_row = 2*N+1;  
gausFilter = fspecial('gaussian',[N_row N_row],sigma);      %matlab 自带高斯模板滤波  
f  = imfilter(f,gausFilter,'conv');  
figure,imshow(f);

[m,n] = size(f);

f = double(f);

f_mean = mean(mean(f));

f_variance = (f - f_mean).^2;
f_variance_mean = mean(mean(f_variance));
f_variance_max = max(max(f_variance));
f_variance_min = min(min(f_variance));


f_D = (f_variance - f_variance_min)./(f_variance_max - f_variance_min);



for i = 1:m
    for j = 1:n
        if f_D(i,j) > 1
            f_D(i,j) = 1;
        end
    end
end 

f_D_mean = mean(mean(f_D));
f_D_variance_standardDeviation = sqrt((f_D - f_D_mean).^2);

threshold = f_D_mean + 0.55 * f_D_variance_standardDeviation;

for i = 1:m
    for j = 1:n
        if f_D(i,j) > threshold
            f_D(i,j) = 255;
        else
            f_D(i,j) = 0;
        end
    end
end 

f_D = uint8(f_D);
figure,imshow(f_D);

%% 
%利用外接矩形获取目标的初始中心
left_index = n;
right_index = 0;
up_index = m;
bottom_index = 0;
for i = 1:m
    for j = 1:n
        if f_D(i,j) == 255
            f_D(i,j) = 255;
            if i < up_index
                up_index = i;
            end
            if i > bottom_index
                bottom_index = i;
            end
            if j < left_index
                left_index = j;
            end
            if j > right_index
                right_index = j;
            end
        end
    end
end 

center_row = uint8((bottom_index + up_index)./2);
center_col = uint8((left_index + right_index)./2);

diagonal_length = sqrt((right_index - left_index).^2 + (bottom_index - up_index).^2);


f_D(center_row,center_col) = 0;
figure,imshow(f_D);

%%
%形态学处理
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