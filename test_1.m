clc;
clear;

%直方图统计
image = rgb2gray(imread('noise.bmp'));

figure;imhist(image);

image = double(image);
image_hist = imread('saliency_noise.jpg');

figure;imhist(image_hist);

[image_hist_row, image_hist_col] = size(image_hist);

fore = zeros(1,256);
back = zeros(1,256);

for i = 1:image_hist_row
    for j = 1:image_hist_col
        if image_hist(i,j) > 110
            fore(image(i,j) + 1) = fore(image(i,j) + 1) + 1;
            image(i,j) = 255;
        else
            back(image(i,j) + 1) = back(image(i,j) + 1) + 1;
            image(i,j) = 0;
        end
    end
end

figure,imshow(image);
figure;plot(fore);
figure;plot(back);
m

% syms m d r
% u = solve('(m*u + 2*(m-2)*u*ln(u))*d.^2 + r=0','u')
%  u = solve('(m*u + 2*(m-2)*u*(u-1))*d*d + r=0','u')
 
% syms d m
% r = solve(' (4*d - d*m + (d^2*m^2 - 8*d^2*m + 16*d^2 - 8*r*m + 16*r)^(1/2))/(4*(2*d - d*m)) = 1', 'r')

% syms d1 d2 m
% u = solve('(m*u + 2*(m-2)*u*(u-1))*d1*d1 - d2*d2*m = 0', 'u')
 
%模糊聚类
% image = double(rgb2gray(imread('FLIR0172.jpg')));
% 
% R = raylrnd(image);
% image = min(R+image,255);
a=rgb2gray(imread('FLIR0359.jpg'));  
image=imnoise(a,'gaussian',0.1,0.008);%对图像加入高斯噪声  

figure,imshow(image);title('nosiyImage');


image = double(image);

% image = double(rgb2gray(imread('FLIR0359.jpg')));

% image = double(rgb2gray(imread('FLIR0371.jpg')));

% image = double(rgb2gray(imread('FLIR0178.jpg')));

[image_row, image_col] = size(image);

Y = image(:);

% Y=mapminmax(Y',0,1);
% 
% Y = Y';

tic;

[center,U,obj_fcn] = FCMClust(Y,3); 

U_col = size(U,2);
A = zeros(1,U_col);

%两类
% for i = 1:U_col
%     if U(1,i) > U(2,i)
% %     if U(1,i) > 0.7
%         A(i) = 0;
%     else
%         A(i) = 1;
%     end
% end

%三类
for i = 1:U_col
    if U(1,i) > U(2,i) && U(1,i) > U(3,i)
%     if U(1,i) > 0.7
        A(i) = 0;
    else if  U(2,i) > U(1,i) && U(2,i) > U(3,i)
        A(i) = 0.5;
        else
            A(i) = 1; 
        end
    end
end

toc;

B = reshape(A,image_row,image_col);

figure;imshow(B,[]);


%imwrite(uint8(B),'FLIR0322_uv_nom.jpg');


