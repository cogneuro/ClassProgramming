% =======================================
% PRAC 1: Introduction to Gabor
% =======================================

% First, let's visit;
% https://www.cogsci.nl/gabor-generator

imSize = 100;   % image size: n X n
lamda = 20;     % wavelength (number of pixels per cycle)
theta = 15;     % grating orientation
sigma = 10;     % gaussian standard deviation in pixels
phase = .25;    % phase (0 -> 1)
trim = .001;    % trim off gaussian values smaller than this

g = makeGabor(imSize, lamda, theta, sigma, phase, trim);

figure;
imagesc( g, [-1 1] );
axis image;
colormap gray;
truesize;

% figure
% surf(g)
% view(20,120)

%%
% =======================================
% PRAC 2: Original image & filtering
% =======================================

figure;
A = imread('shapes.jpg');
A = rgb2gray(A);
imshow(A);

%%
H = ones(13,13);
H = H/sum(H(:));

imshow(imfilter(A, H))

%%
% =======================================
% PRAC 3: Gabor filtering
% =======================================

%% Gabor bank
wavelength = [50,100];
orientation = [0, 45, 90, 135];
g = gabor(wavelength,orientation);

figure;
subplot(length(wavelength),length(orientation),1)
for p = 1:length(g)
    subplot(length(wavelength),length(orientation),p)
    imshow(real(g(p).SpatialKernel),[]);
    title(sprintf('\\lambda = %d, \\theta = %d',g(p).Wavelength,g(p).Orientation));
end

%% Apply the filters to the image
outMag = imgaborfilt(A,g);

figure;
subplot(length(wavelength),length(orientation),1)
for p = 1:length(g)
    subplot(length(wavelength),length(orientation),p)
    imshow(outMag(:,:,p),[]);
    lambda = g(p).Wavelength;
    theta  = g(p).Orientation;
    title(sprintf('\\lambda = %d, \\theta = %d',lambda,theta));
end

%------ EOF.