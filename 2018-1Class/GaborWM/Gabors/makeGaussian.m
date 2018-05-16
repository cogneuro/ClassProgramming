function gauss = makeGaussian (imSize, sigma, trim)
% function gabor = makeGaussian (imSize, sigma, trim)
% 
% http://www.icn.ucl.ac.uk/courses/MATLAB-Tutorials/Elliot_Freeman/html/gabor_tutorial.html
% parameters
% imSize = 100;   % image size: n X n
% sigma = 10;     % gaussian standard deviation in pixels
% trim = .005;    % trim off gaussian values smaller than this
if nargin < 1 || isempty(imSize)
	imSize = 100;
end

if nargin < 2 || isempty(sigma)
	sigma = 10;
end

if nargin < 3 || isempty(trim)
	trim = .005;
end

% linear ramp
X = 1:imSize;                           % X is a vector from 1 to imageSize
X0 = (X / imSize) - .5;                 % rescale X -> -.5 to .5

% make a 2D grating
[Xm, Ym] = meshgrid(X0, X0);             % 2D matrices

% make a Gaussian mask
s = sigma / imSize;                     % gaussian width as fraction of imageSize
gauss = exp( -(((Xm.^2)+(Ym.^2)) ./ (2* s^2)) ); % formula for 2D gaussian

% multply grating and gaussian to get a GABOR
gauss(gauss < trim) = 0;                 % trim around edges (for 8-bit colour displays)
