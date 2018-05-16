function gabor = makeGabor (imSize, lamda, theta, sigma, phase, trim)
% function gabor = makeGabor (imSize, lambda, theta, sigma, phase, trim)
% 
% http://www.icn.ucl.ac.uk/courses/MATLAB-Tutorials/Elliot_Freeman/html/gabor_tutorial.html
% parameters
% imSize = 100;   % image size: n X n
% lamda = 10;     % wavelength (number of pixels per cycle)
% theta = 15;     % grating orientation
% sigma = 10;     % gaussian standard deviation in pixels
% phase = .25;    % phase (0 -> 1)
% trim = .005;    % trim off gaussian values smaller than this
if nargin < 1 || isempty(imSize)
	imSize = 100;
end

if nargin < 2 || isempty(lamda)
	lamda = 10;
end

if nargin < 3 || isempty(theta)
	theta = 15;
end

if nargin < 4 || isempty(sigma)
	sigma = 10;
end

if nargin < 5 || isempty(phase)
	phase = .25;
end

if nargin < 6 || isempty(trim)
	trim = .005;
end

% frequency
freq = imSize/lamda;                    % compute frequency from wavelength
phaseRad = (phase * 2* pi);             % convert to radians: 0 -> 2*pi

% linear ramp
X = 1:imSize;                           % X is a vector from 1 to imageSize
X0 = (X / imSize) - .5;                 % rescale X -> -.5 to .5

% make a 2D grating
[Xm, Ym] = meshgrid(X0, X0);             % 2D matrices

% change orientation by adding Xm and Ym together in different proportions
thetaRad = (theta / 360) * 2*pi;        % convert theta (orientation) to radians
Xt = Xm * cos(thetaRad);                % compute proportion of Xm for given orientation
Yt = Ym * sin(thetaRad);                % compute proportion of Ym for given orientation
XYt = Xt + Yt;							% sum X and Y components
XYf = XYt * freq * 2*pi;                % convert to radians and scale by frequency
grating = sin( XYf + phaseRad);         % make 2D sinewave

% make a Gaussian mask
s = sigma / imSize;                     % gaussian width as fraction of imageSize
gauss = exp( -(((Xm.^2)+(Ym.^2)) ./ (2* s^2)) ); % formula for 2D gaussian

% multply grating and gaussian to get a GABOR
gauss(gauss < trim) = 0;                 % trim around edges (for 8-bit colour displays)
gabor = grating .* gauss;                % use .* dot-product