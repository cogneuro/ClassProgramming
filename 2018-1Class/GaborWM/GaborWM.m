function GaborWM9(varargin)
% 161205. from GaborWM8 class practice

try
    main(varargin{:});
catch me % error handling
    ListenChar(0);
    fclose('all');
    Screen('CloseAll');
    ShowCursor();
    
    fprintf(2, '\n\n???:: %s\n\n', me.message); % main error message
	for k = 1:(length(me.stack) - 1)
        current = me.stack(k);
        fprintf('Error in ==> ');
        fprintf('<a href="matlab: opentoline(''%s'',%d,0)">%s at %d</a>\n\n',...
            current.file, current.line, current.name, current.line);
	end
	fclose('all');
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MAIN EXPERIMENT
function main(varargin)

fclose('all');
format short;
ClockRandSeed;

% Check 'help SyncTrouble' 
Screen('Preference', 'SkipSyncTests', 1); 

global w	% in fact, not necessary in the current script
% ------------------------------------------------- EXPERIMENT INFO
% get user parameters
prompt = {'Enter subject number: ','Initials (e.g. YI)'};
defaults = {'99',''};
answer = inputdlg(prompt, 'Experimental setup information', 1, defaults);
[SN, NM] = deal(answer{:});

% open files and specify file format
baseName = ['GaborWMs' sprintf('%02d', str2double(SN)) NM];
fileName = ['datRaw' baseName '.csv'];	dataFile = fopen(fileName, 'a');
fileName = ['datStm' baseName '.txt'];	stimFile = fopen(fileName, 'a');
paramName =['datVar' baseName '.mat'];

% ------------------------------------------------- KEYBOARD INFO
% We collect mouse responses and use only the 'spacebar' from a keyboard.
Experimenter = -1;			% default keyboard
SPC	= KbName('space');		% spacebar

% -------------------------------------------------- SCREEN SETUP
whichScreen = max(Screen('Screens'));	 % determining a dual monitor setup.
bcolor  = 127;		% gray background

wRect = [0, 0, 800, 600];
w = Screen('OpenWindow', whichScreen, bcolor, wRect);
% [w, wRect] = Screen('OpenWindow', whichScreen);

[cx,cy] = RectCenter(wRect);

% Enable alpha blending with proper blend-function. We need it
% for drawing of smoothed points:
Screen(w,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Removes the blue screen flash and minimize extraneous warnings.
oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);

commandwindow;
% ------------------------------------------------- TIMING SETUP
ifi = Screen('GetFlipInterval', w);
durMEM	 = (round(1.0/ifi)-0.5)*ifi;	% duration of memory sample
durISI	 = (round(1.0/ifi)-0.5)*ifi;	% retention duration
dur500	 = (round(0.5/ifi)-0.5)*ifi;
 
% ------------------------------------------------- STIMULUS SETUP
% Below is for, so called, "Real World Unit."
% See http://www.psychopy.org/general/monitors.html
vWidth	= 27;		% horizontal dimension of viewable screen (cm)
vDist	= 57;		% viewing distance (cm)
ppd = pi * RectWidth(wRect) / atan(vWidth/vDist/2) / 360; % pixels per degree

gSize  = round(ppd * 6);	% gabor size (degree x ppd)
gLamda = round(ppd * 0.5);	% pixels per cycle (degree x ppd)
% gTheta = 25;				% gabor orientation
gPhase = .5;				% gabor phase (0 -> 1)
gSigma = round(ppd * 0.5);	% sd of gaussian aperture (degree x ppd)

gSrcRect = [0,0,1,1] * gSize;

gEccent	= round(ppd * 5);	% eccentricity of gabor (degree x ppd)
gNumLoc	= 6;				% max number of gabors
gXY		= zeros(gNumLoc, 2);% xy coordinates of gabors
for mm=1:gNumLoc
	deg=360/gNumLoc/2 + (360/gNumLoc)*(mm-1);
	gXY(mm,1) = cx + gEccent * sin(deg*pi/180);
	gXY(mm,2) = cy - gEccent * cos(deg*pi/180);
end
gXY = round(gXY);

ThetaRange = 0:.5:179.5;	% range of gabor orientation: IV
RespRange  = circshift(ThetaRange-90, [0, 180]); % target-centered: DV

pSize  = round(ppd * 2.5 );	
pSrcRect = [0,0,1,1] * pSize;
pThick = round(ppd * [.2, .02]);
pColor = 50;

% ------------------------------------------- EXPERIMENTAL DESIGN
SSIZE	= 2;	% small, large
REP		= 50;
nTRIAL  = SSIZE * REP;

xIndex=0:(nTRIAL-1);

data.cSSize = 0:(nTRIAL-1);
data.cSSize(mod(xIndex,SSIZE)==0)=2;
data.cSSize(mod(xIndex,SSIZE)==1)=4;

data.xOrder	= Shuffle(1:nTRIAL);

data.xTargTheta = NaN(1,nTRIAL);
data.xTargLoc	= NaN(1,nTRIAL);

data.dRespTheta0 = NaN(1,nTRIAL);	% raw theta
data.dRespTheta1 = NaN(1,nTRIAL);	% rounded off nearest 0.5
data.dRespTheta2 = NaN(1,nTRIAL);	% target-centered theta (-90 -> 90)
data.dCursorX = NaN(1,nTRIAL);
data.dCursorY = NaN(1,nTRIAL);
	
% -------------------------------------------------- GET STARTED
% Logging
fprintf(		  '\n***** Task begins at %s\n', datestr(now, 0));
fprintf(stimFile, '\n***** Task begins at %s\n', datestr(now, 0));
fprintf(stimFile, 'SN, trial, SS, tTheta, tLoc, rTheta, rounded, centered, x, y\n');
	
xText = sprintf('Press spacebar to begin');
bound = Screen('TextBounds', w, xText);
Screen('DrawText', w, xText, cx - bound(3)/2, cy - bound(4)/2, 0);
Screen('Flip', w);
while 1
 	[keyIsDown, ~, keyCode] = KbCheck(Experimenter);
	if keyIsDown
		if keyCode(SPC)
			break;
		end
	end
end

HideCursor(w);
ListenChar(2);
GetSecs;

Screen('FillRect', w, bcolor);
tflip = Screen('Flip', w);
WaitSecs(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRIAL LOOP
for trial = 1: nTRIAL

 	Screen('FillRect', w, bcolor);
 	Screen('DrawDots', w, [0,0;0,0], [12,10], [255,255,255;0,0,0]', [cx, cy], 1);
 	tflip = Screen('Flip', w, tflip + dur500);
	
	% ------------------------------------------- current condition
	it = data.xOrder(trial);
	
	xSSize = data.cSSize(it);
	clear xTheta xGLoc;
	xTheta = ThetaRange(randi(360, [1,xSSize]));	% 1st = target theta
	data.xTargTheta(it) = xTheta(1);
	
	xGLoc = zeros(1, xSSize);	% 1st = target loc
	if xSSize==2
		xGLoc(1) = randi(gNumLoc);
		xGLoc(2) = mod(xGLoc(1)+2,6)+1; 
	elseif xSSize==3
		xGLoc(1) = randi(gNumLoc);
		xGLoc(2:3) = mod(xGLoc(1)+[1,3],6)+1;
	elseif xSSize==4
		xGLoc(1:2) = randperm(3,2);
		xGLoc(3:4) = randperm(3,2) + 3;
		xGLoc = Shuffle(xGLoc);
	else
		xGLoc  = randperm(gNumLoc, xSSize);
	end
	data.xTargLoc(it) = xGLoc(1);

	xPThick = ones(1, xSSize) + pThick(2);
	xPThick(1) = pThick(1);
	%-------------------------------------------- SHOW SHOW SHOW
	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Memory display
	for mm=1:xSSize
		newgabor = myMakeGabor(gSize, gLamda, xTheta(mm), gSigma, gPhase);
		xtex = Screen('MakeTexture', w, newgabor);	clear newgabor;
		newRect = CenterRectOnPoint(gSrcRect, gXY(xGLoc(mm),1), gXY(xGLoc(mm),2));
		Screen('DrawTexture', w, xtex, gSrcRect, newRect);	Screen('Close', xtex);
	end
	Screen('DrawDots', w, [0,0;0,0], [12,10], [255,255,255;0,0,0]', [cx, cy], 1);
	initTrial = Screen('Flip', w, tflip + dur500);

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Blank interval
% 	for mm=1:xSSize
% 		newgabor = myMakeGabor(gSize, gLamda, xTheta(mm), gSigma, gPhase);
% 		xtex = Screen('MakeTexture', w, newgabor);	clear newgabor;
% 		newRect = CenterRectOnPoint(gSrcRect, gXY(xGLoc(mm),1), gXY(xGLoc(mm),2));
% 		Screen('DrawTexture', w, xtex, gSrcRect, newRect);	Screen('Close', xtex);
% 	end

	Screen('DrawDots', w, [0,0;0,0], [12,10], [255,255,255;0,0,0]', [cx, cy], 1);
	tflip = Screen('Flip', w, initTrial + durMEM);

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Probe display
% 	for mm=1:xSSize
% 		newgabor = myMakeGabor(gSize, gLamda, xTheta(mm), gSigma, gPhase);
% 		xtex = Screen('MakeTexture', w, newgabor);	clear newgabor;
% 		newRect = CenterRectOnPoint(gSrcRect, gXY(xGLoc(mm),1), gXY(xGLoc(mm),2));
% 		Screen('DrawTexture', w, xtex, gSrcRect, newRect);	Screen('Close', xtex);
% 	end
 
	for mm=1:1%xSSize
		newRect = CenterRectOnPoint(pSrcRect, gXY(xGLoc(mm),1), gXY(xGLoc(mm),2));
		Screen('FrameOval', w, pColor, newRect, xPThick(mm), xPThick(mm));
	end
	
	% Draw a probe gabor of random orientation + fixation.
	newRect		 = CenterRectOnPoint(gSrcRect, cx, cy);
	mm = round(rand)*180;
	respTheta	= ThetaRange(randi(360)) + mm;	% random initial orientation
	newgabor	= myMakeGabor(gSize, gLamda, respTheta, gSigma, gPhase);
	xtex = Screen('MakeTexture', w, newgabor);	clear newgabor;
	Screen('DrawTexture', w, xtex, gSrcRect, newRect);	Screen('Close', xtex);

	% Show target orientation, probe's initial orientation, and probe's
	% target-centered orientation.
	RR = circshift(RespRange, [0, xTheta(1)*2]);
% 	xText = sprintf('target: %1.1f, probe: %1.1f, targ-centered: %1.1f', ...
% 				xTheta(1), mod(respTheta, 180), RR(mod(respTheta, 180)*2+1));
% 	bound = Screen('TextBounds', w, xText);
% 	Screen('DrawText', w, xText, cx - bound(3)/2, cy - bound(4)/2 + 120, 0);

	% Show a cursor somewhere between fixation and search array.
	x = cx + gEccent * 1.5 * sin(respTheta*pi/180);
	y = cy - gEccent * 1.5 * cos(respTheta*pi/180);
	SetMouse(x,y,w);
% 	[x,y,~] = GetMouse(w);
% 	r = sqrt((x-cx)^2 + (y-cy)^2);
% 	respTheta = ((cx<x)*2-1) * acos((cy-y)/r) * 180/pi;
% 	respTheta = round(respTheta*2)/2;	% round off nearest 0.5
% 	xText = sprintf('mouse: %d, %d, %1.1f', round(x), round(y), respTheta);
% 	bound = Screen('TextBounds', w, xText);
% 	Screen('DrawText', w, xText, cx - bound(3)/2, cy - bound(4)/2 + 80, 0);

	Screen('Flip', w, tflip + durISI);
	ShowCursor;

	% ---------------------------- Collect response
	FlushEvents('mouseDown');
	while 1
		if (x ~= cx || y ~= cy)	% if cursor moved
			r = sqrt((x-cx)^2 + (y-cy)^2);
			respTheta = ((cx<x)*2-1) * acos((cy-y)/r) * 180/pi;

% 			for mm=1:xSSize
% 				newgabor = myMakeGabor(gSize, gLamda, xTheta(mm), gSigma, gPhase);
% 				xtex = Screen('MakeTexture', w, newgabor);	clear newgabor;
% 				newRect = CenterRectOnPoint(gSrcRect, gXY(xGLoc(mm),1), gXY(xGLoc(mm),2));
% 				Screen('DrawTexture', w, xtex, gSrcRect, newRect);	Screen('Close', xtex);
% 			end

			for mm=1:1%xSSize
				newRect = CenterRectOnPoint(pSrcRect, gXY(xGLoc(mm),1), gXY(xGLoc(mm),2));
				Screen('FrameOval', w, pColor, newRect, xPThick(mm), xPThick(mm));
			end
	
			% Draw the probe gabor with a new orietation + fixation
			newRect	 = CenterRectOnPoint(gSrcRect, cx, cy);
			newgabor = myMakeGabor(gSize, gLamda, respTheta, gSigma, gPhase);
			xtex = Screen('MakeTexture', w, newgabor);	clear newgabor;
			Screen('DrawTexture', w, xtex, gSrcRect, newRect);	Screen('Close', xtex);

% 			xText = sprintf('mouse: %d, %d, %1.1f', round(x), round(y), respTheta);
% 			bound = Screen('TextBounds', w, xText);
% 			Screen('DrawText', w, xText, cx - bound(3)/2, cy - bound(4)/2 + 80, 0);	

			respTheta1 = round(respTheta*2)/2;	% round off nearest 0.5
% 			xText = sprintf('target: %1.1f, probe: %1.1f, targ-centered: %1.1f', ...
% 				xTheta(1), mod(respTheta1, 180), RR(mod(respTheta1, 180)*2+1));
% 			bound = Screen('TextBounds', w, xText);
% 			Screen('DrawText', w, xText, cx - bound(3)/2, cy - bound(4)/2 + 120, 0);
			Screen('Flip', w);
		end
		[x,y,buttons] = GetMouse(w);
		if buttons(1)
			break;
		end		
	end
	tflip = Screen('Flip', w);
	
	% ---------------------------- Logging response
	data.dRespTheta0(it) = respTheta;
	data.dRespTheta1(it) = mod(respTheta1, 180); % rounded off nearest 0.5
	data.dRespTheta2(it) = RR(mod(respTheta1, 180)*2+1); % target-centered
	data.dCursorX = x;
	data.dCursorY = y;

	fprintf('sn: %s, %d, con: %d, %3.1f, %d, rsp: %3.4f, %3.1f, acc: %3.1f, %3.4f, %3.4f\n', ...
		SN, trial, xSSize, xTheta(1), xGLoc(1), ...
		respTheta, mod(respTheta1, 180), RR(mod(respTheta1, 180)*2+1), x, y);
	fprintf(dataFile, '%s, %d, %d, %3.1f, %d, %3.4f, %3.1f, %3.1f, %3.4f, %3.4f\n', ...
		SN, trial, xSSize, xTheta(1), xGLoc(1), ...
		respTheta, mod(respTheta1, 180), RR(mod(respTheta1, 180)*2+1), x, y);
	fprintf(stimFile, '%s, %d, %d, %3.1f, %d, %3.4f, %3.1f, %3.1f, %3.4f, %3.4f\n', ...
		SN, trial, xSSize, xTheta(1), xGLoc(1), ...
		respTheta, mod(respTheta1, 180), RR(mod(respTheta1, 180)*2+1), x, y);

	HideCursor(w);
	Screen('FillRect', w, bcolor);
	tflip = Screen('Flip', w, tflip + dur500);
	
	% ---------------------------- Rest
	if mod(trial,20)==0 && trial<nTRIAL
		xText = sprintf('Take a rest. Press spacebar to continue.');
		bound = Screen('TextBounds', w, xText);
		Screen('DrawText', w, xText, cx - bound(3)/2, cy - bound(4)/2, 0);		
		Screen('Flip', w, tflip + dur500);
		while 1
			[keyIsDown, ~, keyCode] = KbCheck(Experimenter);
			if keyIsDown
				if keyCode(SPC)
					break;
				end
			end
		end
		Screen('Flip', w);
		WaitSecs(2);
 		
		Screen('FillRect', w, bcolor);
		Screen('DrawDots', w, [0,0;0,0], [12,10], [255,255,255;0,0,0]', [cx, cy], 1);
		Screen('Flip', w);
		WaitSecs(1);
	end	
	Screen('DrawDots', w, [0,0;0,0], [12,10], [255,255,255;0,0,0]', [cx, cy], 1);
	tflip = Screen('Flip', w);
end
%-------------------------------------------- Finishing...
fprintf(		  '************** RUN ended at %s\n', datestr(now, 0));	
fprintf(stimFile, '************** RUN ended at %s\n', datestr(now, 0));	

xText = sprintf('Done. Press spacebar to finish.');
bound = Screen('TextBounds', w, xText);
Screen('DrawText', w, xText, cx - bound(3)/2, cy - bound(4)/2, 0);		
Screen('Flip', w);
while 1
	[keyIsDown, ~, keyCode] = KbCheck(Experimenter);
	if keyIsDown
		if keyCode(SPC)
			break;
		end
	end 
end

Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);

% clean up and go home
fclose('all');
ListenChar(0);
ShowCursor;
Screen('CloseAll');
save(paramName);

return

function gabor = myMakeGabor (imSize, lamda, theta, sigma, phase, trim)
% function gabor = myMakeGabor (imSize, lambda, theta, sigma, phase, trim)
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

gabor = zeros(imSize, imSize, 4);
gabor(:,:,1:3) = repmat(Scale(grating)*255,1,1,3);
gabor(:,:,4)   = gauss*255;
%-------------------------------- EOF.
