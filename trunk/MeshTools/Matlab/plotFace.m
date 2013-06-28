function plotFace(varargin)
% plotFace(FP,P,s)
% plotFace(F,FP,P,s)
%
% P:  3d point coordinates
% FP: face-to-point connectivity (aka elem)
% F:  face indeces to plot (within FP), optional
% s:  plot style
%     example: 'b' plots a blue patch
%              'b-' plots a blue line
%

%% SETUP

switch nargin
  case 3
    s  = varargin{3};
    P  = varargin{2};
    FP = varargin{1};
    F  = (1:size(FP,1))';
  case 4
    s  = varargin{4};
    P  = varargin{3};
    FP = varargin{2};
    F  = varargin{1};
end
  
[nF,nFP] = size(FP);


%% COORDS

X = zeros(nF,nFP+1);
Y = zeros(nF,nFP+1);
Z = zeros(nF,nFP+1);

lFP = [1:nFP 1];
for iFP = 1:nFP+1;
  X(:,iFP) = P(FP(F,lFP(iFP)),1);
  Y(:,iFP) = P(FP(F,lFP(iFP)),2);
  Z(:,iFP) = P(FP(F,lFP(iFP)),3);
end


%% PLOT

% plot style
if length(s) == 1
  patch(X',Y',Z',s)
else
  plot3(X',Y',Z',s,'LineWidth',1)
end


%% DONE





