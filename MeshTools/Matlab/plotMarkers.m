function plotMarkers(varargin)
% plotMarkers(meshdata,plotstyle)
% plotMarkers(meshdata,marklist,plotstyle)
%
% plots markers from the meshdata
%
% meshdata: mesh data struct
% marklist: cell array of marker names, optional
% plotstyle: element plot style
%            i.e.: 'b'  plots a blue patch
%                  'b-' plots a blue line
%

%% SETUP

switch nargin
  case 2
    meshdata  = varargin{1};
    marklist  = {meshdata.mark(:).tag};
    plotstyle = varargin{2};
  case 3
    meshdata  = varargin{1};
    marklist  = varargin{2};
    plotstyle = varargin{3};
end

% number of markers
nmark = meshdata.nmark;

% dimension
dim = meshdata.dim;

% point array
poin = meshdata.poin(:,1:dim);
if dim == 2
  poin(:,3) = 0;
end


%% MARKERS

% plot nodes for each marker
for imark = 1:nmark
  
  % pull marker tag
  this_tag = meshdata.mark(imark).tag;

  % check if this tag exists
  if any(strcmp(this_tag,marklist))

    % pull elements and plot
    this_elem = meshdata.mark(imark).elem;
    plotElem(poin,this_elem,plotstyle);

  end
  
end


%% DONE




