function meshout = ReadSU2(meshname)
% meshout = ReadSU2(meshname)
% reads meshname.su2 into struct meshout 
% 
% meshout = 
%       dim: mesh dimension
%     npoin: number of points
%      poin: point array [x y (z) ip]
%     nelem: number of elements
%      elem: element array [ type p1 p2 ... pn ielem ]
%     nmark: number of markers
%      mark: marker struct array
% meshout.mark(imark) = 
%       tag: marker tag
%     nelem: number of marker elements
%      elem: element array [ type p1 p2 ... pn ]
%
% note: element array is right-padded with zeros to
% accomodate multiple element types.  use getElemTypeInfo()
% to get how many points are in an element.
%

%% Open File

fprintf('Read %s \n',meshname);
fid = fopen(meshname,'r');

% make sure file opened
if fid == -1
  error('file open error')
end

%% Setup

% initialize data
meshout.dim   = 0;
meshout.npoin = 0;
meshout.poin  = [0];
meshout.nelem = 0;
meshout.elem  = [0];
meshout.nmark = 0;
meshout.mark(1).tag  = '';
meshout.mark(1).nelem = 0;
meshout.mark(1).elem  = [0];

% marker index
imarkr = 1;


%% Read

% scan file until end of file
while ~feof(fid)
  
  % look for mesh keywords
  % read one line
  C = textscan(fid,'%s %s',1,'CommentStyle','%');
  
  % check for end of file
  if feof(fid)
    break
  end
    
  % tokenize line
  key = C{1}{1};
  val = C{2}{1};
  
  % check if val is a double or a string
  temp = str2double(val);
  if ~isnan(temp)
    val = temp;
  end
  
  % handle the keyword
  switch key
    
    % dimension number
    case 'NDIME='
      % save to data struct
      meshout.dim=val;
    
    % elements
    case 'NELEM='
      % save to data struct
      meshout.nelem = val;
      
      % read next lines for element data
      % builds a potentially oversized matrix
      E = textscan(fid, '%u64 %u64 %u64 %u64 %u64 %u64 %u64 %u64 %u64 %u64', ...
                        meshout.nelem, ...
                        'CollectOutput',true,'CommentStyle','%' );
      E = E{1};
      
      % save to data struct
      meshout.elem = E;
           
    % points
    case 'NPOIN='
      % save to data struct
      meshout.npoin = val;
      
      % read next lines for point data
      E = textscan(fid, '%f64 %f64 %f64 %f64 %f64', ...
                        meshout.npoin   , ...
                        'CollectOutput',true,'CommentStyle','%' );
      E = E{1};
      
      % save to data struct
      meshout.poin = E(:,1:(meshout.dim+1));
      
    % number of markers
    case 'NMARK='
      % save to data struct
      meshout.nmark = val;
      
    % a marker
    case 'MARKER_TAG='
      % save name to data.mark struct
      meshout.mark(imarkr).tag=val;
      
      % read next line to get marker nelems
      C = textscan(fid,'%s %s',1,'CommentStyle','%');
      % tokenize line
      key = C{1}{1};
      val = C{2}{1};
      
      if ~strcmp(key,'MARKER_ELEMS=')
        error('marker specification error')
      end
        
      % save number of elements to data.mark struct
      meshout.mark(imarkr).nelem = str2double(val);
      
      % read next lines to get marker elements
      E = textscan(fid, '%u64 %u64 %u64 %u64 %u64 %u64 %u64 %u64 %u64 %u64', ...
                        meshout.mark(imarkr).nelem, ...
                        'CollectOutput',true,'CommentStyle','%' );
      E = E{1};
      
      % save elements to data.mark struct
      meshout.mark(imarkr).elem = E;
      
      % increment marker index
      imarkr = imarkr+1;
            
  end % switch key
  
  
  
end % while read file line


%% Done

% close
fclose(fid);




