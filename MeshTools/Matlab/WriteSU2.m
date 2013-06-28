function WriteSU2(meshout,meshname)
% WriteSU2(meshout,meshname)
% writes meshname.su2 from struct meshout 
% for some reason it will only write single precision
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

fprintf('Write %s \n',meshname);
fid = fopen(meshname,'w');


%% Dimension

fprintf(fid,'%% \n%% Problem dimension \n%% \n');
fprintf(fid,'NDIME= %i\n%',meshout.dim);


%% Elements

fprintf(fid,'%% \n%% Inner element connectivity \n%% \n');
fprintf(fid,'NELEM= %i\n%',meshout.nelem);
% fprintf(fid,'%-2i %5i %5i %5i %5i \n',data.elem');
for i = 1:meshout.nelem
  this_elem = meshout.elem(i,:);
  elem_type = this_elem(1);
  elem_leng = getElemTypeInfo(elem_type);
  elem_range = 1:(2+elem_leng);
  printElement( fid, this_elem(elem_range) );
end


%% Nodes

fprintf(fid,'%% \n%% Node coordinates  \n%% \n');
fprintf(fid,'NPOIN= %i\n%',meshout.npoin);
if meshout.dim == 2
  fprintf(fid,'%#28.16g %#28.18g  %i \n',meshout.poin');
else
  fprintf(fid,'%#28.18g %#28.18g %#28.18g  %i \n',meshout.poin');
end


%% Markers

fprintf(fid,'%%\n%% Boundary elements  \n%% \n');
fprintf(fid,'NMARK= %i\n%',meshout.nmark);

nM = meshout.nmark;
for i = 1:nM
  fprintf(fid,'MARKER_TAG= %s\n',meshout.mark(i).tag);
  fprintf(fid,'MARKER_ELEMS= %i\n',meshout.mark(i).nelem);
  for j = 1:meshout.mark(i).nelem
    this_elem = meshout.mark(i).elem(j,:);
    elem_type = this_elem(1);
    elem_leng = getElemTypeInfo(elem_type);
    elem_range = 1:(1+elem_leng);
    printElement( fid, this_elem(elem_range) );
  end
end

fclose all;


%% Helper Function
function printElement(fid,elem)
fprintf(fid,'%-2i ',elem(1));
fprintf(fid,'%5i ',elem(2:end));
fprintf(fid,'\n');





