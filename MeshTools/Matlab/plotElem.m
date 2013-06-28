function plotElem(poin,elem,plotstyle)
% plotElem(poin,elem,plotstyle)
% plots points of the element array
%
% poin: 3d point locations
% elem: element connectivity
% plotstyle: element plot style
%            i.e. 'b'  plots a blue patch
%                 'b-' plots a blue line
%

%% SETUP

% element type plot list
etype_list = [3 5 9]; % line, triangle, rectangle


%% ELEMENTS

for elem_type = etype_list

  % find element rows of this type
  i_plotelem = elem(:,1) == elem_type;

  if any(i_plotelem)

    % pull how many points
    elem_leng  = getElemTypeInfo(elem_type);
    elem_range = 2:(1+elem_leng);
    
    % pull plotable element array
    this_elem  = elem(i_plotelem,elem_range);
    this_elem  = this_elem + 1;

    % plot elements
    plotFace(this_elem,poin,plotstyle)

  end

end


%% DONE


