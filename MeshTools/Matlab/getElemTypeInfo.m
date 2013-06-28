function [elem_leng,elem_type] = getElemTypeInfo(elem_type)
% get vtk element type data
% input - 
% elem_type = name or vtk type number
% output - 
% elem_leng = number of nodes in element
% elem_type = vtk type number


if isnumeric(elem_type)
  switch elem_type
    case 3
      elem_leng = 2;
    case 5
      elem_leng = 3;
    case 9
      elem_leng = 4;
    case 10
      elem_leng = 4;
    case 12
      elem_leng = 8;
    case 13
      elem_leng = 6;
    case 14
      elem_leng = 5;   
  end 
  
else
  switch elem_type
    case 'line'
      elem_type = 3;
      elem_leng = 2;
    case 'triangle'
      elem_type = 5; 
      elem_leng = 3;
    case 'rectangle'
      elem_type = 9;
      elem_leng = 4;
    case 'tetrahedral'
      elem_type = 10;
      elem_leng = 4;
    case 'hexahedral'
      elem_type = 12;
      elem_leng = 8;
    case 'wedge'
      elem_type = 13;
      elem_leng = 6;
    case 'pyramid'
      elem_type = 14;
      elem_leng = 5;   
  end  
end