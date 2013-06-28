function [dcm] = angle2dcm(r,seq)
% dcm = angle2dcm(r1,r2,r3,seq)
% builds euler angle rotation matrix
% 
% r = [r1 r2 r3]
% seq = 'ZYX' (default)
%       'ZXZ'
% dcm = direction cosine matrix

if nargin == 1
  seq = 'ZYX';
end

switch seq
  case 'ZYX'
    dcm = Tx(r(3))*Ty(r(2))*Tz(r(1));
  case 'ZXZ'
    dcm = Tz(r(3))*Tx(r(2))*Tz(r(1)); 
end

function A = Tx(a)
A = [1 0 0;0 cosd(a) sind(a);0 -sind(a) cosd(a)];

function A = Ty(a)
A = [cosd(a) 0 -sind(a);0 1 0;sind(a) 0 cosd(a)];

function A = Tz(a)
A = [cosd(a) sind(a) 0;-sind(a) cosd(a) 0;0 0 1];