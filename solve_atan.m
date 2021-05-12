% Solve angle  z(k)=atan(x(1)/x(2))
% Input: x: [x(1) x(2)]
% Output: angle  z according to x
% 
function z=solve_atan(x)
if x(2)>0
    z=atan(x(1)/x(2));
elseif x(2)<0
    if x(1)>0 
       z=pi+atan(x(1)/x(2));
   else  % when x(1) change from + to - ,z has a big change
       z=-pi-atan(x(1)/x(2));
   end
elseif x(1)>0  
    z=.5*pi;
elseif x(1)<0
    z=-.5*pi;
else
    disp('singularity point');
    return;
end