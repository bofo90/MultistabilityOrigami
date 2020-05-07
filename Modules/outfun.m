function stop = outfun(x,optimValues,state)

stop = false;
 
switch state
   case 'init'
%        hold on
   case 'iter'
       setGlobalx(x)
   case 'done'
%        hold off
   otherwise
end


function setGlobalx(val)
global data_min
data_min.x = [data_min.x val];
%Turn the flag on only after making an iteration on the fmincon
data_min.flag = 1;