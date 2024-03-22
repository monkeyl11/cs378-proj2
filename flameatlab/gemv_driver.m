pathstoadd = genpath('../libflameatlab');
addpath(pathstoadd);

m = 4; n = 6;

y = randi ([-3,3], [m 1]);
A = randi( [-3,3], [m n]);
x = randi( [-3,3], [n 1]);

yref = y + A*x; 

myy = gemv_unb(y, A, x);
if (isequal(yref, myy))
    disp('All is well');
else
    disp('Trouble in Paradise');
end
