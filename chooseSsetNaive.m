%% choose S
 
function active  = chooseSsetNaive(t,n)
% input t n neib
active = zeros(1,n);
if mod(t,n) == 0
    active(n) = 1;
else
    active(mod(t,n)) = 1;
end
