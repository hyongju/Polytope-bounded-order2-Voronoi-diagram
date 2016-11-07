%% choose S
 
function active  = chooseSset(t,n,neib1)
% input t n neib
active = zeros(1,n);
for i = 1:n
    if mod(t,n) == i || (mod(t,n)==0 && i == n)
        active(i) = 1;
    end
end
for i = 1:n
    if active(i) ~= 1 
        if active(1,neib1{i}) == 0
%             if i < neib1{i}
            active(i) = 1;
%             end
        end
    end
end
