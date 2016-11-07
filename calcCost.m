function [sum2,indx] = calcCost(neib4,voronoi_rg,pos,p2,eta,n1,adv,type)
% input: active, neib3, voronoi_rg, pos,p2
% output: sum2
sum2 = 0;           % going to be the total cost...
indx = 0;    
cnt = 0;
for i = 1:size(neib4,2)             % for each i 
    for j = 1:length(neib4{i})      % for each j associated with i
%         clear clcv;
        clcv = voronoi_rg{i,neib4{i}(j)};     % for each order-2 voronoi region guarted by (i,j)  
        if ~isempty(clcv)                     % check if the region is empty. If not, proceed...
            in1 = inhull(p2,clcv,[],1e-15);   % obtain all sample points inside the region
            cl = find(in1);                   % find the index
            q1 = p2(cl,:);                    % return the positions of the sample points  
            if type == 1 || type == 3         % 
                if ~ismember(i,adv) && ~ismember(neib4{i}(j),adv)
                    for l = 1:size(q1,1);
        %                 sum2 = sum2 + (eta^2*norm(q1(l,:)-pos(neib4{i}(j),:))^2 * norm(q1(l,:)-pos(i,:))^2)/size(q1,1);
                        sum2 = sum2 + (eta^2*norm(q1(l,:)-pos(neib4{i}(j),:))^2 * norm(q1(l,:)-pos(i,:))^2)/n1;
                        indx(cl(l)) = (eta^2*norm(q1(l,:)-pos(neib4{i}(j),:))^2 * norm(q1(l,:)-pos(i,:))^2)/n1;
                    end
                elseif ismember(i,adv) && ~ismember(neib4{i}(j),adv)
                    for l = 1:size(q1,1);                
                        sum2 = sum2 + (eta*norm(q1(l,:)-pos(neib4{i}(j),:))^2)/n1;
                        indx(cl(l)) = (eta*norm(q1(l,:)-pos(neib4{i}(j),:))^2)/n1;
                    end
                elseif ~ismember(i,adv) && ismember(neib4{i}(j),adv)
                    for l = 1:size(q1,1);                
                        sum2 = sum2 + (eta*norm(q1(l,:)-pos(i,:))^2)/n1;
                        indx(cl(l)) = (eta*norm(q1(l,:)-pos(i,:))^2)/n1;
                    end
                else
                    for l = 1:size(q1,1);                
                        sum2 = sum2 + 1/n1;
                        indx(cl(l)) = 1/n1;
                    end
                end
            else
                for l = 1:size(q1,1);           % for all sample points...
                    sum2 = sum2 + (eta^2*norm(q1(l,:)-pos(neib4{i}(j),:))^2 * norm(q1(l,:)-pos(i,:))^2)/n1;
                    indx(cl(l)) = (eta^2*norm(q1(l,:)-pos(neib4{i}(j),:))^2 * norm(q1(l,:)-pos(i,:))^2)/n1;
                    cnt = cnt + 1;
                end
            end
        end
    end
end
% cnt
