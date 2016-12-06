function [sum2,indx] = calcCostMinWait(neib4,voronoi_rg2,pos,p2,coef,n1,adv,type,prob_int,voronoi_rg_1,neib1)


%% restore full voronoi_rg
for i = 1:size(neib1,2)
    for j = 1:length(neib1{i})
        if ~isempty(voronoi_rg2{i,neib1{i}(j)})
            voronoi_rg{i,neib1{i}(j)} =  voronoi_rg2{i,neib1{i}(j)};
            voronoi_rg{neib1{i}(j),i} =  voronoi_rg2{i,neib1{i}(j)};
        end
    end
end



%%
penalty = 10;
%%
% input: active, neib3, voronoi_rg, pos,p2
% output: sum2
sum2 = 0;           % going to be the total cost...
indx = 0;    
cnt = 0;
for i = 1:size(neib4,2)             % for each i 
    for j = 1:length(neib4{i})      % for each j associated with i
%         clear clcv;
        clear clcv clcv0;
        clcv = voronoi_rg{i,neib4{i}(j)};     % for each order-2 voronoi region guarted by (i,j)  
        clcv0 = voronoi_rg_1{neib4{i}(j)};
        if ~isempty(clcv) && ~isempty(clcv0)
        in1_0 = inhull(p2,clcv,[],1e-15);   % obtain all sample points inside the region
        in1_1 = inhull(p2,clcv0,[],1e-15);   % obtain all sample points inside the region
        in1 = in1_0 .* in1_1;
        cl = find(in1);
        if ~isempty(cl)                     % check if the region is empty. If not, proceed...

%             in1 = inhull(p2,clcv,[],1e-15);   % obtain all sample points inside the region
%             cl = find(in1);                   % find the index
            q1 = p2(cl,:);                    % return the positions of the sample points  
            p_int1 = prob_int(cl,:);
            if type == 1 || type == 3         % 
                if ~ismember(i,adv) && ~ismember(neib4{i}(j),adv)
                    for l = 1:size(q1,1);
        %                 sum2 = sum2 + (eta^2*norm(q1(l,:)-pos(neib4{i}(j),:))^2 * norm(q1(l,:)-pos(i,:))^2)/size(q1,1);
                        
                        sum2 = sum2 + norm(q1(l,:)-pos(i,:)) * p_int1(l,:);
%                         indx(cl(l)) = (1-f_exp(q1(l,:),pos(neib4{i}(j),:),coef)) * (1-f_exp(q1(l,:),pos(i,:),coef)) *p_int1(l,:);
        
%                         sum2 = sum2 + (eta^2*norm(q1(l,:)-pos(neib4{i}(j),:))^2 * norm(q1(l,:)-pos(i,:))^2)/n1;
%                         indx(cl(l)) = (eta^2*norm(q1(l,:)-pos(neib4{i}(j),:))^2 * norm(q1(l,:)-pos(i,:))^2)/n1;
                    end
                elseif ismember(i,adv) && ~ismember(neib4{i}(j),adv)
                    for l = 1:size(q1,1);                
                        sum2 = sum2 + norm(q1(l,:)-pos(neib4{i}(j),:)) * p_int1(l,:);
%                         indx(cl(l)) = (1-f_exp(q1(l,:),pos(neib4{i}(j),:),coef))*p_int1(l,:);
                    end
                elseif ~ismember(i,adv) && ismember(neib4{i}(j),adv)
                    for l = 1:size(q1,1);                
                        sum2 = sum2 + norm(q1(l,:)-pos(i,:)) * p_int1(l,:);
%                         indx(cl(l)) = (1-f_exp(q1(l,:),pos(i,:),coef))*p_int1(l,:);
                    end
                else
                    for l = 1:size(q1,1);                
                        sum2 = sum2 + penalty * p_int1(l,:);
%                         indx(cl(l)) = p_int1(l,:);
                    end
                end
            else
                for l = 1:size(q1,1);           % for all sample points...
                    sum2 = sum2 + norm(q1(l,:)-pos(i,:)) * p_int1(l,:);
%                     indx(cl(l)) = (1-f_exp(q1(l,:),pos(neib4{i}(j),:),coef)) * (1-f_exp(q1(l,:),pos(i,:),coef))*p_int1(l,:);
                    cnt = cnt + 1;
                end
            end
        end
        end
    end
end
% cnt