function reslt = calcGradientMinWait_hotfix(active,neib1,v1,pos,p2,n,coef,prob_int,voronoi_rg_1,neib2,adv,type,bnd_pnts)
%% restore to complete the Voronoi partition
for i  = 1: size(v1,1)
    for j = 1: size(v1,2)
        if ~isempty(v1{i,j}) && (j > i)
            v1{j,i} = v1{i,j};
        end
    end
end
%% compute the critical point------
% for each i
for i = 1:n
    if active(i) == 1
        t_sum1{i} = [0 0];
        t_sum2{i} = 0;
        for j = 1:size(neib1{i},2)
            k = 0;
            sum1{j} = [0 0];
            sum2{j} = 0;
            if ~isempty(v1{i,neib1{i}(j)})
%                 in1{i,j} = inhull(p2,voronoi_rg_1{neib1{i}(j)},v1{i,neib1{i}(j)},[],1e-15);
                in1_0 = inhull(p2,voronoi_rg_1{neib1{i}(j)},[],1e-15);
                in1_1 = inhull(p2,v1{i,neib1{i}(j)},[],1e-15);
                in1{i,j} = in1_0.* in1_1;
                q1{i,j} = p2(find(in1{i,j}),:);
                p_int1{i,j} = prob_int(find(in1{i,j}),:);
                %%%%%%%%%%%
                for l = 1:size(q1{i,j},1)
                    % function f_val = f_exp(q,p,coef)
%                     sum1{j} = sum1{j} + q1{i,j}(l,:)/norm(q1{i,j}(l,:)-pos(i,:))*p_int1{i,j}(l,:);
                    
                    sum2{j} = sum2{j} - 1/norm(q1{i,j}(l,:)-pos(i,:))*(q1{i,j}(l,:)-pos(i,:))*p_int1{i,j}(l,:);
    
                    
%                     sum1{j} = sum1{j} + norm(q1{i,j}(l,:) - pos(neib1{i}(j),:))^2 * q1{i,j}(l,:);
%                     sum2{j} = sum2{j} + norm(q1{i,j}(l,:) - pos(neib1{i}(j),:))^2;
                end
                if ~isempty(q1{i,j}) 
%                     t_sum1{i} = t_sum1{i} + sum1{j};    
                    t_sum2{i} = t_sum2{i} + sum2{j};
                end                
            end
                %%%%%%%%%%%
        end
%         reslt{i} = t_sum1{i} ./ t_sum2{i};
        reslt{i} = pos(i,:)- t_sum2{i};
    end
end
alph = 1;
for i =1:n
    if active(i) == 1
        pos_tmp(i,:) = pos(i,:) - alph*t_sum2{i};
    else
        pos_tmp(i,:) = pos(i,:);
    end
end
cst1 = calcCostMinWait(neib1,v1,pos,p2,coef,[],adv,type,prob_int,voronoi_rg_1,neib2);
[~,v_1_r_1,~,~] = polybnd_voronoi(pos_tmp,bnd_pnts);
[v_1_r,neib1,neib2] = polybnd_order2voronoi(pos_tmp,bnd_pnts);
cst2 = calcCostMinWait(neib1,v_1_r,pos_tmp,p2,coef,[],adv,type,prob_int,v_1_r_1,neib2);

% cst1 = 0;
% cst2 = 1;
% max_Iter = 30;
k = 0;
while(cst2 > cst1)
    k = k+1;
    alph = alph * 0.9;
    for i =1:n
        if active(i) == 1
            pos_tmp(i,:) = pos(i,:) - alph*t_sum2{i};
        else
            pos_tmp(i,:) = pos(i,:);
        end
    end
    [~,v_1_r_1,~,~] = polybnd_voronoi(pos_tmp,bnd_pnts);
    [v_1_r,neib1,neib2] = polybnd_order2voronoi(pos_tmp,bnd_pnts);
    cst2 = calcCostMinWait(neib1,v_1_r,pos_tmp,p2,coef,[],adv,type,prob_int,v_1_r_1,neib2);    
end
for i = 1:n
    reslt{i} = pos_tmp(i,:);
end
cst1
cst2
k

% [cst,~] = calcCostMinWait(neib1,v1,pos,p2,coef,n1,adv,type,prob_int,voronoi_rg_1,neib2);




