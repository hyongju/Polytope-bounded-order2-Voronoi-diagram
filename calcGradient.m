function reslt = calcGradient(active,neib1,v1,pos,p2,n)
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
                in1{i,j} = inhull(p2,v1{i,neib1{i}(j)},[],0.001);
                q1{i,j} = p2(find(in1{i,j}),:);
                %%%%%%%%%%%
                for l = 1:size(q1{i,j},1)
                    sum1{j} = sum1{j} + norm(q1{i,j}(l,:) - pos(neib1{i}(j),:))^2 * q1{i,j}(l,:);
                    sum2{j} = sum2{j} + norm(q1{i,j}(l,:) - pos(neib1{i}(j),:))^2;
                end
                if ~isempty(q1{i,j}) 
                    t_sum1{i} = t_sum1{i} + sum1{j};    
                    t_sum2{i} = t_sum2{i} + sum2{j};
                end                
            end
                %%%%%%%%%%%
        end
        reslt{i} = t_sum1{i} ./ t_sum2{i};
    end
end

