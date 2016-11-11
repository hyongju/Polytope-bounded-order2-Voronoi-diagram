% THIS IS A VERSION WITH NO FAILURES....
clear all;close all;clc
%% generate random samples
n = 20;
% % m = 15;
d = 2;
eta = 1/2;
p1_0 = haltonset(d,'Skip',1e3,'Leap',1e2);
p1_1 = scramble(p1_0,'RR2');
pos = 1/4 *(net(p1_1,n) -0.5 * ones(n,d)) + 0.5 * ones(n,d);
% pos = net(p1_1,n);
figure,plot(pos(:,1),pos(:,2),'x');axis([0 1 0 1]);
% % pos = rand(n,d);
% % pos = [0.15 1-0.15;1-0.15 1-0.15;0.5 0.5;0.15 0.15;1-0.15 0.15];
vmax = 1/eps;       % velocity constraint is RELAXED....
stage = 30;
coef = 1000;

% bnd_pnts = [0 1;1 1;1 0;0 0];
m = 50;             % # of vertices for the outter polytope
p2_0 = net(p1_1,m);
bnd_idx = convhull(p2_0);
bnd_pnts = p2_0(bnd_idx,:);

%% new codes: normal distribution multi-variate (uni-modal)
% mean = [0.75 0.75], covariance = 0.1*eye(2)



mu = [0.75 0.75]; 
SIGMA = 0.05*eye(2);
X_test = mvnrnd(mu,SIGMA,10000); 
p_test = mvnpdf(X_test,mu,SIGMA); 
k = 0;

[A_bnd,b_bnd] = vert2lcon(bnd_pnts);     
in_y = inhull(X_test,bnd_pnts,[],1e-15);

X_int = X_test(in_y,:);
p_int_unnorm = p_test(in_y,:);

% for i = 1:size(X_test,1)
%     if X_test(i,1) >=0 && X_test(i,2) <= 1 && X_test(i,1) <=1 && X_test(i,2) >= 0
%         k = k+1;
%         X_int(k,:) = X_test(i,:);
%         p_int_unnorm(k) = p_test(i);
%     end
% end
p_int = p_int_unnorm / sum(p_int_unnorm);
h_00 = figure('position',[100 100 600 600],'Color',[1 1 1]);
plot3(X_int(:,1),X_int(:,2),p_int(:),'Marker','.','MarkerSize',1,'LineStyle','none');
hold on;
bdp = convhull(bnd_pnts);
plot(bnd_pnts(bdp,1),bnd_pnts(bdp,2),'b-');
hold on;
plot(mu(1),mu(2),'Marker','o','MarkerSize',10,'Color','r','LineWidth',2);
axis('square')
axis([0 1 0 1 0 0.0003])
% axis('off')
set(gca,'xtick',[0 1]);
set(gca,'ytick',[0 1]);

% grid on;
% axis('equal');
xlabel('X');ylabel('Y');zlabel('Target distribution');
view(0,90)
%%
% n1 = 10000;         % number of quasi-random samples to compute the cost (Monte-Carlo)
% n1 should be at least > 1000 for accuracy...                
% p2  = net(p1_1,n1); % generate samples points for Monte-Carlo
p_sav{1} = pos;
% adv = [1 2 3 4 5];
adv = [];           % index set of faulty nodes
type = 1;
p_sav{1} = pos;
% n1 = size(p2,1);
p2 = X_int;
n1 = size(X_int,1);
%% call function
for t = 1:stage
    t
    [voronoi_rg{t},neib1{t},neib2{t}] = polybnd_order2voronoi(pos,bnd_pnts);
    
    %   order2 Voronoi regions - {cell} data structure
    %   neib1: neighbors indice can be repeated
    %   neib2: neighbors indice cannot be repeated

    active{t}  = chooseSset(t,size(pos,1),neib1{t});
    l_min{t} = calcGradientExp(active{t},neib1{t},voronoi_rg{t},pos,p2,size(pos,1),coef,p_int);   % l_min: local minimizer
    [cst(t),~] = calcCostExp(neib2{t},voronoi_rg{t},pos,p2,coef,n1,adv,type,p_int);
    idx{t} = find(active{t});
    if type == 2 || type == 3
        for y = 1:length(idx{t})
            if ~ismember(idx{t}(y),adv)
                if norm(l_min{t}{idx{t}(y)}- pos(idx{t}(y),:)) <= vmax
                    pos(idx{t}(y),:) = l_min{t}{idx{t}(y)};
                else
                    pos(idx{t}(y),:) = pos(idx{t}(y),:) + vmax * (l_min{t}{idx{t}(y)}- pos(idx{t}(y),:))/norm(l_min{t}{idx{t}(y)}- pos(idx{t}(y),:));
                end
            end
        end
    else
        for y = 1:length(idx{t})
            if norm(l_min{t}{idx{t}(y)}- pos(idx{t}(y),:)) <= vmax
                pos(idx{t}(y),:) = l_min{t}{idx{t}(y)};
            else
                pos(idx{t}(y),:) =  pos(idx{t}(y),:) + vmax* (l_min{t}{idx{t}(y)}- pos(idx{t}(y),:))/norm(l_min{t}{idx{t}(y)}- pos(idx{t}(y),:));
            end
       end
    end
    p_sav{t+1} = pos;
end
cst
figure,plot(1:stage,cst,'-s');xlabel('stage');ylabel('cost');
h0 = figure('position',[0 0 700 700],'Color',[1 1 1]);
k = 0;
t = 16;
for i = 1:size(voronoi_rg{t},1)*size(voronoi_rg{t},2)
    col(i,:)= rand(1,3);
%     col(i,:) = [i/(size(voronoi_rg{t},1)*size(voronoi_rg{t},2)) 1 1];
end
col = distinguishable_colors(size(voronoi_rg{t},1)*size(voronoi_rg{t},2));
for i = 1:size(voronoi_rg{t},1)
    for j = 1:size(voronoi_rg{t},2)
        if ~isempty(voronoi_rg{t}{i,j})
            k = k+1;
            if ismember(i,adv) && ismember(j,adv)
                 patch(voronoi_rg{t}{i,j}(:,1),voronoi_rg{t}{i,j}(:,2),[0.9 0.9 0.9]);
                 hold on;
            end
            plot(voronoi_rg{t}{i,j}(:,1),voronoi_rg{t}{i,j}(:,2),'-','Color','b');
            hold on;
             patch(voronoi_rg{t}{i,j}(:,1),voronoi_rg{t}{i,j}(:,2),col(k,:));
        end
    end
end
bdp = convhull(bnd_pnts);
plot(bnd_pnts(bdp,1),bnd_pnts(bdp,2),'b-');
hold on;
plot(p_sav{t}(:,1),p_sav{t}(:,2),'Marker','o','MarkerSize',12,'MarkerFaceColor','r','Color','b','LineStyle','none');hold on;
plot(p_sav{t}(adv,1),p_sav{t}(adv,2),'Marker','o','MarkerSize',24,'MarkerFaceColor','r','Color','b','LineStyle','none'); hold on;
axis('equal')
axis([0 1 0 1]);
axis('off');
set(gca,'xtick',[]);
set(gca,'ytick',[]);  

p_ani{1} = p_sav{1};
p_ani{2} = p_sav{stage+1};

res = 60;
k_p = 5;
M = animationFunc(p_ani',res,n,0,bnd_pnts,k_p,adv);
h2 = figure('position',[100 100 600 600],'Color',[1 1 1]);
axis('square')
axis([0 1 0 1 0 1])
axis('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
movie(h2,M,1,30);
% movie2avi(M, 'vid_exp_50_0_75_0_75.avi', 'compression', 'None','quality',100,'fps',15);
% % 
% % 
% % 
% h0 = figure('position',[0 0 800 700],'Color',[1 1 1]);
% x = p2(:,1);
% y = p2(:,2);
% z = indx{stage};
%  
% h = scatter3(x,y,z,'Marker','o','MarkerFaceColor',[0 .75 .75], 'MarkerEdgeColor','k');
% hChildren = get(h, 'Children');
% set(hChildren, 'MarkerSize', 0.5)
%  
% % axis('equal')
% % axis([0 1 0 1 0 1]);
% set(gca,'xtick',[0 1]);
% set(gca,'ytick',[0 1]);
% set(gca,'ztick',[0 1]);
% view(3)
 
