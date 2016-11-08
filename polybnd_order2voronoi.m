function [voronoi_rg,vornb,vornb2] = polybnd_order2voronoi(pos,bnd_pnts)

%% =======================================================
% Order-2 Voronoi Diagram with set of points in 2D/3D polygon
% ========================================================
% version 1.01
% by Hyongju Park
%---------------------------------------------------------
% inputs: bnd_pnts      boundary points                m x 2
%         pos           points inside the boundary     n x 2
%---------------------------------------------------------
% outputs: voronoi_rg   order-2 Voronoi regions        ? x n  
%          vornb        Voronoi neighbors              1 x n
%          vornb2       Voronoi neighbors              1 x ? (repeating
%          neighbors removed)
% =========================================================================
% This functions works for d = 2, 3
% -------------------------------------------------------------------------
% This function requires:
%       vert2lcon.m (Matt Jacobson / Michael Keder)
%       pbisec.m (by me)
%       con2vert.m (Michael Keder)
%       inhull.m (John D'Errico)
% -------------------------------------------------------------------------
% Written by Hyongju Park, hyongju@gmail.com / park334@illinois.edu
% Change logs:
% 11 Aug 2015: skip error messages (version 1.01) 
% 5  May 2015: initial release (version 1.0)
% =========================================================================
% Known issues:
% Input points must satisfy assumptions such as non co-circularity and
% general position assumptions
% -------------------------------------------------------------------------
%%
[vornb,~,~] = polybnd_voronoi(pos,bnd_pnts);      % 
% obtain set of voronoi neighbors/vertices
[Abnd,bbnd] = vert2lcon(bnd_pnts);              % obtain series of linear inequalities that defined Voronoi regions
%% create list
for i = 1:size(pos,1)
    k = 0;
    for j = 1:size(vornb{i},2)
        if vornb{i}(1,j) > i
            k = k + 1;
            vornb2{i}(1,k) = vornb{i}(1,j);
        end
    end
end
for m1 =1:size(vornb2,2)
    for j = 1:size(vornb2{m1},2)
        c1 = m1;
        c2 = vornb2{m1}(1,j);
        clear Aag1 bag1 Aag2 bag2 Aagmt bagmt pos1 pos2
        % given (c1,c2) where c1< c2
        % remove c2, compute voronoi vertices of c1
        k = 0;
        for i = 1:size(pos,1)
            if i ~= c2
                k = k + 1;
                pos1(k,:) = pos(i,:);
            end
        end
        [~,~,Aag1,bag1] = polybnd_voronoi(pos1,bnd_pnts);
        % remove c1, compute voronoi vertices of c2
        k = 0;
        for i = 1:size(pos,1)
            if i ~= c1
                k = k +1;
                pos2(k,:) = pos(i,:);
            end
        end
        [~,~,Aag2,bag2] = polybnd_voronoi(pos2,bnd_pnts);
        Aagmt = [Aag1{c1};Aag2{c2-1};Abnd];
        bagmt = [bag1{c1};bag2{c2-1};bbnd];
        Vl{c1,c2}= MY_con2vert(Aagmt,bagmt);
        if ~isempty(Vl{c1,c2})
            % remove Voronoi regions that are not in the interior of the polytope
            if inhull(Vl{c1,c2},bnd_pnts,[],0.001)
%             if inhull(Vl{c1,c2},vorvx{c1},[],0.001) & inhull(Vl{c1,c2},vorvx{c2},[],1.e-13*mean(abs(bnd_pnts(:)))) 
                IDl{c1,c2} = convhull(Vl{c1,c2});
                voronoi_rg{c1,c2} = Vl{c1,c2}(IDl{c1,c2},:);
            else
                voronoi_rg{c1,c2} = [];
            end
        else
            voronoi_rg{c1,c2} = [];
        end
    end
end
