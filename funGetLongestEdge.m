%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-07-05(yyyy-mm-dd)
% 深度搜索给定节点序列
% h: header
% t: branch
% e: 枝信息
%--------------------------------------------------------------------------
function [mTree] = funGetLongestEdge(h, t, e, eSource, eLoad)
mTree = [];
if e(eSource).node1 ~= 1
    nst = e(eSource).node1;
else
    nst = e(eSource).node2;
end
if e(eLoad).node1 ~= 1
    nsp = e(eLoad).node1;
else
    nsp = e(eLoad).node2;
end
% 0:V,1:I,2:R,3:L,4:C
np = 1;% 路径数量
nv = zeros(size(e));% 边判重，防止进入循环
[nf, np, nv] = funDFS(h, t, e, nst, 0, nsp, [], np, nv);
[m, n] = size(nf);
mm = zeros(1, m-1);
for jj=1:m-1
    mTree(jj, 1) = nst;
    for ii=2:n
        mTree(jj, ii) = nf(jj, mTree(jj, ii-1));
        if mTree(jj, ii) == nsp
            mm(jj) = ii;
            break;
        end
    end
end
[m_max, i_max] = max(mm);
mTree = mTree(i_max, :);
% if e==0
% end
function [nf, np, nv] = funDFS(h, t, e, node_c, node_pre, nsp, nf, np, nv)
next_node = h{node_c};
m = length(next_node);
for ii=1:m
    cur_edge = next_node(ii);
    if nv(cur_edge) == 0
        nv(cur_edge) = 1;
        cur_node = e(cur_edge).node2;
        if nsp == cur_node
            nf(np, e(cur_edge).node1) = e(cur_edge).node2;
            np = np + 1;
            [a, b] = size(nf);
            for jj=1:b
                nf(np, jj) = nf(np-1, jj);
            end
    %         fprintf('%d,%d,%d,%s\n', e(cur_edge).node1, e(cur_edge).node2, e(cur_edge).iType, e(cur_edge).strName{1});
    %         break;
        elseif cur_node ~= node_pre && cur_node ~= 1
            nf(np, e(cur_edge).node1) = e(cur_edge).node2;
    %         fprintf('%d,%d,%d,%s\n', e(cur_edge).node1, e(cur_edge).node2, e(cur_edge).iType, e(cur_edge).strName{1});
            [nf, np, nv] = funDFS(h, t, e, cur_node, node_c, nsp, nf, np, nv);
        end
    end
end

