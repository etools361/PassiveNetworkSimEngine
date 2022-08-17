%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-07-05(yyyy-mm-dd)
% 由网表生成原理图
% 如何由网表生成原理图？
% 前提，需要给定起点(源)和终点(负载),其他电路均为待定电路网络，也即将整个电路理解为一个二端口网络。
%       0节点为接地点，不用考虑。
% 使用图论方法，建立从起点到终点的最长树，再在这棵树上寻每个枝节，可能存在主树上长出一个连通的小树。
% 数据结构：使用
%--------------------------------------------------------------------------
function [img] = funGenSchematic(iType, Value, cellNode1, CellNode2, cellName, DispEn)
% netlist standard
plot(0,0);
[node1, node2] = funSimNetlistRenode(cellNode1, CellNode2);
iSource = [find(iType==0), find(iType==1)];
nSource = length(iSource);
for ii=1:nSource
    if Value(iSource(ii)) ~= 0
        iSource = iSource(ii);
        break;
    end
end
if isempty(iSource) || length(iSource) > 1
    iSource = 1;
end
[a, b] = ismember('RL', cellName);
if a
    iLoad = b;
else
    iLoad   = find(iType==2);% is R and R to GND
    nLoad   = length(iLoad);
    for ii=1:nLoad
        if node1(iLoad(ii)) == 0 || node2(iLoad(ii)) == 0
            iLoad = iLoad(ii);
            break;
        end
    end
    if isempty(iLoad)
        iLoad = length(iType);
    end
end
% find the main tree
if node1(iSource) ~= 0
    nodeS = node1(iSource);
else
    nodeS = node2(iSource);
end
if node1(iLoad) ~= 0
    nodeL = node1(iLoad);
else
    nodeL = node2(iLoad);
end
node1t  = node1;
node2t  = node2;
nodeCh  = nodeS;
n = length(iType);
nodeLen  = zeros(1, n);
histNode = nodeS;
nodeChT  = nodeS;
index = 0;
% 递归找出从Source到Load的树
[nodeLen, nodeNew, nodeChT, index, histNode] = funFindNodeNext(nodeChT, node1, node2, nodeLen, index, histNode, nodeS);

% 先绘制出主树干的器件
nTree = length(histNode)-1;
isDrawNode = zeros(1, n);
hist_mNode = 0;
nBranchStartPoint    = zeros(1, nTree+1);% 枝节扩展位置(起始位置)
mExtNodeNum     = zeros(1, nTree+1);% 枝节扩展数量（每个节点）
NodeLength = 0.4;% 枝节扩展长度
for ii=1:nTree
    % 两个相邻节点mT1和mT2
    mT1 = histNode(ii);
    mT2 = histNode(ii+1);
    % find devie between node mT1 and mT2,找到mT1和mT2之间的器件
    [iNode, mNode, mNodeBack] = funFindNode(mT1, mT2, node1, node2);
    % 枝节扩展距离mNode0(每个节点)
    if ii == 1
        mExtNodeNum(ii) = 0;
    else
        mTemp = mNode-hist_mNode-length(iNode);
        if mTemp < 2
            mTemp = 0;
        end
        mExtNodeNum(ii) = mExtNodeNum(ii-1)+mTemp;
    end
    nBranchStartPoint(ii) = ii - 1 + mExtNodeNum(ii)*NodeLength;
    hist_mNode = length(iNode);
    nNode = length(iNode);
    for jj=1:nNode
        iiNode = iNode(jj);
        iTypeD = iType(iiNode);
        ValueD = Value(iiNode);
        r  = 0;
        x2 = nBranchStartPoint(ii);
        y2 = 1;
        [x1, y1] = funDrawDevice(iTypeD, ValueD, x2, y2, r);
        if mNode>2
            funGenPoint(x2, y2, r);
        end
        isDrawNode(iiNode) = 1;
    end
end
mTemp = mNodeBack-hist_mNode-length(iNode);
mExtNodeNum(ii+1) = mExtNodeNum(ii)+mTemp;
nBranchStartPoint(ii+1) = ii + mExtNodeNum(ii+1)*NodeLength;
% 绘制枝节点
iBranch = find(isDrawNode==0);
nBranch = length(iBranch);
for ii=1:nBranch
    CurBrach = iBranch(ii);
    mN1 = node1(CurBrach);
    mN2 = node2(CurBrach);
    iCurBrach1 = find(histNode==mN1);
    iCurBrach2 = find(histNode==mN2);
    iCurBrach = [iCurBrach1,iCurBrach2];
    if length(iCurBrach) == 1 % 此枝节在此节点上，可以直接垂直绘制
        iiNode = CurBrach;
        iTypeD = iType(iiNode);
        ValueD = Value(iiNode);
        r  = 1;
        y2 = nBranchStartPoint(iCurBrach);
        x2 = 0;
        [x1, y1] = funDrawDevice(iTypeD, ValueD, x2, y2, r);
        isDrawNode(iiNode) = 1;
        if mN1==0 || mN2==0
            [x0, y0]=funGenGND(y2, x2, r);
        else
        end
    elseif length(iCurBrach) == 2 % 此枝节在两个横着的节点上，需要水平绘制
    else % 此枝节不属于任何节点，需要继续找垂直节点上的枝节
    end
end
axis equal;
hold off;
img = [];
% if DispEn
%     x0 = 0;
%     y0 = 0;
%     plot(x0,y0);
%     n = length(iType);
%     for ii=1:n
%         r = mod(ii+1, 2);
%         if ii==1
%             r = 1;
%         elseif ii==2
%             r = 0;
%         end
%         if r == 1
%             y2 = x0;
%             x2 = 0;
%         else
%             x2 = x0;
%             y2 = y0;
%             [x2, y2] = funGenLine([x2,x2+0.25], [y2,y2], 0);
%         end
%         %------------------
%          x0 = x1;
%          y0 = y1;
%     end
%     if Value(end) == 0
%         text(x0+0.1, y0+0.1, 'I_o', 'FontSize',12, 'FontWeight', 'bold');
%     else
%         text(x0+0.1, y0+0.1, 'V_o', 'FontSize',12, 'FontWeight', 'bold');
%     end
%     axis equal;
% end

function [x1, y1] = funDrawDevice(iType, Value, x2, y2, r)
% 0:V,1:I,2:R,3:L,4:C
switch iType
    case 0 % V
        [x1, y1] = funGenV(x2, y2, r);
        funGenText(x1, y1, r, Value, 'V');
        text(x1-0.15, y1+0.12, 'V_i', 'FontSize',12, 'FontWeight', 'bold');
        y0 = y1;
    case 1 % I
        [x1, y1] = funGenI(x2, y2, r);
        funGenText(x1, y1, r, Value, 'A');
        text(x1-0.15, y1+0.12, 'I_i', 'FontSize',12, 'FontWeight', 'bold');
        y0 = y1;
    case 2 % R
        if Value == 0
            [x1, y1] = funGenLine([x2,x2+1], [y2,y2], r);
        elseif Value > 1e24
            [x1, y1] = funGenROpen(x2, y2, r);
        else
            [x1, y1] = funGenR(x2, y2, r);
            funGenText(x1, y1, r, Value, '\Omega');
        end
    case 3 % L
        [x1, y1] = funGenL(x2, y2, r);
        funGenText(x1, y1, r, Value, 'H');
    case 4 % C
        [x1, y1] = funGenC(x2, y2, r);
        funGenText(x1, y1, r, Value, 'F');
end
% if r == 1
%     funGenGND(x0, 0, r);
%     if ii>1 && ii < n
%         funGenPoint(x0, y0, 0);
%     end
% else
%     if ii < n
%         [x1, y1] = funGenLine([x1,x1+0.25], [y1,y1], 0);
%     end
% end

function [x, y]=funGenText(x, y, r, Value, strUnit)
hold on;
if r == 0
    text(x-0.1, y+0.15, [Data2Suffix(Value, '0.2'), strUnit], 'FontSize',12, 'FontWeight', 'bold');
else
    text(x+0.1, y-0.95, [Data2Suffix(Value, '0.2'), strUnit], 'FontSize',12, 'FontWeight', 'bold');
end
hold off;

function [x, y]=funGenPoint(x, y, r)
hold on;
if r == 0
    plot(x, y, '-ok', 'LineWidth', 6);
else
    plot(y, x, '-ok', 'LineWidth', 6);
end
hold off;

function [x0,y0]=funGenLine(x, y, r)
hold on;
if r == 0
    plot(x, y, '-k', 'LineWidth', 2);
    x0 = x(end);
    y0 = y(end);
else
    plot(y, x, '-k', 'LineWidth', 2);
    x0 = y(end);
    y0 = x(end);
end
hold off;

function [x0, y0]=funGenR(x, y, r)
rb = 0.2;
ll = 1.0;
rh = 0.25;
rv = ll-2*rb;
hold on;
if r == 0
    plot([x, x+rb], [y, y], '-k', 'LineWidth', 2);
    plot([x+ll-rb, x+ll], [y, y], '-k', 'LineWidth', 2);
    plot([x+rb, x+rb, x+rb+rv, x+rb+rv, x+rb, x+rb], [y+rh/2, y-rh/2, y-rh/2, y+rh/2, y+rh/2, y-rh/2], '-k', 'LineWidth', 4);
    x0 = x+ll;
    y0 = y;
else
    plot([y, y], [x, x+rb], '-k', 'LineWidth', 2);
    plot([y, y], [x+ll-rb, x+ll], '-k', 'LineWidth', 2);
    plot([y+rh/2, y-rh/2, y-rh/2, y+rh/2, y+rh/2, y-rh/2], [x+rb, x+rb, x+rb+rv, x+rb+rv, x+rb, x+rb], '-k', 'LineWidth', 4);
    x0 = y;
    y0 = x+ll;
end
hold off;

function [x0, y0]=funGenROpen(x, y, r)
rb = 0.2;
ll = 1.0;
rh = 0.25;
rv = ll-2*rb;
hold on;
if r == 0
    plot([x, x+rb], [y, y], '-k', 'LineWidth', 2);
    plot([x+ll-rb, x+ll], [y, y], '-k', 'LineWidth', 2);
%     plot([x+rb, x+rb, x+rb+rv, x+rb+rv, x+rb, x+rb], [y+rh/2, y-rh/2, y-rh/2, y+rh/2, y+rh/2, y-rh/2], '-k', 'LineWidth', 4);
    x0 = x+ll;
    y0 = y;
else
    plot([y, y], [x, x+rb], '-k', 'LineWidth', 2);
    plot([y, y], [x+ll-rb, x+ll], '-k', 'LineWidth', 2);
%     plot([y+rh/2, y-rh/2, y-rh/2, y+rh/2, y+rh/2, y-rh/2], [x+rb, x+rb, x+rb+rv, x+rb+rv, x+rb, x+rb], '-k', 'LineWidth', 4);
    x0 = y;
    y0 = x+ll;
end
hold off;

function [x0, y0]=funGenC(x, y, r)
rb = 0.45;
ll = 1.0;
rh = 0.5;
hold on;
if r == 0
    plot([x, x+rb], [y, y], '-k', 'LineWidth', 2);
    plot([x+ll-rb, x+ll], [y, y], '-k', 'LineWidth', 2);
    plot([x+rb, x+rb], [y+rh/2, y-rh/2], '-k', 'LineWidth', 4);
    plot([x+ll-rb, x+ll-rb], [y+rh/2, y-rh/2], '-k', 'LineWidth', 4);
    x0 = x+ll;
    y0 = y;
else
    plot([y, y], [x, x+rb], '-k', 'LineWidth', 2);
    plot([y, y], [x+ll-rb, x+ll], '-k', 'LineWidth', 2);
    plot([y+rh/2, y-rh/2], [x+rb, x+rb], '-k', 'LineWidth', 4);
    plot([y+rh/2, y-rh/2], [x+ll-rb, x+ll-rb], '-k', 'LineWidth', 4);
    x0 = y;
    y0 = x+ll;
end
hold off;

function [x0, y0]=funGenL(x, y, r)
rb = 0.18;
ll = 1.0;
r0 = ll-rb*2;
lx = linspace(x+rb, x+ll-rb, 49);
hold on;
if r == 0
    plot([x, x+rb], [y, y], '-k', 'LineWidth', 2);
    plot([x+ll-rb, x+ll], [y, y], '-k', 'LineWidth', 2);
    plot(lx, y+0.2.*abs(sin((lx-rb-x)./r0.*4.*pi)).^0.5-0.005, '-k', 'LineWidth', 4);
    x0 = x+ll;
    y0 = y;
else
    plot([y, y], [x, x+rb], '-k', 'LineWidth', 2);
    plot([y, y], [x+ll-rb, x+ll], '-k', 'LineWidth', 2);
    plot(y+0.2.*abs(sin((lx-rb-x)./r0.*4.*pi)).^0.5-0.005, lx, '-k', 'LineWidth', 4);
    x0 = y;
    y0 = x+ll;
end
hold off;

function [x0, y0]=funGenGND(x, y, r)
rb = 0.18;
dy = 0.1;
dx0 = 0.4;
dx1 = 0.25;
dx2 = 0.1;
hold on;
if r == 0
    plot([y, y-rb], [x, x], '-k', 'LineWidth', 2);
    plot([y-rb, y-rb], [x-dx0/2, x+dx0/2], '-k', 'LineWidth', 4);
    plot([y-rb-dy, y-rb-dy], [x-dx1/2, x+dx1/2], '-k', 'LineWidth', 4);
    plot([y-rb-dy*2, y-rb-dy*2], [x-dx2/2, x+dx2/2], '-k', 'LineWidth', 4);
    x0 = y;
    y0 = x;
else
    plot([x, x], [y, y-rb], '-k', 'LineWidth', 2);
    plot([x-dx0/2, x+dx0/2], [y-rb, y-rb], '-k', 'LineWidth', 4);
    plot([x-dx1/2, x+dx1/2], [y-rb-dy, y-rb-dy], '-k', 'LineWidth', 4);
    plot([x-dx2/2, x+dx2/2], [y-rb-dy*2, y-rb-dy*2], '-k', 'LineWidth', 4);
    x0 = x;
    y0 = y;
end
hold off;

function [x0, y0]=funGenV(x, y, r)
rb = 0.18;
ll = 1.0;
d0 = ll-rb*2;
d1 = 0.2;
d2 = 0.15;
lx = linspace(0, 2*pi, 41);
hold on;
if r == 0
    plot([x, x+rb], [y, y], '-k', 'LineWidth', 2);
    plot([x+ll-rb, x+ll], [y, y], '-k', 'LineWidth', 2);
    plot(d0/2.*cos(lx)+x+ll/2, d0/2.*sin(lx)+y, '-k', 'LineWidth', 4);
    plot([x+ll-rb-d1+d2, x+ll-rb-d1], [y, y], '-k', 'LineWidth', 4);
    plot([x+ll-rb-d1+d2/2, x+ll-rb-d1+d2/2], [y-d2/2, y+d2/2], '-k', 'LineWidth', 4);
    plot([x+rb+d1-d2/2, x+rb+d1-d2/2], [y-d2/2, y+d2/2], '-k', 'LineWidth', 4);
    x0 = x+ll;
    y0 = y;
else
    plot([y, y], [x, x+rb], '-k', 'LineWidth', 2);
    plot([y, y], [x+ll-rb, x+ll], '-k', 'LineWidth', 2);
    plot(d0/2.*sin(lx)+y, d0/2.*cos(lx)+x+ll/2, '-k', 'LineWidth', 4);
    plot([y, y], [x+ll-rb-d1+d2, x+ll-rb-d1], '-k', 'LineWidth', 4);
    plot([y-d2/2, y+d2/2], [x+ll-rb-d1+d2/2, x+ll-rb-d1+d2/2], '-k', 'LineWidth', 4);
    plot([y-d2/2, y+d2/2], [x+rb+d1-d2/2, x+rb+d1-d2/2], '-k', 'LineWidth', 4);
    x0 = x;
    y0 = y+ll;
end
hold off;

function [x0, y0]=funGenI(x, y, r)
rb = 0.18;
ll = 1.0;
d0 = ll-rb*2;
d1 = 0.2;
d2 = 0.15;
lx = linspace(0, 2*pi, 41);
hold on;
if r == 0
    plot([x, x+rb], [y, y], '-k', 'LineWidth', 2);
    plot([x+ll-rb, x+ll], [y, y], '-k', 'LineWidth', 2);
    plot(d0/2.*cos(lx)+x+ll/2, d0/2.*sin(lx)+y, '-k', 'LineWidth', 4);
    plot([x+ll-rb-d1+d2, x+rb+d1-d2], [y, y], '-k', 'LineWidth', 4);
    plot([x+ll-rb-d1, x+ll-rb-d1+d2, x+ll-rb-d1], [y-d2/2, y, y+d2/2], '-k', 'LineWidth', 4);
    x0 = x+ll;
    y0 = y;
else
    plot([y, y], [x, x+rb], '-k', 'LineWidth', 2);
    plot([y, y], [x+ll-rb, x+ll], '-k', 'LineWidth', 2);
    plot(d0/2.*sin(lx)+y, d0/2.*cos(lx)+x+ll/2, '-k', 'LineWidth', 4);
    plot([y, y], [x+ll-rb-d1+d2, x+rb+d1-d2], '-k', 'LineWidth', 4);
    plot([y-d2/2, y, y+d2/2], [x+ll-rb-d1, x+ll-rb-d1+d2, x+ll-rb-d1], '-k', 'LineWidth', 4);
    x0 = x;
    y0 = y+ll;
end
hold off;

function [nodeLen, nodeNew, nodeChT, index, histNode] = funFindNodeNext(nodeChT, node1, node2, nodeLen, index, histNode, nodeS)
    itoNode = [find(nodeChT == node1),find(nodeChT == node2)];
    index = index + 1;
    nodeLen(itoNode) = index;
    nodeAll = [node1(itoNode),node2(itoNode)];
    histNode(index) = nodeChT;
    nodeNew = nodeAll(~ismember(nodeAll, [histNode,0]));

    nN = length(nodeNew);
    for jj=1:nN
        nodeChT = nodeNew(jj);
        if nodeChT ~= nodeS
            [nodeLen, nodeNew, nodeChT, index, histNode] = funFindNodeNext(nodeChT, node1, node2, nodeLen, index, histNode, nodeS);
        else
            histNode(index+1) = nodeChT;
            break;
        end
    end

function    [iNode, mNodeBefore, mNodeBack] = funFindNode(mT1, mT2, node1, node2)
findN1 = find(node1 == mT1);
findN2 = find(node2 == mT2);
if ~isempty(findN2)
    iNode1 = findN1(findN1 == findN2);
else
    iNode1 = [];
end
findN1 = find(node1 == mT2);
findN2 = find(node2 == mT1);
if ~isempty(findN2)
    iNode2 = findN1(findN1 == findN2);
else
    iNode2 = [];
end
iNode = [iNode1, iNode2];
% 前面有多少个节点
NodeAll = [node1, node2];
ix = find(mT1 == NodeAll);
mNodeBefore = length(ix);
% 后面有多少个节点
ix = find(mT2 == NodeAll);
mNodeBack = length(ix);
