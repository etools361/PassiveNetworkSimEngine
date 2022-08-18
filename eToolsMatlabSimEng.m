%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-07-05(yyyy-mm-dd)
% 通用无源网络仿真引擎
%--------------------------------------------------------------------------
% netlist, PI型，两端接载
strNetlist = {
'V0 V 1 0 1';
'RS R 1 2 1';
'R1 R 2 3 0';
'C1 C 3 0 1';
'L2 L 3 4 2';
'C3 C 4 0 1';
'R2 R 4 5 0';
'RL R 5 0 1'
};
% T型，一端接载
strNetlist1 = {
'V0 V 1 GND 1';
'RS R 1 7 1';
'CS C 7 6 1';
'RS2 R 1 8 1';
'LS L 8 6 1';
'RP R 5 6 10';
'C3 C 3 x 1.333';
'C4 C 3 y 0.2';
'R5 R x y 1.1';
'R6 R y z 1.2';
'R7 R z GND 1.4';
'L2 L 6 3 1.5';
'L4 L 3 5 0.5';
'L5 L z 7 0.5';
'RL R 5 GND 1'
};
strNetlist1 = {
'V0 V 1 0 1';
'RS R 1 2 1';
'R1 R 2 3 0';
'C1 C 3 8 1';
'R3 R 8 9 1';
'L5 L 9 0 1';
'C4 C 8 0 0.2';
'L2 L 3 4 2';
'C3 C 4 0 1';
'R2 R 4 5 0.1';
'RL R 5 0 1'
};
strNetlist1 = {
'V0 V 1 0 1';
'RS R 1 2 1';
'C2 C 2 0 0.314631';
'C3 C 2 3 0.02996';
'L3 L 2 3 0.155487';
'C4 C 3 0 0.396551';
'C5 C 4 3 0.082683';
'L5 L 4 3 0.126308';
'C6 C 4 0 0.314631';
'RL R 4 0 1'
};
% Butterworth BPF
strNetlist1 = {
'V0 V 1 0 1';
'RS R 1 2 1';
'C2 C 2 0 0.491582';
'L2 L 2 0 0.052049';
'C3 C 2 3 0.019881';
'L3 L 3 4 1.287';
'C4 C 4 0 1.591';
'L4 L 4 0 0.016084';
'C5 C 4 5 0.019881';
'L5 L 5 6 1.287';
'C6 C 6 0 0.491582';
'L6 L 6 0 0.052049';
'RL R 6 0 1'
};
strNetlist1 = {
'V0 V 1 0 1';
'RS R 1 2 1';
'C2 C 2 3 0.1';
'R2 R 2 3 0.2';
'C3 C 4 3 1';
'R3 R 4 3 0.2';
'C4 C 4 5 10';
'R4 R 4 5 0.2';
'C5 C 6 5 100';
'R5 R 6 5 0.2';
'RL R 6 0 1'
};
strNetlist1 = {
'V0 V 1 0 1';
'RS R 1 2 1';
'R1 R 3 2 1';
'R2 R 2 0 1';
'R3 R 1 3 10';
'RL R 3 0 1'
};
% split netlist
[iType, Value, cellNode1, CellNode2, cellName] = funSimNetlist2Array(strNetlist);
% netlist standard
[node1, node2] = funSimNetlistRenode(cellNode1, CellNode2);
% netlist analysis
[maxNode, nL, nI, nV, nR0] = funSimNetlistAna(iType, Value, node1, node2);
% 标记给定器件的非GND节点用于获取结果位置
strDevice = 'RL'; % 显示结果的器件
[retNode] = funGetDeviceNode(cellName, node1, node2, strDevice);
% generate matrix, (MM*s+MN)*MX = MV
[MM, MN, MV, MX] = funSimNetlist2Matrix(iType, Value, node1, node2, maxNode, nL, nI, nV, nR0, cellName);

%% -----------------Gen Sch-------------------------
figure(1);
f0 = 1e-3/2/pi;
w0 = 2*pi*f0;
[img] = funGenSchematic2(iType, Value, cellNode1, CellNode2, cellName, w0);
tic;
ylimValue = ylim;
ylim([-ceil(-ylimValue(1)), ceil(ylimValue(2))]);
set(gca,'ytick',[]);
set(gca,'xtick',[]);
set(gca,'Box','off');
set(gcf,'color',[1,1,1]);
% set(gca,'looseInset',[0 0 0 0]);
axis off;
% set(gca,'LooseInset',get(gca,'TightInset'));
% set(gca,'PlotBoxAspectRatio', [8 2 1])
% figure(2);
%% -----------------AC-------------------------
f0 = 1e-2;
f1 = 1e1;
% f0 = 1e4;
% f1 = 1e6;
N  = 500;
freq = logspace(log10(f0), log10(f1), N);
Vo = zeros(1, N);
[a, b]   = ismember('RL',cellName);
[a2, b2] = ismember('iRL',MX);
RL  = Value(b);
nRL = max(node1(b), node2(b));
for ii=1:N
    f = freq(ii);
    s = 1i*2*pi*f;
    V = (MM.*s + MN)\MV';
    if RL == 0
        Vo(ii) = V(b2);
    else
        Vo(ii) = V(nRL);
    end
end
dBVo = 20*log10(abs(Vo));
AgVo = angle(Vo)*180/pi;
uWVo = unwrap(AgVo, 179);
% uWVo = AgVo;
toc;
figure(2)
semilogx(freq, dBVo, '-r', 'LineWidth', 2);
grid on;
xlabel('Freq/Hz');
if RL == 0
    ylabel('I_o Mag/dB');
else
    ylabel('V_o Mag/dB');
end
title('FrequencyResponse');
xlim([min(freq),max(freq)]);
ylim([-80,0]);
figure(3);
semilogx(freq, uWVo, '-r', 'LineWidth', 2);
xlim([min(freq),max(freq)]);
grid on;
xlabel('Freq/Hz');
if RL == 0
    ylabel('I_o Angle/deg');
else
    ylabel('V_o Angle/deg');
end
title('Phase VS. Freq');

%% -----------------Tran-------------------------
f1=2e-2;
% f1=1;
% Tmax=30/f1;
Tmax=2/f1;
Nmax=100000; % 
h=Tmax/Nmax;
t=linspace(0,Tmax,Nmax);   % 构造时间序列
E=1+square(2*pi*f1*t-81.915e-3);% 构造脉冲激励信号
X = [];
X(:,1)=zeros(size(MX));
MNT=inv(1/h*MM+MN);
[a, b]   = ismember('iV0',MX);
[a2, b2] = ismember(sprintf('V%dn', retNode),MX);
for n=2:Nmax
    MV(b) = E(n);
    X(:,n)=MNT*(MV'+1/h.*MM*X(:,n-1));
end
figure(4);
plot(t,E,'-b', 'LineWidth', 2);
hold on;
plot(t,X(b2,:),'-r', 'LineWidth', 2);
hold off;
xlim([min(t),max(t)]);
grid on;
xlabel('Time/s');
if RL == 0
    ylabel('I_o/A');
else
    ylabel('V_o/V');
end
title('V_o VS. t');
legend({'Vi', 'Vo'}, 'location', 'northeast');
