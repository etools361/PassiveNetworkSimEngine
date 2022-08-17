%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-17(yyyy-mm-dd)
% 仿真结果和spectre对比
%--------------------------------------------------------------------------
figure(5);
strPath = './spectre_sim_data/V0.txt';
[t1, vi] = funLoadSpectreSimData(strPath);
strPath = './spectre_sim_data/VRL.txt';
[t2, vo] = funLoadSpectreSimData(strPath);
% plot(t1, vi, '-r', 'Linewidth', 2);
plot(t2, vo, '-b', 'Linewidth', 2);
hold on;
plot(t,X(b2,:),'-r', 'LineWidth', 1);
hold off;
grid on;
legend({'Spectre Sim', 'Matlab Sim'}, 'location', 'northwest');
xlabel('Time/s');
ylabel('V_o/V');
title('V_o VS. t');
% xlim([10.5382,10.5384])
% ylim([-2e-4,2e-4])

