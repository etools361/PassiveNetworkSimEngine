%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-16(yyyy-mm-dd)
% 生成点
%--------------------------------------------------------------------------
function [x, y]=funGenPoint(x, y, r)
hold on;
if r == 0
    plot(x, y, '-ok', 'LineWidth', 6);
else
    plot(y, x, '-ok', 'LineWidth', 6);
end
hold off;