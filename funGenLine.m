%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-16(yyyy-mm-dd)
% 生成线段
%--------------------------------------------------------------------------
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