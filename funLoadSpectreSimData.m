%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2022-08-17(yyyy-mm-dd)
% 读取spectre结果
%--------------------------------------------------------------------------
function [t, v] = funLoadSpectreSimData(strPath)
% strPath = './spectre_sim_data/V0.txt';
fId = fopen(strPath, 'r');
strData = fread(fId, '*char')';
fclose(fId);
cellData = regexp(strData, '\r\n', 'split');
cellData0 = cellData(7:end);
m = length(cellData0);
t = [];
v = [];
for ii=1:m
    strData = cellData0{ii};
    cellData1 = regexp(strData, ',', 'split');
    t(ii) = str2double(cellData1{1});
    v(ii) = str2double(cellData1{2});
end