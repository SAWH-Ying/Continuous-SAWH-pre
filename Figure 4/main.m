clear;clc;close all
% 单天，全地图，全吸附剂的吸附/解吸循环图
tic

evr_list = dir('F:\ywj\matlab\map_NASA_gel\evr\*.nc');
rad_list = dir('F:\ywj\matlab\map_NASA_gel\rad\*.nc');
iso_list = dir('F:\ywj\matlab\map_NASA_gel\Isotherm\*.txt');
evr_num = length(evr_list);iso_num = length(iso_list);

file = ['F:\ywj\matlab\map_NASA_gel\evr\',evr_list(1).name];
lats = ncread(file, 'lat');
lons = ncread(file, 'lon');
% 输出初始化
% WD = zeros(361,576,3);

% 并行运算
p = parpool('local',40);
for n = 1:365
    wt = zeros(361,576,3);
    if n == 1
        file0 = ['F:\ywj\matlab\map_NASA_gel\evr\',evr_list(365).name];
        file1 = ['F:\ywj\matlab\map_NASA_gel\evr\',evr_list(n).name];
        file2 = ['F:\ywj\matlab\map_NASA_gel\evr\',evr_list(n+1).name];
    elseif n == evr_num
        file0 = ['F:\ywj\matlab\map_NASA_gel\evr\',evr_list(n-1).name];
        file1 = ['F:\ywj\matlab\map_NASA_gel\evr\',evr_list(n).name];
        file2 = ['F:\ywj\matlab\map_NASA_gel\evr\',evr_list(1).name];
    else
        file0 = ['F:\ywj\matlab\map_NASA_gel\evr\',evr_list(n-1).name];
        file1 = ['F:\ywj\matlab\map_NASA_gel\evr\',evr_list(n).name];
        file2 = ['F:\ywj\matlab\map_NASA_gel\evr\',evr_list(n+1).name];
    end
    % 读取nc4文件内的信息
    % ncdisp(file)
    
    % 2-meter_air_temperature
    T0 = ncread(file0, 'T2M');T1 = ncread(file1, 'T2M');T2 = ncread(file2, 'T2M');
    Y0 = ncread(file0, 'QV2M');Y1 = ncread(file1, 'QV2M');Y2 = ncread(file2, 'QV2M');
    P0 = ncread(file0, 'PS');P1 = ncread(file1, 'PS');P2 = ncread(file2, 'PS');
    
    T0 = rot90(fliplr(T0));T1 = rot90(fliplr(T1));T2 = rot90(fliplr(T2));
    Y0 = rot90(fliplr(Y0));Y1 = rot90(fliplr(Y1));Y2 = rot90(fliplr(Y2));
    P0 = rot90(fliplr(P0));P1 = rot90(fliplr(P1));P2 = rot90(fliplr(P2));
    
    if n == 1
        file0 = ['F:\ywj\matlab\map_NASA_gel\rad\',rad_list(365).name];
        file1 = ['F:\ywj\matlab\map_NASA_gel\rad\',rad_list(n).name];
        file2 = ['F:\ywj\matlab\map_NASA_gel\rad\',rad_list(n+1).name];
    elseif n == evr_num
        file0 = ['F:\ywj\matlab\map_NASA_gel\rad\',rad_list(n-1).name];
        file1 = ['F:\ywj\matlab\map_NASA_gel\rad\',rad_list(n).name];
        file2 = ['F:\ywj\matlab\map_NASA_gel\rad\',rad_list(1).name];
    else
        file0 = ['F:\ywj\matlab\map_NASA_gel\rad\',rad_list(n-1).name];
        file1 = ['F:\ywj\matlab\map_NASA_gel\rad\',rad_list(n).name];
        file2 = ['F:\ywj\matlab\map_NASA_gel\rad\',rad_list(n+1).name];
    end
    
    R0 = ncread(file0, 'SWGDN');R1 = ncread(file1, 'SWGDN');R2 = ncread(file2, 'SWGDN');
    R0 = rot90(fliplr(R0));R1 = rot90(fliplr(R1));R2 = rot90(fliplr(R2));
    
    parfor i = 1:361
        for j = 1:576
            [wt_s] = wt_par(i,j,T0,Y0,P0,R0,T1,Y1,P1,R1,T2,Y2,P2,R2,iso_num,iso_list);
            wt(i,j,:) = wt_s;
        end
    end
    save(['F:\ywj\matlab\map_NASA_gel\wt_gel_conti\evr_data_save\',num2str(n),'.mat'],'wt')
end
delete(p);

toc