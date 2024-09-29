clear;clc;close all
% 数据预处理
tic

evr_list = dir('F:\ywj\matlab\map_NASA_gel\evr\*.nc');
evr_num = length(evr_list);

RH_aver = zeros(361,576);
for n = 1:evr_num
    file = ['F:\ywj\matlab\map_NASA_gel\evr\',evr_list(n).name];
    T_aver = zeros(361,576);
    Y_aver = zeros(361,576);
    P_aver = zeros(361,576);
    
    T = ncread(file, 'T2M');
    Y = ncread(file, 'QV2M');
    P = ncread(file, 'PS');
    
    T = rot90(fliplr(T));
    Y = rot90(fliplr(Y));
    P = rot90(fliplr(P));
    for k = 1:24
        T_aver = T_aver + T(:,:,k)/24;
        Y_aver = Y_aver + Y(:,:,k)/24;
        P_aver = P_aver + P(:,:,k)/24;
    end
    for i = 1:size(T,1)
        for j = 1:size(T,2)
            RH_aver(i,j) = RH_cal(T_aver(i,j),Y_aver(i,j),P_aver(i,j));
        end
    end
    save(['F:\ywj\matlab\map_NASA_gel\SEC_par_gel\evr_data_save\evr',num2str(n),'.mat'],'T_aver','RH_aver','Y_aver','P_aver');
    
end


toc
