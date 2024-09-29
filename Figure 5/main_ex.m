clear;clc;close all
%能耗主函数
tic

iso_list = dir('F:\ywj\matlab\map_NASA_gel\Isotherm\*.txt');
iso_num = length(iso_list);

p = parpool('local',40);
for n = 181:365
    load(['F:\ywj\matlab\map_NASA_gel\SEC_par_gel\evr_data_save\evr',num2str(n),'.mat']);
    Ex_aver = ones(361,576,6) * 1e7;
    parfor i = 1:361
        for j = 1:576
%             Exergy = 0;
            
            if T_aver(i,j)<273.15 || isnan(RH_aver(i,j))
                continue;
            end
            
            for k = 1:iso_num
                
                data_name = ['F:\ywj\matlab\map_NASA_gel\Isotherm\',iso_list(k).name];
                Exergy = cal_SBCool(RH_aver(i,j),T_aver(i,j)-273.15,P_aver(i,j),data_name);
                
                if Exergy == 0
                    continue;
                end
                
                if Exergy < Ex_aver(i,j,k)
                    Ex_aver(i,j,k) = Exergy;

                end
                
            end
            
        end
    end
    save(['F:\ywj\matlab\map_NASA_gel\SEC_par_gel\SEC_data_save\evr',num2str(n),'.mat'],'Ex_aver');
end
delete(p);
toc



