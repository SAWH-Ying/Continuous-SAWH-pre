clear;clc;close all
% SEC地图主程序
% 按照四个挡位区分取水能力
% 对一个月内无法取水的情况进行补偿
tic

evr_list = dir('F:\ywj\matlab\map_NASA_gel\evr\*.nc');
iso_list = dir('F:\ywj\matlab\map_NASA_gel\Isotherm\*.txt');
evr_num = length(evr_list);iso_num = length(iso_list);

file = ['F:\ywj\matlab\map_NASA_gel\evr\',evr_list(1).name];
lats = ncread(file, 'lat');
lons = ncread(file, 'lon');

% OP = zeros(361,576,iso_num);
%%
RSD_threshold = 1e7;
SEC_threshold = 2e3;

lons(577,1) = 180;
%% figure
% 海岸线修正
long = zeros(9865,1);
load coast
i = 1;len = length(long);
while i < len
    if long(i) <= 180 && long(i+1) > 180
        lat0 = (lat(i)*(long(i+1)-180)+lat(i+1)*(180-long(i)))/(long(i+1)-long(i));
        for j = len:-1:(i+1)
            long(j+3) = long(j);
            lat(j+3) = lat(j);
        end
        long(i+2) = 0/0;lat(i+2) = 0/0;
        long(i+1) = 180;lat(i+1) = lat0;
        long(i+3) = -180;lat(i+3) = lat0;
        len = len + 3; i = i + 3;
    elseif long(i) > 180 && long(i+1) <= 180
        lat0 = (lat(i)*(long(i+1)-180)+lat(i+1)*(180-long(i)))/(long(i+1)-long(i));
        for j = len:-1:(i+1)
            long(j+3) = long(j);
            lat(j+3) = lat(j);
        end
        long(i+2) = 0/0;lat(i+2) = 0/0;
        long(i+1) = -180;lat(i+1) = lat0;
        long(i+3) = 180;lat(i+3) = lat0;
        len = len + 3; i = i + 3;
    end
    i = i + 1;
end
for i = 1:len
    if long(i) > 180
        long(i) = long(i) - 360;
    end
end
%% 循环计算各个陆地

test_data1 = zeros(361,577,5);
test_data2 = zeros(361,577,5);
test_data3 = zeros(361,577,5);

% month = [31,28,31,30,31,30,31,31,30,31,30,31];
month = [59,92,92,91,31];
m_num = 1;num_tra = zeros(361,577,6);
% change made
for n = 1:365
    % Reads the previously calculated daily energy consumption map
    load(['F:\ywj\matlab\map_NASA_gel\SEC_par_gel\SEC_data_save\evr',num2str(n),'.mat']);
    Ex_aver(:,577,:) = Ex_aver(:,1,:);

    num_tra = num_tra + (Ex_aver < 1e7);
    
    if n > sum(month(1:m_num))
        m_num = m_num + 1;
    end
    
    Ex_aver(Ex_aver == 1e7) = 0;
    test_data1(:,:,m_num) = test_data1(:,:,m_num) + Ex_aver(:,:,1);
    test_data2(:,:,m_num) = test_data2(:,:,m_num) + Ex_aver(:,:,2);
    test_data3(:,:,m_num) = test_data3(:,:,m_num) + Ex_aver(:,:,3);
    
    if n == sum(month(1:m_num))
        test_data1(:,:,m_num) = test_data1(:,:,m_num)./num_tra(:,:,1);
        test_data2(:,:,m_num) = test_data2(:,:,m_num)./num_tra(:,:,2);
        test_data3(:,:,m_num) = test_data3(:,:,m_num)./num_tra(:,:,3);
        num_tra = zeros(361,577,6);
    end
end

test_data1(:,:,1) = test_data1(:,:,1)*59/90 + test_data1(:,:,5)*31/90;
test_data2(:,:,1) = test_data2(:,:,1)*59/90 + test_data2(:,:,5)*31/90;
test_data3(:,:,1) = test_data3(:,:,1)*59/90 + test_data3(:,:,5)*31/90;

test_data1(:,:,5) = [];
test_data2(:,:,5) = [];
test_data3(:,:,5) = [];

test_data1(isnan(test_data1)) = 0;
test_data2(isnan(test_data2)) = 0;
test_data3(isnan(test_data3)) = 0;

OP_ave1 = zeros(361,577);OP_RSD1 = zeros(361,577);
OP_ave2 = zeros(361,577);OP_RSD2 = zeros(361,577);
OP_ave3 = zeros(361,577);OP_RSD3 = zeros(361,577);

for i = 1:361
    for j = 1:577
        [OP_ave1(i,j),OP_RSD1(i,j)] = RSD(test_data1(i,j,:));
        [OP_ave2(i,j),OP_RSD2(i,j)] = RSD(test_data2(i,j,:));
        [OP_ave3(i,j),OP_RSD3(i,j)] = RSD(test_data3(i,j,:));
    end
end

%%
OP1 = zeros(361,577);OP2 = zeros(361,577);
for i = 1:361
    for j = 1:577
        
        if OP_ave1(i,j) == 0
            OP_ave1(i,j) = 0/0;
            OP_RSD1(i,j) = 0/0;
        elseif min(test_data1(i,j,:)) == 0
            OP_ave1(i,j) = 4e3;
            OP_RSD1(i,j) = 4e3;
        elseif OP_RSD1(i,j) >= RSD_threshold || max(test_data1(i,j,:)) >= SEC_threshold
            OP_ave1(i,j) = 3e3;
            OP_RSD1(i,j) = 3e3;
        end
        
        if OP_ave2(i,j) == 0
            OP_ave2(i,j) = 0/0;
            OP_RSD2(i,j) = 0/0;
        elseif min(test_data2(i,j,:)) == 0
            OP_ave2(i,j) = 4e3;
            OP_RSD1(i,j) = 4e3;
        elseif OP_RSD2(i,j) >= RSD_threshold || max(test_data2(i,j,:)) >= SEC_threshold
            OP_ave2(i,j) = 3e3;
            OP_RSD2(i,j) = 3e3;
        end
        
        if OP_ave3(i,j) == 0
            OP_ave3(i,j) = 0/0;
            OP_RSD3(i,j) = 0/0;
        elseif min(test_data3(i,j,:)) == 0
            OP_ave3(i,j) = 4e3;
            OP_RSD3(i,j) = 4e3;
        elseif OP_RSD3(i,j) >= RSD_threshold || max(test_data3(i,j,:)) >= SEC_threshold
            OP_ave3(i,j) = 3e3;
            OP_RSD3(i,j) = 3e3;
        end
        
                
        test = [OP_ave1(i,j),OP_ave2(i,j),OP_ave3(i,j)];
        [a,b] = min(test);
        OP1(i,j) = a;
        if isnan(a)
            OP2(i,j) = 0/0;
        elseif a == 3e3
            OP2(i,j) = 0;
        elseif a == 4e3
            OP2(i,j) = -1;
        else
            OP2(i,j) = b;
        end
    end
end

%%
nsum = 0;n = 0;

for i = 60:361
    for j = 1:577
        if inpolygon(lats(i),lons(j),lat,long) == 0
            continue;
        end
        
        nsum = nsum + 1;
        
%         if OP1(i,j)/1000 <=1.4
%             n = n + 1;
%         end
        
        if OP1(i,j) == 3000
            n = n + 1;
        end
        
    end
end





toc