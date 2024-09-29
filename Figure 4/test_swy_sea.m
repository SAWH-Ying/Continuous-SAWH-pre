clear;clc;close all
% 取水地图主程序
% 全年不换吸附剂
% CRPS改，根据全年的波动性跟无法极限取水的情况剔除
tic

evr_list = dir('F:\ywj\matlab\map_NASA_gel\evr\*.nc');
iso_list = dir('F:\ywj\matlab\map_NASA_gel\Isotherm\*.txt');
evr_num = length(evr_list);iso_num = length(iso_list);

file = ['F:\ywj\matlab\map_NASA_gel\evr\',evr_list(1).name];
lats = ncread(file, 'lat');
lons = ncread(file, 'lon');

% OP = zeros(361,576,iso_num);
%%
RSD_threshold = 1;
wt_threshold = 0.08;

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
m_num = 1;
for n = 1:365
    load(['F:\ywj\matlab\map_NASA_gel\wt_gel_conti\evr_data_save\',num2str(n),'.mat']);
    wt(:,577,:) = wt(:,1,:);
    
    if n > sum(month(1:m_num))
        m_num = m_num + 1;
    end
    test_data1(:,:,m_num) = test_data1(:,:,m_num) + wt(:,:,1)/month(m_num);
    test_data2(:,:,m_num) = test_data2(:,:,m_num) + wt(:,:,2)/month(m_num);
    test_data3(:,:,m_num) = test_data3(:,:,m_num) + wt(:,:,3)/month(m_num);
end
test_data1(:,:,1) = test_data1(:,:,1)*59/90 + test_data1(:,:,5)*31/90;
test_data2(:,:,1) = test_data2(:,:,1)*59/90 + test_data2(:,:,5)*31/90;
test_data3(:,:,1) = test_data3(:,:,1)*59/90 + test_data3(:,:,5)*31/90;

test_data1(:,:,5) = [];
test_data2(:,:,5) = [];
test_data3(:,:,5) = [];

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
OP1 = zeros(361,577);OP2 = zeros(361,577);sss = 3;
for i = 1:361
    for j = 1:577
        test = [test_data1(i,j,sss),test_data2(i,j,sss),test_data3(i,j,sss)];
        [a,b] = max(test);
        OP1(i,j) = a;
        if a == 0
            OP1(i,j) = 0/0;
            OP2(i,j) = 0/0;
        elseif a < 0.08
            OP1(i,j) = 0;
            OP2(i,j) = 0;
        else
            OP2(i,j) = b;
        end
    end
end


% % 平均取水
mycol = [128 128 126;
    206,204,199;
    206,204,199;
    210,227,243;
    170,207,229;
    104,172,213;
    56,136,192;
    16,92,164;
    8,51,110]/255;

mycolor = interp1(linspace(0,1,size(mycol,1)),mycol,linspace(0,1,256),'cubic');
mycolor(mycolor<0) = 0;
mycolor(mycolor>1) = 1;

figure
set(gca,'LooseInset',[0 0 0 0])
pcolor(lons,lats,OP1);
shading flat
c = colorbar;
colormap(mycolor)
caxis([-1 2.5])
hold on
title(' ')
plot(long,lat,'black')
set(gcf,'position',[150,150,1200,600]);
set(gca,'Xtick',(-180:90:180));
set(gca,'Ytick',(-90:45:90));

%对应最多次吸附剂
mycol = [
    128 128 126;
    206,204,199;
    230,111,81;
    243,162,97;
    232,197,107;
    138,176,125;
    41,157,143;
    40,114,113;]/255;
figure
set(gca,'LooseInset',[0 0 0 0])
pcolor(lons,lats,OP2);
shading flat
c = colorbar;
colormap(mycol)
caxis([-1 7])
hold on
title(' ')
plot(long,lat,'black')
set(gcf,'position',[150,150,1200,600]);
set(gca,'Xtick',(-180:90:180));
set(gca,'Ytick',(-90:45:90));
