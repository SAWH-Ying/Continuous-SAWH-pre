clear;clc;close all
% nano energy验证
% 工况恒定的连续式系统
% 输入tAD和tDE，输出连续图像
% OP - wt随时间变化
% OP2 - 总取水量
% OP3 - 每分钟累计取水量
tic
%% Input

Iso = load('nano.txt');
% 吸附速率常数[min -1]
kAD = 0.0018;
kDE = 0.0819;


% 转一圈所需的时间,min
tAD = 67;
tDE = 23;
t_total = 4*60;
% 线性升温 [oC]
T_delta = 10;
% 初始含湿量
WT0 = 0.7;
% 温度[oC]
Tamb = 22;
TAD = Tamb + 3;
TDE = Tamb + sqrt((3.555/0.0580)^2+(1.85*0.5*600)/(1.85*0.029))-3.555/0.058-9;
TCD = Tamb + 18;
% Humidity ratio [kg/kga]
YAD = 9.9/1000;
% RH [%]
RHAD = RH_cal(TAD + 273.15,YAD,101325);
for YCD = 0:1e-4:1e4
    RHCD = RH_cal(TCD + 273.15,YCD/1000,101325);
    if RHCD == 100
        break;
    end
end
YCD = YCD/1000;
RHDE = RH_cal(TDE + 273.15,YCD,101325);
% 温度[oC]
T = ones(tAD+tDE,1)*TAD;
%independent variable
t = 0:1:t_total;t = t';
% 记录每个单元每个时间段的含湿量
WD = ones(t_total+1,tAD+tDE)*WT0;
% WD_flag 表示吸附剂的状态0 - AD, 1 - DE
WD_f = zeros(tAD+tDE,1);
% 计算累计的取水量
WT = zeros(t_total+1,1);
% 假设unit1到unittDE一开始处于解吸状态
WD_f(1:tDE) = 1;

%% 双循环更新状态
for k = 2:t_total+1
    % 更新WD_f
    if  mod(k+tDE-1,tAD+tDE) == 0
        WD_f(tAD+tDE) = 1;
    else
        WD_f(mod(k+tDE-1,tAD+tDE)) = 1;
    end
    if mod(k-1,tAD+tDE) == 0
        WD_f(tAD+tDE) = 0;
    else
        WD_f(mod(k-1,tAD+tDE)) = 0;
    end
    % 更新T
    for n = 1:tAD+tDE
        if WD_f(n) == 0 && T(n) > TAD
            T(n) = T(n) - T_delta;
            if T(n) < TAD
                T(n) = TAD;
            end
        end
        if WD_f(n) == 1 && T(n) <= TDE
            T(n) = T(n) + T_delta;
            if T(n) > TDE
                T(n) = TDE;
            end
        end
    end

    
    for s_num = 1:(tAD+tDE)
        if WD_f(s_num) == 0
            RH = RH_cal(T(s_num)+273.15,YAD,101325);
            WD_eq = interp1(Iso(:,1),Iso(:,2),RH,'linear');
            WD(k,s_num) = (1-exp(-kAD)) * (WD_eq - WD(k-1,s_num)) + WD(k-1,s_num);
        else
            if T(s_num) < 55
                WD(k,s_num) = WD(k-1,s_num);
            elseif T(s_num) >= 55
                RH = RH_cal(T(s_num)+273.15,YCD,101325);
                WD_eq = interp1(Iso(:,1),Iso(:,2),RH,'linear');
                if WD(k-1,s_num) > WD_eq
                    WD(k,s_num) = (1-exp(-kDE)) * (WD_eq - WD(k-1,s_num)) + WD(k-1,s_num);
                    WT(k) = WT(k) + (WD(k-1,s_num) - WD(k,s_num))/(tAD+tDE);
                else
                    WD(k,s_num) = WD(k-1,s_num);
                end
            end
        end
    end
    
end
OP = WD;
OP2 = sum(WT);
OP3 = zeros(t_total+1,1);
for i = 1:t_total+1
    OP3(i) = sum(WT(1:i));
end
%% Output
mesh(OP)
toc