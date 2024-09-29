function [wt] = wt_par(i,j,T0,Y0,P0,R0,T1,Y1,P1,R1,T2,Y2,P2,R2,iso_num,iso_list)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
wt = zeros(1,1,3);

T(1:24,1) = T0(i,j,:);T(25:48,1) = T1(i,j,:);T(49:72,1) = T2(i,j,:);
Y(1:24,1) = Y0(i,j,:);Y(25:48,1) = Y1(i,j,:);Y(49:72,1) = Y2(i,j,:);
P(1:24,1) = P0(i,j,:);P(25:48,1) = P1(i,j,:);P(49:72,1) = P2(i,j,:);
R(1:24,1) = R0(i,j,:);R(25:48,1) = R1(i,j,:);R(49:72,1) = R2(i,j,:);

% 总运行时间 8*60/5=96
max_t = 192/2;

% 时差纠正，原数据为0：30-71：30，每0.625度有两分半时差
% 选择当地时间第二天的8：00-16：00，共8小时，8*60/2.5=192
x = 1:1/24:72;
T_amb = interp1(1:1:72,T,x','spline');
Y_amb = interp1(1:1:72,Y,x','spline');
P_amb = interp1(1:1:72,P,x','spline');
R_amb = interp1(1:1:72,R,x','spline');

T_amb = T_amb(769-(j-289)-12 :2: 961-(j-289)-12,1);
Y_amb = Y_amb(769-(j-289)-12 :2: 961-(j-289)-12,1);
P_amb = P_amb(769-(j-289)-12 :2: 961-(j-289)-12,1);
R_amb = R_amb(769-(j-289)-12 :2: 961-(j-289)-12,1);

ad_temp = zeros(max_t+1,1);
ad_hum = zeros(max_t+1,1);
de_temp = zeros(max_t+1,1);
de_hum = zeros(max_t+1,1);

for k = 1:max_t+1
    ad_temp(k,1) = T_amb(k,1)+3;
    ad_hum(k,1) = RH_cal(ad_temp(k,1),Y_amb(k,1),P_amb(k,1));
    ad_hum(k,1) = 100*exp(log(ad_hum(k,1)/100)*(ad_temp(k,1)/298.15));
    if ad_hum(k,1) > 100
        ad_hum(k,1) = 100;
    elseif ad_hum(k,1) < 0
        ad_hum(k,1) = 0;
    end
    
    Tc = T_amb(k,1)+18-273.15;
    de_temp(k,1) = T_amb(k,1)+sqrt((3.555/0.0580)^2+(1.85*0.5*R_amb(k))/(1.85*0.029))-3.555/0.058-9;
    Psat = exp(-1.0440397e4/(Tc*9/5+32+459.67) ...
        -1.129465e1 ...
        -2.7022352e-2*(Tc*9/5+32+459.67) ...
        +1.289036e-5*(Tc*9/5+32+459.67)^2 ...
        -2.4780681e-9*(Tc*9/5+32+459.67)^3 ...
        +6.5459673*log(Tc*9/5+32+459.67))*6890;
    Y_de = 621.945/1000 * 100/100 / (P_amb(k,1)/Psat - 100/100);
    de_hum(k,1) = RH_cal(de_temp(k,1),Y_de,P_amb(k,1));
    de_hum(k,1) = 100*exp(log(de_hum(k,1)/100)*(de_temp(k,1)/343.15));
end

%     figure
%     hold on
%     plot(1:289,ad_hum)
%     plot(1:145,de_hum)

% 动态吸附计算
for u = 1:iso_num
    Iso = load(['F:\ywj\matlab\map_NASA_gel\Isotherm\',iso_list(u).name]);
    if u == 1
        k1 = 0.0055;k2 = 0.094;
        tAD = 105/5; tDE = 35/5;
    elseif u == 2
        k1 = 0.0152;k2 = 0.049;
        tAD = 80/5; tDE = 50/5;
    elseif u ==3
        k1 = 0.0108;k2 = 0.047;
        tAD = 90/5; tDE = 45/5;
    else
        disp('error');
        break;
    end
    
    % 温度储存
    Td = ones(tAD+tDE,1)*ad_temp(1,1);
    % 数据储存，行；时间单元，列：空间单元
    WD0 = zeros(tAD+tDE,max_t+1);
    % WD_flag 表示单位吸附剂的状态0 - AD, 1 - DE
    WD_f = zeros(tAD+tDE,1);
    WD_f(1:tDE) = 1;
    
    % 开始旋转
    for k = 2:size(WD0,2)
        % 状态翻转
        if mod(k+tDE-1,tAD+tDE) == 0
            WD_f(tAD+tDE) = 1;
        else
            WD_f(mod(k+tDE-1,tAD+tDE)) = 1;
        end
        if mod(k-1,tAD+tDE) == 0
            WD_f(tAD+tDE) = 0;
        else
            WD_f(mod(k-1,tAD+tDE)) = 0;
        end
        
        % 更新Td
        for ss = 1:tAD+tDE
            if WD_f(ss) == 0 && Td(ss) > ad_temp(k,1)
                Td(ss) = Td(ss) + ad_temp(k,1) - ad_temp(k-1,1) - 20;
                if Td(ss) < ad_temp(k,1)
                    Td(ss) = ad_temp(k,1);
                end
            end
            if WD_f(ss) == 1 && Td(ss) <= de_temp(k,1)
                Td(ss) = Td(ss) + de_temp(k,1) - de_temp(k-1,1) + 20;
                if Td(ss) > de_temp(k,1)
                    Td(ss) = de_temp(k,1);
                end
            end
        end
        
        for s_num = 1:size(WD_f,1)
            % 对于吸附，低于0，不变，低于55，按照吸附自动变化，高于55，热响应但不冷凝
            if WD_f(s_num) == 0
                if ad_temp < 273.15
                    WD0(s_num,k) = WD0(s_num,k-1);
                elseif ad_temp < 328.15
                    WD_eq = interp1(Iso(:,1),Iso(:,2),ad_hum(k),'linear');
                    WD0(s_num,k) = (1-exp(-k1 * 2.5)) * (WD_eq - WD0(s_num,k-1)) + WD0(s_num,k-1);
                else
                    WD_eq = interp1(Iso(:,3),Iso(:,4),ad_hum(k),'linear');
                    if WD0(s_num,k-1) > WD_eq
                        WD0(s_num,k) = (1-exp(-k2 * 2.5*2)) * (WD_eq - WD0(s_num,k-1)) + WD0(s_num,k-1);
                    else
                        WD0(s_num,k) = WD0(s_num,k-1);
                    end
                end
            else
                % 对于解吸，低于0，不变，低于55，不变，高于55，热响应且冷凝
                if de_temp < 328.15
                    WD0(s_num,k) = WD0(s_num,k-1);
                else
                    WD_eq = interp1(Iso(:,3),Iso(:,4),de_hum(k),'linear');
                    if WD0(s_num,k-1) > WD_eq
                        WD0(s_num,k) = (1-exp(-k2 * 2.5*2)) * (WD_eq - WD0(s_num,k-1)) + WD0(s_num,k-1);
                        wt(1,1,u) = wt(1,1,u) + (WD0(s_num,k-1) - WD0(s_num,k))/size(WD_f,1);
                    else
                        WD0(s_num,k) = WD0(s_num,k-1);
                    end
                end
            end
        end
    end
end
% end u

