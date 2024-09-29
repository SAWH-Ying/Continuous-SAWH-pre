function [exergy_OP] = cal_SBCool(RHG_ADin,TG_ADin,Patm,data_name)
% this function aims to fit into parfor
% input:
%       RHG_ADin入口空气湿度
%       TG_ADin入口空气温度
%       Patm surface pressre pa
%       data_name材料isotherm，四列，分别为吸附和解吸
% output:exergy_OP:能耗 kJ/kg

%% Parameters
% create a SatPressure according to Patm
data_Sat = zeros(101,2);
for i = 1:101
    data_Sat(i,1) = i-1;
    
    data_Sat(i,2) = exp(-1.0440397e4/(data_Sat(i,1)*9/5+32+459.67) ...
        -1.129465e1 ...
        -2.7022352e-2*(data_Sat(i,1)*9/5+32+459.67) ...
        +1.289036e-5*(data_Sat(i,1)*9/5+32+459.67)^2 ...
        -2.4780681e-9*(data_Sat(i,1)*9/5+32+459.67)^3 ...
        +6.5459673*log(data_Sat(i,1)*9/5+32+459.67))*6890/1000;
end

data_Iso = load(data_name);
hfg = 2500;

%计算调节量
%吸附剂循环吸附量
WDcycle = 1.5;
%吸附剂物性
MassRatio = 1; Qst = 2700; cps = 1;

%与装置不可逆度有关
% TytG_AD = 0.75;
TytG_DE = 0.8;

TL_AD_min = -5;TL_AD_dlt = 1;
TL_CD_min = -5;TL_CD_dlt = 1;
WDnum = 100;

%%
PG_ADin = exp(-1.0440397e4/(TG_ADin*9/5+32+459.67) ...
    -1.129465e1 ...
    -2.7022352e-2*(TG_ADin*9/5+32+459.67) ...
    +1.289036e-5*(TG_ADin*9/5+32+459.67)^2 ...
    -2.4780681e-9*(TG_ADin*9/5+32+459.67)^3 ...
    +6.5459673*log(TG_ADin*9/5+32+459.67))*6890;
YG_ADin = 621.945 * RHG_ADin/100 / (Patm/PG_ADin - RHG_ADin/100);
exergy_OP = -10;

%迭代最佳吸附冷源温度
for TL_AD = TL_AD_min:TL_AD_dlt:TG_ADin
    exergy_OP1 = -10;
    for YGrate = 8:1:8

        DLYS_DE = 0.4 * YGrate;
        DLYS_ADend = 0.2 * YGrate;
        DLTS_DE = 6 + 0.2 * YGrate;
        DLTS_ADend = 4 + 0.2 * YGrate;
        DLTG_CDout = 4 + 0.2 * YGrate;

        % 吸附结束状态
        TS_ADend = TL_AD + DLTS_ADend;
        YS_ADend = YG_ADin - DLYS_ADend;
        PS_ADend = exp(-1.0440397e4/(TS_ADend*9/5+32+459.67) ...
            -1.129465e1 ...
            -2.7022352e-2*(TS_ADend*9/5+32+459.67) ...
            +1.289036e-5*(TS_ADend*9/5+32+459.67)^2 ...
            -2.4780681e-9*(TS_ADend*9/5+32+459.67)^3 ...
            +6.5459673*log(TS_ADend*9/5+32+459.67))*6890;
        RHS_ADend = Patm * 100 / ...
            (621.945 * PS_ADend / YS_ADend + PS_ADend);
        if RHS_ADend <= 0
            continue;
        end
        
        % 计算WD2
        WD_ADend = interp1(data_Iso(:,1),data_Iso(:,2),RHS_ADend,'linear');
        
        % 解吸计算
        WD_DEend = WD_ADend - WDcycle;
        if WD_DEend <= 0
            continue;
        end
        RHS_DEend = interp1(data_Iso(:,4),data_Iso(:,3),WD_DEend,'linear');

        % 计算吸附剂平均吸附/解吸相对湿度
        RHS_ADaver = 0;RHS_DEaver = 0;
        WD_DEtemp = WD_DEend;
        for k = 1:WDnum
            RHS_ADaver = RHS_ADaver+interp1(data_Iso(:,2),data_Iso(:,1),WD_DEtemp,'linear')/WDnum;
            RHS_DEaver = RHS_DEaver+interp1(data_Iso(:,4),data_Iso(:,3),WD_DEtemp,'linear')/WDnum;
            WD_DEtemp = WD_DEtemp + WDcycle/WDnum;
        end
        

        % 吸附出口温度
        YS_ADaver = 621.945 * RHS_ADaver/(Patm/PS_ADend-RHS_ADaver/100)/100;
        FlowRatio = max((YGrate+DLYS_DE)/(YG_ADin-YS_ADaver),2);
        TytG_AD = 0.8/(FlowRatio^0.2);
        
        if TG_ADin <= TS_ADend
            TG_ADout = TG_ADin;
        else
            TG_ADout = TytG_AD * TS_ADend + (1-TytG_AD) * TG_ADin;
        end
        
        exergy_OP2 = -10;
        
        % 冷凝过程
        for TL_CD = TL_CD_min:TL_CD_dlt:TG_ADin
            % 冷凝出口状态
            TG_CDout = TL_CD + DLTG_CDout;
            PG_CDout = exp(-1.0440397e4/(TG_CDout*9/5+32+459.67) ...
                -1.129465e1 ...
                -2.7022352e-2*(TG_CDout*9/5+32+459.67) ...
                +1.289036e-5*(TG_CDout*9/5+32+459.67)^2 ...
                -2.4780681e-9*(TG_CDout*9/5+32+459.67)^3 ...
                +6.5459673*log(TG_CDout*9/5+32+459.67))*6890;
            
            YG_CDout = 621.945/(Patm/PG_CDout - 1);
            
            % 解吸过程
            YG_DEout = YG_CDout + YGrate;
            YS_DEaver = YG_DEout + DLYS_DE;
            PS_DEaver = Patm * 100/(621.945 * RHS_DEaver / YS_DEaver + RHS_DEaver);
            TS_DE = interp1(data_Sat(:,2),data_Sat(:,1),PS_DEaver/1000,'linear');
            
            PS_DEend = Patm * 100/(621.945 * RHS_DEend / YG_CDout + RHS_DEend);

            if PS_DEend/1000 > data_Sat(end,2)
                continue;
            end
            TS_DEcheck = interp1(data_Sat(:,2),data_Sat(:,1),PS_DEend/1000,'linear');
            
            TL_DE = max(TS_DE,TS_DEcheck) + DLTS_DE;
            
            if TL_DE >= 100
                continue;
            end
            
            %解吸出口温度
            if TL_DE - DLTS_DE >= TG_CDout
                TG_DEout = TytG_DE * (TL_DE - DLTS_DE) + (1-TytG_DE)*TG_CDout;
            else
                TG_DEout = TG_CDout;
            end
            
            DLI1 = TG_DEout * 1.01 + YG_DEout / 1000 * ...
                (Qst + 1.84 * TG_DEout) - TG_CDout * 1.01 - ...
                YG_CDout / 1000 * (Qst + 1.84 * TG_CDout) + ...
                YGrate / WDcycle / 1000 * (TS_DE - TS_ADend) * cps + ...
                YGrate / WDcycle / 1000 * (TL_DE - TL_AD) * cps * MassRatio;
            
            DLI2 = FlowRatio * (TG_ADin * 1.01 + YG_ADin / 1000 * ...
                (Qst + 1.84 * TG_ADin) - TG_ADout * 1.01 - ...
                (YG_ADin - YGrate / FlowRatio) / 1000 * (Qst + 1.84 * TG_ADout)) + ...
                YGrate / WDcycle / 1000 * (TS_DE - TS_ADend) * cps + ...
                YGrate / WDcycle / 1000 * (TL_DE - TL_AD) * cps * MassRatio;
            
            DLI3 = TG_DEout * 1.01 + YG_DEout / 1000 * (hfg + 1.84 * TG_DEout) - ...
                TG_CDout * 1.01 - YG_CDout / 1000 * (hfg + 1.84 * TG_CDout);
            
            exergy = (DLI1 * (TL_DE - TG_ADin) / (TL_DE + 273.15) + ...
                DLI2 * (TG_ADin - TL_AD) / (TL_AD + 273.15) + ...
                DLI3 * (TG_ADin - TL_CD) / (TL_CD + 273.15)) * 1000 / YGrate;
            
            if exergy_OP2 < 0 || exergy <= exergy_OP2
                exergy_OP2 = exergy;
%                 TL_CD_OP2 = TL_CD;
%                 TL_DE_OP2 = TL_DE;
            end
            
        end
        
        if exergy_OP2 > 0 && (exergy_OP1 < 0 || exergy_OP2 <= exergy_OP1)
            exergy_OP1 = exergy_OP2;
%             TL_CD_OP1 = TL_CD_OP2;
%             TL_DE_OP1 = TL_DE_OP2;
%             YGrate_OP1 = YGrate;
%             FlowRatio_OP1 = FlowRatio;
        end
        
    end
    
    if exergy_OP1 > 0 && (exergy_OP < 0 || exergy_OP1 <= exergy_OP)
        exergy_OP = exergy_OP1;
%         TL_CD_OP(i,j) = TL_CD_OP1;
%         TL_DE_OP(i,j) = TL_DE_OP1;
%         YGrate_OP = YGrate_OP1;
%         FlowRatio_OP(i,j) = FlowRatio_OP1;
%         TL_AD_OP(i,j) = TL_AD;
%         YGin_OP = YG_ADin;
    end
    
end

% Output Matrix 修正
if exergy_OP < 0
    exergy_OP = 0;
%     TL_CD_OP(i,j) = 0;
%     TL_DE_OP(i,j) = 0;
%     YGrate_OP = 0;
%     FlowRatio_OP(i,j) = 0;
%     TL_AD_OP(i,j) = 0;
%     YGin_OP = 0;
end








