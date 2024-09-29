function [OP1,OP2] = RSD(x)
% 已知一组一维数据，求RSD
ave = mean(x);
OP = 0;
for i = 1:length(x)
    OP = OP + (x(i)-ave)^2;
end
OP = sqrt(OP/(length(x) - 1))/ave;

OP1 = ave;
OP2 = OP;
end

