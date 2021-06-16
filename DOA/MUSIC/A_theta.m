%%%%%%%%%%%%
% 功能：产生导向矢量
% 输入参数：
%    AtnX,AtnY：阵列的位置
%    NumAtn：阵元的数目
%    Doa：   预设估计的角度
%    w0：载频角频率
%%%%%%%%%%%%
function A = A_theta(AtnX,AtnY,NumAtn,w0,Doa)
    c = 3e8;
    Td = (AtnX'.*cos(Doa*pi/180)+AtnY'.*sin(Doa*pi/180))/c;%时延差
    A = exp(-1j*w0*Td);%导向矢量???不取-
end