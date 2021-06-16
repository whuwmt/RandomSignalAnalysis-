%%%%%%%%%%%%%%%
% 功能：传统的MUSIC估计
% 参数输入：
%     flag：    为1表示使用特征分解，为2表示使用SVD奇异值分解
%     X：       输入快扫数据矢量
%    NumSource：信源数
%    NumAtn：   阵元数=天线数=通道数
%    Snap：     快拍数   
%    AtnX,AtnY：天线位置
%    w0：       载频角频率
% 参数输出：
%    P：        谱估计
%%%%%%%%%%%%%%%
function [P,theta] = MUSIC_Tranditon(X,NumSource,NumAtn,Snap,AtnX,AtnY,flag,w0,Research)
    c = 3e8;
    % 最大似然估计
    Rxx = X*X'/Snap;
    
    % 提取信号子空间与噪声子空间，flag=1为特征值分解，flag=2为SVD分解
    if flag == 1
        [V,D] = eig(Rxx);
        EVA = diag(D); %提取D中的对角元素，即特征值
        [EVA , I] = sort(EVA);%对特征值从小到大排序
        Us = V(:,I(NumAtn-NumSource+1:NumAtn,:));
        Un = V(:,I(1:NumAtn-NumSource,:));
    else if flag == 2
            [S,V,D] = svd(Rxx);
            EVA = diag(V); %提取V中的对角元素，即特征值
            [EVA , I] = sort(EVA);%对特征值从小到大排序
            Us = S(:,I(NumAtn-NumSource+1:NumAtn,:));
            Un = S(:,I(1:NumAtn-NumSource,:));
        end
    end
    P = [];
    theta = [Research(1):Research(3):Research(2)]*pi/180;
    for ii = 1:length(theta)
        td = (AtnX'.*cos(theta(ii))+AtnY'.*sin(theta(ii)))/c;%时延差
        A_Research = exp(-1j*w0*td);
        P(ii) = 1/abs(A_Research'*Un*Un'*A_Research);
    end
end