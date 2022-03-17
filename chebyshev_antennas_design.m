function I = chebyshev_antennas_design(N, dlamr, str, value)
% CHEBYSHEV_ANTENNAS_DESIGN 切比雪夫天线电流设计
%   I = CHEBYSHEV_ANTENNAS_DESIGN(N, dlamr, str, value)
%       N       天线单元个数
%       dlamr   天线单元间距与波长之比 (distance and lambda ratio)
%       str     指定已知条件，只能为 SLL 或 PNBW
%       value   str 的值
%   输出 I      归一化的电流
%
%   例如：
%         clear, clc;
%         N = 8;
%         dlamr = 0.5;
%         str = 'SLL';
%         value = -25;
%         I = chebyshev_antennas_current(N, dlamr, str, value)
%   结果为：
%         I =
%               0.37783485957707
%              0.584272242824945
%              0.842415295145896
%                              1
%                              1
%              0.842415295145896
%              0.584272242824945
%               0.37783485957707


    % 根据不同给定条件，使用不同的方法求解系数 a0
    % 一般给定旁瓣电平（SLL），或零功率波瓣宽度（PNBW）
    if strcmp(str, 'SLL') || strcmp(str, 'sll')
        a0 = sll2a0(N, value);
    elseif strcmp(str, 'PNBW') || strcmp(str, 'pnbw')
        a0 = fnbw2a0(N, dlamr, value);
    else
        error('请正确指定 SLL 或 PNBW！');
    end

    % 求出 N-1 阶切比雪夫多项式的系数
    [p, T] = chebyshev(N-1,1);
    T = T(:, end:-2:1);     % 取出需要的部分
    
    % 构造等式右边，使得其满足切比雪夫多项式
    A0 = a0*ones(N,1);
    A0 = vander(A0);
    A0 = A0(1, :)';
    Tp = p .* A0;
    
    % 去除全零行
    index = Tp ~= 0;
    T = T(index, :);
    Tp = Tp(index, :);

    % 判断奇偶性
    isNodd = mod(N,2);
    if isNodd
        T(end) = T(end)/2;
    end

    % 求解方程系数，即电流值
    I = T\Tp;
    
    % 得到每个天线单元阵的电流比例
    I = I/max(I);
    if isNodd
        I = [I; flipud(I(1:end-1))];
    else
        I = [I; flipud(I)];
    end
    
    % 方向图
    theta = linspace(-pi, pi, 360);
    x = cos(2*pi*dlamr * sin(theta) /2);    
    fn = polyval(p, a0*x);
    f = abs(fn);
    f = f/max(f);
%     plot(theta/pi*180, f);       % 直角坐标系的方向图
    figure;
    polar(theta,f);
    view(90, -90);
    title([num2str(N), ' 元切比雪夫天线的方向图']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% a0 的求解
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a0 = sll2a0(N,SLL)
    R0 = 10^(-SLL/20);
    
    Rt = sqrt(R0^2 - 1);
    a0 = 0.5 * ( power(R0 + Rt, 1/(N-1)) + power(R0 + Rt, -1/(N-1)) );
end
function a0 = fnbw2a0(N, dlamr, theta01)
    x01 = cos(pi*dlamr * sind(theta01/2));
    a0 = cos(pi/2/(N-1))./x01;
end
