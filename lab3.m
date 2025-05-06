%% 常微分方程数值解法的MATLAB实现
% 求解初值问题: y'=(2x)/(3y^2)，y(0)=1，在区间[0,1]上取h=0.1
% 定义微分方程右端函数 f(x, y)
fprintf('于智同23013214@%s\n',datetime);
f = @(x, y) (2*x)/(3*y^2);
% 定义精确解
exact = @(x) (1 + x.^2).^(1/3);
% 参数设置
x0 = 0;      % 初始x值
y0 = 1;      % 初始y值 (y(0) = 1)
h = 0.1;     % 步长
xn = 1;      % 终止x值
n = (xn-x0)/h;  % 迭代次数
% 创建x网格
x = x0:h:xn;

%% 1. 欧拉方法
y_euler = zeros(1, n+1);
y_euler(1) = y0;

for i = 1:n
    y_euler(i+1) = y_euler(i) + h * f(x(i), y_euler(i));
end

%% 2. 改进的欧拉方法（也称为梯形法或Heun法）
y_improved = zeros(1, n+1);
y_improved(1) = y0;

for i = 1:n
    k1 = f(x(i), y_improved(i));
    k2 = f(x(i+1), y_improved(i) + h*k1);
    y_improved(i+1) = y_improved(i) + h/2 * (k1 + k2);
end

%% 3. 四阶龙格-库塔方法
y_rk4 = zeros(1, n+1);
y_rk4(1) = y0;

for i = 1:n
    k1 = f(x(i), y_rk4(i));
    k2 = f(x(i) + h/2, y_rk4(i) + h/2*k1);
    k3 = f(x(i) + h/2, y_rk4(i) + h/2*k2);
    k4 = f(x(i) + h, y_rk4(i) + h*k3);
    y_rk4(i+1) = y_rk4(i) + h/6 * (k1 + 2*k2 + 2*k3 + k4);
end

%% 4. 四阶亚当姆斯方法（也称为亚当姆斯-巴什福斯法）
y_adams = zeros(1, n+1);
y_adams(1) = y0;
% 使用四阶龙格-库塔方法计算前四个值
for i = 1:3
    k1 = f(x(i), y_adams(i));
    k2 = f(x(i) + h/2, y_adams(i) + h/2*k1);
    k3 = f(x(i) + h/2, y_adams(i) + h/2*k2);
    k4 = f(x(i) + h, y_adams(i) + h*k3);
    y_adams(i+1) = y_adams(i) + h/6 * (k1 + 2*k2 + 2*k3 + k4);
end

% 使用四阶亚当姆斯方法计算后续值
for i = 4:n
    y_adams(i+1) = y_adams(i) + h/24 * (55*f(x(i), y_adams(i)) - 59*f(x(i-1), y_adams(i-1)) + 37*f(x(i-2), y_adams(i-2)) - 9*f(x(i-3), y_adams(i-3)));
end

%% 5. 改进的四阶亚当姆斯预估校正系统
y_adams_pc = zeros(1, n+1);
y_adams_pc(1) = y0;
% 使用四阶龙格-库塔方法计算前四个值
for i = 1:3
    k1 = f(x(i), y_adams_pc(i));
    k2 = f(x(i) + h/2, y_adams_pc(i) + h/2*k1);
    k3 = f(x(i) + h/2, y_adams_pc(i) + h/2*k2);
    k4 = f(x(i) + h, y_adams_pc(i) + h*k3);
    y_adams_pc(i+1) = y_adams_pc(i) + h/6 * (k1 + 2*k2 + 2*k3 + k4);
end

% 使用四阶亚当姆斯预估校正系统计算后续值
for i = 4:n
    % 预估（使用亚当姆斯-巴什福斯法）
    y_pred = y_adams_pc(i) + h/24 * (55*f(x(i), y_adams_pc(i)) - 59*f(x(i-1), y_adams_pc(i-1)) + 37*f(x(i-2), y_adams_pc(i-2)) - 9*f(x(i-3), y_adams_pc(i-3)));
    
    % 校正（使用亚当姆斯-莫尔顿法）
    y_adams_pc(i+1) = y_adams_pc(i) + h/24 * (9*f(x(i+1), y_pred) + 19*f(x(i), y_adams_pc(i)) - 5*f(x(i-1), y_adams_pc(i-1)) + f(x(i-2), y_adams_pc(i-2)));
end

%% 计算精确解并比较误差
y_exact = exact(x);

% 计算每种方法的误差
error_euler = abs(y_euler - y_exact);
error_improved = abs(y_improved - y_exact);
error_rk4 = abs(y_rk4 - y_exact);
error_adams = abs(y_adams - y_exact);
error_adams_pc = abs(y_adams_pc - y_exact);

%% 结果显示
% 绘制各方法的数值解与精确解的比较图
figure(1);
plot(x, y_exact, 'k-', 'LineWidth', 2, 'DisplayName', '精确解');
hold on;
plot(x, y_euler, 'r--o', 'DisplayName', '欧拉方法');
plot(x, y_improved, 'g--*', 'DisplayName', '改进的欧拉方法');
plot(x, y_rk4, 'b--d', 'DisplayName', '四阶龙格-库塔方法');
plot(x, y_adams, 'c--s', 'DisplayName', '四阶亚当姆斯方法');
plot(x, y_adams_pc, 'm--^', 'DisplayName', '改进的四阶亚当姆斯方法');
hold off;
xlabel('x');
ylabel('y');
title('各种数值方法与精确解的比较');
legend('Location', 'best');
grid on;

% 绘制各方法的误差
figure(2);
semilogy(x, error_euler, 'r-o', 'DisplayName', '欧拉方法');
hold on;
semilogy(x, error_improved, 'g-*', 'DisplayName', '改进的欧拉方法');
semilogy(x, error_rk4, 'b-d', 'DisplayName', '四阶龙格-库塔方法');
semilogy(x, error_adams, 'c-s', 'DisplayName', '四阶亚当姆斯方法');
semilogy(x, error_adams_pc, 'm-^', 'DisplayName', '改进的四阶亚当姆斯方法');
hold off;
xlabel('x');
ylabel('绝对误差 (对数尺度)');
title('各种数值方法的误差比较');
legend('Location', 'best');
grid on;

% 打印结果表格
fprintf('数值解与精确解的比较:\n');
fprintf('%-2s %-10s %-10s %-10s %-10s %-10s %-10s\n', 'x', '精确解', '欧拉法', '改进欧拉法', 'RK4法', '亚当姆斯法', '改进亚当姆斯法');
for i = 1:n+1
    fprintf('%.1f %-10.6f %-10.6f %-15.6f %-15.6f %-15.6f %-15.6f\n', x(i), y_exact(i), y_euler(i), y_improved(i), y_rk4(i), y_adams(i), y_adams_pc(i));
end

fprintf('\n最大绝对误差:\n');
fprintf('欧拉方法: %.6e\n', max(error_euler));
fprintf('改进的欧拉方法: %.6e\n', max(error_improved));
fprintf('四阶龙格-库塔方法: %.6e\n', max(error_rk4));
fprintf('四阶亚当姆斯方法: %.6e\n', max(error_adams));
fprintf('改进的四阶亚当姆斯方法: %.6e\n', max(error_adams_pc));