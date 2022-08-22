clear;clc;

debug = false;
%% 读入数据并画图
if debug % 正常模式：从文件读入散点
    x = load(); % 读入散点，n*2矩阵
else % 否则自动生成散点
    theta = (-pi:0.1:pi)';
    a0 = 1.5; b0 = 1; cx0 = 1; cy0 = 0.5; phi00 = pi/6;
    x(:,1) = a0*cos(theta) + 0.1*rand(length(theta),1) + cx0;
    x(:,2) = b0*sin(theta) + 0.1*rand(length(theta),1) + cy0;
    R = [cos(phi00), -sin(phi00); sin(phi00), cos(phi00)];
    x = (R*x')';
end

xmin = min(x(:,1));
xmax = max(x(:,1));
ymin = min(x(:,2));
ymax = max(x(:,2));
n = size(x,1);

figure('Name', '椭圆拟合并求面积');
plot(x(:,1), x(:,2), 'ro'); hold on;
xtitle('x');
ytitle('y');

%% 椭圆拟合：
% p = fitellipse(x(:,1),x(:,2));
% Fplot = @(x) p(1)*x(:,1).^2 + p(2)*x(:,1).*x(:,2) + p(3)*x(:,2).^2 + p(4)*x(:,1) + p(5)*x(:,2) + p(6); % 椭圆一般方程
% % fimplicit(Fplot);
% ezplot(@(x, y) Fplot(p, [x, y]), [-1 + xmin, 1 + xmax, -1  +ymin, 1 + ymax]);

%% 涉及的函数
% 输出的是列向量，表示椭圆一般方程的各个参数
function W = fitellipse(x,y)
% 构造矩阵
D = [x.*x, x.*y, y.*y, x, y, ones(size(x))];
S = D'*D;
G = zeros(6);
G(1,3) = 2; G(3,1) = 2; G(2,2) = -1;

% 求解
[vec, val] = eig(S\G);
[~, idx] = find(val>0&~isinf(val));
W = vec(:,idx);
W = sqrt(1/(W'*S*W))*W;
end