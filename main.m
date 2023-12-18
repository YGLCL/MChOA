%% 清除环境变量
clear
clc

%% 参数设置
Max_iteration = 500;        % 最大迭代次数
N = 30;         % 种群规模
Function_name = 'F1';  % 从F1到F23的测试函数的名称(Please adjust epsilon to 2 when using F8.)
% 加载所选基准函数的详细信息
[lb, ub, dim, fobj] = Get_Functions_details(Function_name);

%% 优化
X = initialization(N, dim, ub, lb);
[Best_score, Best_pos, DSF_ChoA_Curve] = MChOA(X, N, Max_iteration, lb, ub, dim, fobj);
mean_MChOA = mean(Best_score);

%% 绘图
% 1、画出所选基准函数的三维立体图形
figure;
func_plot(Function_name);
title(Function_name)
xlabel('x_1');
ylabel('x_2');
zlabel([Function_name,'( x_1 , x_2 )'])

% 2、画出目标函数值变化曲线图
figure;
t = 1:Max_iteration;
semilogy(t, DSF_ChoA_Curve, 'r-', 'linewidth', 2); 
title(Function_name); 
xlabel('Iteration'); 
ylabel('Best score');
axis fill 
grid on
box on
legend('DSF-ChoA');

%% 显示结果
display(['The optimal solution of ', Function_name, ' is: ', num2str(Best_pos)]);
display(['The optimal value of ', Function_name, ' is : ', num2str(Best_score)]);

