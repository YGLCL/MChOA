%% �����������
clear
clc

%% ��������
Max_iteration = 500;        % ����������
N = 30;         % ��Ⱥ��ģ
Function_name = 'F1';  % ��F1��F23�Ĳ��Ժ���������(Please adjust epsilon to 2 when using F8.)
% ������ѡ��׼��������ϸ��Ϣ
[lb, ub, dim, fobj] = Get_Functions_details(Function_name);

%% �Ż�
X = initialization(N, dim, ub, lb);
[Best_score, Best_pos, DSF_ChoA_Curve] = MChOA(X, N, Max_iteration, lb, ub, dim, fobj);
mean_MChOA = mean(Best_score);

%% ��ͼ
% 1��������ѡ��׼��������ά����ͼ��
figure;
func_plot(Function_name);
title(Function_name)
xlabel('x_1');
ylabel('x_2');
zlabel([Function_name,'( x_1 , x_2 )'])

% 2������Ŀ�꺯��ֵ�仯����ͼ
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

%% ��ʾ���
display(['The optimal solution of ', Function_name, ' is: ', num2str(Best_pos)]);
display(['The optimal value of ', Function_name, ' is : ', num2str(Best_score)]);

