function [fitnessgbest, gbest, zz] = MChOA(X, N, maxgen, lb, ub, dim, fobj)


w = 2.5;
S = 2;     % 空翻因子
aini = 2; afin = 0;         % 控制参数a的初值和终值
epsilon = 6; % 计算适应度值
for i = 1:N
    fitness(i) = fobj(X(i, :));
end
X = GPSinitialization(N, dim, lb, ub);
%% 排序
[bestfitness, bestindex] = sort(fitness);
gbest = X(bestindex(1), :);          % 群体最优极值
fitnessgbest = bestfitness(1);     % 种群最优适应度值
% 计算alpha, beta和delta_pos
Attacker_pos = gbest;
Barrier_pos = X(bestindex(2), :);
Chaser_pos = X(bestindex(3), :);
Driver_pos=X(bestindex(4), :);

%% 迭代
for gen = 1:maxgen
    l=gen;
       f = aini-(aini-afin)*tan(l*pi/(epsilon*maxgen)); 
     C1G1=1.95-((2*l^(1/3))/(maxgen^(1/3)));
    C2G1=(2*l^(1/3))/(maxgen^(1/3))+0.5;
        
    %Group 2
    C1G2= 1.95-((2*l^(1/3))/(maxgen^(1/3)));
    C2G2=(2*(l^3)/(maxgen^3))+0.5;
    
    %Group 3
    C1G3=(-2*(l^3)/(maxgen^3))+2.5;
    C2G3=(2*l^(1/3))/(maxgen^(1/3))+0.5;
    
    %Group 4
    C1G4=(-2*(l^3)/(maxgen^3))+2.5;
    C2G4=(2*(l^3)/(maxgen^3))+0.5;  
    for i = 1:N
        X_old(i, :) = X(i, :);
        for j = 1:dim  
            %% 动态扰动因子策略

            E = randn*((sin(pi*gen/(2*maxgen)))^w+cos(pi*gen/(2*maxgen))-1);
             r11=C1G1*rand(); % r1 is a random number in [0,1]
            r12=C2G1*rand(); % r2 is a random number in [0,1]
            
            r21=C1G2*rand(); % r1 is a random number in [0,1]
            r22=C2G2*rand(); % r2 is a random number in [0,1]
            
            r31=C1G3*rand(); % r1 is a random number in [0,1]
            r32=C2G3*rand(); % r2 is a random number in [0,1]
            
            r41=C1G4*rand(); % r1 is a random number in [0,1]
            r42=C2G4*rand(); % r2 is a random number in [0,1]
            
            A1=2*f*r11-f+E; % Equation (3)
            C1=2*r12; % Equation (4)
           
%% % Please note that to choose various Chaotic maps you should use the related Chaotic maps strategies
            m=chaos(3); % Equation (5)
            D_Attacker=abs(C1*Attacker_pos(j)-m*X(i,j)); % Equation (6)
            X1=Attacker_pos(j)-A1*D_Attacker; % Equation (7)
                       
            A2=2*f*r21-f+E; % Equation (3)
            C2=2*r22; % Equation (4)
            
                   
            D_Barrier=abs(C2*Barrier_pos(j)-m*X(i,j)); % Equation (6)
            X2=Barrier_pos(j)-A2*D_Barrier; % Equation (7)     
            
        
            
            A3=2*f*r31-f+E; % Equation (3)
            C3=2*r32; % Equation (4)
            
            D_Driver=abs(C3*Chaser_pos(j)-m*X(i,j)); % Equation (6)
            X3=Chaser_pos(j)-A3*D_Driver; % Equation (7)      
            
            A4=2*f*r41-f+E; % Equation (3)
            C4=2*r42; % Equation (4)
            
            D_Driver=abs(C4*Driver_pos(j)-m*X(i,j)); % Equation (6)
            X4=Chaser_pos(j)-A4*D_Driver; % Equation (7)
            W1=abs(X1)/(abs(X1)+abs(X2)+abs(X3)+abs(X4));
            W2=abs(X2)/(abs(X1)+abs(X2)+abs(X3)+abs(X4));
            W3=abs(X3)/(abs(X1)+abs(X2)+abs(X3)+abs(X4));
            W4=abs(X4)/(abs(X1)+abs(X2)+abs(X3)+abs(X4));
            
            X(i,j)=((W1*X1+W2*X2+W3*X3+W4*X4)/4)*(1-l/maxgen)+X1*l/maxgen;% Equation (8  
        end
        % 边界处理
      X(i, :) = BoundCheck(X(i, :), lb, ub);
        % 判断
        if fobj(X(i, :)) > fitness(i)
            X(i, :) = X_old(i, :);
        else
            fitness(i) = fobj(X(i, :));
        end
        %% 翻筋斗觅食策略
        if rand < gen/maxgen
            r1 = rand(); r2 = rand();
            X_new(i, :) = X(i, :)+S*(r1*Attacker_pos-r2*X(i, :));
            % 边界处理
               X_new(i, :) = BoundCheck(X_new(i, :), lb, ub);
            % 贪婪选择
            if fobj(X_new(i, :)) < fitness(i)
                X(i, :) = X_new(i, :);
                fitness(i) = fobj(X_new(i, :));
            end
        end
    end
    
    [bestfitness, bestindex] = sort(fitness);
   Attacker_pos = X(bestindex(1), :);
    Attacker_score = bestfitness(1);
    Barrier_pos = X(bestindex(2), :);
    Chaser_pos = X(bestindex(3), :);
    Driver_pos=X(bestindex(4), :);
    
    %% 每一代群体最优值存入zz数组
    zz(gen) = Attacker_score;
    gbest = Attacker_pos;
    fitnessgbest = Attacker_score; 
    
   
end
end
%% 边界处理函数
function s = BoundCheck(s, Lb, Ub)
Flag4ub = s > Ub;
Flag4lb = s <Lb;
s = (s.*(~(Flag4ub+Flag4lb)))+Ub.*Flag4ub+Lb.*Flag4lb;
end