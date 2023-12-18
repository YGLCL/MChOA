%% logistic ªÏ„Á”≥…‰
function O=chaos(max_iter)

O=zeros(1,max_iter);
x(1)=0.7;
%Logistic map
a=4;
for i=1:max_iter
    x(i+1)=a*x(i)*(1-x(i));
    G(i)=x(i);
end
O = G(max_iter);
end
