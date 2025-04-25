fprintf('于智同23013214@%s\n',datetime);
f = @(x) 4./(1+x.^2);
real_int_f = integral(f,0,1);
fprintf('f的实际值:%.14f\n',real_int_f);
plot_x = linspace(0,1,100);
plot(plot_x,f(plot_x));
function y = T_n(f,a,b,n)
    y = 0;
    h = (b-a)/n;
    for i = 0:n-1
        y = y + h/2*(f(a+i*h)+f(a+(i+1)*h));
    end
    fprintf('复化求积公式:%.14f\n误差%.14f\n',y,pi-y);
end
function y = S_n(f,a,b,n)
    y = 0;
    h = (b-a)/n;
    for i = 0:n-1
        y = y+ h/6*(f(a+i*h)+4*f(a+(i+0.5)*h)+f(a+(i+1)*h));
    end
    fprintf('Simpson求积公式:%.14f\n误差%.14f\n',y,pi-y);
end
function y = Romberg(f,a,b,epsilon)
    %初始化
    h = b-a;n=1;
    R=zeros(20,4);
    %生成表格前五行
    R(1,1)=h/2*(f(a)+f(b));
    for i = 2:10
        tmp = 0;
        for j = 0:n-1
            tmp = tmp+f(a+(j+0.5)*h);
        end
        R(i,1) = 0.5*R(i-1,1)+h/2*tmp;
        h = h/2;n = n*2;
    end
    for i = 2:10
        R(i,2)= 4/3*R(i,1)-1/3*R(i-1,1);
    end
    for i = 3:10
        R(i,3) = 16/15*R(i,2)-1/15*R(i-1,2);
    end
    for i = 4:10
        R(i,4) = 64/63*R(i,2)-1/63*R(i-1,3);
    end
    for i = 5:10
        y = R(i,4);
        diff = abs(R(i,4)-R(i,4));
        if diff<epsilon
            fprintf('Romberg加速算法:%.14f\n误差%.14f\n',y,pi-y);
            return;
        end
    end
end
function y=Gauss(f,a,b)
    m = (b-a)/2;
    n = (b+a)/2;
    y = m*(5/9*f(n-sqrt(3/5)*m)+8/9*f(n)+5/9*f(n+sqrt(3/5)*m));
    fprintf('三点Gauss算法:%.14f\n误差%.14f\n',y,pi-y);
end
y = T_n(f,0,1,16);
y = S_n(f,0,1,32);
y = Romberg(f,0,1,0.5*10^(-7));
y = Gauss(f,0,1);
