clear all % Clears all variables.
clc % Clears the command window.
close all % Closes all figures/plots.
format long 

%grid size
n50 = 50;
k50 = n50^2;

%temperature distribution 
T50 = zeros(n50,n50);

%temperature entries cauculation grid by grid
T_050 = zeros (k50,1);
T_new50 = zeros (k50,1);

%norm of the difference
difference50 = 1;
e50=1;

%error tolerance
tol=10e-5;

tic
while difference50 > tol
  
for j = 1:n50 %makes the temperature at the boundaries equal to zero
    for i = 1:n50
      t = i+(j-1)*n50;
      if t <= n50 || rem(t,n50) == 0 || rem(t,n50) == 1 || t > (n50)*((n50)-1)
        T_new50(t) = 0;
      else
        
        K50 =  k50value(i,j,n50); %Code for Newtons Method
        k50 = K50*(n50)^2;
        T_new50(t) = f50value(i,j,n50)-(T_new50(t))^4;
        T_new50(t) = T_new50(t)+ k50*T_050(t-1)+ k50*T_050(t+1) + k50*T_050(t+n50)+ k50*T_050(t-n50);
        if i==j
            T_new50(t) = T_new50(t)/((4+4*(T_new50(t))^3)*(K50)*(n50)^2);
        else
        T_new50(t) = T_new50(t)/(4*(K50)*(n50)^2);  
        end

      end
    end  
end

difference50 = norm(T_new50 - T_050);
T_050 = T_new50;

err50(e50) = log10(difference50);
time50(e50)=toc;
e50=e50+1;
end


t = 1;
for j = 1:n50
    for i = 1:n50
        T50 (i,j) = T_new50(t,1);
        t = t + 1;
    end
end

j=1;

for t = (((n50^2)/2)+1) :(((n50^2)/2)+n50)
    X50(j)= T_new50(t,1);
    j=j+1;
end

figure(1)
surf(T50)
%shading interp
colormap(jet)  
xlabel('x')
ylabel('y')
zlabel('Temperature')
axis tight
title(['Temperature Field using Newtons Method for ' num2str(n50) 'x' num2str(n50)])



%grid size
n100 = 100;
k100 = n100^2;

%temperature distribution 
T100 = zeros(n100,n100);

%temperature entries cauculation grid by grid
T_100 = zeros (k100,1);
T_new100 = zeros (k100,1);

%norm of the difference
difference100 = 1;
e100=1;

%error tolerance
tol=10e-5;

tic
while difference100 > tol
  
for j = 1:n100 %makes the temperature at the boundaries equal to zero
    for i = 1:n100
      t = i+(j-1)*n100;
      if t <= n100 || rem(t,n100) == 0 || rem(t,n100) == 1 || t > (n100)*((n100)-1)
        T_new100(t) = 0;
        
      else
    
        K100 =  k100value(i,j,n100); %Code for Newtons Method
        k100 = K100*(n100)^2;
        T_new100(t) = f100value(i,j,n100)-(T_new100(t))^4;
        T_new100(t) = T_new100(t)+ k100*T_100(t-1)+ k100*T_100(t+1) + k100*T_100(t+n100)+ k100*T_100(t-n100);
        if i==j
            T_new100(t) = T_new100(t)/((4+4*(T_new100(t))^3)*(K100)*(n100)^2);
        else
        T_new100(t) = T_new100(t)/(4*(K100)*(n100)^2);
        end

      end
    end  
end

difference100 = norm(T_new100 - T_100);
T_100 = T_new100;

err100(e100) = log10(difference100);
time100(e100)=toc;
e100=e100+1;
end


t = 1;
for j = 1:n100
    for i = 1:n100
        T100 (i,j) = T_new100(t,1);
        t = t + 1;
    end
end

j=1;

for t = (((n100^2)/2)+1) :(((n100^2)/2)+n100)
    X100(j)= T_new100(t,1);
    j=j+1;
end

figure(2)
surf(T100)
%shading interp
colormap(jet)  
xlabel('x')
ylabel('y')
zlabel('Temperature')
axis tight
title(['Temperature Field using Newtons Method for ' num2str(n100) 'x' num2str(n100)])



%grid size
n200 = 200;
k200 = n200^2;

%temperature distribution 
T200 = zeros(n200,n200);

%temperature entries cauculation grid by grid
T_200 = zeros (k200,1);
T_new200 = zeros (k200,1);

%norm of the difference
difference200 = 1;
e200=1;

%error tolerance
tol=10e-5;

tic
while difference200 > tol
  
for j = 1:n200 %makes the temperature at the boundaries equal to zero
    for i = 1:n200
      t = i+(j-1)*n200;
      if t <= n200 || rem(t,n200) == 0 || rem(t,n200) == 1 || t > (n200)*((n200)-1)
        T_new200(t) = 0;
      else
         
        K200 =  k200value(i,j,n200); %Code for Newtons Method
        k200 = K200*(n200)^2;
        T_new200(t) = f200value(i,j,n200)-(T_new200(t))^4;
        T_new200(t) = T_new200(t)+ k200*T_200(t-1)+ k200*T_200(t+1) + k200*T_200(t+n200)+ k200*T_200(t-n200);
        if i==j
            T_new200(t) = T_new200(t)/((4+4*(T_new200(t))^3)*(K200)*(n200)^2);
        else
        T_new200(t) = T_new200(t)/(4*(K200)*(n200)^2); 
        end

      end
    end  
end

difference200 = norm(T_new200 - T_200);
T_200 = T_new200;

err200(e200) = log10(difference200);
time200(e200)=toc;
e200=e200+1;
end


t = 1;
for j = 1:n200
    for i = 1:n200
        T200 (i,j) = T_new200(t,1);
        t = t + 1;
    end
end

j=1;

for t = (((n200^2)/2)+1) :(((n200^2)/2)+n200)
    X200(j)= T_new200(t,1);
    j=j+1;
end

figure(3)
surf(T200)
%shading interp
colormap(jet)  
xlabel('x')
ylabel('y')
zlabel('Temperature')
axis tight
title(['Temperature Field using Newtons Method for ' num2str(n200) 'x' num2str(n200)])

p50 = linspace(0,1,50);
p100 = linspace(0,1,100);
p200 = linspace(0,1,200);

figure(4)
hold on
plot(p50, X50)
plot(p100, X100)
plot(p200, X200)
xlabel('x')
ylabel('Temperature')
legend({'50x50','100x100', '200x200'},'Location','northeast', 'Fontsize',15)
%axis tight
grid on
hold off
title('Centre-line Temperature using Newtons Method for all cases', 'Fontsize',15)


figure(5)
hold on
plot(err50, 'Linewidth',2)
plot(err100, 'Linewidth',2)
plot(err200, 'Linewidth',2)
xlabel('Iterations')
ylabel('Temperature Difference')
%axis tight
set(gca,'FontSize',20)
legend({'50x50','100x100', '200x200'},'Location','northeast')
grid on
hold off
title('Temperature Difference vs. Iterations using Newtons Method for all cases')


figure(6)
hold on
plot(time50,err50, 'Linewidth',2)
plot(time100,err100,'Linewidth',2)
plot(time200,err200, 'Linewidth',2)
xlabel('Time (sec)')
ylabel('Temperature Difference')
%axis tight
legend({'50x50','100x100', '200x200'},'Location','northeast')
set(gca,'FontSize',20)
grid on
hold off
title('Temperature Difference vs. CPU Time using Newtons Method for all cases')


function [k50] = k50value(i,j,n50)

k50 = 2+ cos(i/n50+ j/n50);

end 


function[F50] = f50value(i,j,n50)

F50= exp(-(i/n50)^2-(j/n50)^2);
  
end

function [k100] = k100value(i,j,n100)

k100 = 2+ cos(i/n100+ j/n100);

end 


function[F100] = f100value(i,j,n100)

F100= exp(-(i/n100)^2-(j/n100)^2);
  
end

function [k200] = k200value(i,j,n200)

k200 = 2+ cos(i/n200+ j/n200);

end 


function[F200] = f200value(i,j,n200)

F200= exp(-(i/n200)^2-(j/n200)^2);
  
end