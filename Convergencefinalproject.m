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
%res = zeros (k50,1);

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
          
        K50 =  k50value(i,j,n50); %Code for Jacobi
        k50 = K50*(n50)^2;
        T_new50(t) = f50value(i,j,n50);
        T_new50(t) = T_new50(t)+ k50*T_050(t-1)+ k50*T_050(t+1) + k50*T_050(t+n50)+ k50*T_050(t-n50);
        T_new50(t) = T_new50(t)/(4*(K50)*(n50)^2);

      end
    end  
end

difference50 = norm(T_new50 - T_050);
T_050 = T_new50;

err50(e50) = log(difference50);
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

n100 = 100;
k100 = n100^2;

%temperature distribution 
T100 = zeros(n100,n100);

%temperature entries cauculation grid by grid
T_100 = zeros (k100,1);
T_new100 = zeros (k100,1);
r = zeros (k100,1);
%res = zeros (k100,1);

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
          
        K100 =  k100value(i,j,n100); %Code for Jacobi
        k100 = K100*(n100)^2;
        T_new100(t) = f100value(i,j,n100);
        T_new100(t) = T_new100(t)+ k100*T_100(t-1)+ k100*T_100(t+1) + k100*T_100(t+n100)+ k100*T_100(t-n100);
        T_new100(t) = T_new100(t)/(4*(K100)*(n100)^2);

      end
    end  
end

difference100 = norm(T_new100 - T_100);
T_100 = T_new100;

err100(e100) = log(difference100);
time100(e100)=toc;
e100=e100+1;
end


n150 = 150;
k150 = n150^2;

%temperature distribution 
T150 = zeros(n150,n150);

%temperature entries cauculation grid by grid
T_150 = zeros (k150,1);
T_new150 = zeros (k150,1);

%norm of the difference
difference150 = 1;
e150=1;

%error tolerance
tol=10e-5;

tic
while difference150 > tol
  
for j = 1:n150 %makes the temperature at the boundaries equal to zero
    for i = 1:n150
      t = i+(j-1)*n150;
      if t <= n150 || rem(t,n150) == 0 || rem(t,n150) == 1 || t > (n150)*((n150)-1)
        T_new150(t) = 0;
      else
          
        K150 =  k150value(i,j,n150); %Code for Jacobi
        k150 = K150*(n150)^2;
        T_new150(t) = f150value(i,j,n150);
        T_new150(t) = T_new150(t)+ k150*T_150(t-1)+ k150*T_150(t+1) + k150*T_150(t+n150)+ k150*T_150(t-n150);
        T_new150(t) = T_new150(t)/(4*(K150)*(n150)^2);      

      end
    end  
end

difference150 = norm(T_new150 - T_150);
T_150 = T_new150;

err150(e150) = log(difference150);
time150(e150)=toc;
e150=e150+1;
end


%Grid Size
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
          
        K200 =  k200value(i,j,n200); %Code for Jacobi
        k200 = K200*(n200)^2;
        T_new200(t) = f200value(i,j,n200);
        T_new200(t) = T_new200(t)+ k200*T_200(t-1)+ k200*T_200(t+1) + k200*T_200(t+n200)+ k200*T_200(t-n200);
        T_new200(t) = T_new200(t)/(4*(K200)*(n200)^2);      

      end
    end  
end

difference200 = norm(T_new200 - T_200);
T_200 = T_new200;

err200(e200) = log(difference200);
time200(e200)=toc;
e200=e200+1;
end


%Grid Size
n250 = 250;
k250 = n250^2;

%temperature distribution 
T250 = zeros(n250,n250);

%temperature entries cauculation grid by grid
T_250 = zeros (k250,1);
T_new250 = zeros (k250,1);

%norm of the difference
difference250 = 1;
e250=1;

%error tolerance
tol=10e-5;

tic
while difference250 > tol
  
for j = 1:n250 %makes the temperature at the boundaries equal to zero
    for i = 1:n250
      t = i+(j-1)*n250;
      if t <= n250 || rem(t,n250) == 0 || rem(t,n250) == 1 || t > (n250)*((n250)-1)
        T_new250(t) = 0;
      else
          
        K250 =  k250value(i,j,n250); %Code for Jacobi
        k250 = K250*(n250)^2;
        T_new250(t) = f250value(i,j,n250);
        T_new250(t) = T_new250(t)+ k250*T_250(t-1)+ k250*T_250(t+1) + k250*T_250(t+n250)+ k250*T_250(t-n250);
        T_new250(t) = T_new250(t)/(4*(K250)*(n250)^2);      

      end
    end  
end

difference250 = norm(T_new250 - T_250);
T_250 = T_new250;

err250(e250) = log(difference250);
time250(e250)=toc;
e250=e250+1;
end


x50=linspace(50,250,3);

Error(1) = norm(T_new250)-norm(T_new50);
Error(2) = norm(T_new250)-norm(T_new100);
%Error(2) = norm(T_new250)-norm(T_new150);
%Error(2) = norm(T_new250)-norm(T_new200);
Error(3) = norm(T_new250)-norm(T_new250);

figure(1)
plot(log10(x50),log10(Error), 'linewidth', 2)
set(gca,'FontSize',20)
xlabel('log of Grid Size')
ylabel('log of Error')
axis tight
grid on
title('Order of Convergence using Jacobi')


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

function [k150] = k150value(i,j,n150)

k150 = 2+ cos(i/n150+ j/n150);

end 


function[F150] = f150value(i,j,n150)

F150= exp(-(i/n150)^2-(j/n150)^2);
  
end

function [k200] = k200value(i,j,n200)

k200 = 2+ cos(i/n200+ j/n200);

end 


function[F200] = f200value(i,j,n200)

F200= exp(-(i/n200)^2-(j/n200)^2);
  
end

function [k250] = k250value(i,j,n250)

k250 = 2+ cos(i/n250+ j/n250);

end 


function[F250] = f250value(i,j,n250)

F250= exp(-(i/n250)^2-(j/n250)^2);
  
end
