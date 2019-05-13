clear all % Clears all variables.
clc % Clears the command window.
close all % Closes all figures/plots.
format long 

%grid size
nj100 = 100;
kj100 = nj100^2;

%temperature distribution 
TJ100 = zeros(nj100,nj100);

%temperature entries cauculation grid by grid
Tj_100 = zeros (kj100,1);
Tj_new100 = zeros (kj100,1);

%norm of the difference
differencej100 = 1;
ej100=1;

%error tolerance
tol=10e-5;

tic
while differencej100 > tol
  
for j = 1:nj100 %makes the temperature at the boundaries equal to zero
    for i = 1:nj100
      t = i+(j-1)*nj100;
      if t <= nj100 || rem(t,nj100) == 0 || rem(t,nj100) == 1 || t > (nj100)*((nj100)-1)
        Tj_new100(t) = 0;
      else
          
        Kj100 =  kj100value(i,j,nj100); %Code for Jacobi
        kj100 = Kj100*(nj100)^2;
        Tj_new100(t) = fj100value(i,j,nj100);
        Tj_new100(t) = Tj_new100(t)+ kj100*Tj_100(t-1)+ kj100*Tj_100(t+1) + kj100*Tj_100(t+nj100)+ kj100*Tj_100(t-nj100);
        Tj_new100(t) = Tj_new100(t)/(4*(Kj100)*(nj100)^2);      

      end
    end  
end

differencej100 = norm(Tj_new100 - Tj_100);
Tj_100 = Tj_new100;

errj100(ej100) = log10(differencej100);
timej100(ej100)=toc;
ej100=ej100+1;
end


t = 1;
for j = 1:nj100
    for i = 1:nj100
        TJ100 (i,j) = Tj_new100(t,1);
        t = t + 1;
    end
end

j=1;

for t = (((nj100^2)/2)+1) :(((nj100^2)/2)+nj100)
    Xj100(j)= Tj_new100(t,1);
    j=j+1;
end


%grid size
ngs100 = 100;
kgs100 = ngs100^2;

%temperature distribution 
Tgs100 = zeros(ngs100,ngs100);

%temperature entries cauculation grid by grid
Tgs_100 = zeros (kgs100,1);
Tgs_new100 = zeros (kgs100,1);

%norm of the difference
differencegs100 = 1;
egs100=1;

%error tolerance
tol=10e-5;

tic
while differencegs100 > tol
  
for j = 1:ngs100 %makes the temperature at the boundaries equal to zero
    for i = 1:ngs100
      t = i+(j-1)*ngs100;
      if t <= ngs100 || rem(t,ngs100) == 0 || rem(t,ngs100) == 1 || t > (ngs100)*((ngs100)-1)
        Tgs_new100(t) = 0;
      else
          
        Kgs100 =  kgs100value(i,j,ngs100); %Code for Gauss Seidel
        kgs100 = Kgs100*(ngs100)^2;
        Tgs_new100(t) = fgs100value(i,j,ngs100);
        Tgs_new100(t) = Tgs_new100(t)+ kgs100*Tgs_new100(t-1)+ kgs100*Tgs_new100(t+1) + kgs100*Tgs_new100(t+ngs100)+ kgs100*Tgs_new100(t-ngs100);
        Tgs_new100(t) = Tgs_new100(t)/(4*(Kgs100)*(ngs100)^2);      

      end
    end  
end

differencegs100 = norm(Tgs_new100 - Tgs_100);
Tgs_100 = Tgs_new100;

errgs100(egs100) = log10(differencegs100);
timegs100(egs100)=toc;
egs100=egs100+1;
end


t = 1;
for j = 1:ngs100
    for i = 1:ngs100
        Tgs100 (i,j) = Tgs_new100(t,1);
        t = t + 1;
    end
end

j=1;

for t = (((ngs100^2)/2)+1) :(((ngs100^2)/2)+ngs100)
    Xgs100(j)= Tgs_new100(t,1);
    j=j+1;
end

%grid size
nsor100 = 100;
ksor100 = nsor100^2;
w100 = (4/(2+((sqrt(4-((2*cos(pi/nsor100))^2))))));


%temperature distribution 
Tsor100 = zeros(nsor100,nsor100);

%temperature entries cauculation grid by grid
Tsor_100 = zeros (ksor100,1);
Tsor_new100 = zeros (ksor100,1);

%norm of the difference
differencesor100 = 1;
esor100=1;

%error tolerance
tol=10e-5;

tic
while differencesor100 > tol
  
for j = 1:nsor100 %makes the temperature at the boundaries equal to zero
    for i = 1:nsor100
      t = i+(j-1)*nsor100;
      if t <= nsor100 || rem(t,nsor100) == 0 || rem(t,nsor100) == 1 || t > (nsor100)*((nsor100)-1)
        Tsor_new100(t) = 0;
      else
          
        Ksor100 =  ksor100value(i,j,nsor100); %Code for SOR
        ksor100 = (w100)*Ksor100*(nsor100)^2;
        Tsor_new100(t) = (w100)*fsor100value(i,j,nsor100);
        Tsor_new100(t) = Tsor_new100(t)+ ksor100*Tsor_new100(t-1)+ ksor100*Tsor_new100(t+1) + ksor100*Tsor_new100(t+nsor100)+ ksor100*Tsor_new100(t-nsor100);
        Tsor_new100(t) = Tsor_new100(t)/(4*(w100)*(Ksor100)*(nsor100)^2);
        Tsor_new100(t) = (w100)*Tsor_new100(t) + (1-w100)*Tsor_100(t);

      end
    end  
end

differencesor100 = norm(Tsor_new100 - Tsor_100);
Tsor_100 = Tsor_new100;

errsor100(esor100) = log10(differencesor100);
timesor100(esor100)=toc;
esor100=esor100+1;
end


t = 1;
for j = 1:nsor100
    for i = 1:nsor100
        Tsor100 (i,j) = Tsor_new100(t,1);
        t = t + 1;
    end
end

j=1;

for t = (((nsor100^2)/2)+1) :(((nsor100^2)/2)+nsor100)
    Xsor100(j)= Tsor_new100(t,1);
    j=j+1;
end


figure(5)
hold on
plot(errj100, 'Linewidth',2)
plot(errgs100, 'Linewidth',2)
plot(errsor100, 'Linewidth',2)
set(gca,'FontSize',20)
xlabel('Iterations')
ylabel('log (Residual)')
legend({'Jacobi','Gauss-Seidel', 'SOR'},'Location','northeast')
% axis tight
grid on
hold off
title('Log of Residual vs. Iterations for 100x100 grid using all methods')


figure(6)
hold on
plot(timej100,errj100, 'Linewidth',2)
plot(timegs100,errgs100, 'Linewidth',2)
plot(timesor100,errsor100, 'Linewidth',2)
set(gca,'FontSize',20)
legend({'Jacobi','Gauss-Seidel', 'SOR'},'Location','northeast')
xlabel('Time (sec)')
ylabel('log (Residual)')
% axis tight
grid on
hold off
title('Log of Residual vs. CPU Time for 100x100 grid using all methods')


function[Fj100] = fj100value(i,j,nj100)

Fj100= exp(-(i/nj100)^2-(j/nj100)^2);
  
end

function [kj100] = kj100value(i,j,nj100)

kj100 = 2+ cos(i/nj100+ j/nj100);

end 

function[Fgs100] = fgs100value(i,j,ngs100)

Fgs100= exp(-(i/ngs100)^2-(j/ngs100)^2);
  
end

function [kgs100] = kgs100value(i,j,ngs100)

kgs100 = 2+ cos(i/ngs100+ j/ngs100);

end 

function[Fsor100] = fsor100value(i,j,nsor100)

Fsor100= exp(-(i/nsor100)^2-(j/nsor100)^2);
  
end

function [ksor100] = ksor100value(i,j,nsor100)

ksor100 = 2+ cos(i/nsor100+ j/nsor100);

end 