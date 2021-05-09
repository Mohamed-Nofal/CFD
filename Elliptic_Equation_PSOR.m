%% Solve Elliptic Equation by PSOR 

clc;clear all;close all

%% Given or Arbitrary data

nmax=2000;   imax=51;     jmax=51;  Lx=1;        Ly=1;
a_omega=[0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0];

%% Calculated data

dx          = Lx/(imax-1);
dy          = Ly/(jmax-1);
dx_Over_dy  = dx/dy;

for ii=1:8
    omega=a_omega(ii);
    
    % initial values   
    for i=1:imax
        for j=1:jmax
            u(i,j)=0 ;
        end
    end
    %Boundary values
    for j=1:jmax
        u(1,j)=-(dy*(j-1))^2;
        u(imax,j)=(dx*(imax-1))^2-(dy*(j-1))^2;
        up(1,j)=u(1,j);
        up(imax,j)=u(imax,j);
    end
    for i=1:imax
        u(i,1)=(dx*(i-1))^2;
        u(i,jmax)=(dx*(i-1))^2-(dy*(jmax-1))^2;
        up(i,1)=u(i,1);
        up(i,jmax)=u(i,jmax);
    end
    % Iteration procedure
    for n=1:nmax
        % General points calculations
        for j=2:jmax-1
            for i=2:imax-1
                up(i,j)=.5*(up(i-1,j)+u(i+1,j)+dx_Over_dy^2*(up(i,j-1)+u(i,j+1)))/(1+dx_Over_dy^2);
                up(i,j)=u(i,j)+omega*(up(i,j)-u(i,j));
            end
        end
       	% Error calculations
        rmser=0   ;  maxer=0;   imaxer=1;    jmaxer=1;
        for j=2:jmax-1
            for i=2:imax-1
                uexat=(dx*(i-1))^2-(dy*(j-1))^2;
                er=abs(up(i,j)-uexat);
                if(er>maxer)
                    maxer=er;
                    imaxer=i;
                    jmaxer=j ;
                end
                rmser=rmser+er*er;
                u(i,j)=up(i,j);
            end
        end
        rmser=sqrt(rmser)/((imax-2)*(jmax-2));
        RMS_ER(n)=log10(rmser);
        Max_ER(n)=log10(maxer);
        an(n)=n;
    end
    
    if omega==0.25
        RMS_ER_025=RMS_ER;
        Max_ER_025=Max_ER;
    elseif omega==0.5
        RMS_ER_050=RMS_ER;
        Max_ER_050=Max_ER;
    elseif omega==0.75
        RMS_ER_075=RMS_ER;
        Max_ER_075=Max_ER;
    elseif omega==1.0
        RMS_ER_100=RMS_ER;
        Max_ER_100=Max_ER;
    elseif omega==1.25
        RMS_ER_125=RMS_ER;
        Max_ER_125=Max_ER;
    elseif omega==1.5
        RMS_ER_150=RMS_ER;
        Max_ER_150=Max_ER;
    elseif omega==1.75
        RMS_ER_175=RMS_ER;
        Max_ER_175=Max_ER;
    elseif omega==2.0
        RMS_ER_200=RMS_ER;
        Max_ER_200=Max_ER;
    end
end
% Results output

figure
plot(an,RMS_ER_025,an,RMS_ER_050,an,RMS_ER_075,an,RMS_ER_100,an,RMS_ER_125,an,RMS_ER_150,an,RMS_ER_175,an,RMS_ER_200)
xlabel('Iteration No', 'fontsize',12)
ylabel('Log_1_0 RMS(Error)', 'fontsize',12)
title('Convergence history using PSOR & 51X51 grid points','fontsize',12)
legend('\omega=0.25','\omega=0.50','\omega=0.75','\omega=1.00','\omega=1.25','\omega=1.50','\omega=1.75','\omega=2.00','Location','best')
grid on
set(findall(gcf,'type','line'),'linewidth',2.6)

figure
plot(an,Max_ER_025,an,Max_ER_050,an,Max_ER_075,an,Max_ER_100,an,Max_ER_125,an,Max_ER_150,an,Max_ER_175,an,Max_ER_200)
xlabel('Iteration No', 'fontsize',12)
ylabel('Log_1_0 Max(Error)', 'fontsize',12)
title('Convergence history using PSOR & 51X51 grid points','fontsize',12)
legend('\omega=0.25','\omega=0.50','\omega=0.75','\omega=1.00','\omega=1.25','\omega=1.50','\omega=1.75','\omega=2.00','Location','best')
grid on
set(findall(gcf,'type','line'),'linewidth',2.6)

