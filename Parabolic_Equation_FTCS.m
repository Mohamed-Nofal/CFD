clc;clear all;close all;
figure(1)
hold on
grid on
legend('on','Location','northwest')
xlabel('u');ylabel('y')
%% Geometric DATA
h=4/100;
U=40;
v=0.0002;
plot([0 U], [0 h],'--','LineWidth',1.5,'DisplayName','Exact')
%% Mesh Data
j_max=5;
n_max=5;
d_t_vec=[0.05,0.08,0.1,0.2,0.4,1];
ii=0;
for t_max=n_max.*d_t_vec
ii=ii+1;
d_t=d_t_vec(ii);
d_y=h/(j_max-1);
d=v*d_t/d_y^2  % diffusion Number
u=nan(j_max,n_max) ;
if d>0.5
%     error("Solution Fails")
end
%% Initial Conditions
u(:,1)=0*ones(j_max,1);

%% Boundary Conditions
u(1,:)=0*ones(1,n_max);
u(j_max,:)=U*ones(1,n_max);

%% Solution FTCS - explicit
for n=1:n_max-1
    for j=2:j_max-1
        u(j,n+1)=(1-2*d)*u(j,n)+d*(u(j+1,n)+u(j-1,n));
    end
end

%% Result
u_f=array2table(flipud(u),'VariableNames',strcat('T',cellstr(string(0:n_max-1))),'RowNames',strcat('j=', cellstr(string(j_max:-1:1))))

%% Plot
figure(1)
set(gca,'ytick',linspace(0,h,j_max))
plot(u(:,n_max),(0:d_y:h),'-o','LineWidth',1.5,'DisplayName',num2str(d))
end