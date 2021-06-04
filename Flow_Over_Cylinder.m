%% Flow Over Clinder Using PSOR scheme

clc;clear all;close all;    
imax=51;jmax=51;
a     = 1;
R     = 4*a;
omega = 1;  
V_inf = 100 ;
% etta_1   = theta/pi
etta_2   = log(R/a)/pi;
dzeta    = 1/(imax-1);
deta     = etta_2/(jmax-1);
beta     = dzeta/deta;

% Inner circle generation 
for i = 1:imax
    j = 1;    
    zeta  = (i-1)/(imax-1); 
    eta   = (j-1)*etta_2/(jmax-1);
    x_c(i)= a*exp(pi*eta)*cos(pi*zeta);
    y_c(i)= a*exp(pi*eta)*sin(pi*zeta);
end

% Outer circle generation 
for i = 1:imax 
    j = jmax;
    zeta = (i-1)/(imax-1);
    eta  = (j-1)*etta_2/(jmax-1);
    x_o(i)= a*exp(pi*eta)*cos(pi*zeta);
    y_o(i)= a*exp(pi*eta)*sin(pi*zeta);
end

%Grid Generation :
for i=1:imax 
    for j = 1:jmax
      zeta   = (i-1)/(imax-1);
      eta    = (j-1)*etta_2/(jmax-1);
      x(i,j) = a*exp(pi*eta)*cos(pi*zeta);
      y(i,j) = a*exp(pi*eta)*sin(pi*zeta);
    end
end

for i=1:imax
    for j = 1:jmax
       xm(j,i) = x(i,j);
       ym(j,i) = y(i,j);
    end
end

% initial and boundary conditions

epsii=zeros(imax,jmax);

for i=1:imax
    for  j=1:jmax
        epsi=0;
        if i==imax && i==1
            epsi=0;
        end
        if j==1
            epsi=0;
        end
        if j==jmax
              zeta=(i-1)/(imax-1);
              etta=(etta_2*(j-1))/(jmax-1);
              r=a*exp(pi*etta);
              theta=pi*zeta;
              epsi=V_inf*r*sin(theta)*(1-(a^2/r^2));
        end
        epsii(i,j)=epsi;
    end 
end
% linear interpolations
for i=2:imax-1
    for j=2:jmax-1
        epsi=epsii(i,1)+(((j-1)/(jmax-1))*(epsii(i,jmax)-epsii(i,1)));
        epsii(i,j)=epsi;
    end
end
Epsi_initial=epsii;  adder=0;

% solving using PSOR scheme

n=1;   rms=1;
while rms>1e-9
   for i=1:imax;
      for j=1:jmax ;
        eps_new=epsii(i,j);
        epsi_ex=epsii(i,j);
          if i>1 && i<imax
                if j>1 && j<jmax
                     epsi_i_j=epsii(i,j);
                     epsi_ip1_j=epsii(i+1,j);
                     epsi_im1_j_np1=epsii_new(i-1,j);
                     epsi_i_jp1=epsii(i,j+1);
                     epsi_i_jm1_np1=epsii_new(i,j-1);
                     epsi_dash=(1/(2*(1+(beta^2))))*(epsi_ip1_j+epsi_im1_j_np1+((beta^2)*(epsi_i_jp1+ epsi_i_jm1_np1)));
                     eps_new=epsi_i_j+(omega*((epsi_dash)-epsi_i_j));
                 end
          end
       epsii_new(i,j)=eps_new;    epsiii_new    =eps_new;
       s1=(epsiii_new-epsi_ex)^2;   s2=adder+s1;     adder=s2;
      end
  end
   rms=(adder^(0.5))/((imax-1)*(jmax-1));   epsii=epsii_new;
   adder=0;  
   n=n+1  ; 
   log_rms=log10(rms);
   a_n(n)=n;a_log_rms(n)=log_rms;
end

% Convergence History
figure(1)
plot (a_n,a_log_rms)
xlabel('iteration number','fontsize',14);ylabel ('log_1_0(RMS.Error)','fontsize',14) ;
title('Convergence History Using PSOR scheme','fontsize',14)
grid
set(findall(gcf,'type','line'),'linewidth',2.6)
axis tight 

% Calculation  of u ,v ,V ,Cp
Epsi=epsii;
for i=1:imax
    for j=1:jmax
       zeta=(i-1)/(imax-1);
       eta=(j-1)*etta_2/(jmax-1);
        x_zeta=-pi*a*exp(pi*eta)*(sin(pi*zeta));
        y_zeta=pi*a*exp(pi*eta)*(cos(pi*zeta));
        x_eta=pi*a*exp(pi*eta)*(cos(pi*zeta));
        y_eta=pi*a*exp(pi*eta)*(sin(pi*zeta));
        Jac=x_zeta*y_eta-x_eta*y_zeta;
        zeta_x(i,j)=(y_eta)/Jac;
        zeta_y(i,j)=(-x_eta)/Jac;
        eta_x(i,j)=(-y_zeta)/Jac;
        eta_y(i,j)=(x_zeta)/Jac;
    end
end
 
%1) v   DISTRIBUTION    &     %2) cp  dISTRIBUTION
for j=1:jmax
    for i=1:imax
         if i==1
             depsi_dzeta=(-3*Epsi(i,j)+4*Epsi(i+1,j)-Epsi(i+2,j))/2/dzeta;
           elseif (i>1)&&(i<imax)
             depsi_dzeta=(Epsi(i+1,j)-Epsi(i-1,j))/2/dzeta;
           elseif  i==imax
           depsi_dzeta=(Epsi(i-2,j)-4*Epsi(i-1,j)+3*Epsi(i,j))/2/dzeta;
         end
        if j==1
           depsi_deta=(-3*Epsi(i,j)+4*Epsi(i,j+1)-Epsi(i,j+2))/2/deta;
        elseif (j>1)&&(j<jmax)
            depsi_deta=(Epsi(i,j+1)-Epsi(i,j-1))/2/deta;
        elseif  j==jmax
           depsi_deta=(Epsi(i,j-2)-4*Epsi(i,j-1)+3*Epsi(i,j))/2/deta;
        end
        
        u(i,j)=depsi_dzeta*zeta_y(i,j)+depsi_deta*eta_y(i,j);
        v(i,j)=-(depsi_dzeta*zeta_x(i,j)+depsi_deta*eta_x(i,j));  
        V(i,j)=(u(i,j)^2+v(i,j)^2)^0.5;       % THE VELOCIT
        Cp(i,j)=1-(V(i,j)/V_inf)^2;          % CP DISTRIBTUION
    end
end

%*************  RESULTS ************************

%  coordinates
 figure(2)
 plot(x,y,'k')
 hold on
 plot(xm,ym,'k'); hold on
 plot(x_o,y_o,'k','linewidth',2);hold on
 plot(x_c,y_c,'k','linewidth',2) ; xlabel('x','fontsize',14);ylabel ('y','fontsize',14) ;
 title('Grid over a circular cylinder','fontsize',14)
 axis tight
% Initial stream lines
% figure(3)
% contour(x,y,Epsi_initial,25)
% hold on
% plot(x_c,y_c); hold on
% plot(x_o,y_o) ; xlabel('x','fontsize',18);ylabel ('y','fontsize',18) ; title(' INITIAL STREAM LINES DISTRIBUTIONS','fontsize',18)
%  stream lines
figure(4)
contour(x,y,Epsi,40,'b')
hold on
area(x_c,y_c,'FaceColor','r');hold on
plot(x_o,y_o,'k','linewidth',2) ; xlabel('x','fontsize',14);ylabel ('y','fontsize',14) ; 
title(' Stream lines contour','fontsize',14)

% Non dimension velocity distribution and velocity contour
for i=1:imax  ;  Vel(i)=V(i,1)/V_inf;  end

figure(5)
plot(x_c,Vel); xlabel('x','fontsize',14);ylabel ('V/V_o_o','fontsize',14) ; 
title('Non dimensional velocity distribution over the circular cylinder','fontsize',12)
grid
set(findall(gcf,'type','line'),'linewidth',2.6)

figure(6)
contourf(x,y ,V,40)
hold on
plot(x_o,y_o);hold on
plot(x_c,y_c) ; xlabel('x','fontsize',14);ylabel ('y','fontsize',14) ; 
title('Velocity contour','fontsize',14)
colormap(jet)
colorbar

% Cp distribution and pressure contour
for i=1:imax ;  cpp(i)=Cp(i,1);  end  
figure(7)
plot(x_c,cpp); xlabel('x','fontsize',14);ylabel ('C_p','fontsize',14) ; 
title('Pressure distribution over the circular cylinder')
grid
set(findall(gcf,'type','line'),'linewidth',2.6)

figure(8)
contourf(x,y ,Cp,20)
hold on
plot(x_c,y_c);hold on
plot(x_o,y_o);
xlabel('x','fontsize',14);ylabel ('y','fontsize',14) ;
title('Pressure coefficient contour','fontsize',14)
colormap(jet)
colorbar
