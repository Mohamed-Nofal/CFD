clear all;   close all;   clc;
tic
% '                 SOLUTIION OF LAPLACE EQUATION
% '                 -----------------------------
% '
% '                    Two dimensional problem
% '                    -----------------------
% '   [ for an arbitrary shape with Dirichlet boundary conditions]
% '   ------------------------------------------------------------
% '    Definitions:
% '            ps(i,j)  = the dependent variable at iteration number "n"
% '            psp(i,j) = the dependent variable at iteration number "n+1"
% '            il      = the position number in i-direction of the leading  edge point
% '            it      = the position number in i-direction of the trailing edge point
% '            imax    = the number of points in x-direction
% '            jmax    = the number of points in y-direction
% '            lx      = the maximum length in x-direction
% '            ly      = the maximum length in y-direction
% '            nx      = number of interval in x-direction
% '            ny      = number of interval in y-direction
% '            dx      = step size in x-direction
% '            dy      = step size in y-direction
% '            x(i,j)  = the x-coordinate of the grid points
% '            y(i,j)  = the y-coordinate of the grid points
% '**********************************************************************


global x y imax jmax jair il it cord yal yau ps psp dx dy r d1 d2 omega Vinf cosa sina

% Data

             Vinf = 100; alfad = 10; cord = 1; nitd = 0;
             il = 31; it = 71; imax = 101; jair = 26; jmax = 51;
             omega = 1; per = .000001; nmax = 20;
% ' ----------------------------------------------------------------
             alfa = alfad * pi / 180;
             cosa = cos(alfa); sina = sin(alfa);
             uxinf = Vinf * cosa; uyinf = Vinf * sina;
             iimax = 2 * imax - 1; jjmax = 2 * jmax - 1; jjair = 2 * jair - 1;
             lx = 3 * cord; ly = cord; M = imax;
             nx = imax - 1; ny = jmax - 1; d1 = 1 / (it - il); d2 = 1 / ny;
             dx = cord / (it - il); dy = ly / ny;
             r = dx / dy; t1 = omega / (2 * (1 + r * r));
             t2 = t1 * r * r;
geom

            xxx = 0 ; % "if you want the SOR by Point method"
 %           xxx = 1 ; % "if you want the SOR by Line method"


% Boundary values  & Initialization
       
                 ps(1, 1) = 0;
       i = 1;   for j = 1 : jmax - 1
                 ii=2*i-1; jj=2*j-1;
                 ps(i, j + 1) = ps(i, j) + uxinf * (y(ii, jj + 2) - y(ii, jj));
                 end 
       j = 1;    for i = 1 : imax - 1
                 ii=2*i-1; jj=2*j-1;
                 ps(i + 1, j) = ps(i, j) - uyinf * (x(ii + 2, jj) - x(ii, jj));
                 end 
       i = imax; for j = 1 : jmax - 1
                 ii=2*i-1; jj=2*j-1;
                 ps(i, j + 1) = ps(i, j) + uxinf * (y(ii, jj + 2) - y(ii, jj));
                 end
       j = jmax; for i = 1 : imax - 1
                 ii=2*i-1; jj=2*j-1;
                 ps(i + 1, j) = ps(i, j) - uyinf * (x(ii + 2, jj) - x(ii, jj));
                 end 
                
                 for i = 2 : imax - 1
                 for j = 2 : jmax - 1
                     
        ps(i, j) = ps(i, 1) + (ps(i, jmax) - ps(i, 1)) * (j - 1) / (jmax - 1);
                 end 
                 end
                 psp = ps; 
                 
          
          for n=1:nmax
         %while erps > .001 ;   %n = nitd;   %n = n + 1;

% Calculation of new values of ps(i,j) = psp(i,j)

            if xxx == 0 ;  P_SOR ; end
            if xxx == 1 ;  L_SOR  ; end
% Errors calculation

            mder = 0;
            for i = 2 : nx
              for j = 2 : ny
              der = abs(psp(i, j) - ps(i, j));
              if der > mder ;  mder = der; ier = i; jer = j  ; end
              ps(i, j) = psp(i, j);
              end
            end
           
            if mder > 0 ; lmder = log10(mder) ;  end
            
% Updating the value of psi on the airfoil

       psold = psp(it, jair);
       psnew = psp(it + 1, jair);
       erps = abs(psnew - psold);
      
       j = jair;
       for i = il : it
       ps(i, j) = psnew;
       psp(i, j) = psnew;
       end 
       a_n(n)=n ; a_lmder(n)=lmder;

% Check convergence

%'            if ((mder > per) AND (n <= nmax)) THEN GO: iter
          end
                figure
                plot(a_n,a_lmder)
                grid on
                xlabel('Iteration number', 'fontsize',18)
                ylabel('Log_1_0 (Error)', 'fontsize',18)
 title(['Convergence history using Line-SOR for the flow past NACA-0012 airfoil with angle of attack =10^o'],'fontsize',8)

                results
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c11,c12,c22]=coef(ip,jp)
global x y d1 d2
% calculate the metric terms and the jacobian.
d1x = (x(ip + 1, jp) - x(ip - 1, jp)) / d1;
d1y = (y(ip + 1, jp) - y(ip - 1, jp)) / d1;
d2x = (x(ip, jp + 1) - x(ip, jp - 1)) / d2;
d2y = (y(ip, jp + 1) - y(ip, jp - 1)) / d2;
jaco = d1x * d2y - d1y * d2x;
c11 = (d2x * d2x + d2y * d2y) / jaco;
c12 = -(d1x * d2x + d1y * d2y) / jaco;
c22 = (d1x * d1x + d1y * d1y) / jaco;
end
% '____________________________________________________________________
function P_SOR
% '______________________________________________________
% PointSOR;   'Subprogram for determining psp(i,j)
% '              using the SOR by points method
% '              ------------------------------

global x y imax jmax jair il it cord yal yau ps psp dx dy r d1 d2 omega Vinf cosa sina
iimax = 2*imax-1 ; jjmax = 2*jmax-1;jjair = 2*jair-1;

 for j = 2 : jmax - 1
    if j == jair - 1 ; for ii = 1 : iimax; y(ii, jjair) = yal(ii); end; end
    if j == jair + 1 ; for ii = 1 : iimax; y(ii, jjair) = yau(ii); end; end

         for i = 2 : imax - 1
         ii = 2 * i - 1;
         jj = 2 * j - 1;

         ip = ii + 1; jp = jj; [c11ip c12ip c22ip]=coef(ip,jp);
                               
         ip = ii - 1; jp = jj; [c11im c12im c22im]=coef(ip,jp);
                               
         ip = ii; jp = jj + 1; [c11jp c12jp c22jp]=coef(ip,jp);
                               
         ip = ii; jp = jj - 1; [c11jm c12jm c22jm]=coef(ip,jp);
                               
%   ' The difference equation is written in the following form
%   ' sij * ps(i,j) = sim * ps(i-1,j) + sip * ps(i+1,j) + sjm * ps(i,j-1)
%   '              + sjp * ps(i,j+1) + smm * ps(i-1,j-1) + smp * ps(i-1,j+1)
%   '              + spm * ps(i+1,j-1) + spp * ps(i+1,j+1)

         sij = c11ip + c11im + r * r * (c22jp + c22jm);
         sim = c11im - r * r * (c12jp - c12jm) / 4;
         sip = c11ip + r * r * (c12jp - c12jm) / 4;
         sjm = r * r * c22jm - r * (c12ip - c12im) / 4;
         sjp = r * r * c22jp + r * (c12ip - c12im) / 4;
         smm = r * (c12im + c12jm) / 4;
         smp = -r * (c12im + c12jp) / 4;
         spm = -r * (c12ip + c12jm) / 4;
         spp = r * (c12ip + c12jp) / 4;

 psp(i, j) = sim * ps(i - 1, j) + sip * ps(i + 1, j) + sjm * ps(i, j - 1) + sjp * ps(i, j + 1);
 psp(i, j) = psp(i, j) + smm * ps(i - 1, j - 1) + smp * ps(i - 1, j + 1);
 psp(i, j) = psp(i, j) + spm * ps(i + 1, j - 1) + spp * ps(i + 1, j + 1);
 psp(i, j) = psp(i, j) / sij;
           
            psp(i, j) = ps(i, j) + omega * (psp(i, j) - ps(i, j));
         end
if j == jair; for i=il:it ;psp(i,j)=ps(i,j);end; end
                     
            
 end 
end       
function L_SOR

% '____________________________________________________________________________
% LineSOR;    'Subprogram for determining psp(i,j)
% '              using the SOR by lines method
% '              -----------------------------

global x y imax jmax jair il it cord yal yau ps psp dx dy r d1 d2 omega Vinf cosa sina

iimax = 2*imax-1 ; jjmax = 2*jmax-1;jjair = 2*jair-1;

for j = 2 : jmax-1
     if (j == jair - 1) ; for ii = 1 : iimax; y(ii, jjair) = yal(ii); end ; end
     if (j == jair + 1) ; for ii = 1 : iimax; y(ii, jjair) = yau(ii); end ; end
          
%        ' Calculation of the coefficients b(i),d(i),a(i) and c(i)
%        '    for   i=1
            b(1) = 0; d(1) = 1; a(1) = 0; c(1) = ps(1, j);
%        '    for   i=imax
            b(imax) = 0; d(imax) = 1; a(imax) = 0; c(imax) = ps(imax, j);
        

         for i = 2 : imax-1
 
         ii = 2 * i - 1;
         jj = 2 * j - 1;

         ip = ii + 1; jp = jj; [c11ip c12ip c22ip]=coef(ip,jp);
                               
         ip = ii - 1; jp = jj; [c11im c12im c22im]=coef(ip,jp);
                               
         ip = ii; jp = jj + 1; [c11jp c12jp c22jp]=coef(ip,jp);
                               
         ip = ii; jp = jj - 1; [c11jm c12jm c22jm]=coef(ip,jp);
                 
%   ' The difference equation is written in the following form
%   ' sij * ps(i,j) = sim * ps(i-1,j) + sip * ps(i+1,j) + sjm * ps(i,j-1)
%   '              + sjp * ps(i,j+1) + smm * ps(i-1,j-1) + smp * ps(i-1,j+1)
%   '              + spm * ps(i+1,j-1) + spp * ps(i+1,j+1)

         sij = c11ip + c11im + r * r * (c22jp + c22jm);
         sim = c11im - r * r * (c12jp - c12jm) / 4;
         sip = c11ip + r * r * (c12jp - c12jm) / 4;
         sjm = r * r * c22jm - r * (c12ip - c12im) / 4;
         sjp = r * r * c22jp + r * (c12ip - c12im) / 4;
         smm = r * (c12im + c12jm) / 4;
         smp = -r * (c12im + c12jp) / 4;
         spm = -r * (c12ip + c12jm) / 4;
         spp = r * (c12ip + c12jp) / 4;

         b(i) = -sim; d(i) = sij; a(i) = -sip;
    c(i) = sjm * psp(i, j - 1) + sjp * ps(i, j + 1) + smm * psp(i - 1, j - 1);
    c(i) = c(i) + smp * ps(i - 1, j + 1) + spm * psp(i + 1, j - 1) + spp * ps(i + 1, j + 1);

      if ((j == jair) && (i >= il) && (i <= it)) ;
          b(i) = 0; d(i) = 1; a(i) = 0; c(i) = psp(i, j);
      end
    
         end 

            ps_p=tri_sol(a,b,c,d,imax); 
            for i=1:imax ; psp(i,j)=ps_p(i); end
 end 
end
function geom

global x y imax jmax jair il it cord yal yau 
% function determines the x(i,j) and y(i,j) for all points
%  of H-grid for NACA-0012
 
%  il = i of the leading edge
%  it = i of the trailing edge
%  cord = chord length
   figure
   iil = 2 * il - 1;
   iit = 2 * it - 1;
   iimax = 2 * imax - 1;
   jjmax = 2 * jmax - 1;
   jjair = 2 * jair - 1;
  
   for jj = 1 : jjmax
   for ii = 1 : iimax
   x(ii, jj) = (cord / (iit - iil)) * (ii - iil);
   end  
   end
   
 %toc is the thickness to chord ratio for NACA 0012  toc=.12
   toc = .12;
   for ii = 1 : iil; y(ii, jjair) = 0; end 
   for ii = iil : iit; xp = x(ii, jjair);
 y(ii, jjair) = -5 * toc * (.2969 * sqrt(xp) - .126 * xp - .3537 * xp ^ 2 + .2843 * xp ^ 3 - .1015 * xp ^ 4);
   end 
   for ii = iit : iimax; y(ii, jjair) = 0; end
  
   for ii = 1 : iimax; yal(ii) = y(ii, jjair); end

   for ii = 1 : iimax; y(ii, 1) = -cord; end 

   for ii = 1 : iimax
   for jj = 2 : jjair - 1
 y(ii, jj) = y(ii, 1) + (jj - 1) * (y(ii, jjair) - y(ii, 1)) / (jjair - 1);
   end
   end 
  
%'toc is the thickness to chord ratio for NACA 0012  toc=.12
   toc = .12;
   for ii = 1 : iil; y(ii, jjair) = 0; end
   for ii = iil : iit; xp = x(ii, jjair);
 y(ii, jjair) = 5 * toc * (.2969 * sqrt(xp) - .126 * xp - .3537 * xp ^ 2 + .2843 * xp ^ 3 - .1015 * xp ^ 4);
   end
   
   for ii = iit : iimax; y(ii, jjair) = 0; end

   for ii = 1 : iimax; yau(ii) = y(ii, jjair); end 
  
   for ii = 1 : iimax; y(ii, jjmax) = cord; end 

   for ii = 1 : iimax
   for jj = jjair + 1 : jjmax - 1
 y(ii, jj) = y(ii, jjair) + (jj - jjair) * (y(ii, jjmax) - y(ii, jjair)) / (jjmax - jjair);
   end ; 
   end 


 % Plot the H-Grid 
 % ---------------
for j = 1 : jair - 1;    jj = 2 * j - 1;
         x1 = x(:,jj); y1=y(:,jj);plot (x1,y1); hold on; 
end
         jj = jjair ;x1 = x(:,jj); y1=yal(:);plot (x1,y1); hold on
         x1 = x(:,jj); y1=yau(:);plot (x1,y1); hold on
for j = jair + 1 : jmax;  jj = 2 * j - 1;
         x1 = x(:, jj); y1 = y(:, jj);plot (x1,y1); hold on; 
end

         y(:, jjair) = yal(:);
         for i=1:il; ii=2*i-1; x1 = x(ii,:);
         y1 = y(ii,:);plot (x1,y1); hold on; end 
         
         y(:, jjair) = yal(:);
         for i=il+1:it-1; ii=2*i-1;
             for j=1:jair; jj=2*j-1; 
                 x2(j)= x(ii,jj);y2(j)= y(ii,jj);
             end
             plot (x2,y2); hold on; 
         end 
         
         y(:, jjair) = yau(:);
         for i=il+1:it-1;ii=2*i-1;
              for j=jair:jmax; jj=2*j-1;k=j-jair+1; 
                  x2(k)= x(ii,jj);y2(k)= y(ii,jj);
              end
                  plot (x2,y2); hold on; 
         end 
     
     for i=it:imax;ii=2*i-1;x1 = x(ii,:);
         y1 = y(ii,:);plot (x1,y1); hold on; 
     end 
         
         xlabel('X-axis', 'fontsize',18)
         ylabel('Y-axis', 'fontsize',18)
         title(['H-Grid for NACA-0012 airfoil '],'fontsize',18)
end
function results

global x y imax jmax jair il it cord yal yau ps psp dx dy r d1 d2 omega Vinf cosa sina
iimax = 2*imax-1 ; jjmax = 2*jmax-1;jjair = 2*jair-1;
% '____________________________________________________________________
%'____________________________________________________________________________
% results;   'function of calculation of the velocity and the pressure coefficients
% '-------    ---------------------------------------------------------------------
          uxinf = Vinf * cosa; uyinf = Vinf * sina;    
          i=1    ; for j=1:jmax ; a_vx(i,j)=uxinf; a_vy(i,j)=uyinf; end
          i=imax ; for j=1:jmax ; a_vx(i,j)=uxinf; a_vy(i,j)=uyinf; end
          j=1    ; for i=1:imax ; a_vx(i,j)=uxinf; a_vy(i,j)=uyinf; end
          j=jmax ; for i=1:imax ; a_vx(i,j)=uxinf; a_vy(i,j)=uyinf; end
          
           
           for i=2:imax-1
                    for j=2:jmax-1
                        ii=2*i-1;jj=2*j-1;
               
               d1x = (x(ii + 1, jj) - x(ii - 1, jj)) / d1;
               d1y = (y(ii + 1, jj) - y(ii - 1, jj)) / d1;
               d2x = (x(ii, jj + 2) - x(ii, jj)) / d2;
               d2y = (y(ii, jj + 2) - y(ii, jj)) / d2;
               jaco = d1x * d2y - d1y * d2x;
               et1x = d2y / jaco; et1y = -d2x / jaco;
               et2x = -d1y / jaco; et2y = d1x / jaco;
           
               d1u = (ps(i + 1, j) - ps(i - 1, j)) / 2 / d1;
              
               d2u = (ps(i, j + 1) - ps(i, j - 1)) / 2 / d2;

               a_vx(i,j) = d1u * et1y + d2u * et2y;
               a_vy(i,j) = -(d1u * et1x + d2u * et2x);
                    end
           end
          
%      
            j = jair;
% ' for upper surface      
% ' -----------------

  for ii = 1 : iimax; y(ii, jjair) = yau(ii); end 
            for i = il : it
            ii = 2 * i - 1; jj = 2 * j - 1;
              
               d1x = (x(ii + 1, jj) - x(ii - 1, jj)) / d1;
               d1y = (y(ii + 1, jj) - y(ii - 1, jj)) / d1;
               d2x = (x(ii, jj + 2) - x(ii, jj)) / d2;
               d2y = (y(ii, jj + 2) - y(ii, jj)) / d2;
               jaco = d1x * d2y - d1y * d2x;
               et1x = d2y / jaco; et1y = -d2x / jaco;
               et2x = -d1y / jaco; et2y = d1x / jaco;
      
           d1u = (ps(i + 1, j) - ps(i - 1, j)) / 2 / d1;
           d2u = (4 * ps(i, j + 1) - 3 * ps(i, j) - ps(i, j + 2)) / 2 / d2;

            vx = d1u * et1y + d2u * et2y;
            vy = -(d1u * et1x + d2u * et2x);
            a_vx_up(i)=vx ; a_vy_up(i)=vy;
           % a_vx(i,j)=vx; a_vy(i,j)=vy;
            Vel = sqrt(vx * vx + vy * vy);
            Vru(i) = Vel / Vinf;
            Cpu(i) = 1 - Vru(i) * Vru(i);
            xup(i) = x(ii, jj);
            yup(i) = y(ii, jj);

            end 
             
             
% ' for lower surface      
% ' -----------------

            j = jair;
  for ii = 1 : iimax; y(ii, jjair) = yal(ii); end 
            for i = il : it
            ii = 2 * i - 1; jj = 2 * j - 1;
             
               d1x = (x(ii + 1, jj) - x(ii - 1, jj)) / d1;
               d1y = (y(ii + 1, jj) - (y(ii - 1, jj))) / d1;
               d2x = (x(ii, jj) - x(ii, jj - 2)) / d2;
               d2y = (y(ii, jj) - y(ii, jj - 2)) / d2;
               jaco = d1x * d2y - d1y * d2x;
               et1x = d2y / jaco; et1y = -d2x / jaco;
               et2x = -d1y / jaco; et2y = d1x / jaco;
              
               d1u = (ps(i + 1, j) - ps(i - 1, j)) / 2 / d1;
               d2u = (-4 * ps(i, j - 1) + 3 * ps(i, j) + ps(i, j - 2)) / 2 / d2;

            vx = d1u * et1y + d2u * et2y;
            vy = -(d1u * et1x + d2u * et2x);
            a_vx_lo(i)=vx ; a_vy_lo(i)=vy;
            %a_vx(i,j)=vx; a_vy(i,j)=vy;
            Vel = sqrt(vx * vx + vy * vy);
            Vrl(i) = Vel / Vinf;
            Cpl(i) = 1 - Vrl(i) * Vrl(i);
            xlo(i) = x(ii, jj);
            ylo(i) = y(ii, jj);

            end
            
            figure
            hold on
            grid on
            plot(xup,yup,xlo,ylo)
            plot(xup,Vru,xlo,Vrl)
            xlabel('Chord line', 'fontsize',18)
            ylabel('Non-dimensional velocity', 'fontsize',18)
 title(['Non-dimensional velocity over NACA-0012 airfoil surface(angle of attack =10^o)'],'fontsize',12)
 legend('upper surface','lower surface','Location','best'),grid  
            
            figure
            hold on
            grid on
            plot(xup,yup,xlo,ylo)
            plot(xup,Cpu,xlo,Cpl)
            xlabel('Chord line', 'fontsize',18)
            ylabel('Pressure coefficient', 'fontsize',18)
 title(['Pressure coefficient over NACA-0012 airfoil surface(angle of attack =10^o)'],'fontsize',12)
 legend('upper surface','lower surface','Location','best'),grid  
 
 figure
 i=1; ii=2*i-1; for j=1:jmax ; jj=2*j-1; x1(j)=x(ii,jj); y1(j)=y(ii,jj);end
 i=imax; ii=2*i-1; for j=1:jmax ; jj=2*j-1; x2(j)=x(ii,jj); y2(j)=y(ii,jj);end
 j=1; jj=2*j-1; for i=1:imax ; ii=2*i-1; x3(i)=x(ii,jj); y3(i)=y(ii,jj);end
 j=jmax; jj=2*j-1; for i=1:imax ; ii=2*i-1; x4(i)=x(ii,jj); y4(i)=y(ii,jj);end
 
 plot (x1,y1,x2,y2,x3,y3,x4,y4)
 hold on
  xlabel('X-axis', 'fontsize',18)
                ylabel('Y-axis', 'fontsize',18)
 title(['Stream lines for the flow past NACA-0012 airfoil with angle of attack =10^o'],'fontsize',12)         
      
     j=jair; jj=2*j-1;
    for i=il:it; ii=2*i-1;k=i-il+1;x5(k)= x(ii,jj);y5(k)=yal(ii); end
      plot (x5,y5); hold on;           
    for i=il:it; ii=2*i-1;k=i-il+1;x6(k)= x(ii,jj);y6(k)=yau(ii); end
      plot (x6,y6); hold on;    
               
               for ii = 1 : iimax; y(ii, jjair) = yau(ii); end 
               for i=1:imax
                   ii=2*i-1;
                    for j=jair:jmax
                        jj=2*j-1;
                         k=j-jair+1; x7(i,k)=x(ii,jj); y7(i,k)=y(ii,jj);p7(i,k)=ps(i,j);
                    end
               end
               
               contour(x7,y7,p7,50)
               hold on;
               
               for ii = 1 : iimax; y(ii, jjair) = yal(ii); end 
               for i=1:imax
                   ii=2*i-1;
                    for j=1:jair
                        jj=2*j-1;
                          x8(i,j)=x(ii,jj); y8(i,j)=y(ii,jj);p8(i,j)=ps(i,j);
                    end
               end
 %              figure
               contour(x8,y8,p8,50)
               hold on;
               

                figure
 i=1; ii=2*i-1; for j=1:jmax ; jj=2*j-1; x1(j)=x(ii,jj); y1(j)=y(ii,jj);end
 i=imax; ii=2*i-1; for j=1:jmax ; jj=2*j-1; x2(j)=x(ii,jj); y2(j)=y(ii,jj);end
 j=1; jj=2*j-1; for i=1:imax ; ii=2*i-1; x3(i)=x(ii,jj); y3(i)=y(ii,jj);end
 j=jmax; jj=2*j-1; for i=1:imax ; ii=2*i-1; x4(i)=x(ii,jj); y4(i)=y(ii,jj);end
 
 plot (x1,y1,x2,y2,x3,y3,x4,y4)
 hold on
  xlabel('X-axis', 'fontsize',18)
                ylabel('Y-axis', 'fontsize',18)
 title(['Velocity vector for the flow past NACA-0012 airfoil with angle of attack =10^o'],'fontsize',12)         
      
     j=jair; jj=2*j-1;
    for i=il:it; ii=2*i-1;k=i-il+1;x5(k)= x(ii,jj);y5(k)=yal(ii); end
      plot (x5,y5); hold on;           
    for i=il:it; ii=2*i-1;k=i-il+1;x6(k)= x(ii,jj);y6(k)=yau(ii); end
      plot (x6,y6); hold on;    
               
               for ii = 1 : iimax; y(ii, jjair) = yau(ii);end 
               for i=il:it ; a_vx(i,jair)=a_vx_up(i);a_vy(i,jair)=a_vy_up(i); end 
               for i=1:imax
                   ii=2*i-1;
                    for j=jair:jmax
                        jj=2*j-1;
                         k=j-jair+1; x7(i,k)=x(ii,jj); y7(i,k)=y(ii,jj);a_vx_7(i,k)=a_vx(i,j);a_vy_7(i,k)=a_vy(i,j);
                    end
               end
               quiver(x7,y7,a_vx_7,a_vy_7)
               hold on;
               
               for ii = 1 : iimax; y(ii, jjair) = yal(ii); end 
               for i=il:it ;a_vx(i,jair)=a_vx_lo(i);a_vy(i,jair)=a_vy_lo(i); end
               for i=1:imax
                   ii=2*i-1;
                    for j=1:jair
                        jj=2*j-1;
                          x8(i,j)=x(ii,jj); y8(i,j)=y(ii,jj);a_vx_8(i,j)=a_vx(i,j);a_vy_8(i,j)=a_vy(i,j);
                    end
               end
               quiver(x8,y8,a_vx_8,a_vy_8)
               hold on;
               xlabel('X-axis', 'fontsize',18)
               ylabel('Y-axis', 'fontsize',18)
  title(['Velocity vector for the flow past NACA-0012 airfoil with angle of attack =10^o'],'fontsize',12)
%             

               
               
%              j=jair; jj=2*j-1;
%     for i=il:it; ii=2*i-1;k=i-il+1;x2(k)= x(ii,jj);y2(k)=yal(ii); end
%       plot (x2,y2); hold on;           
%     for i=il:it; ii=2*i-1;k=i-il+1;x2(k)= x(ii,jj);y2(k)=yau(ii); end
%       plot (x2,y2); hold on;    
%                 xlabel('X-axis', 'fontsize',18)
%                 ylabel('Y-axis', 'fontsize',18)
%  title(['Stream lines for the flow past NACA-0012 airfoil with angle of attack =10^o'],'fontsize',12)         
  
     
%                for i=1:imax
%                    ii=2*i-1;
%                     for j=1:jair
%                         jj=2*j-1;
%                           if (j == jair ) ; for ii = 1 : iimax; y(ii, jjair) = yal(ii); end ;end
%                           x1(i,j)=x(ii,jj); y1(i,j)=y(ii,jj);p1(i,j)=ps(i,j);
%                     end
%                end
%                figure
%                contour(x1,y1,p1,30)
%                hold on;
%                
%                for i=1:imax
%                    ii=2*i-1;
%                     for j=jair:jmax
%                         jj=2*j-1;
%                           if (j == jair ) ; for ii = 1 : iimax; y(ii, jjair) = yau(ii); end ; end
%                           k=j-jair+1;x1(k,j)=x(ii,jj); y1(k,j)=y(ii,jj);p1(k,j)=ps(i,j);
%                     end
%                end
%                figure
%                contour(x1,y1,p1,30)
%                hold on;
%                
%                for i=1:imax
%                     for j=1:jmax
%                         ii=2*i-1;jj=2*j-1; 
%                         x1(i,j)=x(ii,jj);y1(i,j)=y(ii,jj);
%                     end
%                end
%                 figure
%                 contour(x1,y1,ps,50)
%                 hold on;
%                 
%                 j=jair; jj=2*j-1;
%     for i=il:it; ii=2*i-1;k=i-il+1;x2(k)= x(ii,jj);y2(k)=yal(ii); end
%       plot (x2,y2); hold on;           
%     for i=il:it; ii=2*i-1;k=i-il+1;x2(k)= x(ii,jj);y2(k)=yau(ii); end
%       plot (x2,y2); hold on;    
%                 xlabel('X-axis', 'fontsize',18)
%                 ylabel('Y-axis', 'fontsize',18)
%  title(['Stream lines for the flow past NACA-0012 airfoil with angle of attack =10^o'],'fontsize',12)         
% 
%                 figure
%                 quiver(x1,y1,a_vx,a_vy,2)
%                 hold on
%                 for i=il:it; ii=2*i-1;k=i-il+1;x2(k)= x(ii,jj);y2(k)=yal(ii); end
%       plot (x2,y2); hold on;           
%                 for i=il:it; ii=2*i-1;k=i-il+1;x2(k)= x(ii,jj);y2(k)=yau(ii); end
%       plot (x2,y2); hold on; 
%                 xlabel('X-axis', 'fontsize',18)
%                 ylabel('Y-axis', 'fontsize',18)
%  title(['Velocity vector for the flow past NACA-0012 airfoil with angle of attack =10^o'],'fontsize',12)
%             
% 'Calculation of the lift and drag coefficients    
% '---------------------------------------------   
      cx = 0; cy = 0;
      j = jair; jj = jjair;
      for i = il : it - 1
      ii = 2 * i - 1;
 x1 = x(ii, jj); y1 = yau(ii); x2 = x(ii + 2, jj); y2 = yau(ii + 2);
 cx = cx + .5 * (Cpu(i) + Cpu(i + 1)) * (y2 - y1) / cord;
 cy = cy - .5 * (Cpu(i) + Cpu(i + 1)) * (x2 - x1) / cord;
 x1 = x(ii, jj); y1 = yal(ii); x2 = x(ii + 2, jj); y2 = yal(ii + 2);
 cx = cx + .5 * (Cpl(i) + Cpl(i + 1)) * (y2 - y1) / cord;
 cy = cy + .5 * (Cpl(i) + Cpl(i + 1)) * (x2 - x1) / cord;
      end 
    cl = cy * cosa - cx * sina
    cd = cy * sina + cx * cosa

end
function e=tri_sol(a,b,c,d,M)
% '____________________________________________________________
% tri:     ' Subroutine determine the
% '          solution of a Tridiagonal matrix given by
% '          b(i)*e(i-1) + d(i)*e(i) + a(i)*e(i+1) = c(i)
% '          **********************************************
           
            for i=2:M
            t = b(i) / d(i - 1);
            d(i) = d(i) - t * a(i - 1);
            c(i) = c(i) - t * c(i - 1);
            end
                e(M) = c(M) / d(M);
            for k=2:M
            i = M - k + 1;
            e(i) = (c(i) - a(i) * e(i + 1)) / d(i);
            end
           % RETURN
% '____________________________________________________________________
end