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
