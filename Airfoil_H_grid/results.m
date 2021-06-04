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

