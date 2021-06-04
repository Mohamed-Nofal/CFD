function geom

global x y imax jmax jair il it cord yal yau ps psp dx dy r d1 d2 omega Vinf cosa sina


%        ' function determines the x(i,j) and y(i,j) for all points
%          of H-grid for NACA-0012
%        ' ----------------------------------------------------------
%        '  il = i of the leading edge
%        '  it = i of the trailing edge
%        '  cord = chord length
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
           
 %'toc is the thickness to chord ratio for NACA 0012  toc=.12
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
