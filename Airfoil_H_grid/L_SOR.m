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

