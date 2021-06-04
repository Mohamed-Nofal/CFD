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
            
