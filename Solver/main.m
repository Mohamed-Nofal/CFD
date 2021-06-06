                        %% Created by Mo7aMeD Adel %%
                     %% Computitional Fluid Dynamics %%
                             %% 22 / 3 / 2016 %%
clc
clear all
close all
% Notes:
% 1) That this code SOLVES the flow over NACA 4-digit, Joukowski Airfoils and Cylinders.
% 2) This code solves the PDE: [[Exx + Eyy = 0]]       "where E is the Stream Function"
% 3) Number of nodes in x and y should be EQUAL.
% 4) The Grid used in this solution is a H-Grid .
% 5) The solver uses PSOR (Non-Updating) Scheme.
%% Body Equation
AirfoilTypeAsk = 'Enter "J" for Jokowski Airfoil , "N" for NACA Airfoils or "C" for Cylinder:  ';      % Body Type
Type = input(AirfoilTypeAsk);
if Type == 'N'
AirfoilAsk = 'Enter The Airfoil Number in Raw Vector as [x x x x]:  ';                                 % Airfoil number
Airfoil = input(AirfoilAsk);
NACA = Airfoil;
ChordAsk = 'Enter The Airfoil Chord Length:  ';                                                        % Airfoil Chord
Chord = input(ChordAsk);
[x_upper,y_upper,x_lower,y_lower] = NACAFoil(NACA,Chord);
elseif Type == 'J'
AirfoilAsk = 'Enter The Airfoil max. Thickness/Chord and max. Camber/Chord as [t/c C/c]:  ';           % Airfoil Thickness and Camber
Airfoil = input(AirfoilAsk);
ChordAsk = 'Enter The Airfoil Chord Length:  ';                                                        % Airfoil Chord
Chord = input(ChordAsk);
t_c = Airfoil(1);
C_c = Airfoil(2);
[x_upper,y_upper,x_lower,y_lower] = JFoil(t_c,C_c,Chord);    
elseif Type == 'C'
CylinderAsk = 'Enter The Radius:  '; 
Raduis = input(CylinderAsk);
Chord = 2*Raduis;
[x_upper,y_upper,x_lower,y_lower] = cylinder(Chord);    
end
%% System Data
V_inf = 100;                               % Free Stream Velocity
Alpha = 0;                                % Angle of Attack
i_max = 100;                               % Number of points on x-axis (Should be positive even)
j_max = 100;                               % Number of points on y-axis (Should be positive even)
n_max = 51;                               % Maximum number of iterations
error = 1e-4;                              % Max. Error Allowed for Convergance Checks
%% Grid initialization
x1 = -1*Chord;                              % starting x
x2 = 2*Chord;                               % ending x
y1 = -1*Chord;                              % starting y
y2 = 1*Chord;                               % ending y
I_max = i_max+2;                            % Because of some removal of repeated points will occure
J_max = j_max+2;                            % Because of some removal of repeated points will occure
x  = linspace(x1,0,ceil(I_max/3));
xx  = linspace(0,Chord,floor(I_max/3));     % To ensure that the leading and trailing edges are on the grid
xxx  = linspace(Chord,x2,floor(I_max/3));
x = [x xx xxx];
ind1 = find (x == 0);
ind2 = find (x == Chord);
x = [x(1:ind1(1)) x(ind1(end)+1:ind2(1)) x(ind2(end)+1:end)];
ind3 = find (x == 0);
ind4 = find (x == Chord);
yu_airfoil = polyval(polyfit(x_upper,y_upper,10),x(ind3:ind4));
yl_airfoil = polyval(polyfit(x_lower,y_lower,10),x(ind3:ind4));
yu_airfoil(1) = 0;
yl_airfoil(1) = 0;
yu_airfoil(end) = 0;
yl_airfoil(end) = 0;
for m = 1:length(x)
   if m <= ind3 || m >= ind4
         y{m} = [linspace(y1,0,j_max/2) linspace(0,y2,j_max/2)];
   else
         y{m} = [linspace(y1,yl_airfoil(m-ind3+1),j_max/2) linspace(yu_airfoil(m-ind3+1),y2,j_max/2)];
   end   
end
X = meshgrid(x);
for m=1:length(y)
    for n = 1:j_max
        Y(n,m) = y{m}(n);
    end
end
XX = X';
YY = Y';
ETA1 = (XX-x1)/(x2-x1);
ETA2 = (YY(1,:)-y1)/(y2-y1);
ETA2 = meshgrid(ETA2);
deta1 = ETA1(2,1)-ETA1(1,1);
deta2 = ETA2(1,2)-ETA2(1,1);
x_eta1 = x2-x1;
x_eta2 = 0;
y_eta1 = 0;
y_eta2 = y2-y1;
J = x_eta1*y_eta2-y_eta1*x_eta2;
eta1_x = y_eta2/J;
eta2_x = -y_eta1/J;
eta1_y = -x_eta2/J;
eta2_y = x_eta1/J;
C11 = (x_eta2^2+y_eta2^2)/J;
C22 = (x_eta1^2+y_eta1^2)/J;
C12 = -(x_eta1*x_eta2+y_eta1*y_eta2)/J;
%% B.Cs
%%% N.B. dE = Ex*dx  +Ey*dy  = -V_inf * sind(Alpha) dx   + V_inf * cosd(Alpha) dy
E = zeros(i_max,j_max);
i_le = find(x == 0);                        % i corresponding to the L.E
i_te = find(x == Chord);                    % i corresponding to the T.E
j_airfoil = find(y{end} == 0);              % j corresponding to the Airfoil surface
for i = 2:i_max
    E(i,1) = E(i-1,1)+ -V_inf * sind(Alpha)*(x(i)-x(i-1));             % at j = 1
end
for j = 2:j_max
    E(1,j) = E(1,j-1)+ V_inf * cosd(Alpha)*(y2-y1)/j_max;              % at i = 1
    if j == j_airfoil(2) 
                    E(1,j) = E(1,j-1);
    end
end
for i = 2:i_max
    E(i,j_max) = E(i-1,j_max)+ -V_inf * sind(Alpha)*(x(i)-x(i-1));     % at j = j_max
end
for j = 2:j_max
    E(i_max,j) = E(i_max,j-1)+ V_inf * cosd(Alpha)*(y2-y1)/j_max;      % at i = i_max
    if j == j_airfoil(2) 
                    E(i_max,j) = E(i_max,j-1);
    end
end
% Initial Values (Linear Interpolation)
%----------------------------------------------------------%
for i = 2:i_max-1
    for j = 2:j_max-1
        E(i,j) = E(i,j-1)+((E(i,j_max)-E(i,1))/(j_max-1));
    end
end
%----------------------------------------------------------%
ind6 = find(y{end}==0);
E_airfoil = E(1,ind6(1));
E(i_le:i_te,j_airfoil) = E_airfoil;
%% PSOR Method
%%% N.B. For a spesific point(i,j) coordinate  (ETA1(j,i),ETA2(j,i))
%%% N.B. Make a cond. that same points have same Epsi (if there are repeated points)
Omega = 1;
cutta = 0;
E_old = E;
for n = 2:n_max             % Iteration Level
        E = E_old;
    for i = 2:i_max-1
        for j = 2:j_max-1
% constants Calculations
%------------------------------------------------------------------------%
r = deta1/deta2;%abs((XX(2,1)-XX(1,1))/(YY(1,2)-YY(1,1))); %
if j == j_airfoil(1)
%-------------%
XX_i1j2 = (XX(i+1,j+2)+XX(i+1,j))/2;
XX_i_1j2 = (XX(i-1,j+2)+XX(i-1,j))/2;
YY_i1j2 = (YY(i+1,j+2)+YY(i+1,j))/2;
YY_i_1j2 = (YY(i-1,j+2)+YY(i-1,j))/2;
XX_i1j_2 = (XX(i+1,j)+XX(i+1,j-1))/2;
XX_i_1j_2 = (XX(i-1,j)+XX(i-1,j-1))/2;
YY_i1j_2 = (YY(i+1,j)+YY(i+1,j-1))/2;
YY_i_1j_2 = (YY(i-1,j)+YY(i-1,j-1))/2;
%-------------%
x_eta1_ip = (XX(i+1,j)-XX(i,j))/deta1;
y_eta1_ip = (YY(i+1,j)-YY(i,j))/deta1;
x_eta1_im = (XX(i,j)-XX(i-1,j))/deta1;
y_eta1_im = (YY(i,j)-YY(i-1,j))/deta1;
x_eta1_jp = (XX_i1j2-XX_i_1j2)/2/deta1;
y_eta1_jp = (YY_i1j2-YY_i_1j2)/2/deta1;
x_eta1_jm = (XX_i1j_2-XX_i_1j_2)/2/deta1;
y_eta1_jm = (YY_i1j_2-YY_i_1j_2)/2/deta1;
%-------------%
XX_i2j1 = (XX(i+1,j+2)+XX(i,j+2))/2;
XX_i2j_1 = (XX(i+1,j-1)+XX(i,j-1))/2;
YY_i2j1 = (YY(i+1,j+2)+YY(i,j+2))/2;
YY_i2j_1 = (YY(i+1,j-1)+YY(i,j-1))/2;
XX_i_2j1 = (XX(i,j+2)+XX(i-1,j+2))/2;
XX_i_2j_1 = (XX(i,j-1)+XX(i-1,j-1))/2;
YY_i_2j1 = (YY(i,j+2)+YY(i-1,j+2))/2;
YY_i_2j_1 = (YY(i,j-1)+YY(i-1,j-1))/2;
%-------------%
x_eta2_ip = (XX_i2j1-XX_i2j_1)/2/deta2;
y_eta2_ip = (YY_i2j1-YY_i2j_1)/2/deta2;
x_eta2_im = (XX_i_2j1-XX_i_2j_1)/2/deta2;
y_eta2_im = (YY_i_2j1-YY_i_2j_1)/2/deta2;
x_eta2_jp = (XX(i,j+2)-XX(i,j))/deta2;
y_eta2_jp = (YY(i,j+2)-YY(i,j))/deta2;
x_eta2_jm = (XX(i,j)-XX(i,j-1))/deta2;
y_eta2_jm = (YY(i,j)-YY(i,j-1))/deta2;
elseif j == j_airfoil(2)
%-------------%
XX_i1j2 = (XX(i+1,j+1)+XX(i+1,j))/2;
XX_i_1j2 = (XX(i-1,j+1)+XX(i-1,j))/2;
YY_i1j2 = (YY(i+1,j+1)+YY(i+1,j))/2;
YY_i_1j2 = (YY(i-1,j+1)+YY(i-1,j))/2;
XX_i1j_2 = (XX(i+1,j)+XX(i+1,j-2))/2;
XX_i_1j_2 = (XX(i-1,j)+XX(i-1,j-2))/2;
YY_i1j_2 = (YY(i+1,j)+YY(i+1,j-2))/2;
YY_i_1j_2 = (YY(i-1,j)+YY(i-1,j-2))/2;
%-------------%
x_eta1_ip = (XX(i+1,j)-XX(i,j))/deta1;
y_eta1_ip = (YY(i+1,j)-YY(i,j))/deta1;
x_eta1_im = (XX(i,j)-XX(i-1,j))/deta1;
y_eta1_im = (YY(i,j)-YY(i-1,j))/deta1;
x_eta1_jp = (XX_i1j2-XX_i_1j2)/2/deta1;
y_eta1_jp = (YY_i1j2-YY_i_1j2)/2/deta1;
x_eta1_jm = (XX_i1j_2-XX_i_1j_2)/2/deta1;
y_eta1_jm = (YY_i1j_2-YY_i_1j_2)/2/deta1;
%-------------%
XX_i2j1 = (XX(i+1,j+1)+XX(i,j+1))/2;
XX_i2j_1 = (XX(i+1,j-2)+XX(i,j-2))/2;
YY_i2j1 = (YY(i+1,j+1)+YY(i,j+1))/2;
YY_i2j_1 = (YY(i+1,j-2)+YY(i,j-2))/2;
XX_i_2j1 = (XX(i,j+1)+XX(i-1,j+1))/2;
XX_i_2j_1 = (XX(i,j-2)+XX(i-1,j-2))/2;
YY_i_2j1 = (YY(i,j+1)+YY(i-1,j+1))/2;
YY_i_2j_1 = (YY(i,j-2)+YY(i-1,j-2))/2;
%-------------%
x_eta2_ip = (XX_i2j1-XX_i2j_1)/2/deta2;
y_eta2_ip = (YY_i2j1-YY_i2j_1)/2/deta2;
x_eta2_im = (XX_i_2j1-XX_i_2j_1)/2/deta2;
y_eta2_im = (YY_i_2j1-YY_i_2j_1)/2/deta2;
x_eta2_jp = (XX(i,j+1)-XX(i,j))/deta2;
y_eta2_jp = (YY(i,j+1)-YY(i,j))/deta2;
x_eta2_jm = (XX(i,j)-XX(i,j-2))/deta2;
y_eta2_jm = (YY(i,j)-YY(i,j-2))/deta2;
else
%-------------%
XX_i1j2 = (XX(i+1,j+1)+XX(i+1,j))/2;
XX_i_1j2 = (XX(i-1,j+1)+XX(i-1,j))/2;
YY_i1j2 = (YY(i+1,j+1)+YY(i+1,j))/2;
YY_i_1j2 = (YY(i-1,j+1)+YY(i-1,j))/2;
XX_i1j_2 = (XX(i+1,j)+XX(i+1,j-1))/2;
XX_i_1j_2 = (XX(i-1,j)+XX(i-1,j-1))/2;
YY_i1j_2 = (YY(i+1,j)+YY(i+1,j-1))/2;
YY_i_1j_2 = (YY(i-1,j)+YY(i-1,j-1))/2;
%-------------%
x_eta1_ip = (XX(i+1,j)-XX(i,j))/deta1;
y_eta1_ip = (YY(i+1,j)-YY(i,j))/deta1;
x_eta1_im = (XX(i,j)-XX(i-1,j))/deta1;
y_eta1_im = (YY(i,j)-YY(i-1,j))/deta1;
x_eta1_jp = (XX_i1j2-XX_i_1j2)/2/deta1;
y_eta1_jp = (YY_i1j2-YY_i_1j2)/2/deta1;
x_eta1_jm = (XX_i1j_2-XX_i_1j_2)/2/deta1;
y_eta1_jm = (YY_i1j_2-YY_i_1j_2)/2/deta1;
%-------------%
XX_i2j1 = (XX(i+1,j+1)+XX(i,j+1))/2;
XX_i2j_1 = (XX(i+1,j-1)+XX(i,j-1))/2;
YY_i2j1 = (YY(i+1,j+1)+YY(i,j+1))/2;
YY_i2j_1 = (YY(i+1,j-1)+YY(i,j-1))/2;
XX_i_2j1 = (XX(i,j+1)+XX(i-1,j+1))/2;
XX_i_2j_1 = (XX(i,j-1)+XX(i-1,j-1))/2;
YY_i_2j1 = (YY(i,j+1)+YY(i-1,j+1))/2;
YY_i_2j_1 = (YY(i,j-1)+YY(i-1,j-1))/2;
%-------------%
x_eta2_ip = (XX_i2j1-XX_i2j_1)/2/deta2;
y_eta2_ip = (YY_i2j1-YY_i2j_1)/2/deta2;
x_eta2_im = (XX_i_2j1-XX_i_2j_1)/2/deta2;
y_eta2_im = (YY_i_2j1-YY_i_2j_1)/2/deta2;
x_eta2_jp = (XX(i,j+1)-XX(i,j))/deta2;
y_eta2_jp = (YY(i,j+1)-YY(i,j))/deta2;
x_eta2_jm = (XX(i,j)-XX(i,j-1))/deta2;
y_eta2_jm = (YY(i,j)-YY(i,j-1))/deta2;
end
J_ip = x_eta1_ip*y_eta2_ip-y_eta1_ip*x_eta2_ip;
C11_ip = (x_eta2_ip^2+y_eta2_ip^2)/J_ip;
C22_ip = (x_eta1_ip^2+y_eta1_ip^2)/J_ip;
C12_ip = -(x_eta1_ip*x_eta2_ip+y_eta1_ip*y_eta2_ip)/J_ip;
J_im = x_eta1_im*y_eta2_im-y_eta1_im*x_eta2_im;
C11_im = (x_eta2_im^2+y_eta2_im^2)/J_im;
C22_im = (x_eta1_im^2+y_eta1_im^2)/J_im;
C12_im = -(x_eta1_im*x_eta2_im+y_eta1_im*y_eta2_im)/J_im;
J_jp = x_eta1_jp*y_eta2_jp-y_eta1_jp*x_eta2_jp;
C11_jp = (x_eta2_jp^2+y_eta2_jp^2)/J_jp;
C22_jp = (x_eta1_jp^2+y_eta1_jp^2)/J_jp;
C12_jp = -(x_eta1_jp*x_eta2_jp+y_eta1_jp*y_eta2_jp)/J_jp;
J_jm = x_eta1_jm*y_eta2_jm-y_eta1_jm*x_eta2_jm;
C11_jm = (x_eta2_jm^2+y_eta2_jm^2)/J_jm;
C22_jm = (x_eta1_jm^2+y_eta1_jm^2)/J_jm;
C12_jm = -(x_eta1_jm*x_eta2_jm+y_eta1_jm*y_eta2_jm)/J_jm;
Sij = C11_ip + C11_im + r^2 * (C22_jp + C22_jm);
Sim = C11_im - r^2 * (C12_jp - C12_jm) / 4;
Sip = C11_ip + r^2 * (C12_jp - C12_jm) / 4;
Sjm = r^2 * C22_jm - r * (C12_ip - C12_im) / 4;
Sjp = r^2 * C22_jp + r * (C12_ip - C12_im) / 4;
Smm = r * (C12_im + C12_jm) / 4;
Smp = -r * (C12_im + C12_jp) / 4;
Spm = -r * (C12_ip + C12_jm) / 4;
Spp = r * (C12_ip + C12_jp) / 4;
%------------------------------------------------------------------------%
% The following condition is to make the E on airfoil constant as a B.C.
%------------------------------------------------------------------------%
                if i >= i_le && i <= i_te && (j == j_airfoil(1) || j == j_airfoil(2))
                        E(i,j) = E_airfoil;
%------------------------------------------------------------------------% 
% The following condition is to omit the repeated line before the airfoil
%------------------------------------------------------------------------% 
                elseif j == j_airfoil(1) && (i<i_le || i>i_te)
                        E_telta = (Sim*E(i-1,j)+Sip*E(i+1,j)+Sjm*E(i,j-1)+Sjp*E(i,j+2)+Smm*E(i-1,j-1)+Smp*E(i-1,j+2)+Spm*E(i+1,j-1)+Spp*E(i+1,j+2))/Sij;
                        E(i,j) = E(i,j) + Omega*(E_telta-E(i,j));
                elseif j == j_airfoil(2) && (i<i_le || i>i_te) 
                        E(i,j) = E(i,j-1);
%------------------------------------------------------------------------%
                else
                        E_telta = (Sim*E(i-1,j)+Sip*E(i+1,j)+Sjm*E(i,j-1)+Sjp*E(i,j+1)+Smm*E(i-1,j-1)+Smp*E(i-1,j+1)+Spm*E(i+1,j-1)+Spp*E(i+1,j+1))/Sij;
                        E(i,j) = E(i,j) + Omega*(E_telta-E(i,j));
                end
         end               
    end
% Check on Cutta Condition
%------------------------------------------------%
% if Type == 'N'   || Type =='C' 
    if  abs(E(i_te+1,j_airfoil(1))-E_airfoil) <= error
        cutta = 1;
    else
        E_airfoil = E(i_te+1,j_airfoil(1));
        cutta = 0;
    end
% else
%     if  abs(E(i_te+1,j_airfoil(1)-1)-E_airfoil) <= error
%         cutta = 1;
%     else
%         E_airfoil = E(i_te+1,j_airfoil(1)-1);
%         cutta = 0;
%     end
% end
%------------------------------------------------%    
        Conv = sqrt((E - E_old).^2./i_max./j_max);            % Convergence Check
        Conv_max(n-1) = max(max(Conv));
        n_it = n;
        if Conv_max(n-1) <= error && cutta == 1             % psi is contant and cutta condition satisfied
            disp('Steady-State Solution Obtained')
            disp(['Max. Error = ', num2str(error)])
            disp(['Number of Iterations = ', num2str(n)])
            break
        elseif n == n_max
            disp('Solution is too far')
        end
        E_old = E;
end
%% Velocities
%%% Ex = u , Ey = -v
u = zeros(i_max,j_max);
v = zeros(i_max,j_max);
for j = 1:j_max
    for i = 1:i_max
% constants Calculations
%------------------------------------------------------------------------%
if i == 1
    if j == j_airfoil(1)
    x_eta1 = (-3*XX(i,j)+4 *XX(i+1,j)-XX(i+2,j))/2/deta1;
    y_eta1 = (-3*YY(i,j)+4 *YY(i+1,j)-YY(i+2,j))/2/deta1;
    x_eta2 = (XX(i,j+2)-XX(i,j-1))/2/deta2;
    y_eta2 = (YY(i,j+2)-YY(i,j-1))/2/deta2;
    elseif j == j_airfoil(2)
    x_eta1 = (-3*XX(i,j)+4 *XX(i+1,j)-XX(i+2,j))/2/deta1;
    y_eta1 = (-3*YY(i,j)+4 *YY(i+1,j)-YY(i+2,j))/2/deta1;
    x_eta2 = (XX(i,j+1)-XX(i,j-2))/2/deta2;
    y_eta2 = (YY(i,j+1)-YY(i,j-2))/2/deta2;
    elseif j == 1
    x_eta1 = (-3*XX(i,j)+4 *XX(i+1,j)-XX(i+2,j))/2/deta1;
    y_eta1 = (-3*YY(i,j)+4 *YY(i+1,j)-YY(i+2,j))/2/deta1;
    x_eta2 = (-3*XX(i,j)+4 *XX(i,j+1)-XX(i,j+2))/2/deta2;
    y_eta2 = (-3*YY(i,j)+4 *YY(i,j+1)-YY(i,j+2))/2/deta2;
    elseif j == j_max
    x_eta1 = (-3*XX(i,j)+4 *XX(i+1,j)-XX(i+2,j))/2/deta1;
    y_eta1 = (-3*YY(i,j)+4 *YY(i+1,j)-YY(i+2,j))/2/deta1;
    x_eta2 = (3*XX(i,j)-4 *XX(i,j-1)+XX(i,j-2))/2/deta2;
    y_eta2 = (3*YY(i,j)-4 *YY(i,j-1)+YY(i,j-2))/2/deta2;
    else
    x_eta1 = (-3*XX(i,j)+4 *XX(i+1,j)-XX(i+2,j))/2/deta1;
    y_eta1 = (-3*YY(i,j)+4 *YY(i+1,j)-YY(i+2,j))/2/deta1;
    x_eta2 = (XX(i,j+1)-XX(i,j-1))/2/deta2;
    y_eta2 = (YY(i,j+1)-YY(i,j-1))/2/deta2;
    end
elseif i == i_max
    if j == j_airfoil(1)
    x_eta1 = (3*XX(i,j)-4 *XX(i-1,j)+XX(i-2,j))/2/deta1;
    y_eta1 = (3*YY(i,j)-4 *YY(i-1,j)+YY(i-2,j))/2/deta1;
    x_eta2 = (XX(i,j+2)-XX(i,j-1))/2/deta2;
    y_eta2 = (YY(i,j+2)-YY(i,j-1))/2/deta2;
    elseif j == j_airfoil(2)
    x_eta1 = (3*XX(i,j)-4 *XX(i-1,j)+XX(i-2,j))/2/deta1;
    y_eta1 = (3*YY(i,j)-4 *YY(i-1,j)+YY(i-2,j))/2/deta1;
    x_eta2 = (XX(i,j+1)-XX(i,j-2))/2/deta2;
    y_eta2 = (YY(i,j+1)-YY(i,j-2))/2/deta2;
    elseif j == 1
    x_eta1 = (3*XX(i,j)-4 *XX(i-1,j)+XX(i-2,j))/2/deta1;
    y_eta1 = (3*YY(i,j)-4 *YY(i-1,j)+YY(i-2,j))/2/deta1;
    x_eta2 = (-3*XX(i,j)+4 *XX(i,j+1)-XX(i,j+2))/2/deta2;
    y_eta2 = (-3*YY(i,j)+4 *YY(i,j+1)-YY(i,j+2))/2/deta2;
    elseif j == j_max
    x_eta1 = (3*XX(i,j)-4 *XX(i-1,j)+XX(i-2,j))/2/deta1;
    y_eta1 = (3*YY(i,j)-4 *YY(i-1,j)+YY(i-2,j))/2/deta1;
    x_eta2 = (3*XX(i,j)-4 *XX(i,j-1)+XX(i,j-2))/2/deta2;
    y_eta2 = (3*YY(i,j)-4 *YY(i,j-1)+YY(i,j-2))/2/deta2;
    else
    x_eta1 = (3*XX(i,j)-4 *XX(i-1,j)+XX(i-2,j))/2/deta1;
    y_eta1 = (3*YY(i,j)-4 *YY(i-1,j)+YY(i-2,j))/2/deta1;
    x_eta2 = (XX(i,j+1)-XX(i,j-1))/2/deta2;
    y_eta2 = (YY(i,j+1)-YY(i,j-1))/2/deta2;
    end
elseif j == 1
    x_eta1 = (XX(i+1,j)-XX(i-1,j))/2/deta1;
    y_eta1 = (YY(i+1,j)-YY(i-1,j))/2/deta1;
    x_eta2 = (-3*XX(i,j)+4 *XX(i,j+1)-XX(i,j+2))/2/deta2;
    y_eta2 = (-3*YY(i,j)+4 *YY(i,j+1)-YY(i,j+2))/2/deta2;
elseif j == j_max
    x_eta1 = (XX(i+1,j)-XX(i-1,j))/2/deta1;
    y_eta1 = (YY(i+1,j)-YY(i-1,j))/2/deta1;
    x_eta2 = (3*XX(i,j)-4 *XX(i,j-1)+XX(i,j-2))/2/deta2;
    y_eta2 = (3*YY(i,j)-4 *YY(i,j-1)+YY(i,j-2))/2/deta2;
else
    if j == j_airfoil(1) && (i<i_le || i>i_te)
    x_eta1 = (XX(i+1,j)-XX(i-1,j))/2/deta1;
    y_eta1 = (YY(i+1,j)-YY(i-1,j))/2/deta1;
    x_eta2 = (XX(i,j+2)-XX(i,j-1))/2/deta2;
    y_eta2 = (YY(i,j+2)-YY(i,j-1))/2/deta2;
    elseif j == j_airfoil(2) && (i<i_le || i>i_te)
    x_eta1 = (XX(i+1,j)-XX(i-1,j))/2/deta1;
    y_eta1 = (YY(i+1,j)-YY(i-1,j))/2/deta1;
    x_eta2 = (XX(i,j+1)-XX(i,j-2))/2/deta2;
    y_eta2 = (YY(i,j+1)-YY(i,j-2))/2/deta2;
    elseif j == j_airfoil(1)
    x_eta1 = (XX(i+1,j)-XX(i-1,j))/2/deta1;
    y_eta1 = (YY(i+1,j)-YY(i-1,j))/2/deta1;
    x_eta2 = (3*XX(i,j)-4 *XX(i,j-1)+XX(i,j-2))/2/deta2;
    y_eta2 = (3*YY(i,j)-4 *YY(i,j-1)+YY(i,j-2))/2/deta2;
    elseif j == j_airfoil(2)
    x_eta1 = (XX(i+1,j)-XX(i-1,j))/2/deta1;
    y_eta1 = (YY(i+1,j)-YY(i-1,j))/2/deta1;
    x_eta2 = (-3*XX(i,j)+4 *XX(i,j+1)-XX(i,j+2))/2/deta2;
    y_eta2 = (-3*YY(i,j)+4 *YY(i,j+1)-YY(i,j+2))/2/deta2;
    else
    x_eta1 = (XX(i+1,j)-XX(i-1,j))/2/deta1;
    y_eta1 = (YY(i+1,j)-YY(i-1,j))/2/deta1;
    x_eta2 = (XX(i,j+1)-XX(i,j-1))/2/deta2;
    y_eta2 = (YY(i,j+1)-YY(i,j-1))/2/deta2;
    end
end
J = x_eta1*y_eta2-y_eta1*x_eta2;
eta1_x = y_eta2/J;
eta2_x = -y_eta1/J;
eta1_y = -x_eta2/J;
eta2_y = x_eta1/J;
%------------------------------------------------------------------------%
        if i > 1 && i ~= i_max && j > 1 && j ~= j_max && (i<i_le || i>i_te)
            if j == j_airfoil(1) 
            u(i,j) = (E(i+1,j)-E(i-1,j))/2/deta1*eta1_y + (E(i,j+2)-E(i,j-1))/2/deta2*eta2_y;
            v(i,j) = -(E(i+1,j)-E(i-1,j))/2/deta1*eta1_x - (E(i,j+2)-E(i,j-1))/2/deta2*eta2_x;
            elseif j == j_airfoil(2)
            u(i,j) = (E(i+1,j)-E(i-1,j))/2/deta1*eta1_y + (E(i,j+1)-E(i,j-2))/2/deta2*eta2_y;
            v(i,j) = -(E(i+1,j)-E(i-1,j))/2/deta1*eta1_x - (E(i,j+1)-E(i,j-2))/2/deta2*eta2_x;
            else
            u(i,j) = (E(i+1,j)-E(i-1,j))/2/deta1*eta1_y + (E(i,j+1)-E(i,j-1))/2/deta2*eta2_y;
            v(i,j) = -(E(i+1,j)-E(i-1,j))/2/deta1*eta1_x - (E(i,j+1)-E(i,j-1))/2/deta2*eta2_x;
            end
        elseif i == 1 && j > 1 && j ~= j_max
            if j == j_airfoil(1) 
            u(i,j) = (-3*E(i,j)+4 *E(i+1,j)-E(i+2,j))/2/deta1*eta1_y + (E(i,j+2)-E(i,j-1))/2/deta2*eta2_y;
            v(i,j) = -(-3*E(i,j)+4 *E(i+1,j)-E(i+2,j))/2/deta1*eta1_x - (E(i,j+2)-E(i,j-1))/2/deta2*eta2_x;
            elseif j == j_airfoil(2) 
            u(i,j) = (-3*E(i,j)+4 *E(i+1,j)-E(i+2,j))/2/deta1*eta1_y + (E(i,j+1)-E(i,j-2))/2/deta2*eta2_y;
            v(i,j) = -(-3*E(i,j)+4 *E(i+1,j)-E(i+2,j))/2/deta1*eta1_x - (E(i,j+1)-E(i,j-2))/2/deta2*eta2_x;
            else
            u(i,j) = (-3*E(i,j)+4 *E(i+1,j)-E(i+2,j))/2/deta1*eta1_y + (E(i,j+1)-E(i,j-1))/2/deta2*eta2_y;
            v(i,j) = -(-3*E(i,j)+4 *E(i+1,j)-E(i+2,j))/2/deta1*eta1_x - (E(i,j+1)-E(i,j-1))/2/deta2*eta2_x;    
            end
        elseif i > 1 && i ~= i_max && j == 1 
            u(i,j) = (E(i+1,j)-E(i-1,j))/2/deta1*eta1_y + (-3*E(i,j)+4 *E(i,j+1)-E(i,j+2))/2/deta2*eta2_y;
            v(i,j) = -(E(i+1,j)-E(i-1,j))/2/deta1*eta1_x - (-3*E(i,j)+4 *E(i,j+1)-E(i,j+2))/2/deta2*eta2_x;
        elseif i == i_max && j > 1 && j ~= j_max
            if j == j_airfoil(1)
            u(i,j) = (3*E(i,j)-4 *E(i-1,j)+E(i-2,j))/2/deta1*eta1_y + (E(i,j+2)-E(i,j-1))/2/deta2*eta2_y;
            v(i,j) = -(3*E(i,j)-4 *E(i-1,j)+E(i-2,j))/2/deta1*eta1_x - (E(i,j+2)-E(i,j-1))/2/deta2*eta2_x;
            elseif j == j_airfoil(2)
            u(i,j) = (3*E(i,j)-4 *E(i-1,j)+E(i-2,j))/2/deta1*eta1_y + (E(i,j+1)-E(i,j-2))/2/deta2*eta2_y;
            v(i,j) = -(3*E(i,j)-4 *E(i-1,j)+E(i-2,j))/2/deta1*eta1_x - (E(i,j+1)-E(i,j-2))/2/deta2*eta2_x;
            else
            u(i,j) = (3*E(i,j)-4 *E(i-1,j)+E(i-2,j))/2/deta1*eta1_y + (E(i,j+1)-E(i,j-1))/2/deta2*eta2_y;
            v(i,j) = -(3*E(i,j)-4 *E(i-1,j)+E(i-2,j))/2/deta1*eta1_x - (E(i,j+1)-E(i,j-1))/2/deta2*eta2_x;    
            end
        elseif i > 1 && i ~= i_max && j == j_max
            u(i,j) = (E(i+1,j)-E(i-1,j))/2/deta1*eta1_y + (3*E(i,j)-4 *E(i,j-1)+E(i,j-2))/2/deta2*eta2_y;
            v(i,j) = -(E(i+1,j)-E(i-1,j))/2/deta1*eta1_x - (3*E(i,j)-4 *E(i,j-1)+E(i,j-2))/2/deta2*eta2_x;
        elseif i == 1 && j == 1
            u(i,j) = (-3*E(i,j)+4 *E(i+1,j)-E(i+2,j))/2/deta1*eta1_y + (-3*E(i,j)+4 *E(i,j+1)-E(i,j+2))/2/deta2*eta2_y;
            v(i,j) = -(-3*E(i,j)+4 *E(i+1,j)-E(i+2,j))/2/deta1*eta1_x - (-3*E(i,j)+4 *E(i,j+1)-E(i,j+2))/2/deta2*eta2_x;
        elseif i == 1 && j == j_max
            u(i,j) = (-3*E(i,j)+4 *E(i+1,j)-E(i+2,j))/2/deta1*eta1_y + (3*E(i,j)-4 *E(i,j-1)+E(i,j-2))/2/deta2*eta2_y;
            v(i,j) = -(-3*E(i,j)+4 *E(i+1,j)-E(i+2,j))/2/deta1*eta1_x - (3*E(i,j)-4 *E(i,j-1)+E(i,j-2))/2/deta2*eta2_x;
        elseif i == i_max && j == 1
            u(i,j) = (3*E(i,j)-4 *E(i-1,j)+E(i-2,j))/2/deta1*eta1_y + (-3*E(i,j)+4 *E(i,j+1)-E(i,j+2))/2/deta2*eta2_y;
            v(i,j) = -(3*E(i,j)-4 *E(i-1,j)+E(i-2,j))/2/deta1*eta1_x - (-3*E(i,j)+4 *E(i,j+1)-E(i,j+2))/2/deta2*eta2_x;
        elseif i == i_max && j == j_max
            u(i,j) = (3*E(i,j)-4 *E(i-1,j)+E(i-2,j))/2/deta1*eta1_y + (3*E(i,j)-4 *E(i,j-1)+E(i,j-2))/2/deta2*eta2_y;
            v(i,j) = -(3*E(i,j)-4 *E(i-1,j)+E(i-2,j))/2/deta1*eta1_x - (3*E(i,j)-4 *E(i,j-1)+E(i,j-2))/2/deta2*eta2_x;
        elseif j == j_airfoil(2) && (i>i_le && i<i_te) % Upper Surface
            u(i,j) = (E(i+1,j)-E(i-1,j))/2/deta1*eta1_y + (-3*E(i,j)+4 *E(i,j+1)-E(i,j+2))/2/deta2*eta2_y;
            v(i,j) = -(E(i+1,j)-E(i-1,j))/2/deta1*eta1_x - (-3*E(i,j)+4 *E(i,j+1)-E(i,j+2))/2/deta2*eta2_x;
        elseif j == j_airfoil(1) && (i>i_le && i<i_te) % lower Surface
            u(i,j) = (E(i+1,j)-E(i-1,j))/2/deta1*eta1_y + (3*E(i,j)-4 *E(i,j-1)+E(i,j-2))/2/deta2*eta2_y;
            v(i,j) = -(E(i+1,j)-E(i-1,j))/2/deta1*eta1_x - (3*E(i,j)-4 *E(i,j-1)+E(i,j-2))/2/deta2*eta2_x;
        elseif j == j_airfoil(2) && i == i_le % L.E
            u(i,j) = (3*E(i,j)-4 *E(i-1,j)+E(i-2,j))/2/deta1*eta1_y + (E(i,j+1)-E(i,j-2))/2/deta2*eta2_y;
            v(i,j) = -(3*E(i,j)-4 *E(i-1,j)+E(i-2,j))/2/deta1*eta1_x - (E(i,j+1)-E(i,j-2))/2/deta2*eta2_x;
        elseif j == j_airfoil(1) && i == i_le % L.E
            u(i,j) = (3*E(i,j)-4 *E(i-1,j)+E(i-2,j))/2/deta1*eta1_y + (E(i,j+2)-E(i,j-1))/2/deta2*eta2_y;
            v(i,j) = -(3*E(i,j)-4 *E(i-1,j)+E(i-2,j))/2/deta1*eta1_x - (E(i,j+2)-E(i,j-1))/2/deta2*eta2_x;
        elseif j == j_airfoil(2) && i == i_te % T.E
            u(i,j) = (-3*E(i,j)+4 *E(i+1,j)-E(i+2,j))/2/deta1*eta1_y + (E(i,j+1)-E(i,j-2))/2/deta2*eta2_y;
            v(i,j) = -(-3*E(i,j)+4 *E(i+1,j)-E(i+2,j))/2/deta1*eta1_x - (E(i,j+1)-E(i,j-2))/2/deta2*eta2_x;
        elseif j == j_airfoil(1) && i == i_te % T.E
            u(i,j) = (-3*E(i,j)+4 *E(i+1,j)-E(i+2,j))/2/deta1*eta1_y + (E(i,j+2)-E(i,j-1))/2/deta2*eta2_y;
            v(i,j) = -(-3*E(i,j)+4 *E(i+1,j)-E(i+2,j))/2/deta1*eta1_x - (E(i,j+2)-E(i,j-1))/2/deta2*eta2_x;
        else
            u(i,j) = (E(i+1,j)-E(i-1,j))/2/deta1*eta1_y + (E(i,j+1)-E(i,j-1))/2/deta2*eta2_y;
            v(i,j) = -(E(i+1,j)-E(i-1,j))/2/deta1*eta1_x - (E(i,j+1)-E(i,j-1))/2/deta2*eta2_x;
        end
    end
end
V = sqrt(u.^2+v.^2);
V_airfoil_l = V(i_le:i_te,j_airfoil(1));
V_airfoil_u = V(i_le:i_te,j_airfoil(2));
%% Pressure Coefficient
Cp = 1-(V./V_inf).^2;
Cp_u = 1-(V_airfoil_u./V_inf).^2;
Cp_l = 1-(V_airfoil_l./V_inf).^2;
%% Plots
figure      % Airfoil
hold on
grid on
axis equal
plot(0:0.0001:Chord,polyval(polyfit(x_upper,y_upper,10),0:0.0001:Chord),0:0.0001:Chord,polyval(polyfit(x_lower,y_lower,10),0:0.0001:Chord),'LineWidth',1.5,'color','b')
plot(0:0.0001:Chord,(polyval(polyfit(x_upper,y_upper,10),0:0.0001:Chord)+polyval(polyfit(x_lower,y_lower,10),0:0.0001:Chord))/2,'--','LineWidth',1.5,'color','r')
if Type == 'N'
title(['NACA ',  num2str(NACA(1)) num2str(NACA(2)) num2str(NACA(3)) num2str(NACA(4)),'  Airfoil'])
elseif Type == 'J'
title(['Flow Stream Lines around  Joukowski Airfoil:    t_{max}/c = ',  num2str(t_c*100),'%    C_{max}/c = ' ,num2str(C_c*100),'%'])    
end
xlabel('x')
ylabel('y')
legend('Upper Surface','Lower Surface','Camber Line')
figure      % Physical Grid
hold on
grid on
axis equal
plot(XX,YY)
for m =1:length(X)          %Virtical Lines
    if m > ind3 && m < ind4
    ind = find(Y(:,m) == yl_airfoil(m-ind3+1));
    plot(X(1:ind,m),Y(1:ind,m))            
    plot(X(ind+1:end,m),Y(ind+1:end,m))
    else
    plot(X(:,m),Y(:,m))            
    end
end
title('Physical Grid')
xlabel('x')
ylabel('y')
figure      % Computitional Grid
hold on
grid on
axis equal
plot(ETA1,ETA2,'color','b')
plot(ETA2,ETA1,'color','b')
title('Computitional Grid')
xlabel('\eta_1')
ylabel('\eta_2')
figure      % Convergance
hold on
grid on
% axis equal
plot(2:1:n_it,log10(Conv_max))
title('Convergance of \psi')
xlabel('Num. of Iteratios (n)')
ylabel('1og_{10}(Max. RMS)')
figure      % Stream Lines
hold on
grid on
axis equal
contour(XX,YY,E,50)
fill([x_upper x_lower],[y_upper y_lower],'b')
if Type == 'N'
title(['Flow Stream Lines around  NACA ',  num2str(NACA(1)) num2str(NACA(2)) num2str(NACA(3)) num2str(NACA(4)),' , \alpha = ',num2str(Alpha),'^o'])
elseif Type == 'J'
title(['Flow Stream Lines around  Joukowski Airfoil:    t_{max}/c = ',  num2str(t_c*100),'%    C_{max}/c = ' ,num2str(C_c*100),'%'])    
end
figure      % Vector Field
hold on
grid on
axis equal
quiver(XX,YY,u,v)
plot(x_upper,y_upper,'b')
plot(x_lower,y_lower,'b')
if Type == 'N'
title(['Vector Field around  NACA ',  num2str(NACA(1)) num2str(NACA(2)) num2str(NACA(3)) num2str(NACA(4)),' , \alpha = ',num2str(Alpha),'^o'])
else
end
figure      % Velocity Distribution on airfoil
hold on
grid on
plot(x_upper,y_upper,x_lower,y_lower,'LineWidth',1.5,'color','b')
plot(xx,V_airfoil_l/V_inf,'LineWidth',1.5,'color','g')
plot(xx,V_airfoil_u/V_inf,'LineWidth',1.5,'color','r')
title('Velocity Distribution on the Airfoil Surface')
legend('Airfoil surface','Airfoil surface','Lower Surface','Upper Surface')
figure      % Pressure Coefficient Distribution on airfoil
hold on
grid on
plot(x_upper,y_upper,x_lower,y_lower,'LineWidth',1.5,'color','b')
plot(xx,Cp_l,'LineWidth',1.5,'color','g')
plot(xx,Cp_u,'LineWidth',1.5,'color','r')
title('Pressure Coefficient Distribution on the Airfoil Surface')
legend('Airfoil surface','Airfoil surface','Lower Surface','Upper Surface')
figure      % Velocity contures Lines
hold on
grid on
axis equal
contour(XX,YY,V./V_inf,50)
fill([x_upper x_lower],[y_upper y_lower],'b')
if Type == 'N'
title(['Velocity Conture Lines around  NACA ',  num2str(NACA(1)) num2str(NACA(2)) num2str(NACA(3)) num2str(NACA(4)),' , \alpha = ',num2str(Alpha),'^o'])
elseif Type == 'J'
title(['Velocity Conture Lines around  Joukowski Airfoil:    t_{max}/c = ',  num2str(t_c*100),'%    C_{max}/c = ' ,num2str(C_c*100),'%'])    
end
figure
surf(XX,YY,V./V_inf);
view(2)
shading interp
colormap jet
colorbar
hold on
fill([x_upper x_lower],[y_upper y_lower],'b')
if Type == 'N'
title(['Velocity around  NACA ',  num2str(NACA(1)) num2str(NACA(2)) num2str(NACA(3)) num2str(NACA(4)),' , \alpha = ',num2str(Alpha),'^o'])
elseif Type == 'J'
title(['Velocity around  Joukowski Airfoil:    t_{max}/c = ',  num2str(t_c*100),'%    C_{max}/c = ' ,num2str(C_c*100),'%'])    
end
figure      % Pressure contures Lines
hold on
grid on
axis equal
contour(XX,YY,Cp,50)
fill([x_upper x_lower],[y_upper y_lower],'b')
if Type == 'N'
title(['Pressure Conture Lines around  NACA ',  num2str(NACA(1)) num2str(NACA(2)) num2str(NACA(3)) num2str(NACA(4)),' , \alpha = ',num2str(Alpha),'^o'])
elseif Type == 'J'
title(['Pressure Conture Lines around  Joukowski Airfoil:    t_{max}/c = ',  num2str(t_c*100),'%    C_{max}/c = ' ,num2str(C_c*100),'%'])    
end
figure
surf(XX,YY,Cp);
view(2)
shading interp
colormap jet
colorbar
hold on
fill([x_upper x_lower],[y_upper y_lower],'k')
if Type == 'N'
title(['Pressure around  NACA ',  num2str(NACA(1)) num2str(NACA(2)) num2str(NACA(3)) num2str(NACA(4)),' , \alpha = ',num2str(Alpha),'^o'])
elseif Type == 'J'
title(['Pressure around around  Joukowski Airfoil:    t_{max}/c = ',  num2str(t_c*100),'%    C_{max}/c = ' ,num2str(C_c*100),'%'])    
end