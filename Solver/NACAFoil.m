function [x_upper,y_upper,x_lower,y_lower] = NACAFoil(NACA,Chord)
x = 0:0.00001:Chord;
if length(NACA) == 4
    disp(['NACA 4 Digit Series:  NACA ',  num2str(NACA(1)) num2str(NACA(2)) num2str(NACA(3)) num2str(NACA(4))])
    if NACA(1) == 0 && NACA(2) == 0
        Symm = 1;
        disp('Symmetric Airfoil')
    else
        Symm = 0;
        disp('Cambered Airfoil')
    end
end
if Symm == 1
    t = str2num([num2str(NACA(3)),num2str(NACA(4))])/100;
    y_upper = 5*t*Chord*(0.2969*sqrt(x/Chord)-0.126*(x/Chord)-0.3516*(x/Chord).^2+0.2843*(x/Chord).^3-0.1015*(x/Chord).^4); 
    y_lower = -y_upper;
    x_upper = x;
    x_lower = x;
else
    m = NACA(1)/100;
    p = NACA(2)*Chord/10;
    t = str2num([num2str(NACA(3)),num2str(NACA(4))])/100;
    for i = 1:length(x)
        if x(i)/Chord<=p
            y_camber(i) = m*x(i)/p^2*(2*p-x(i)/Chord);
            dy_camber(i) = 2*m/p^2*(p-x(i)/Chord);
        else
            y_camber(i) = m*(Chord-x(i))/(1-p)^2*(1+x(i)/Chord-2*p);
            dy_camber(i) = 2*m/(1-p)^2*(p-x(i)/Chord);
        end
    end
    y_t = 5*t*Chord*(0.2969*sqrt(x/Chord)-0.126*(x/Chord)-0.3516*(x/Chord).^2+0.2843*(x/Chord).^3-0.1015*(x/Chord).^4); 
    theta = atan(dy_camber);
    x_upper = x-y_t.*sin(theta);
    x_lower = x+y_t.*sin(theta);
    y_upper = y_camber+y_t.*cos(theta);
    y_lower = y_camber-y_t.*cos(theta);
    y_upper(end) = 0;           % Closed T.E
    y_lower(end) = 0;           % Closed T.E
end
end