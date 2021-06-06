function [x_upper,y_upper,x_lower,y_lower] = cylinder(Chord)
x = 0:0.00001:Chord;
a = Chord/2;
% (x-a)^2+y^2=a^2
y_upper = sqrt(a^2-(x-a).^2);
y_lower = -sqrt(a^2-(x-a).^2);
x_upper = x;
x_lower = x;
end