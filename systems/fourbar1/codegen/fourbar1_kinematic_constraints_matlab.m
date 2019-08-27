% Verwende Cosinussatz für Minimalwinkel
l1=4.5; l2=2; l3=3; l4=4;

a = l2+l3;
b = l1;
c = l4;
q_min = [acos((a^2+b^2-c^2)/(2*a*b))];
q_max = [pi];
