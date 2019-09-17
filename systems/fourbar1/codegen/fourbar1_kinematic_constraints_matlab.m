% Verwende Cosinussatz für Minimalwinkel
l1=4.5; l2=2; l3=3; l4=4;

la = l2+l3;
lb = l1;
lc = l4;
q_min = [acos((la^2+lb^2-lc^2)/(2*la*lb))];
q_max = [pi];
