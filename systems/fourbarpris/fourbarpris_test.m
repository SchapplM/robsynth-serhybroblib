% Teste Klassendefinition für Viergelenkkette mit Schubgelenk
% Dieser Mechanismus entspricht einem Robotergelenk mit angetriebenem
% Schubgelenk zur Umsetzung in einer Drehbewegung.
% Implementierungen
%
% Verfügbare Beispiele: Siehe Tabelle models.csv

% Quelle:
% https://de.wikipedia.org/wiki/Schubkurbel
% (Hinweis: Die Kinematik entspricht nicht direkt der üblichen Schubkurbel)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2020-04
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clear
clc

%% Definition der Roboterklassen
RS_TE = hybroblib_create_robot_class('fourbarpris', 'TE', 'fourbarprisBsp1');
RS_DE1 = hybroblib_create_robot_class('fourbarpris', 'DE1', 'fourbarprisBsp1');
RS_DE2 = hybroblib_create_robot_class('fourbarpris', 'DE2', 'fourbarprisBsp1');

%% Vergleich der Implementierungen

TSS = RS_TE.gen_testsettings();
I_io = false(size(TSS.Q,1),1);
for i = 1:size(TSS.Q,1)
  q = TSS.Q(i);
  try
    T_DE1 = RS_DE1.fkine(q);
    T_DE2 = RS_DE2.fkine(q);
    T_TE = RS_TE.fkine(q);
  catch
    continue
  end
  I_io(i) = true;
  test1 = T_DE1-T_DE2;
  test2 = T_DE1-T_TE;
  if any(abs([test1(:); test2(:)]) > 1e-10)
    error('Methoden DE1, DE2 und TE stimmen nicht überein');
  end
end
fprintf('Kinematik der Schubkurbel für %d Kombinationen in unterschiedlichen Implementierungen getestet. %d nicht möglich\n', ...
  sum(I_io), size(TSS.Q,1)-sum(I_io));

% Ab jetzt nur noch mit einem Modell weiterrechnen und dieses Parametrieren
RS = copy(RS_DE1);

figure(10); hold on;
plot(TSS.Q(I_io), 1, 'gv');
if sum(~I_io) > 0, plot(TSS.Q(~I_io), 1, 'rx'); end
title('Übersicht über zulässige Antriebspositionen');
xlabel('q1 in m'); ylabel('');
%% End-Effektor definieren
% Liegt auf beweglicher Schwinge
l = RS.pkin(1);
r_N_E = [3/4*l;-l/3;0];
RS.update_EE(r_N_E);


%% Gelenk-Trajektorie mit einem Hub
QE = RS.qlim';

[Q,QD,QDD,t] = traj_trapez2_multipoint(QE, 1, 1e-1, 1e-2, 1e-3, 0.25);
X = NaN(size(Q, 1), 6);
for ii = 1:size(Q, 1)
  T_0_Ei = RS.fkineEE(Q(ii,:)');
  X(ii,:) = RS.t2x(T_0_Ei);
end

%% Plotten
s_plot = struct( 'ks', [1:RS.NJ, RS.NJ+2], 'straight', 0);
q = mean(RS.qlim);
figure(2);clf;
hold on;
grid on;
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
view(3);
title(RS.descr);
RS.plot( q, s_plot );
RS.jointvar(q)

s_anim = struct( 'gif_name', '');
figure(3);clf;
hold on;
plot3(X(:,1), X(:,2), X(:,3));
grid on;
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
view(3);
title(RS.descr);
RS.anim( Q(1:50:end,:), [], s_anim, s_plot);