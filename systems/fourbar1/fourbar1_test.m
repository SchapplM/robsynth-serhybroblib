% Teste Klassendefinition für Viergelenkkette in unterschiedlichen
% Implementierungen

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

clear
clc

repopath = fileparts(which('hybrid_mdl_test_path_init.m'));


%% Definition der Roboterklassen
mdlprefix = 'fourbar1';
for variante2 = {'DE1', 'DE2', 'TE'}
  varkn = variante2{1};
  addpath(fullfile(repopath,'systems','fourbar1',sprintf('matlabfcn_fourbar1%s',varkn)));
  addpath(fullfile(repopath,'systems','fourbar1',sprintf('testfcn_fourbar1%s',varkn)));
  eval(sprintf('TSS = fourbar1%s_varpar_testfunctions_parameter();', varkn));
  Par_struct = struct('beta', TSS.beta, 'b', TSS.b, ...
                      'alpha', TSS.alpha, 'a', TSS.a, ...
                      'theta', TSS.theta, 'd', TSS.d, ...
                      'sigma', TSS.sigma, 'offset', TSS.q_offset, ...
                      'v', TSS.v, ...
                      'pkin', TSS.pkin, ...
                      'mu', TSS.mu, ...
                      'm', TSS.m, 'mrSges', TSS.mrSges, 'Ifges', TSS.Ifges, ...
                      'NJ', TSS.NJ, 'NL', TSS.NL, 'NQJ', TSS.NQJ);   
  eval(sprintf('Par_struct_%s = Par_struct;', varkn));
end

RS_DE1 = SerRob(Par_struct_DE1, sprintf('%s%s',mdlprefix,'DE1'));
RS_DE1 = RS_DE1.fill_fcn_handles(false);
RS_DE2 = SerRob(Par_struct_DE2, sprintf('%s%s',mdlprefix,'DE2'));
RS_DE2 = RS_DE2.fill_fcn_handles(false);
RS_TE  = SerRob(Par_struct_TE, sprintf('%s%s',mdlprefix,'TE'));
RS_TE = RS_TE.fill_fcn_handles(false);

%% Vergleich der Implementierungen
for q = linspace(0, 2*pi, 1000)
  T_DE1 = RS_DE1.fkine(q);
  T_DE2 = RS_DE2.fkine(q);
  T_TE = RS_TE.fkine(q);
  test1 = T_DE1-T_DE2;
  test2 = T_DE1-T_TE;
  if any(abs([test1(:); test2(:)]) > 1e-10)
    error('Methoden DE1, DE2 und TE stimmen nicht überein');
  end
end
fprintf('Kinematik der Viergelenkkette für 1000 Kombinationen in unterschiedlichen Implementierungen getestet\n');

%% Gelenk-Trajektorie mit einem Umlauf 
QE = [0, 2*pi]';

[Q,QD,QDD,t] = traj_trapez2_multipoint(QE, 1, 1e-1, 1e-2, 1e-3, 0.25);
X = NaN(size(Q, 1), 6);
for ii = 1:size(Q, 1)
  T_0_Ei = RS_DE1.fkineEE(Q(ii,:)');
  X(ii,:) = RS_DE1.t2x(T_0_Ei);
end

%% Gelenkmomentenverlauf berechnen

nt = size(Q,1);
TAU = NaN(nt, RS_DE1.NQJ);
for i = 1:nt
  q_i = Q(i,:)';
  qD_i = QD(i,:)';
  qDD_i = QDD(i,:)';
  
  tau_i = RS_DE1.invdyn(q_i, qD_i, qDD_i);
  TAU(i,:) = tau_i;
end

%% Plotten
figure(1);clf;
subplot(4,1,1);
plot(t, 180/pi*Q);
xlabel('t [s]');
ylabel('q [deg]');
grid on;
subplot(4,1,2);
plot(t, 180/pi*QD);
xlabel('t [s]');
ylabel('qD [deg/s]');
grid on;
subplot(4,1,3);
plot(t, 180/pi*QDD);
xlabel('t [s]');
ylabel('qDD [deg/s^2]');
grid on;
subplot(4,1,4);
plot(t, TAU);
xlabel('t [s]');
ylabel('tau [Nm]');
grid on;


s_plot = struct( 'ks', [1:RS_DE1.NJ, RS_DE1.NJ+2], 'straight', 0);
q = 0;
figure(2);clf;
hold on;
grid on;
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
view(3);
RS_DE1.plot( q, s_plot );


s_anim = struct( 'gif_name', '');
figure(3);clf;
hold on;
plot3(X(:,1), X(:,2), X(:,3));
grid on;
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
view(3);
RS_DE1.anim( Q(1:50:end,:), s_anim, s_plot);