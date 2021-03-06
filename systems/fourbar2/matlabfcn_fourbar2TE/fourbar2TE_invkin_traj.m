% Inverse Kinematik für komplette Trajektorie für
% fourbar2TE
% Allgemeine, stark parametrierbare Funktion zum Aufruf mit allen möglichen
% Einstellungen
% Iterative Lösung der inversen Kinematik mit inverser Jacobi-Matrix
% Zusätzlich Nutzung der differentiellen Kinematik für schnellere Konvergenz
% (Es wird vorausgesetzt und ausgenutzt, dass die EE-Trajektorie stetig ist)
%
% Eingabe:
% X
%   Trajektorie von EE-Lagen (Sollwerte)
% XD
%   Trajektorie von EE-Geschwindigkeiten (Sollwerte)
%   (Die Orientierung wird durch Euler-Winkel-Zeitableitung dargestellt)
% XDD
%   Trajektorie von EE-Beschleunigungen (Sollwerte)
%   Orientierung bezogen auf Euler-Winkel
% PHI
%   Kinematische Zwangsbedingungen über die Trajektorie
% T
%   Zeitbasis der Trajektorie (Alle Zeit-Stützstellen)
% q0
%   Anfangs-Gelenkwinkel für Algorithmus
% s
%   Struktur mit Eingabedaten. Felder, siehe Quelltext dieser Funktion und
%   von `fourbar2TE_invkin_eulangresidual`
%
% Ausgabe:
% Q
%   Trajektorie von Gelenkpositionen (Lösung der IK)
% QD
%   Trajektorie von Gelenkgeschwindigkeiten
% QDD
%   Trajektorie von Gelenkbeschleunigungen
%
% Siehe auch: ParRob/invkin_traj

% TODO: Nullraumbewegung mit Nebenbedingung
% TODO: Erfolg der IK prüfen
% TODO: EE-Trajektorie auch als Winkelgeschwindigkeit und zusätzlicher
%       Schalter für Umrechnung auf analytische Jacobi-Matrix
% TODO: Schalter für Weglassen der Beschleunigung

% Quelle:
% [1] Aufzeichnungen Schappler vom 28.11.2018
% [2] Aufzeichnungen Schappler vom 11.12.2018

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-11 13:28
% Revision: cdfd09ebe5e500df1cf15b0f4129ff9f592109b7 (2019-10-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-02
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Q, QD, QDD, PHI] = fourbar2TE_invkin_traj(X, XD, XDD, T, q0, s)

%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,6]),coder.newtype('double',[inf,6]),
%$cgargs  coder.newtype('double',[inf,6]),coder.newtype('double',[inf,1]),
%$cgargs zeros(1,1), struct(
%$cgargs            'pkin', zeros(2,1),
%$cgargs          'sigmaJ', zeros(1,1),
%$cgargs             'NQJ', 0,
%$cgargs            'qlim', zeros(1,2),
%$cgargs            'I_EE', true(1,6),
%$cgargs     'phiconv_W_E', uint8(2),
%$cgargs        'I_EElink', uint8(0),
%$cgargs            'reci', true,
%$cgargs           'T_N_E', zeros(4,4),
%$cgargs               'K', zeros(1,1),
%$cgargs              'Kn', zeros(1,1),
%$cgargs              'wn', zeros(2,1),
%$cgargs       'scale_lim', 0,
%$cgargs      'maxrelstep', 0.1,
%$cgargs       'normalize', false,
%$cgargs           'n_min', 0,
%$cgargs           'n_max', 1000,
%$cgargs        'Phit_tol', 1.0000e-10,
%$cgargs        'Phir_tol', 1.0000e-10,
%$cgargs     'retry_limit', 100)}

%% Initialisierung
% Vorbelegung der Ausgabe
Q = NaN(length(T), s.NQJ);
QD = Q;
QDD = Q;
PHI = NaN( length(T), sum(s.I_EE) );

% Einstellungsvariablen aus Struktur herausholen
I_EE = s.I_EE;
nt = length(T);
qk0 = q0;
link_index = s.I_EElink;
pkin = s.pkin;
r_i_i_C = s.T_N_E(1:3,4);


%% Iterative Berechnung der gesamten Trajektorie
for k = 1:nt
  %% Gelenk-Position berechnen
  % Inverse Kinematik für aktuellen Bahnpunkt. Nutze Anfangswert aus der
  % differentiellen Kinematik hiernach von der letzten Iteration (k-1)
  [q_k, Phi_k] = fourbar2TE_invkin_eulangresidual(X(k,:)', qk0, s);

  %% Gelenk-Geschwindigkeit berechnen (Siehe [1]).
  % Geometrische Jacobi-Matrix in analytische Jacobi umrechnen
  Jg = fourbar2TE_jacobig_sym_varpar(q_k, link_index, r_i_i_C, pkin);
  Tw = euljac(X(k,4:6)', s.phiconv_W_E);
  J_x = [Jg(1:3,:); Tw \ Jg(4:6,:)];
  % Gelenk-Geschwindigkeit mit inverser Jacobi
  qD_k = J_x(I_EE,:) \ XD(k,I_EE)';

  %% Gelenk-Beschleunigung berechnen
  % Zeitableitung der geometrischen Jacobi-Matrix
  JgD = fourbar2TE_jacobigD_sym_varpar(q_k, qD_k, link_index, r_i_i_C, pkin);
  % Zeitableitung der Euler-Transformationsmatrix
  TDw = euljacD(X(k,4:6)', XD(k,4:6)', s.phiconv_W_E);
  % Zeitableitung der inversen Euler-Transformationsmatrix
  TwD_inv = -Tw\TDw/Tw;
  % Zeitableitung der analytischen Jacobi (Rotationsteil)
  JeD = Tw\JgD(4:6,:) + TwD_inv *Jg(4:6,:);
  % Zeitableitung analytische Jacobi komplett
  JD_x = [JgD(1:3,:); JeD];
  % Gelenk-Beschleunigung mit inverser Jacobi berechnen
  qDD_k = J_x(I_EE,:) \ (XDD(k,I_EE)' - JD_x(I_EE,:)*qD_k);

  %% Verbesserter Anfangswert für Positionsberechnung in nächster Iteration
  % Aus Geschwindigkeit berechneter neuer Winkel für den nächsten Zeitschritt
  % Taylor-Reihe bis 2. Ordnung für Position (Siehe [2])
  if k < nt
    dt = T(k+1)-T(k);
    qk0 = q_k + qD_k*dt + 0.5*qDD_k*dt^2;
  end

  %% Ergebnisse speichern
  Q(k,:) = q_k;
  QD(k,:) = qD_k;
  QDD(k,:) = qDD_k;
  PHI(k,:) = Phi_k;
end


