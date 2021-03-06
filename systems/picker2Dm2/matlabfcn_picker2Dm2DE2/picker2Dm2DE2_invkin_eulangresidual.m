% Inverse Kinematik basierend auf Euler-Winkel-Residuum für
% picker2Dm2DE2
%
% Eingabe:
% xE_soll
%   EE-Lage (Sollwert)
% q0
%   Anfangs-Gelenkwinkel für Algorithmus
% s
%   Struktur mit Eingabedaten. Felder, siehe Quelltext.
%
% Ausgabe:
% q
%   Lösung der IK
% Phi
%   Restfehler mit der IK-Lösung

% Quelle:
% [SchapplerTapOrt2019] Schappler, M. et al.: Resolution of Functional
% Redundancy for 3T2R Robot Tasks using Two Sets of Reciprocal Euler
% Angles, Proc. of the 15th IFToMM World Congress, 2019
% [1] Aufzeichnungen Schappler vom 3.8.2018

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-17 01:41
% Revision: 3973e870210db52a56606103c46cf40226abfb99 (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-02
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [q, Phi] = picker2Dm2DE2_invkin_eulangresidual(xE_soll, q0, s)

%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(2,1),struct(
%$cgargs            'pkin', zeros(9,1),
%$cgargs          'sigmaJ', zeros(2,1),
%$cgargs             'NQJ', 0,
%$cgargs            'qlim', zeros(2,2),
%$cgargs            'I_EE', true(1,6),
%$cgargs     'phiconv_W_E', uint8(2),
%$cgargs        'I_EElink', uint8(0),
%$cgargs            'reci', true,
%$cgargs           'T_N_E', zeros(4,4),
%$cgargs               'K', zeros(2,1),
%$cgargs              'Kn', zeros(2,1),
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

% Einstellungsvariablen aus Struktur herausholen
phiconv_W_E = s.phiconv_W_E;
K = s.K;
Kn = s.Kn;
n_min = s.n_min;
n_max = s.n_max;
wn = s.wn;
Phit_tol = s.Phit_tol;
Phir_tol = s.Phir_tol;
retry_limit = s.retry_limit;
reci = s.reci;
I_EElink = s.I_EElink;
T_N_E = s.T_N_E;
pkin = s.pkin;
NQJ = s.NQJ;
qmin = s.qlim(:,1);
qmax = s.qlim(:,2);
if all(~isnan(s.qlim(:)))
  limits_set = true;
else
  % Grenzen sind nicht wirksam
  qmin(:) = -Inf;
  qmax(:) =  Inf;
  limits_set = false;
end
delta_qlim = qmax - qmin;

scale_lim = s.scale_lim;
maxrelstep = s.maxrelstep;
n_Phi_t = sum(s.I_EE(1:3)); % Anzahl der translatorischen Zwangsbedingungen
q = q0; %#ok<NASGU> % Zuweisung notwendig für mex

if any(wn ~= 0)
  nsoptim = true;
else
  % Keine zusätzlichen Optimierungskriterien
  nsoptim = false;
end

% Indizes für kinematische Zwangsbedingungen festlegen
if reci
  % Reihenfolge hergeleitet in [SchapplerTapOrt2019]; siehe z.B. Gl. (21)
  I_IK2 = [1 2 3 6 5 4];
else
  I_IK2 = [1 2 3 4 5 6];
end
I_IK = I_IK2(s.I_EE);

% Damit der Roboter einen Nullraum für Nebenoptimierungen hat, muss er min.
% 7FG für 6FG-Aufgaben und 6FG für 5FG-Aufgaben haben.
if nsoptim && NQJ <= length(I_IK)
  nsoptim = false;
end

% Ausgabe belegen
Phi = NaN(length(I_IK),1);

success = false;
q1 = q0;
%% Iterative Berechnung der inversen Kinematik
for rr = 0:retry_limit % Schleife über Neu-Anfänge der Berechnung

  % Variablen zum Speichern der Zwischenergebnisse
  q1 = q0;
  for jj = 2:n_max % Schleife über iteratives Annähern mit Newton-Verfahren

    % Gradienten-Matrix, siehe [SchapplerTapOrt2019]/(23)
    Jdk_voll = picker2Dm2DE2_constr2grad(q1, xE_soll, pkin, T_N_E, phiconv_W_E, I_EElink, reci);
    % Zwangsbedingungen
    Phi_voll = picker2Dm2DE2_constr2(q1, xE_soll, pkin, T_N_E, phiconv_W_E, I_EElink, reci);

    %% Reduktion auf betrachtete FG
    Jdk = Jdk_voll(I_IK,:);
    Phi = Phi_voll(I_IK);
    %% Nullstellensuche für Positions- und Orientierungsfehler
    % (Optimierung der Aufgabe)
    % Normale Invertierung der Jacobi-Matrix der seriellen Kette
    delta_q_T = Jdk \ (-Phi);
    %% Optimierung der Nebenbedingungen (Nullraum)
    delta_q_N = zeros(size(delta_q_T));
    if nsoptim && jj < n_max-10 % die letzten Iterationen sind zum Ausgleich des Positionsfehlers (ohne Nullraum)
      % Berechne Gradienten der zusätzlichen Optimierungskriterien
      v = zeros(NQJ, 1);
      if wn(1) ~= 0
        [~, hdq] = invkin_optimcrit_limits1(q1, [qmin, qmax]);
        v = v - wn(1)*hdq';
      end
      if wn(2) ~= 0
        [~, hdq] = invkin_optimcrit_limits2(q1, [qmin, qmax]);
        v = v - wn(2)*hdq';
      end
      % [SchapplerTapOrt2019]/(35); [1], Gl. (24)
      delta_q_N(:) = (eye(NQJ) - pinv(Jdk)* Jdk) * v; % "(:)" notwendig für mex
    end

    % [SchapplerTapOrt2019]/(35); [1], Gl. (23)
    delta_q = K.*delta_q_T + Kn.*delta_q_N;
    
    % Reduziere Schrittweite auf einen Maximalwert bezogen auf
    % Winkelgrenzen
    if limits_set && ~isnan(maxrelstep)
      % Bestimme Inkrement relativ zur Spannbreite der Grenzen
      abs_delta_q_rel = abs(delta_q ./ delta_qlim);
      if any(abs_delta_q_rel > maxrelstep)
        % Ein Element hat ein zu großes Inkrement. Normiere den
        % Inkrement-Vektor damit
        delta_q = delta_q .* maxrelstep / max(abs_delta_q_rel);
      end
    end
    
    q2 = q1 + delta_q;
    
    % Prüfe, ob die Gelenkwinkel ihre Grenzen überschreiten und reduziere
    % die Schrittweite, falls das der Fall ist
    if scale_lim
      delta_ul_rel = (qmax - q2)./(qmax-q1); % Überschreitung der Maximalwerte: <0
      delta_ll_rel = (-qmin + q2)./(q1-qmin); % Unterschreitung Minimalwerte: <0
      if any([delta_ul_rel;delta_ll_rel] < 0)
        % Berechne die prozentual stärkste Überschreitung
        % und nutze diese als Skalierung für die Winkeländerung
        % Reduziere die Winkeländerung so, dass die gröte Überschreitung auf
        % das Erreichen der Grenzen reduziert wird.
        if min(delta_ul_rel)<min(delta_ll_rel)
          % Verletzung nach oben ist die größere
          [~,I_max] = min(delta_ul_rel);
          scale = (qmax(I_max)-q1(I_max))./(delta_q(I_max));
        else
          % Verletzung nach unten ist maßgeblich
          [~,I_min] = min(delta_ll_rel);
          scale = (qmin(I_min)-q1(I_min))./(delta_q(I_min));
        end
        % Mit `scale` werden die Grenzen direkt für ein Gelenk erreicht.
        % Durch `scale_lim` wird dieses Erreichen weiter nach "innen" gezogen
        q2 = q1 + scale_lim*scale*delta_q;
      end
    end

    if any(isnan(q2)) || any(isinf(q2))
      break; % ab hier kann das Ergebnis nicht mehr besser werden wegen NaN/Inf
    end

    % Nachverarbeitung der Ergebnisse der Iteration
    q1 = q2;

    % Abbruchbedingungen prüfen
    if jj > n_min ... % Mindestzahl Iterationen erfüllt
        && all(abs(Phi(1:n_Phi_t)) < Phit_tol) && all(abs(Phi(n_Phi_t+1:end)) < Phir_tol) && ... % Haupt-Bedingung ist erfüllt
        ( ~nsoptim || ... %  und keine Nebenoptimierung läuft
        nsoptim && all(abs(delta_q_N) < 1e-10) ) % oder die Nullraumoptimierung läuft noch
      success = true;
      break;
    end
  end
  if success
    break;
  end
  % Beim vorherigen Durchlauf kein Erfolg. Generiere neue Anfangswerte
  q0(:) = qmin + rand(NQJ,1).*(qmax-qmin); % "(:)" notwendig für mex
end
q = q1;
if s.normalize
  sigmaJ = s.sigmaJ; % Marker für Dreh-/Schubgelenk (in den Minimalkoordinaten)
  q(sigmaJ==0) = normalize_angle(q(sigmaJ==0)); % nur Winkel normalisieren
end

