% Calculate kinetic energy for
% palh3m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [9x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m2DE2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE2_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_energykin_fixb_slag_vp2: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE2_energykin_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2DE2_energykin_fixb_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2DE2_energykin_fixb_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 02:14:31
% EndTime: 2020-05-07 02:14:32
% DurationCPUTime: 0.81s
% Computational Cost: add. (788->117), mult. (1471->208), div. (0->0), fcn. (1636->20), ass. (0->61)
t210 = cos(pkin(16));
t220 = cos(pkin(15));
t233 = sin(pkin(16));
t236 = sin(pkin(15));
t190 = t210 * t236 + t233 * t220;
t191 = t210 * t220 - t233 * t236;
t208 = pkin(17) + pkin(18);
t203 = sin(t208);
t204 = cos(t208);
t180 = (t190 * t203 - t191 * t204) * qJD(1);
t207 = qJD(2) + qJD(3);
t218 = cos(qJ(3));
t234 = pkin(1) * qJD(2);
t226 = pkin(4) * t207 - t218 * t234;
t214 = sin(qJ(3));
t181 = t190 * t218 + t191 * t214;
t182 = -t190 * t214 + t191 * t218;
t215 = sin(qJ(2));
t219 = cos(qJ(2));
t175 = t181 * t219 + t215 * t182;
t176 = -t181 * t215 + t182 * t219;
t228 = t175 * t203 - t176 * t204;
t229 = t175 * t204 + t176 * t203;
t230 = t214 * t234;
t170 = t228 * t226 + t229 * t230;
t239 = t170 ^ 2;
t238 = m(9) / 0.2e1;
t179 = (t190 * t204 + t191 * t203) * qJD(1);
t237 = t179 / 0.2e1;
t235 = sin(pkin(18));
t232 = qJD(2) ^ 2 * pkin(1) ^ 2;
t231 = (t214 ^ 2 + t218 ^ 2) * t232;
t200 = -t219 * pkin(1) - pkin(12);
t188 = (t214 * t215 - t218 * t219) * qJD(1);
t187 = -pkin(4) * t188 + t200 * qJD(1);
t212 = cos(pkin(18));
t225 = -t236 * t212 - t220 * t235;
t223 = qJD(1) ^ 2;
t221 = cos(pkin(14));
t217 = cos(qJ(4));
t216 = sin(pkin(14));
t213 = sin(qJ(4));
t209 = qJ(2) + qJ(3);
t206 = cos(t209);
t205 = sin(t209);
t196 = t200 ^ 2 * t223;
t195 = t236 * t216 + t220 * t221;
t194 = -t220 * t216 + t236 * t221;
t193 = t212 * t220 - t236 * t235;
t189 = (-t214 * t219 - t215 * t218) * qJD(1);
t186 = -t194 * t215 + t195 * t219;
t185 = t194 * t219 + t195 * t215;
t184 = -t193 * t219 - t215 * t225;
t183 = t193 * t215 - t219 * t225;
t178 = qJD(4) - t180;
t177 = (cos(pkin(17)) * t193 + t225 * sin(pkin(17))) * pkin(3) + t200;
t174 = -pkin(8) * t180 - pkin(10) * t179 + t187;
t172 = -t229 * t226 + t228 * t230;
t169 = t172 * t217 + t174 * t213;
t168 = -t172 * t213 + t174 * t217;
t1 = m(5) * (t172 ^ 2 + t187 ^ 2 + t239) / 0.2e1 + m(6) * (t168 ^ 2 + t169 ^ 2 + t239) / 0.2e1 + m(4) * (t196 + t231) / 0.2e1 + Ifges(4,1) * t189 ^ 2 / 0.2e1 + t231 * t238 + m(8) * (t196 + (t183 ^ 2 + t184 ^ 2) * t232) / 0.2e1 + (Ifges(2,3) / 0.2e1 + t177 ^ 2 * t238 + m(3) * pkin(12) ^ 2 / 0.2e1 + m(7) * pkin(6) ^ 2 / 0.2e1) * t223 + (Ifges(4,4) * t189 + Ifges(4,2) * t188 / 0.2e1) * t188 + (-t187 * mrSges(5,1) + t172 * mrSges(5,3) + Ifges(5,2) * t180 / 0.2e1) * t180 + (t168 * mrSges(6,1) - t169 * mrSges(6,2) + Ifges(6,3) * t178 / 0.2e1) * t178 + (Ifges(4,5) * t189 + Ifges(4,6) * t188 + (Ifges(9,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t207) * t207 + ((Ifges(7,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * qJD(2) + ((-t214 * t188 + t218 * t189) * mrSges(4,3) + ((-mrSges(4,1) - mrSges(9,1)) * t218 + (mrSges(4,2) + mrSges(9,2)) * t214) * t207) * pkin(1)) * qJD(2) + (Ifges(5,4) * t180 + t187 * mrSges(5,2) + t170 * mrSges(5,3) + Ifges(5,1) * t237 + (Ifges(6,1) * t217 * t237 + t170 * mrSges(6,2) - t168 * mrSges(6,3) + Ifges(6,5) * t178) * t217 + (t170 * mrSges(6,1) - t169 * mrSges(6,3) - Ifges(6,6) * t178 + (-Ifges(6,4) * t217 + Ifges(6,2) * t213 / 0.2e1) * t179) * t213) * t179 + (t200 * (-mrSges(4,1) * t188 + mrSges(4,2) * t189) + ((pkin(12) * mrSges(3,1) + Ifges(3,2) * t219 / 0.2e1) * t219 + (t177 * mrSges(9,1) + Ifges(9,2) * t206 / 0.2e1) * t206 + (t200 * mrSges(8,1) + Ifges(8,2) * t193 / 0.2e1) * t193 + (-pkin(6) * mrSges(7,1) + Ifges(7,2) * t186 / 0.2e1) * t186 + (-pkin(12) * mrSges(3,2) + Ifges(3,4) * t219 + Ifges(3,1) * t215 / 0.2e1) * t215 + (-t177 * mrSges(9,2) + Ifges(9,4) * t206 + Ifges(9,1) * t205 / 0.2e1) * t205 - (t200 * mrSges(8,2) - Ifges(8,4) * t193 - Ifges(8,1) * t225 / 0.2e1) * t225 + (pkin(6) * mrSges(7,2) + Ifges(7,4) * t186 + Ifges(7,1) * t185 / 0.2e1) * t185) * qJD(1) + t207 * (-Ifges(9,5) * t205 - Ifges(9,6) * t206) + (Ifges(3,5) * t215 + Ifges(7,5) * t185 + Ifges(3,6) * t219 + Ifges(7,6) * t186 + ((-t205 * t218 + t206 * t214) * mrSges(9,3) + (-t183 * t193 + t184 * t225) * mrSges(8,3)) * pkin(1)) * qJD(2)) * qJD(1);
T = t1;
