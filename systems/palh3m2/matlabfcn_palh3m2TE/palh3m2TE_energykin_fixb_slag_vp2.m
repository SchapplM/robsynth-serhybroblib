% Calculate kinetic energy for
% palh3m2TE
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
% Datum: 2020-05-07 01:49
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m2TE_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2TE_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2TE_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2TE_energykin_fixb_slag_vp2: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2TE_energykin_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2TE_energykin_fixb_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2TE_energykin_fixb_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:42:08
% EndTime: 2020-05-07 01:42:08
% DurationCPUTime: 0.74s
% Computational Cost: add. (788->117), mult. (1471->208), div. (0->0), fcn. (1636->20), ass. (0->61)
t209 = cos(pkin(16));
t219 = cos(pkin(15));
t231 = sin(pkin(16));
t234 = sin(pkin(15));
t189 = t209 * t234 + t219 * t231;
t190 = t209 * t219 - t231 * t234;
t207 = pkin(17) + pkin(18);
t202 = sin(t207);
t203 = cos(t207);
t179 = (t189 * t202 - t190 * t203) * qJD(1);
t206 = qJD(2) + qJD(3);
t217 = cos(qJ(3));
t232 = pkin(1) * qJD(2);
t224 = pkin(4) * t206 - t217 * t232;
t213 = sin(qJ(3));
t180 = t189 * t217 + t190 * t213;
t181 = -t189 * t213 + t190 * t217;
t214 = sin(qJ(2));
t218 = cos(qJ(2));
t174 = t180 * t218 + t181 * t214;
t175 = -t180 * t214 + t181 * t218;
t226 = t174 * t202 - t175 * t203;
t227 = t174 * t203 + t175 * t202;
t228 = t213 * t232;
t169 = t224 * t226 + t227 * t228;
t237 = t169 ^ 2;
t236 = m(9) / 0.2e1;
t178 = (t189 * t203 + t190 * t202) * qJD(1);
t235 = t178 / 0.2e1;
t233 = sin(pkin(18));
t230 = qJD(2) ^ 2 * pkin(1) ^ 2;
t229 = (t213 ^ 2 + t217 ^ 2) * t230;
t199 = -pkin(1) * t218 - pkin(12);
t187 = (t213 * t214 - t217 * t218) * qJD(1);
t186 = -pkin(4) * t187 + t199 * qJD(1);
t211 = cos(pkin(18));
t191 = t211 * t234 + t219 * t233;
t222 = qJD(1) ^ 2;
t220 = cos(pkin(14));
t216 = cos(qJ(4));
t215 = sin(pkin(14));
t212 = sin(qJ(4));
t208 = qJ(3) + qJ(2);
t205 = cos(t208);
t204 = sin(t208);
t195 = t199 ^ 2 * t222;
t194 = t215 * t234 + t220 * t219;
t193 = -t219 * t215 + t220 * t234;
t192 = t211 * t219 - t233 * t234;
t188 = (-t213 * t218 - t214 * t217) * qJD(1);
t185 = -t193 * t214 + t194 * t218;
t184 = t193 * t218 + t194 * t214;
t183 = t191 * t214 - t192 * t218;
t182 = t191 * t218 + t192 * t214;
t177 = qJD(4) - t179;
t176 = (cos(pkin(17)) * t192 - t191 * sin(pkin(17))) * pkin(3) + t199;
t173 = -pkin(8) * t179 - pkin(10) * t178 + t186;
t171 = -t224 * t227 + t226 * t228;
t168 = t171 * t216 + t173 * t212;
t167 = -t171 * t212 + t173 * t216;
t1 = m(5) * (t171 ^ 2 + t186 ^ 2 + t237) / 0.2e1 + m(6) * (t167 ^ 2 + t168 ^ 2 + t237) / 0.2e1 + m(4) * (t195 + t229) / 0.2e1 + Ifges(4,1) * t188 ^ 2 / 0.2e1 + t229 * t236 + m(8) * (t195 + (t182 ^ 2 + t183 ^ 2) * t230) / 0.2e1 + (m(7) * pkin(6) ^ 2 / 0.2e1 + m(3) * pkin(12) ^ 2 / 0.2e1 + t176 ^ 2 * t236 + Ifges(2,3) / 0.2e1) * t222 + (Ifges(4,4) * t188 + Ifges(4,2) * t187 / 0.2e1) * t187 + (-t186 * mrSges(5,1) + t171 * mrSges(5,3) + Ifges(5,2) * t179 / 0.2e1) * t179 + (t167 * mrSges(6,1) - t168 * mrSges(6,2) + Ifges(6,3) * t177 / 0.2e1) * t177 + (Ifges(4,5) * t188 + Ifges(4,6) * t187 + (Ifges(9,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t206) * t206 + ((Ifges(3,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * qJD(2) + ((-t213 * t187 + t217 * t188) * mrSges(4,3) + ((-mrSges(4,1) - mrSges(9,1)) * t217 + (mrSges(4,2) + mrSges(9,2)) * t213) * t206) * pkin(1)) * qJD(2) + (t186 * mrSges(5,2) + Ifges(5,4) * t179 + t169 * mrSges(5,3) + Ifges(5,1) * t235 + (Ifges(6,1) * t216 * t235 + t169 * mrSges(6,2) - t167 * mrSges(6,3) + t177 * Ifges(6,5)) * t216 + (t169 * mrSges(6,1) - t168 * mrSges(6,3) - t177 * Ifges(6,6) + (-Ifges(6,4) * t216 + Ifges(6,2) * t212 / 0.2e1) * t178) * t212) * t178 + (t199 * (-mrSges(4,1) * t187 + mrSges(4,2) * t188) + ((pkin(12) * mrSges(3,1) + Ifges(3,2) * t218 / 0.2e1) * t218 + (t176 * mrSges(9,1) + Ifges(9,2) * t205 / 0.2e1) * t205 + (t199 * mrSges(8,1) + Ifges(8,2) * t192 / 0.2e1) * t192 + (-pkin(6) * mrSges(7,1) + Ifges(7,2) * t185 / 0.2e1) * t185 + (-pkin(12) * mrSges(3,2) + Ifges(3,4) * t218 + Ifges(3,1) * t214 / 0.2e1) * t214 + (-t176 * mrSges(9,2) + Ifges(9,4) * t205 + Ifges(9,1) * t204 / 0.2e1) * t204 + (t199 * mrSges(8,2) - Ifges(8,4) * t192 + Ifges(8,1) * t191 / 0.2e1) * t191 + (pkin(6) * mrSges(7,2) + Ifges(7,4) * t185 + Ifges(7,1) * t184 / 0.2e1) * t184) * qJD(1) + t206 * (-Ifges(9,5) * t204 - Ifges(9,6) * t205) + (Ifges(3,5) * t214 + Ifges(7,5) * t184 + Ifges(3,6) * t218 + Ifges(7,6) * t185 + ((-t204 * t217 + t205 * t213) * mrSges(9,3) + (-t182 * t192 - t183 * t191) * mrSges(8,3)) * pkin(1)) * qJD(2)) * qJD(1);
T = t1;
