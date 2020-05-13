% Calculate kinetic energy for
% palh1m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% qJD [13x1]
%   Generalized joint velocities
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [11x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:30
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m2OL_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_energykin_fixb_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m2OL_energykin_fixb_slag_vp2: qJD has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_energykin_fixb_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2OL_energykin_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2OL_energykin_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2OL_energykin_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 21:17:24
% EndTime: 2020-05-02 21:17:25
% DurationCPUTime: 1.72s
% Computational Cost: add. (610->162), mult. (1362->290), div. (0->0), fcn. (1048->20), ass. (0->69)
t233 = m(5) / 0.2e1;
t219 = cos(qJ(3));
t232 = pkin(1) * t219;
t200 = qJD(2) + qJD(7);
t231 = pkin(4) * t200;
t230 = sin(qJ(10));
t229 = pkin(1) * qJD(2);
t228 = pkin(5) * qJD(3);
t220 = cos(qJ(2));
t227 = t219 * t220;
t226 = qJD(2) ^ 2 * pkin(1) ^ 2;
t211 = sin(qJ(3));
t193 = pkin(1) * t211 + pkin(5);
t210 = sin(qJ(4));
t218 = cos(qJ(4));
t174 = (t193 * t210 - t218 * t232) * qJD(2) + t210 * t228;
t201 = qJD(2) + qJD(3);
t199 = qJD(2) + qJD(8);
t212 = sin(qJ(2));
t182 = (t211 * t220 + t212 * t219) * qJD(1);
t185 = (-t211 * t212 + t227) * qJD(1);
t169 = -t182 * t210 + t185 * t218;
t175 = qJD(2) * (t193 * t218 + t210 * t232) + t218 * t228;
t222 = qJD(1) ^ 2;
t217 = cos(qJ(5));
t216 = cos(qJ(6));
t215 = cos(qJ(7));
t214 = cos(qJ(8));
t213 = cos(qJ(9));
t209 = sin(qJ(5));
t208 = sin(qJ(6));
t207 = sin(qJ(7));
t206 = sin(qJ(8));
t205 = sin(qJ(9));
t204 = cos(qJ(10));
t203 = cos(pkin(19));
t202 = sin(pkin(19));
t197 = qJD(4) + t201;
t196 = qJD(9) + t199;
t195 = qJD(10) + t200;
t192 = t212 * pkin(1) - pkin(15);
t190 = t192 ^ 2 * t222;
t189 = -t202 * t230 - t203 * t204;
t188 = t202 * t204 - t203 * t230;
t187 = t203 * t231 + t215 * t229;
t186 = t202 * t231 + t207 * t229;
t184 = (-t207 * t212 + t215 * t220) * qJD(1);
t183 = (-t206 * t212 + t214 * t220) * qJD(1);
t181 = (-t207 * t220 - t212 * t215) * qJD(1);
t180 = (-t206 * t220 - t212 * t214) * qJD(1);
t178 = -pkin(5) * t227 + (pkin(5) * t211 + pkin(1)) * t212 - pkin(15);
t176 = -pkin(2) * t180 - qJD(1) * pkin(15);
t173 = -pkin(9) * t197 - t175;
t172 = pkin(11) * t197 + t174;
t171 = t182 * t218 + t185 * t210;
t170 = -t180 * t205 - t183 * t213;
t168 = -t180 * t213 + t183 * t205;
t167 = qJD(5) - t169;
t166 = qJD(1) * t192 + (-t181 * t203 - t184 * t202) * pkin(4);
t165 = t171 * t217 + t197 * t209;
t164 = -t171 * t209 + t197 * t217;
t163 = t186 * t189 + t187 * t188;
t162 = -t186 * t188 + t187 * t189;
t161 = t181 * t188 + t184 * t189;
t160 = t181 * t189 - t184 * t188;
t159 = -pkin(9) * t169 - pkin(11) * t171 + qJD(1) * t178;
t158 = t159 * t209 + t172 * t217;
t157 = t159 * t217 - t172 * t209;
t1 = (-t166 * mrSges(11,1) + t163 * mrSges(11,3) + Ifges(11,4) * t161 + Ifges(11,6) * t195 + Ifges(11,2) * t160 / 0.2e1) * t160 + (t166 * mrSges(11,2) - t162 * mrSges(11,3) + Ifges(11,5) * t195 + Ifges(11,1) * t161 / 0.2e1) * t161 + Ifges(10,3) * t196 ^ 2 / 0.2e1 + (t174 * mrSges(5,3) + Ifges(5,4) * t171 + Ifges(5,6) * t197 + Ifges(5,2) * t169 / 0.2e1) * t169 + (Ifges(4,4) * t185 + Ifges(4,5) * t201 + Ifges(4,1) * t182 / 0.2e1) * t182 + (Ifges(2,3) / 0.2e1 + t178 ^ 2 * t233 + m(7) * pkin(14) ^ 2 / 0.2e1 + (m(3) / 0.2e1 + m(9) / 0.2e1) * pkin(15) ^ 2) * t222 + Ifges(9,1) * t183 ^ 2 / 0.2e1 + (t176 * mrSges(10,2) + Ifges(10,5) * t196 + Ifges(10,1) * t170 / 0.2e1) * t170 + (Ifges(3,3) * qJD(2) / 0.2e1 + (t207 * (-mrSges(8,2) * t200 + mrSges(8,3) * t181) + t211 * (mrSges(4,1) * t201 - mrSges(4,3) * t182) + t215 * (mrSges(8,1) * t200 - mrSges(8,3) * t184) - t219 * (-mrSges(4,2) * t201 + mrSges(4,3) * t185)) * pkin(1)) * qJD(2) + (Ifges(8,4) * t184 + Ifges(8,6) * t200 + Ifges(8,2) * t181 / 0.2e1) * t181 + m(11) * (t162 ^ 2 + t163 ^ 2 + t166 ^ 2) / 0.2e1 + Ifges(8,3) * t200 ^ 2 / 0.2e1 + Ifges(4,3) * t201 ^ 2 / 0.2e1 + (Ifges(8,5) * t200 + Ifges(8,1) * t184 / 0.2e1) * t184 + m(6) * (t157 ^ 2 + t158 ^ 2 + t173 ^ 2) / 0.2e1 + (-t176 * mrSges(10,1) + Ifges(10,4) * t170 + Ifges(10,6) * t196 + Ifges(10,2) * t168 / 0.2e1) * t168 + (Ifges(9,5) * t183 + Ifges(9,6) * t180 + Ifges(9,3) * t199 / 0.2e1 + (-t205 * (-mrSges(10,2) * t196 + mrSges(10,3) * t168) - t213 * (mrSges(10,1) * t196 - mrSges(10,3) * t170)) * pkin(2)) * t199 + (Ifges(9,4) * t183 + Ifges(9,2) * t180 / 0.2e1) * t180 + m(10) * (t176 ^ 2 + (t205 ^ 2 + t213 ^ 2) * pkin(2) ^ 2 * t199 ^ 2) / 0.2e1 + (t174 ^ 2 + t175 ^ 2) * t233 + (Ifges(4,6) * t201 + Ifges(4,2) * t185 / 0.2e1) * t185 + (t157 * mrSges(6,1) - t158 * mrSges(6,2) + Ifges(6,3) * t167 / 0.2e1) * t167 + (-t173 * mrSges(6,1) + t158 * mrSges(6,3) + Ifges(6,4) * t165 + Ifges(6,6) * t167 + Ifges(6,2) * t164 / 0.2e1) * t164 + (t162 * mrSges(11,1) - t163 * mrSges(11,2) + Ifges(11,3) * t195 / 0.2e1) * t195 + m(4) * (t190 + (t211 ^ 2 + t219 ^ 2) * t226) / 0.2e1 + m(8) * (t190 + (t207 ^ 2 + t215 ^ 2) * t226) / 0.2e1 + (-t175 * mrSges(5,3) + Ifges(5,1) * t171 / 0.2e1) * t171 + (t173 * mrSges(6,2) - t157 * mrSges(6,3) + Ifges(6,5) * t167 + Ifges(6,1) * t165 / 0.2e1) * t165 + (t175 * mrSges(5,1) - t174 * mrSges(5,2) + Ifges(5,5) * t171 + Ifges(5,3) * t197 / 0.2e1) * t197 + (t178 * (-mrSges(5,1) * t169 + mrSges(5,2) * t171) - pkin(15) * (-mrSges(9,1) * t180 + mrSges(9,2) * t183) + ((-pkin(15) * mrSges(3,2) + Ifges(3,1) * t220 / 0.2e1) * t220 + (-pkin(14) * mrSges(7,1) + Ifges(7,2) * t216 / 0.2e1) * t216 + (-mrSges(3,1) * pkin(15) - Ifges(3,4) * t220 + Ifges(3,2) * t212 / 0.2e1) * t212 + (pkin(14) * mrSges(7,2) + Ifges(7,4) * t216 + Ifges(7,1) * t208 / 0.2e1) * t208) * qJD(1) + (-mrSges(4,1) * t185 - mrSges(8,1) * t181 + mrSges(4,2) * t182 + mrSges(8,2) * t184) * t192 + qJD(6) * (Ifges(7,5) * t208 + Ifges(7,6) * t216) + qJD(2) * (Ifges(3,5) * t220 - Ifges(3,6) * t212)) * qJD(1) + Ifges(7,3) * qJD(6) ^ 2 / 0.2e1;
T = t1;
