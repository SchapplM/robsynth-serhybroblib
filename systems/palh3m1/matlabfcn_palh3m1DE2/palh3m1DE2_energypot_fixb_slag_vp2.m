% Calculate potential energy for
% palh3m1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 16:51
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m1DE2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(19,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE2_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1DE2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE2_energypot_fixb_slag_vp2: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1DE2_energypot_fixb_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1DE2_energypot_fixb_slag_vp2: mrSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-19 19:34:30
% EndTime: 2020-04-19 19:34:37
% DurationCPUTime: 4.09s
% Computational Cost: add. (69661->120), mult. (104514->142), div. (4968->6), fcn. (66105->38), ass. (0->89)
t251 = pkin(3) * m(9);
t250 = -m(6) - m(5);
t187 = sin(qJ(4));
t193 = cos(qJ(4));
t249 = m(6) * pkin(9) + t193 * mrSges(6,1) - t187 * mrSges(6,2) + mrSges(5,1);
t248 = -m(9) - m(4) - m(8);
t247 = m(6) * pkin(11) - mrSges(5,2) + mrSges(6,3);
t246 = -t187 * mrSges(6,1) - t193 * mrSges(6,2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3);
t189 = sin(qJ(2));
t191 = sin(pkin(16));
t195 = cos(qJ(2));
t197 = cos(pkin(16));
t170 = t189 * t191 - t195 * t197;
t233 = pkin(5) * t170;
t226 = (-0.2e1 * t233 + pkin(1)) * pkin(1);
t243 = -pkin(6) + pkin(2);
t244 = -pkin(6) - pkin(2);
t163 = sqrt(-((pkin(5) - t243) * (pkin(5) + t243) + t226) * ((pkin(5) - t244) * (pkin(5) + t244) + t226));
t222 = pkin(5) ^ 2 + t226;
t224 = pkin(2) ^ 2 - pkin(6) ^ 2;
t165 = t222 + t224;
t167 = pkin(1) - t233;
t171 = t189 * t197 + t191 * t195;
t161 = pkin(5) * t165 * t171 + t163 * t167;
t194 = cos(qJ(3));
t166 = 0.1e1 / t222;
t228 = t166 / pkin(2);
t236 = sin(qJ(3)) / 0.2e1;
t230 = t163 * t171;
t160 = -pkin(5) * t230 + t165 * t167;
t240 = -t160 / 0.2e1;
t157 = (t161 * t236 + t194 * t240) * t228;
t239 = t161 / 0.2e1;
t158 = (t160 * t236 + t194 * t239) * t228;
t181 = pkin(18) + pkin(19);
t175 = sin(t181);
t176 = cos(t181);
t148 = -t157 * t176 - t158 * t175;
t234 = pkin(4) * t148;
t227 = -0.2e1 * pkin(3) * t234 + pkin(4) ^ 2;
t223 = pkin(3) ^ 2 + t227;
t225 = pkin(8) ^ 2 - pkin(10) ^ 2;
t141 = t223 - t225;
t144 = -pkin(3) * t148 + pkin(4);
t241 = -pkin(8) + pkin(10);
t242 = -pkin(8) - pkin(10);
t139 = sqrt(-((pkin(3) - t241) * (pkin(3) + t241) + t227) * ((pkin(3) - t242) * (pkin(3) + t242) + t227));
t147 = t157 * t175 - t158 * t176;
t232 = t139 * t147;
t137 = -pkin(3) * t232 + t141 * t144;
t138 = pkin(3) * t141 * t147 + t139 * t144;
t182 = qJ(2) + qJ(3);
t184 = cos(pkin(17));
t142 = 0.1e1 / t223;
t231 = t142 / pkin(10);
t238 = sin(pkin(17)) / 0.2e1;
t133 = atan2((t138 * t184 / 0.2e1 + t137 * t238) * t231, (-t137 * t184 / 0.2e1 + t138 * t238) * t231) + t182;
t131 = sin(t133);
t132 = cos(t133);
t140 = t223 + t225;
t143 = -pkin(3) + t234;
t186 = cos(pkin(19));
t237 = sin(pkin(19)) / 0.2e1;
t152 = qJ(2) + atan2((t160 * t237 + t186 * t239) * t228, (t161 * t237 + t186 * t240) * t228);
t151 = pkin(18) - t152;
t221 = t142 / pkin(8) / 0.2e1;
t136 = -atan2((pkin(4) * t140 * t147 - t139 * t143) * t221, (-pkin(4) * t232 - t140 * t143) * t221) + t151;
t134 = sin(t136);
t135 = cos(t136);
t149 = sin(t152);
t150 = cos(t152);
t164 = t222 - t224;
t168 = pkin(1) * t170 - pkin(5);
t159 = -pkin(1) * t230 - t164 * t168;
t162 = pkin(1) * t164 * t171 - t163 * t168;
t192 = sin(pkin(15));
t229 = t166 / pkin(6);
t235 = cos(pkin(15)) / 0.2e1;
t156 = atan2((t162 * t235 - t159 * t192 / 0.2e1) * t229, (t159 * t235 + t162 * t192 / 0.2e1) * t229);
t153 = sin(t156);
t154 = cos(t156);
t174 = t195 * pkin(1) + pkin(13);
t177 = sin(t182);
t178 = cos(t182);
t245 = -cos(t151) * t251 + m(7) * pkin(7) - m(3) * pkin(13) - mrSges(3,1) * t195 + mrSges(4,1) * t178 - mrSges(7,1) * t154 - mrSges(8,1) * t150 + mrSges(9,1) * t135 + mrSges(3,2) * t189 - mrSges(4,2) * t177 + mrSges(7,2) * t153 + mrSges(8,2) * t149 + mrSges(9,2) * t134 - mrSges(2,1) + t248 * t174 + t250 * (-pkin(4) * t178 + t174) + t249 * t132 + t247 * t131;
t173 = t189 * pkin(1) + pkin(12);
t196 = cos(qJ(1));
t190 = sin(qJ(1));
t1 = (sin(t151) * t251 - m(7) * pkin(14) - t189 * mrSges(3,1) + t177 * mrSges(4,1) - t153 * mrSges(7,1) - t149 * mrSges(8,1) - t134 * mrSges(9,1) - t195 * mrSges(3,2) + t178 * mrSges(4,2) - t154 * mrSges(7,2) - t150 * mrSges(8,2) + t135 * mrSges(9,2) - mrSges(1,3) - mrSges(2,3) + t250 * (-pkin(4) * t177 + t173) + t248 * t173 - t247 * t132 + t249 * t131 + (-m(7) - m(3) - m(2)) * pkin(12)) * g(3) + (t245 * t190 - t246 * t196 - mrSges(1,2)) * g(2) + (t246 * t190 + t245 * t196 - mrSges(1,1)) * g(1);
U = t1;
