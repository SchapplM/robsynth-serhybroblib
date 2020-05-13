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
% rSges [9x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 16:51
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m1DE2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(19,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE2_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1DE2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE2_energypot_fixb_slag_vp1: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1DE2_energypot_fixb_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m1DE2_energypot_fixb_slag_vp1: rSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-19 19:34:33
% EndTime: 2020-04-19 19:34:39
% DurationCPUTime: 4.03s
% Computational Cost: add. (69661->150), mult. (104501->211), div. (4968->6), fcn. (66105->38), ass. (0->100)
t250 = pkin(11) + rSges(6,3);
t249 = -pkin(6) - pkin(2);
t248 = -pkin(6) + pkin(2);
t247 = -pkin(8) - pkin(10);
t246 = -pkin(8) + pkin(10);
t190 = sin(qJ(2));
t192 = sin(pkin(16));
t196 = cos(qJ(2));
t198 = cos(pkin(16));
t170 = t190 * t192 - t196 * t198;
t237 = pkin(5) * t170;
t224 = (-0.2e1 * t237 + pkin(1)) * pkin(1);
t220 = pkin(5) ^ 2 + t224;
t222 = pkin(2) ^ 2 - pkin(6) ^ 2;
t163 = t220 + t222;
t165 = pkin(1) - t237;
t161 = sqrt(-((pkin(5) - t248) * (pkin(5) + t248) + t224) * ((pkin(5) - t249) * (pkin(5) + t249) + t224));
t171 = t190 * t198 + t192 * t196;
t232 = t161 * t171;
t158 = -pkin(5) * t232 + t163 * t165;
t245 = -t158 / 0.2e1;
t159 = pkin(5) * t163 * t171 + t161 * t165;
t244 = t159 / 0.2e1;
t243 = sin(pkin(17)) / 0.2e1;
t242 = sin(pkin(19)) / 0.2e1;
t241 = sin(qJ(3)) / 0.2e1;
t240 = cos(pkin(15)) / 0.2e1;
t195 = cos(qJ(3));
t164 = 0.1e1 / t220;
t230 = t164 / pkin(2);
t155 = (t159 * t241 + t195 * t245) * t230;
t156 = (t158 * t241 + t195 * t244) * t230;
t182 = pkin(18) + pkin(19);
t176 = sin(t182);
t177 = cos(t182);
t146 = -t155 * t177 - t156 * t176;
t238 = pkin(4) * t146;
t225 = -0.2e1 * pkin(3) * t238 + pkin(4) ^ 2;
t221 = pkin(3) ^ 2 + t225;
t223 = pkin(8) ^ 2 - pkin(10) ^ 2;
t139 = t221 - t223;
t142 = -pkin(3) * t146 + pkin(4);
t137 = sqrt(-((pkin(3) - t246) * (pkin(3) + t246) + t225) * ((pkin(3) - t247) * (pkin(3) + t247) + t225));
t145 = t155 * t176 - t156 * t177;
t234 = t137 * t145;
t135 = -pkin(3) * t234 + t139 * t142;
t136 = pkin(3) * t139 * t145 + t137 * t142;
t183 = qJ(2) + qJ(3);
t185 = cos(pkin(17));
t140 = 0.1e1 / t221;
t233 = t140 / pkin(10);
t131 = atan2((t136 * t185 / 0.2e1 + t135 * t243) * t233, (-t135 * t185 / 0.2e1 + t136 * t243) * t233) + t183;
t130 = cos(t131);
t236 = t130 * pkin(9);
t235 = t190 * pkin(1) + pkin(12);
t175 = t196 * pkin(1) + pkin(13);
t231 = t164 / pkin(6);
t188 = sin(qJ(4));
t191 = sin(qJ(1));
t229 = t188 * t191;
t197 = cos(qJ(1));
t228 = t188 * t197;
t194 = cos(qJ(4));
t227 = t191 * t194;
t226 = t194 * t197;
t187 = cos(pkin(19));
t150 = qJ(2) + atan2((t158 * t242 + t187 * t244) * t230, (t159 * t242 + t187 * t245) * t230);
t219 = t140 / pkin(8) / 0.2e1;
t149 = pkin(18) - t150;
t178 = sin(t183);
t218 = -pkin(4) * t178 + t235;
t179 = cos(t183);
t217 = -rSges(4,1) * t179 + rSges(4,2) * t178;
t129 = sin(t131);
t216 = -rSges(5,1) * t130 + rSges(5,2) * t129;
t147 = sin(t150);
t148 = cos(t150);
t215 = rSges(8,1) * t148 - rSges(8,2) * t147;
t162 = t220 - t222;
t166 = pkin(1) * t170 - pkin(5);
t157 = -pkin(1) * t232 - t162 * t166;
t160 = pkin(1) * t162 * t171 - t161 * t166;
t193 = sin(pkin(15));
t154 = atan2((t160 * t240 - t157 * t193 / 0.2e1) * t231, (t157 * t240 + t160 * t193 / 0.2e1) * t231);
t151 = sin(t154);
t152 = cos(t154);
t214 = rSges(7,1) * t152 - rSges(7,2) * t151 - pkin(7);
t213 = rSges(3,1) * t196 - rSges(3,2) * t190 + pkin(13);
t138 = t221 + t223;
t141 = -pkin(3) + t238;
t134 = -atan2((pkin(4) * t138 * t145 - t137 * t141) * t219, (-pkin(4) * t234 - t138 * t141) * t219) + t149;
t132 = sin(t134);
t133 = cos(t134);
t212 = -rSges(9,1) * t133 - rSges(9,2) * t132 + pkin(3) * cos(t149) + t175;
t174 = t197 * t175;
t173 = t191 * t175;
t172 = -pkin(4) * t179 + t175;
t169 = t197 * t172;
t168 = t191 * t172;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t197 - rSges(2,2) * t191) + g(2) * (rSges(2,1) * t191 + rSges(2,2) * t197) + g(3) * (pkin(12) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t190 + rSges(3,2) * t196 + pkin(12)) + (-g(2) * rSges(3,3) + g(1) * t213) * t197 + (g(1) * rSges(3,3) + g(2) * t213) * t191) - m(4) * (g(1) * (rSges(4,3) * t191 + t217 * t197 + t174) + g(2) * (-rSges(4,3) * t197 + t217 * t191 + t173) + g(3) * (-rSges(4,1) * t178 - rSges(4,2) * t179 + t235)) - m(5) * (g(1) * (rSges(5,3) * t191 + t216 * t197 + t169) + g(2) * (-rSges(5,3) * t197 + t216 * t191 + t168) + g(3) * (-rSges(5,1) * t129 - rSges(5,2) * t130 + t218)) - m(6) * (g(1) * (-t197 * t236 + t169 + (-t130 * t226 + t229) * rSges(6,1) + (t130 * t228 + t227) * rSges(6,2)) + g(2) * (-t191 * t236 + t168 + (-t130 * t227 - t228) * rSges(6,1) + (t130 * t229 - t226) * rSges(6,2)) + g(3) * (t250 * t130 + t218) + (g(3) * (-rSges(6,1) * t194 + rSges(6,2) * t188 - pkin(9)) - (g(1) * t197 + g(2) * t191) * t250) * t129) - m(7) * (g(3) * (rSges(7,1) * t151 + rSges(7,2) * t152 + pkin(12) + pkin(14)) + (-g(2) * rSges(7,3) + g(1) * t214) * t197 + (g(1) * rSges(7,3) + g(2) * t214) * t191) - m(8) * (g(1) * (rSges(8,3) * t191 + t215 * t197 + t174) + g(2) * (-rSges(8,3) * t197 + t215 * t191 + t173) + g(3) * (rSges(8,1) * t147 + rSges(8,2) * t148 + t235)) - m(9) * (g(3) * (-pkin(3) * sin(t149) + t132 * rSges(9,1) - t133 * rSges(9,2) + t235) + (-g(2) * rSges(9,3) + g(1) * t212) * t197 + (g(1) * rSges(9,3) + g(2) * t212) * t191);
U = t1;
