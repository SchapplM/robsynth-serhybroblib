% Calculate potential energy for
% palh3m1DE1
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
% Datum: 2020-04-19 19:20
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m1DE1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(19,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE1_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1DE1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE1_energypot_fixb_slag_vp1: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1DE1_energypot_fixb_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m1DE1_energypot_fixb_slag_vp1: rSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-18 10:25:42
% EndTime: 2020-04-18 10:25:52
% DurationCPUTime: 6.86s
% Computational Cost: add. (157249->175), mult. (236340->261), div. (11232->6), fcn. (149713->36), ass. (0->111)
t269 = pkin(11) + rSges(6,3);
t279 = -pkin(6) - pkin(2);
t278 = -pkin(6) + pkin(2);
t277 = -pkin(8) - pkin(10);
t276 = -pkin(8) + pkin(10);
t221 = sin(qJ(2));
t223 = sin(pkin(16));
t227 = cos(qJ(2));
t229 = cos(pkin(16));
t202 = t221 * t223 - t227 * t229;
t266 = pkin(5) * t202;
t258 = (-0.2e1 * t266 + pkin(1)) * pkin(1);
t249 = pkin(5) ^ 2 + t258;
t254 = pkin(2) ^ 2 - pkin(6) ^ 2;
t189 = t249 + t254;
t197 = pkin(1) - t266;
t187 = sqrt(-((pkin(5) - t278) * (pkin(5) + t278) + t258) * ((pkin(5) - t279) * (pkin(5) + t279) + t258));
t204 = t221 * t229 + t223 * t227;
t262 = t187 * t204;
t184 = -pkin(5) * t262 + t189 * t197;
t275 = -t184 / 0.2e1;
t185 = pkin(5) * t189 * t204 + t187 * t197;
t274 = t185 / 0.2e1;
t273 = sin(pkin(17)) / 0.2e1;
t272 = sin(pkin(19)) / 0.2e1;
t220 = sin(qJ(3));
t271 = t220 / 0.2e1;
t270 = cos(pkin(15)) / 0.2e1;
t268 = pkin(1) * t227;
t226 = cos(qJ(3));
t190 = 0.1e1 / t249;
t260 = t190 / pkin(2);
t181 = (t185 * t271 + t226 * t275) * t260;
t182 = (t184 * t271 + t226 * t274) * t260;
t212 = pkin(18) + pkin(19);
t207 = sin(t212);
t208 = cos(t212);
t174 = -t181 * t208 - t182 * t207;
t267 = pkin(4) * t174;
t265 = t221 * pkin(1) + pkin(12);
t259 = -0.2e1 * pkin(3) * t267 + pkin(4) ^ 2;
t160 = sqrt(-((pkin(3) - t276) * (pkin(3) + t276) + t259) * ((pkin(3) - t277) * (pkin(3) + t277) + t259));
t173 = t181 * t207 - t182 * t208;
t264 = t160 * t173;
t252 = pkin(3) ^ 2 + t259;
t169 = 0.1e1 / t252;
t263 = t169 / pkin(10);
t261 = t190 / pkin(6);
t222 = sin(qJ(1));
t209 = t222 * pkin(13);
t257 = t222 * t268 + t209;
t228 = cos(qJ(1));
t210 = t228 * pkin(13);
t256 = t228 * t268 + t210;
t255 = pkin(8) ^ 2 - pkin(10) ^ 2;
t246 = t220 * t227 + t221 * t226;
t253 = -pkin(4) * t246 + t265;
t201 = t220 * t221 - t226 * t227;
t194 = t201 * t222;
t251 = t194 * pkin(4) + t257;
t196 = t201 * t228;
t250 = t196 * pkin(4) + t256;
t248 = t169 / pkin(8) / 0.2e1;
t247 = rSges(3,1) * t227 - rSges(3,2) * t221;
t217 = cos(pkin(19));
t179 = atan2((t184 * t272 + t217 * t274) * t260, (t185 * t272 + t217 * t275) * t260);
t175 = sin(t179);
t176 = cos(t179);
t166 = t175 * t227 + t176 * t221;
t165 = -t175 * t221 + t176 * t227;
t188 = t249 - t254;
t198 = pkin(1) * t202 - pkin(5);
t183 = -pkin(1) * t262 - t188 * t198;
t186 = pkin(1) * t188 * t204 - t187 * t198;
t224 = sin(pkin(15));
t180 = atan2((t186 * t270 - t183 * t224 / 0.2e1) * t261, (t183 * t270 + t186 * t224 / 0.2e1) * t261);
t177 = sin(t180);
t178 = cos(t180);
t245 = rSges(7,1) * t178 - rSges(7,2) * t177 - pkin(7);
t167 = t252 + t255;
t170 = -pkin(3) + t267;
t244 = atan2((pkin(4) * t167 * t173 - t160 * t170) * t248, (-pkin(4) * t264 - t167 * t170) * t248);
t243 = sin(t244);
t225 = cos(qJ(4));
t219 = sin(qJ(4));
t218 = cos(pkin(18));
t216 = sin(pkin(18));
t214 = cos(pkin(17));
t195 = t246 * t228;
t193 = t246 * t222;
t171 = -pkin(3) * t174 + pkin(4);
t168 = t252 - t255;
t164 = t165 * t228;
t163 = t166 * t228;
t162 = t165 * t222;
t161 = t166 * t222;
t159 = pkin(3) * t168 * t173 + t160 * t171;
t158 = -pkin(3) * t264 + t168 * t171;
t157 = cos(t244);
t155 = atan2((t159 * t214 / 0.2e1 + t158 * t273) * t263, (-t158 * t214 / 0.2e1 + t159 * t273) * t263);
t154 = cos(t155);
t153 = sin(t155);
t152 = -t218 * t157 - t216 * t243;
t151 = t157 * t216 - t218 * t243;
t147 = t153 * t201 - t154 * t246;
t146 = -t153 * t246 - t154 * t201;
t145 = t153 * t195 + t154 * t196;
t144 = t153 * t196 - t154 * t195;
t143 = t153 * t193 + t154 * t194;
t142 = t153 * t194 - t154 * t193;
t1 = -m(1) * (rSges(1,1) * g(1) + g(2) * rSges(1,2) + rSges(1,3) * g(3)) - m(2) * (g(1) * (rSges(2,1) * t228 - rSges(2,2) * t222) + g(2) * (rSges(2,1) * t222 + rSges(2,2) * t228) + g(3) * (pkin(12) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,3) * t222 + t228 * t247 + t210) + g(2) * (-rSges(3,3) * t228 + t222 * t247 + t209) + g(3) * (rSges(3,1) * t221 + rSges(3,2) * t227 + pkin(12))) - m(4) * (g(1) * (rSges(4,1) * t196 + rSges(4,2) * t195 + rSges(4,3) * t222 + t256) + g(2) * (rSges(4,1) * t194 + rSges(4,2) * t193 - rSges(4,3) * t228 + t257) + g(3) * (-rSges(4,1) * t246 + rSges(4,2) * t201 + t265)) - m(5) * (g(1) * (rSges(5,1) * t145 - rSges(5,2) * t144 + rSges(5,3) * t222 + t250) + g(2) * (rSges(5,1) * t143 - rSges(5,2) * t142 - rSges(5,3) * t228 + t251) + g(3) * (rSges(5,1) * t147 - rSges(5,2) * t146 + t253)) - m(6) * (g(1) * (t145 * pkin(9) + (t145 * t225 + t219 * t222) * rSges(6,1) + (-t145 * t219 + t222 * t225) * rSges(6,2) + t269 * t144 + t250) + g(2) * (t143 * pkin(9) + (t143 * t225 - t219 * t228) * rSges(6,1) + (-t143 * t219 - t225 * t228) * rSges(6,2) + t269 * t142 + t251) + (t253 + (rSges(6,1) * t225 - rSges(6,2) * t219 + pkin(9)) * t147 + t269 * t146) * g(3)) - m(7) * (g(3) * (rSges(7,1) * t177 + rSges(7,2) * t178 + pkin(12) + pkin(14)) + (-g(2) * rSges(7,3) + g(1) * t245) * t228 + (g(1) * rSges(7,3) + g(2) * t245) * t222) - m(8) * (g(1) * (rSges(8,1) * t164 - rSges(8,2) * t163 + rSges(8,3) * t222 + t256) + g(2) * (rSges(8,1) * t162 - rSges(8,2) * t161 - rSges(8,3) * t228 + t257) + g(3) * (rSges(8,1) * t166 + rSges(8,2) * t165 + t265)) - m(9) * (g(1) * ((-t151 * t163 + t152 * t164) * rSges(9,1) + (-t151 * t164 - t152 * t163) * rSges(9,2) + t222 * rSges(9,3) + t256) + g(2) * ((-t151 * t161 + t152 * t162) * rSges(9,1) + (-t151 * t162 - t152 * t161) * rSges(9,2) - t228 * rSges(9,3) + t257) + g(3) * ((t151 * t165 + t152 * t166) * rSges(9,1) + (-t151 * t166 + t152 * t165) * rSges(9,2) + t265) + (g(1) * (t163 * t216 + t164 * t218) + g(2) * (t161 * t216 + t162 * t218) + g(3) * (-t165 * t216 + t166 * t218)) * pkin(3));
U = t1;
