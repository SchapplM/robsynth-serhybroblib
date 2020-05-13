% Calculate potential energy for
% palh3m1TE
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
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m1TE_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(19,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1TE_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_energypot_fixb_slag_vp1: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1TE_energypot_fixb_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m1TE_energypot_fixb_slag_vp1: rSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-17 15:18:57
% EndTime: 2020-04-17 15:19:07
% DurationCPUTime: 5.20s
% Computational Cost: add. (78745->175), mult. (118428->264), div. (5616->6), fcn. (74977->24), ass. (0->106)
t254 = -pkin(11) - rSges(6,3);
t259 = 0.1e1 / pkin(2);
t258 = pkin(2) - pkin(6);
t257 = -pkin(6) - pkin(2);
t256 = -pkin(8) - pkin(10);
t255 = -pkin(8) + pkin(10);
t253 = sin(pkin(15));
t190 = cos(qJ(2));
t252 = pkin(1) * t190;
t185 = sin(qJ(2));
t187 = sin(pkin(16));
t192 = cos(pkin(16));
t173 = t185 * t187 - t190 * t192;
t251 = pkin(5) * t173;
t250 = cos(pkin(19));
t249 = sin(pkin(19));
t248 = t185 * pkin(1) + pkin(12);
t247 = cos(pkin(17));
t246 = sin(pkin(17));
t243 = (-0.2e1 * t251 + pkin(1)) * pkin(1);
t233 = pkin(5) ^ 2 + t243;
t162 = 0.1e1 / t233;
t245 = t162 / pkin(6);
t160 = sqrt(-((pkin(5) - t258) * (pkin(5) + t258) + t243) * ((pkin(5) - t257) * (pkin(5) + t257) + t243));
t175 = t185 * t192 + t187 * t190;
t244 = t175 * t160;
t186 = sin(qJ(1));
t178 = t186 * pkin(13);
t242 = t186 * t252 + t178;
t191 = cos(qJ(1));
t179 = t191 * pkin(13);
t241 = t191 * t252 + t179;
t240 = -pkin(2) ^ 2 + pkin(6) ^ 2;
t239 = -pkin(8) ^ 2 + pkin(10) ^ 2;
t238 = pkin(18) + pkin(19);
t237 = cos(pkin(15)) / 0.2e1;
t184 = sin(qJ(3));
t189 = cos(qJ(3));
t227 = t184 * t190 + t185 * t189;
t236 = -pkin(4) * t227 + t248;
t172 = t184 * t185 - t189 * t190;
t166 = t172 * t186;
t235 = t166 * pkin(4) + t242;
t168 = t172 * t191;
t234 = t168 * pkin(4) + t241;
t232 = pkin(1) - t251;
t231 = cos(t238);
t230 = sin(t238);
t229 = rSges(3,1) * t190 - rSges(3,2) * t185;
t228 = t233 - t240;
t223 = t259 * (-pkin(5) * t244 + t228 * t232);
t221 = -t223 / 0.2e1;
t222 = t259 * (pkin(5) * t175 * t228 + t160 * t232) / 0.2e1;
t152 = (t221 * t250 + t222 * t249) * t162;
t220 = t223 / 0.2e1;
t153 = (t220 * t249 + t222 * t250) * t162;
t151 = t152 * t190 - t153 * t185;
t150 = t152 * t185 + t153 * t190;
t161 = t233 + t240;
t169 = pkin(1) * t173 - pkin(5);
t224 = -pkin(1) * t244 - t161 * t169;
t225 = pkin(1) * t161 * t175 - t160 * t169;
t154 = (t224 * t237 + t225 * t253 / 0.2e1) * t245;
t155 = (t225 * t237 - t224 * t253 / 0.2e1) * t245;
t226 = rSges(7,1) * t154 - rSges(7,2) * t155 - pkin(7);
t219 = t184 * t220 + t189 * t222;
t218 = t184 * t222 + t189 * t221;
t217 = t162 * (-t218 * t231 - t219 * t230);
t216 = t162 * (t218 * t230 - t219 * t231);
t215 = pkin(4) * t217;
t214 = pkin(3) - t215;
t213 = -pkin(3) * t217 + pkin(4);
t212 = -0.2e1 * pkin(3) * t215 + pkin(4) ^ 2;
t211 = pkin(3) ^ 2 + t212;
t210 = 0.1e1 / t211;
t209 = 0.1e1 / pkin(8) * t210;
t208 = 0.1e1 / pkin(10) * t210;
t207 = t211 + t239;
t206 = t211 - t239;
t205 = sqrt(-((pkin(3) - t255) * (pkin(3) + t255) + t212) * ((pkin(3) - t256) * (pkin(3) + t256) + t212));
t204 = t205 * t216;
t203 = (-pkin(4) * t204 + t206 * t214) * t209;
t202 = (-pkin(3) * t204 + t207 * t213) * t208;
t201 = -(pkin(4) * t206 * t216 + t205 * t214) * t209 / 0.2e1;
t200 = (pkin(3) * t207 * t216 + t205 * t213) * t208 / 0.2e1;
t188 = cos(qJ(4));
t183 = sin(qJ(4));
t182 = cos(pkin(18));
t181 = sin(pkin(18));
t167 = t227 * t191;
t165 = t227 * t186;
t149 = t150 * t191;
t148 = t151 * t191;
t147 = t150 * t186;
t146 = t151 * t186;
t141 = -t182 * t203 / 0.2e1 + t181 * t201;
t140 = t181 * t203 / 0.2e1 + t182 * t201;
t139 = t247 * t200 + t246 * t202 / 0.2e1;
t138 = -t247 * t202 / 0.2e1 + t246 * t200;
t137 = t138 * t172 + t139 * t227;
t136 = -t138 * t227 + t139 * t172;
t135 = t138 * t167 - t139 * t168;
t134 = t138 * t168 + t139 * t167;
t133 = t138 * t165 - t139 * t166;
t132 = t138 * t166 + t139 * t165;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t191 - rSges(2,2) * t186) + g(2) * (rSges(2,1) * t186 + rSges(2,2) * t191) + g(3) * (pkin(12) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,3) * t186 + t191 * t229 + t179) + g(2) * (-rSges(3,3) * t191 + t186 * t229 + t178) + g(3) * (rSges(3,1) * t185 + rSges(3,2) * t190 + pkin(12))) - m(4) * (g(1) * (rSges(4,1) * t168 + rSges(4,2) * t167 + rSges(4,3) * t186 + t241) + g(2) * (rSges(4,1) * t166 + rSges(4,2) * t165 - rSges(4,3) * t191 + t242) + g(3) * (-rSges(4,1) * t227 + rSges(4,2) * t172 + t248)) - m(5) * (g(1) * (rSges(5,1) * t134 + rSges(5,2) * t135 + rSges(5,3) * t186 + t234) + g(2) * (rSges(5,1) * t132 + rSges(5,2) * t133 - rSges(5,3) * t191 + t235) + g(3) * (rSges(5,1) * t136 + rSges(5,2) * t137 + t236)) - m(6) * (g(1) * (t134 * pkin(9) + (t134 * t188 + t183 * t186) * rSges(6,1) + (-t134 * t183 + t186 * t188) * rSges(6,2) + t254 * t135 + t234) + g(2) * (t132 * pkin(9) + (t132 * t188 - t183 * t191) * rSges(6,1) + (-t132 * t183 - t188 * t191) * rSges(6,2) + t254 * t133 + t235) + (t236 + (rSges(6,1) * t188 - rSges(6,2) * t183 + pkin(9)) * t136 + t254 * t137) * g(3)) - m(7) * (g(3) * (rSges(7,1) * t155 + rSges(7,2) * t154 + pkin(12) + pkin(14)) + (-g(2) * rSges(7,3) + g(1) * t226) * t191 + (g(1) * rSges(7,3) + g(2) * t226) * t186) - m(8) * (g(1) * (rSges(8,1) * t148 - rSges(8,2) * t149 + rSges(8,3) * t186 + t241) + g(2) * (rSges(8,1) * t146 - rSges(8,2) * t147 - rSges(8,3) * t191 + t242) + g(3) * (rSges(8,1) * t150 + rSges(8,2) * t151 + t248)) - m(9) * (g(1) * ((-t140 * t149 + t141 * t148) * rSges(9,1) + (-t140 * t148 - t141 * t149) * rSges(9,2) + t186 * rSges(9,3) + t241) + g(2) * ((-t140 * t147 + t141 * t146) * rSges(9,1) + (-t140 * t146 - t141 * t147) * rSges(9,2) - t191 * rSges(9,3) + t242) + g(3) * ((t140 * t151 + t141 * t150) * rSges(9,1) + (-t140 * t150 + t141 * t151) * rSges(9,2) + t248) + (g(1) * (t148 * t182 + t149 * t181) + g(2) * (t146 * t182 + t147 * t181) + g(3) * (t150 * t182 - t151 * t181)) * pkin(3));
U = t1;
