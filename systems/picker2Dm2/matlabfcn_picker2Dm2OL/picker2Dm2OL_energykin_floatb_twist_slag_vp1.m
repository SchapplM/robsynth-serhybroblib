% Calculate kinetic energy for
% picker2Dm2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% qJD [12x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [11x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 23:20
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = picker2Dm2OL_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(6,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2OL_energykin_floatb_twist_slag_vp1: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm2OL_energykin_floatb_twist_slag_vp1: qJD has to be [12x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'picker2Dm2OL_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2OL_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2OL_energykin_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm2OL_energykin_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'picker2Dm2OL_energykin_floatb_twist_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 23:18:48
% EndTime: 2020-05-09 23:18:51
% DurationCPUTime: 2.50s
% Computational Cost: add. (1273->344), mult. (800->446), div. (0->0), fcn. (456->22), ass. (0->185)
t227 = qJ(1) + qJ(2);
t223 = qJ(4) + t227;
t210 = qJ(10) + t223;
t196 = cos(t210);
t279 = t196 / 0.2e1;
t224 = qJ(3) + t227;
t211 = qJ(9) + t224;
t198 = cos(t211);
t278 = -t198 / 0.2e1;
t222 = qJ(6) + t227;
t207 = cos(t222);
t277 = t207 / 0.2e1;
t208 = cos(t223);
t276 = -t208 / 0.2e1;
t209 = cos(t224);
t275 = t209 / 0.2e1;
t225 = pkin(8) + qJ(5);
t213 = cos(t225);
t274 = t213 / 0.2e1;
t226 = qJ(1) + qJ(8);
t219 = cos(t226);
t273 = t219 / 0.2e1;
t220 = cos(t227);
t272 = -t220 / 0.2e1;
t230 = sin(qJ(7));
t271 = t230 / 0.2e1;
t233 = cos(qJ(1));
t270 = -t233 / 0.2e1;
t231 = sin(qJ(1));
t269 = pkin(1) * t231;
t268 = pkin(1) * t233;
t205 = sin(t223);
t267 = pkin(4) * t205;
t266 = pkin(4) * t208;
t206 = sin(t224);
t265 = pkin(6) * t206;
t264 = pkin(6) * t209;
t263 = Icges(2,4) * t231;
t262 = Icges(2,4) * t233;
t218 = sin(t227);
t261 = Icges(3,4) * t218;
t260 = Icges(3,4) * t220;
t259 = Icges(4,4) * t206;
t258 = Icges(5,4) * t205;
t257 = Icges(5,4) * t208;
t212 = sin(t225);
t256 = Icges(6,4) * t212;
t204 = sin(t222);
t255 = Icges(7,4) * t204;
t232 = cos(qJ(7));
t254 = Icges(8,4) * t232;
t217 = sin(t226);
t253 = Icges(9,4) * t217;
t197 = sin(t211);
t252 = Icges(10,4) * t197;
t251 = Icges(10,4) * t198;
t216 = V_base(6) + qJD(1);
t203 = qJD(2) + t216;
t250 = t203 * t218;
t249 = t203 * t220;
t195 = sin(t210);
t248 = Icges(11,4) * t195;
t247 = t216 * t269 + V_base(1);
t246 = V_base(5) * t268 + V_base(3);
t244 = pkin(5) * V_base(6);
t243 = V_base(5) * t220;
t242 = pkin(3) * t250 + t247;
t241 = pkin(2) * t250 + t247;
t240 = pkin(3) * t243 + t246;
t239 = pkin(2) * t243 + t246;
t238 = -pkin(2) * t218 - t269;
t237 = -pkin(3) * t218 - t269;
t236 = -t216 * t268 + V_base(2);
t194 = qJD(3) + t203;
t193 = qJD(4) + t203;
t235 = -pkin(2) * t249 + t236;
t234 = -pkin(3) * t249 + t236;
t229 = cos(pkin(8));
t228 = sin(pkin(8));
t221 = Icges(8,4) * t230;
t215 = V_base(6) + qJD(5);
t214 = V_base(6) + qJD(7);
t202 = qJD(8) + t216;
t201 = Icges(9,4) * t219;
t200 = Icges(6,4) * t213;
t192 = qJD(6) + t203;
t191 = Icges(4,4) * t209;
t190 = Icges(7,4) * t207;
t186 = qJD(9) + t194;
t185 = qJD(10) + t193;
t184 = Icges(11,4) * t196;
t183 = -rSges(2,1) * t233 + t231 * rSges(2,2);
t182 = -rSges(8,1) * t232 + rSges(8,2) * t230;
t181 = -t231 * rSges(2,1) - rSges(2,2) * t233;
t180 = rSges(8,1) * t230 + rSges(8,2) * t232;
t179 = -Icges(2,1) * t233 + t263;
t178 = -Icges(2,1) * t231 - t262;
t177 = -Icges(8,1) * t232 + t221;
t176 = Icges(8,1) * t230 + t254;
t175 = Icges(2,2) * t231 - t262;
t174 = -Icges(2,2) * t233 - t263;
t173 = Icges(8,2) * t230 - t254;
t172 = Icges(8,2) * t232 + t221;
t167 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t166 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t165 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t162 = -rSges(3,1) * t220 + rSges(3,2) * t218;
t161 = rSges(9,1) * t219 - rSges(9,2) * t217;
t160 = -rSges(3,1) * t218 - rSges(3,2) * t220;
t159 = rSges(9,1) * t217 + rSges(9,2) * t219;
t158 = -Icges(3,1) * t220 + t261;
t157 = -Icges(3,1) * t218 - t260;
t156 = Icges(9,1) * t219 - t253;
t155 = Icges(9,1) * t217 + t201;
t154 = Icges(3,2) * t218 - t260;
t153 = -Icges(3,2) * t220 - t261;
t152 = -Icges(9,2) * t217 + t201;
t151 = Icges(9,2) * t219 + t253;
t146 = rSges(6,1) * t213 - rSges(6,2) * t212;
t145 = rSges(6,1) * t212 + rSges(6,2) * t213;
t144 = Icges(6,1) * t213 - t256;
t143 = Icges(6,1) * t212 + t200;
t142 = -Icges(6,2) * t212 + t200;
t141 = Icges(6,2) * t213 + t256;
t138 = rSges(4,1) * t209 - rSges(4,2) * t206;
t137 = -rSges(5,1) * t208 + rSges(5,2) * t205;
t136 = rSges(7,1) * t207 - rSges(7,2) * t204;
t135 = rSges(4,1) * t206 + rSges(4,2) * t209;
t134 = -rSges(5,1) * t205 - rSges(5,2) * t208;
t133 = rSges(7,1) * t204 + rSges(7,2) * t207;
t132 = Icges(4,1) * t209 - t259;
t131 = Icges(4,1) * t206 + t191;
t130 = -Icges(5,1) * t208 + t258;
t129 = -Icges(5,1) * t205 - t257;
t128 = Icges(7,1) * t207 - t255;
t127 = Icges(7,1) * t204 + t190;
t126 = -Icges(4,2) * t206 + t191;
t125 = Icges(4,2) * t209 + t259;
t124 = Icges(5,2) * t205 - t257;
t123 = -Icges(5,2) * t208 - t258;
t122 = -Icges(7,2) * t204 + t190;
t121 = Icges(7,2) * t207 + t255;
t114 = -rSges(10,1) * t198 + rSges(10,2) * t197;
t113 = -rSges(10,1) * t197 - rSges(10,2) * t198;
t112 = -Icges(10,1) * t198 + t252;
t111 = -Icges(10,1) * t197 - t251;
t110 = Icges(10,2) * t197 - t251;
t109 = -Icges(10,2) * t198 - t252;
t106 = rSges(11,1) * t196 - rSges(11,2) * t195;
t105 = rSges(11,1) * t195 + rSges(11,2) * t196;
t104 = Icges(11,1) * t196 - t248;
t103 = Icges(11,1) * t195 + t184;
t102 = -Icges(11,2) * t195 + t184;
t101 = Icges(11,2) * t196 + t248;
t98 = V_base(5) * rSges(2,3) - t181 * t216 + V_base(1);
t97 = V_base(5) * rSges(8,3) - t182 * t214 + V_base(1);
t96 = -V_base(4) * rSges(2,3) + t183 * t216 + V_base(2);
t95 = -V_base(4) * rSges(8,3) + V_base(6) * pkin(7) + t180 * t214 + V_base(2);
t94 = t181 * V_base(4) - t183 * V_base(5) + V_base(3);
t93 = t182 * V_base(4) + V_base(3) + (-pkin(7) - t180) * V_base(5);
t92 = V_base(5) * rSges(6,3) - t145 * t215 - t228 * t244 + V_base(1);
t91 = -V_base(4) * rSges(6,3) + t146 * t215 + t229 * t244 + V_base(2);
t90 = V_base(5) * rSges(3,3) - t160 * t203 + t247;
t89 = V_base(5) * rSges(9,3) - t159 * t202 + t247;
t88 = -V_base(4) * rSges(3,3) + t203 * t162 + t236;
t87 = -V_base(4) * rSges(9,3) + t202 * t161 + t236;
t86 = V_base(5) * rSges(7,3) - t133 * t192 + t247;
t85 = -V_base(4) * rSges(7,3) + t192 * t136 + t236;
t84 = -t162 * V_base(5) + (t160 - t269) * V_base(4) + t246;
t83 = -t161 * V_base(5) + (t159 - t269) * V_base(4) + t246;
t82 = t145 * V_base(4) - t146 * V_base(5) + V_base(3) + (t228 * V_base(4) - t229 * V_base(5)) * pkin(5);
t81 = V_base(5) * rSges(4,3) - t135 * t194 + t241;
t80 = V_base(5) * rSges(5,3) - t134 * t193 + t242;
t79 = -V_base(4) * rSges(4,3) + t194 * t138 + t235;
t78 = -V_base(4) * rSges(5,3) + t193 * t137 + t234;
t77 = -t136 * V_base(5) + (t133 - t269) * V_base(4) + t246;
t76 = V_base(5) * rSges(10,3) - t113 * t186 - t194 * t265 + t241;
t75 = -V_base(4) * rSges(10,3) + t186 * t114 + t194 * t264 + t235;
t74 = V_base(5) * rSges(11,3) - t105 * t185 + t193 * t267 + t242;
t73 = -V_base(4) * rSges(11,3) + t185 * t106 - t193 * t266 + t234;
t72 = -t138 * V_base(5) + (t135 + t238) * V_base(4) + t239;
t71 = -t137 * V_base(5) + (t134 + t237) * V_base(4) + t240;
t70 = (-t114 - t264) * V_base(5) + (t113 + t238 + t265) * V_base(4) + t239;
t69 = (-t106 + t266) * V_base(5) + (t105 + t237 - t267) * V_base(4) + t240;
t1 = (Icges(1,6) * V_base(6) + (-Icges(10,5) * t197 - Icges(10,6) * t198) * t186 + (Icges(7,5) * t204 + Icges(7,6) * t207) * t192 + (-Icges(5,5) * t205 - Icges(5,6) * t208) * t193 + (Icges(4,5) * t206 + Icges(4,6) * t209) * t194 + (Icges(6,5) * t212 + Icges(6,6) * t213) * t215 + (Icges(9,5) * t217 + Icges(9,6) * t219) * t202 + (-Icges(3,5) * t218 - Icges(3,6) * t220) * t203 + (-Icges(8,5) * t232 + Icges(8,6) * t230) * t214 + (-Icges(2,5) * t231 - Icges(2,6) * t233) * t216 + (Icges(11,5) * t195 + Icges(11,6) * t196) * t185 + (Icges(1,2) / 0.2e1 + t151 * t273 + t155 * t217 / 0.2e1 + t153 * t272 - t157 * t218 / 0.2e1 + t141 * t274 + t143 * t212 / 0.2e1 + t125 * t275 + t131 * t206 / 0.2e1 + t121 * t277 + t127 * t204 / 0.2e1 + t123 * t276 - t129 * t205 / 0.2e1 + t101 * t279 + t103 * t195 / 0.2e1 + t109 * t278 - t111 * t197 / 0.2e1 + t174 * t270 - t231 * t178 / 0.2e1 + t173 * t271 - t177 * t232 / 0.2e1) * V_base(5)) * V_base(5) + Icges(11,3) * t185 ^ 2 / 0.2e1 + m(11) * (t69 ^ 2 + t73 ^ 2 + t74 ^ 2) / 0.2e1 + m(10) * (t70 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(5) * (t71 ^ 2 + t78 ^ 2 + t80 ^ 2) / 0.2e1 + m(4) * (t72 ^ 2 + t79 ^ 2 + t81 ^ 2) / 0.2e1 + Icges(5,3) * t193 ^ 2 / 0.2e1 + Icges(6,3) * t215 ^ 2 / 0.2e1 + Icges(3,3) * t203 ^ 2 / 0.2e1 + m(1) * (t165 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + m(8) * (t93 ^ 2 + t95 ^ 2 + t97 ^ 2) / 0.2e1 + m(2) * (t94 ^ 2 + t96 ^ 2 + t98 ^ 2) / 0.2e1 + m(7) * (t77 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(9) * (t83 ^ 2 + t87 ^ 2 + t89 ^ 2) / 0.2e1 + m(3) * (t84 ^ 2 + t88 ^ 2 + t90 ^ 2) / 0.2e1 + m(6) * (t82 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + Icges(9,3) * t202 ^ 2 / 0.2e1 + Icges(7,3) * t192 ^ 2 / 0.2e1 + Icges(10,3) * t186 ^ 2 / 0.2e1 + Icges(2,3) * t216 ^ 2 / 0.2e1 + Icges(8,3) * t214 ^ 2 / 0.2e1 + Icges(4,3) * t194 ^ 2 / 0.2e1 + Icges(1,3) * V_base(6) ^ 2 / 0.2e1 + (((-t175 - t178) * t233 + (t173 - t176) * t232 + (-t179 + t174) * t231 + (t177 + t172) * t230 + (-t157 - t154) * t220 + (t155 + t152) * t219 + (t153 - t158) * t218 + (-t151 + t156) * t217 + (t143 + t142) * t213 + (-t141 + t144) * t212 + (t131 + t126) * t209 + (-t124 - t129) * t208 + (t122 + t127) * t207 + (-t125 + t132) * t206 + (-t130 + t123) * t205 + (t128 - t121) * t204 + (-t111 - t110) * t198 + (t109 - t112) * t197 + (t102 + t103) * t196 + (t104 - t101) * t195) * V_base(5) / 0.2e1 + Icges(1,5) * V_base(6) + Icges(1,4) * V_base(5) + (Icges(11,5) * t196 - Icges(11,6) * t195) * t185 + (-Icges(10,5) * t198 + Icges(10,6) * t197) * t186 + (Icges(7,5) * t207 - Icges(7,6) * t204) * t192 + (-Icges(5,5) * t208 + Icges(5,6) * t205) * t193 + (Icges(4,5) * t209 - Icges(4,6) * t206) * t194 + (Icges(6,5) * t213 - Icges(6,6) * t212) * t215 + (Icges(9,5) * t219 - Icges(9,6) * t217) * t202 + (-Icges(3,5) * t220 + Icges(3,6) * t218) * t203 + (Icges(8,5) * t230 + Icges(8,6) * t232) * t214 + (-Icges(2,5) * t233 + Icges(2,6) * t231) * t216 + (Icges(1,1) / 0.2e1 + t154 * t218 / 0.2e1 + t158 * t272 - t142 * t212 / 0.2e1 + t144 * t274 - t152 * t217 / 0.2e1 + t156 * t273 + t124 * t205 / 0.2e1 + t130 * t276 - t126 * t206 / 0.2e1 + t132 * t275 - t122 * t204 / 0.2e1 + t128 * t277 - t102 * t195 / 0.2e1 + t104 * t279 + t110 * t197 / 0.2e1 + t112 * t278 + t231 * t175 / 0.2e1 + t179 * t270 + t172 * t232 / 0.2e1 + t176 * t271) * V_base(4)) * V_base(4);
T = t1;
