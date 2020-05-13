% Calculate joint inertia matrix for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m2OL_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_inertiaJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2OL_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'palh2m2OL_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'palh2m2OL_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:43:54
% EndTime: 2020-05-03 01:43:55
% DurationCPUTime: 1.45s
% Computational Cost: add. (1323->300), mult. (1827->373), div. (0->0), fcn. (491->58), ass. (0->211)
t251 = -qJ(6) + qJ(3);
t233 = qJ(4) + t251;
t128 = qJ(2) + t233;
t100 = qJ(5) + t128;
t260 = pkin(3) * rSges(7,2);
t279 = rSges(7,1) * rSges(7,3);
t314 = ((t260 - t279) * m(7) + Icges(7,5)) * sin(t100) / 0.2e1;
t164 = qJ(5) - qJ(6);
t148 = qJ(4) + t164;
t263 = pkin(2) * rSges(7,1);
t238 = m(7) * t263;
t63 = cos(t148) * t238;
t262 = pkin(2) * rSges(7,2);
t64 = m(7) * sin(t148) * t262;
t313 = t63 / 0.2e1 + t64 / 0.2e1;
t179 = cos(qJ(4));
t249 = pkin(2) * t179;
t283 = rSges(5,1) * m(5);
t82 = (m(6) + m(7)) * pkin(5) + t283;
t38 = t82 * t249;
t287 = m(6) * rSges(6,1);
t105 = pkin(3) * m(7) + t287;
t165 = qJ(4) + qJ(5);
t259 = t105 * cos(t165);
t39 = pkin(2) * t259;
t311 = 0.2e1 * t38 + 0.2e1 * t39;
t151 = qJ(2) + t251;
t270 = pkin(2) * m(7);
t234 = t270 / 0.2e1;
t223 = rSges(7,2) * t234;
t310 = sin(t151) * t223 + rSges(7,1) * cos(t151) * t234;
t161 = qJ(6) + qJ(2);
t150 = qJ(3) + t161;
t127 = qJ(4) + t150;
t99 = qJ(5) + t127;
t309 = ((t260 + t279) * m(7) - Icges(7,5)) * sin(t99);
t166 = qJ(3) + qJ(4);
t152 = qJ(2) + t166;
t113 = sin(t152);
t118 = cos(t152);
t284 = m(7) * rSges(7,2);
t236 = t284 / 0.2e1;
t227 = pkin(5) * t236;
t243 = pkin(5) * rSges(7,1) * cos(t127);
t285 = m(7) * rSges(7,1);
t237 = t285 / 0.2e1;
t272 = sin(t128) * t227 + pkin(5) * cos(t128) * t237;
t286 = m(6) * rSges(6,3);
t288 = m(5) * rSges(5,3);
t123 = rSges(7,1) * t284 - Icges(7,4);
t192 = -2 * Icges(7,6);
t170 = rSges(7,2) * rSges(7,3);
t261 = pkin(3) * rSges(7,1);
t276 = (-t170 + t261) * m(7);
t277 = (t170 + t261) * m(7);
t194 = rSges(7,2) ^ 2;
t196 = rSges(7,1) ^ 2;
t76 = (-t194 + t196) * m(7) - Icges(7,1) + Icges(7,2);
t80 = cos(t99);
t81 = cos(t100);
t149 = qJ(3) + t165;
t131 = qJ(2) + t149;
t92 = sin(t131);
t96 = cos(t131);
t191 = 2 * qJ(6);
t97 = t191 + t131;
t98 = -(2 * qJ(6)) + t131;
t307 = (t192 + 0.2e1 * t277) * t81 / 0.4e1 + (t192 - 0.2e1 * t276) * t80 / 0.4e1 + (-rSges(6,2) * t286 + Icges(6,6)) * t96 + (cos(t98) / 0.4e1 - cos(t97) / 0.4e1) * t76 + (sin(t97) + sin(t98)) * t123 / 0.2e1 + t309 / 0.2e1 + t314 - t92 * (rSges(6,1) * t286 - Icges(6,5));
t89 = sin(t127);
t308 = t89 * t227 + (-pkin(5) * t286 - rSges(5,3) * t283 + Icges(5,5)) * t113 - m(7) * t243 / 0.2e1 + t272 + t307 - t118 * (rSges(5,2) * t288 - Icges(5,6));
t306 = t39 + t313;
t246 = pkin(5) * t284;
t83 = sin(t164) * t246;
t245 = pkin(5) * t285;
t84 = cos(t164) * t245;
t163 = qJ(5) + qJ(6);
t85 = cos(t163) * t245;
t305 = t83 + t84 + t85;
t304 = t83 / 0.2e1 + t84 / 0.2e1;
t108 = m(6) * rSges(6,2) + rSges(7,3) * m(7);
t172 = sin(qJ(5));
t256 = t108 * t172;
t297 = -0.2e1 * pkin(5);
t303 = t256 * t297 + Icges(5,3) + (rSges(5,1) ^ 2 + rSges(5,2) ^ 2) * m(5);
t171 = sin(qJ(6));
t177 = cos(qJ(6));
t68 = t177 * rSges(7,1) - rSges(7,2) * t171;
t252 = pkin(3) + t68;
t302 = -rSges(7,3) * t92 + t252 * t96;
t226 = pkin(4) * t237;
t147 = qJ(4) + t163;
t129 = qJ(3) + t147;
t90 = sin(t129);
t94 = cos(t129);
t273 = -pkin(4) * t90 * t284 / 0.2e1 + t94 * t226;
t34 = t105 * pkin(4) * cos(t149);
t130 = qJ(5) + t233;
t91 = sin(t130);
t95 = cos(t130);
t301 = pkin(4) * t91 * t236 + t95 * t226 + t273 + t34;
t31 = pkin(4) * t82 * cos(t166);
t300 = t31 + t301;
t230 = 0.2e1 * m(7) * t261;
t296 = -0.2e1 * t171;
t264 = m(7) * t260 * t296 + Icges(6,3);
t299 = t177 * t230 + Icges(7,2) / 0.2e1 + Icges(7,1) / 0.2e1 + t264 - t123 * sin(t191) + t76 * cos(t191) / 0.2e1;
t185 = m(5) + m(6);
t190 = pkin(2) ^ 2;
t173 = sin(qJ(4));
t289 = m(5) * rSges(5,2);
t240 = t173 * t289;
t136 = sin(t165);
t250 = pkin(2) * t136;
t65 = cos(t147) * t238;
t298 = -0.2e1 * pkin(2) * t240 - 0.2e1 * t108 * t250 + Icges(4,3) + (rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * m(4) + t185 * t190 + t63 + t64 + t65;
t295 = m(7) / 0.2e1;
t294 = pkin(4) * m(7);
t198 = pkin(5) ^ 2;
t188 = 0.2e1 * t198;
t153 = t194 + t196;
t189 = pkin(3) ^ 2;
t168 = 0.2e1 * t189;
t193 = rSges(7,3) ^ 2;
t186 = 0.2e1 * t193;
t222 = t168 + t186 + t153;
t216 = t188 + t222;
t293 = t190 + t216 / 0.2e1;
t174 = sin(qJ(3));
t175 = sin(qJ(2));
t180 = cos(qJ(3));
t181 = cos(qJ(2));
t292 = pkin(5) * (t181 * (-t174 * t173 + t180 * t179) - t175 * (t180 * t173 + t174 * t179));
t291 = m(3) * rSges(3,3);
t290 = m(4) * rSges(4,3);
t182 = cos(qJ(1));
t21 = (t180 * pkin(2) + pkin(4)) * t181 - t175 * t174 * pkin(2) + pkin(1);
t17 = t21 * t182;
t281 = t182 * t292 + t17;
t280 = pkin(5) * t172;
t155 = rSges(6,1) ^ 2 + rSges(6,2) ^ 2;
t275 = t155 * m(6);
t109 = sin(t147);
t235 = -t270 / 0.2e1;
t271 = rSges(7,2) * t109 * t235 + t65 / 0.2e1;
t269 = t171 * rSges(7,1);
t268 = t177 * rSges(7,2);
t178 = cos(qJ(5));
t267 = t178 * (-t252 * m(7) - t287);
t266 = (m(4) * rSges(4,1) + (m(7) + t185) * pkin(2)) * t180;
t258 = t105 * t178;
t111 = sin(t149);
t257 = t108 * t111;
t122 = m(7) * t170 - Icges(7,6);
t254 = t122 * t177;
t244 = 0.2e1 * t96;
t242 = m(4) * rSges(4,2) * t174;
t241 = sin(t166) * t289;
t134 = sin(t163);
t239 = t134 * t284;
t232 = -t21 - t292;
t231 = pkin(5) * t239;
t228 = t153 * m(7) + Icges(7,3);
t225 = -t231 / 0.2e1 + t85 / 0.2e1 + t271;
t224 = t171 * (m(7) * t279 - Icges(7,5)) + t254;
t221 = -t286 - t288;
t220 = rSges(6,1) * t96 - rSges(6,2) * t92;
t217 = rSges(5,1) * t118 - rSges(5,2) * t113;
t215 = -t250 - t280;
t214 = -t241 - t257;
t167 = qJ(2) + qJ(3);
t138 = sin(t167);
t145 = cos(t167);
t213 = t181 * pkin(4) + rSges(4,1) * t145 - rSges(4,2) * t138 + pkin(1);
t212 = rSges(3,1) * t181 - rSges(3,2) * t175 + pkin(1);
t211 = -pkin(2) * t173 * t178 - (pkin(5) + t249) * t172;
t210 = (-pkin(2) * t109 - pkin(5) * t134) * rSges(7,2);
t209 = t275 + t299;
t207 = t222 * t295 + t209;
t187 = m(6) * t198;
t206 = t187 + t209 + t303;
t62 = pkin(5) * t258;
t58 = 0.2e1 * t62;
t204 = t206 + t58 + t305;
t203 = t62 + t207 + t225 + t304 + t306;
t202 = t298 + t299 + t303 + t305;
t112 = sin(t150);
t117 = cos(t150);
t201 = -t145 * (rSges(4,2) * t290 - Icges(4,6)) + (t221 * pkin(2) - rSges(4,1) * t290 + Icges(4,5)) * t138 + t112 * t223 + rSges(7,1) * t117 * t235 + t308 + t310;
t40 = t216 * t295;
t200 = t40 + t38 - t231 + (-t108 * t136 - t240) * pkin(2) + t204 + t271 + t306;
t199 = pkin(4) ^ 2;
t176 = sin(qJ(1));
t162 = -qJ(6) + qJ(2);
t140 = cos(t162);
t139 = cos(t161);
t133 = sin(t162);
t132 = sin(t161);
t69 = t182 * rSges(2,1) - t176 * rSges(2,2);
t67 = t268 + t269;
t66 = -t176 * rSges(2,1) - t182 * rSges(2,2);
t12 = t176 * rSges(3,3) + t212 * t182;
t11 = t182 * rSges(3,3) - t212 * t176;
t9 = t176 * rSges(4,3) + t213 * t182;
t8 = t182 * rSges(4,3) - t213 * t176;
t6 = t176 * rSges(5,3) + t217 * t182 + t17;
t5 = t182 * rSges(5,3) + (-t21 - t217) * t176;
t4 = t176 * rSges(6,3) + t220 * t182 + t281;
t3 = t182 * rSges(6,3) + (-t220 + t232) * t176;
t2 = -t176 * t67 + t302 * t182 + t281;
t1 = -t182 * t67 + (t232 - t302) * t176;
t7 = [Icges(3,2) * t181 ^ 2 + Icges(4,2) * t145 ^ 2 + Icges(5,2) * t118 ^ 2 + Icges(2,3) + (Icges(3,1) * t175 + 0.2e1 * Icges(3,4) * t181) * t175 + (Icges(4,1) * t138 + 0.2e1 * Icges(4,4) * t145) * t138 + (Icges(5,1) * t113 + 0.2e1 * Icges(5,4) * t118) * t113 + (Icges(7,3) + Icges(6,2)) * t96 ^ 2 + ((Icges(7,1) * t177 ^ 2 + Icges(6,1) + (-0.2e1 * Icges(7,4) * t177 + Icges(7,2) * t171) * t171) * t92 + (Icges(7,5) * t177 - Icges(7,6) * t171 + Icges(6,4)) * t244) * t92 + m(7) * (t1 ^ 2 + t2 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2) + m(3) * (t11 ^ 2 + t12 ^ 2) + m(2) * (t66 ^ 2 + t69 ^ 2) + m(5) * (t5 ^ 2 + t6 ^ 2) + m(4) * (t8 ^ 2 + t9 ^ 2); t201 + ((t221 - t290) * pkin(4) - rSges(3,1) * t291 + Icges(3,5)) * t175 + ((t132 / 0.2e1 + t133 / 0.2e1) * rSges(7,2) + (-t139 / 0.2e1 + t140 / 0.2e1) * rSges(7,1)) * t294 - t181 * (rSges(3,2) * t291 - Icges(3,6)); t202 + (m(4) + t185) * t199 + 0.2e1 * t34 + (-0.2e1 * t242 - 0.2e1 * t241 - 0.2e1 * t257 + 0.2e1 * t266 + ((-t90 + t91) * rSges(7,2) + (t94 + t95) * rSges(7,1)) * m(7)) * pkin(4) + 0.2e1 * t31 + t58 + (t168 / 0.2e1 + t190 + t186 / 0.2e1 + t188 / 0.2e1 + t194 / 0.2e1 + t196 / 0.2e1 + t199 + t210) * m(7) + (t198 + t155) * m(6) + Icges(3,3) + (rSges(3,1) ^ 2 + rSges(3,2) ^ 2) * m(3) + t311; t201; (t214 - t242 + t266) * pkin(4) + t298 + (t293 + t210) * m(7) + t204 + t300 + t311; t187 + t202 + m(7) * t293 + (-t109 * t284 + 0.2e1 * t82 * t179 + 0.2e1 * t259) * pkin(2) + (-t239 + 0.2e1 * t258) * pkin(5) + t275; t308; t214 * pkin(4) + t200 + t300; t200; t267 * t297 + t206 + t40; t307; (-pkin(4) * t111 + t215) * t108 + t203 + t301; t215 * t108 + t203; (-t256 - t267) * pkin(5) + t207; (t189 + t193 + t194) * m(7) + Icges(7,1) + t275 + (t123 * t296 + t76 * t177 + t230) * t177 + t264; (-Icges(7,6) + t277) * t81 / 0.2e1 + t314 + (Icges(7,6) + t276) * t80 / 0.2e1 - t309 / 0.2e1 + t228 * t244 / 0.2e1 + (-pkin(5) * rSges(7,2) * t89 + 0.2e1 * pkin(1) * t68 - t112 * t262 + t117 * t263 + t243) * t295 + ((-t132 + t133) * rSges(7,2) + (t139 + t140) * rSges(7,1)) * t294 / 0.2e1 + t272 + t310; (-rSges(7,1) * t95 / 0.2e1 - rSges(7,2) * t91 / 0.2e1) * t294 + t224 + t225 + t273 - t304 - t313; t254 - t171 * Icges(7,5) + (t211 * t268 - (-rSges(7,3) - t211) * t269) * m(7); (-t172 * t246 + t122) * t177 - t171 * (Icges(7,5) + (-rSges(7,3) + t280) * t285); t224; t228;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11), t7(16); t7(2), t7(3), t7(5), t7(8), t7(12), t7(17); t7(4), t7(5), t7(6), t7(9), t7(13), t7(18); t7(7), t7(8), t7(9), t7(10), t7(14), t7(19); t7(11), t7(12), t7(13), t7(14), t7(15), t7(20); t7(16), t7(17), t7(18), t7(19), t7(20), t7(21);];
Mq = res;
