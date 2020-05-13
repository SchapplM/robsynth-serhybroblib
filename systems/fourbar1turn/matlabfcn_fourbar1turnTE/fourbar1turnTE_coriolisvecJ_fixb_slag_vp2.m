% Calculate vector of centrifugal and Coriolis load on the joints for
% fourbar1turnTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [2x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = fourbar1turnTE_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_coriolisvecJ_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnTE_coriolisvecJ_fixb_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_coriolisvecJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnTE_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnTE_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnTE_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:18:35
% EndTime: 2020-04-12 19:18:59
% DurationCPUTime: 9.98s
% Computational Cost: add. (90917->345), mult. (128782->778), div. (3868->15), fcn. (34391->4), ass. (0->280)
t121 = pkin(2) ^ 2;
t122 = pkin(1) ^ 2;
t112 = cos(qJ(2));
t296 = pkin(2) * t112;
t253 = -0.2e1 * pkin(1) * t296 + t122;
t105 = t121 + t253;
t340 = pkin(3) ^ 2;
t252 = -pkin(4) ^ 2 + t340;
t100 = t105 + t252;
t107 = pkin(1) * t112 - pkin(2);
t111 = sin(qJ(2));
t318 = -pkin(3) - pkin(4);
t98 = (pkin(2) - t318) * (pkin(2) + t318) + t253;
t317 = pkin(4) - pkin(3);
t99 = (pkin(2) - t317) * (pkin(2) + t317) + t253;
t278 = t98 * t99;
t123 = sqrt(-t278);
t257 = t111 * t123;
t78 = -pkin(1) * t257 - t100 * t107;
t128 = t78 ^ 2;
t71 = 0.1e1 / t128;
t302 = pkin(1) * t100;
t95 = t111 * t302;
t81 = -t107 * t123 + t95;
t285 = t71 * t81;
t102 = 0.1e1 / t105;
t118 = 0.1e1 / pkin(3);
t103 = 0.1e1 / t105 ^ 2;
t262 = t103 * t111;
t227 = pkin(1) * t262;
t188 = qJD(2) * t227;
t168 = pkin(2) * t188;
t160 = -0.2e1 * t168;
t328 = -0.2e1 * t107;
t241 = pkin(2) * t328;
t198 = t100 + t241;
t191 = pkin(1) * pkin(2) * (-t98 - t99);
t86 = t111 * t191;
t85 = qJD(2) * t86;
t90 = 0.1e1 / t123;
t280 = t85 * t90;
t224 = t111 * t280;
t255 = t112 * t123;
t50 = (-t224 + (t111 * t198 - t255) * qJD(2)) * pkin(1);
t38 = (t102 * t50 + t160 * t78) * t118;
t110 = t111 ^ 2;
t260 = t110 * t122;
t222 = pkin(2) * t260;
t185 = qJD(2) * t222;
t270 = t107 * t90;
t246 = qJD(2) * t111;
t211 = t123 * t246;
t245 = qJD(2) * t112;
t275 = pkin(1) * t211 + t245 * t302;
t53 = -t270 * t85 + 0.2e1 * t185 + t275;
t41 = (t102 * t53 + t160 * t81) * t118;
t70 = 0.1e1 / t78;
t334 = -t38 * t285 + t41 * t70;
t77 = t81 ^ 2;
t67 = t71 * t77 + 0.1e1;
t64 = 0.1e1 / t67;
t144 = t64 * t334;
t342 = 0.2e1 * t144;
t101 = t105 - t252;
t106 = pkin(1) - t296;
t79 = -pkin(2) * t257 + t101 * t106;
t74 = 0.1e1 / t79 ^ 2;
t297 = pkin(2) * t111;
t94 = t101 * t297;
t80 = t106 * t123 + t94;
t283 = t74 * t80;
t116 = 0.1e1 / pkin(4);
t313 = t102 / 0.2e1;
t239 = 0.2e1 * t106 * pkin(1);
t51 = (-t224 + (-t255 + (t101 + t239) * t111) * qJD(2)) * pkin(2);
t134 = -t168 * t79 + t313 * t51;
t36 = t134 * t116;
t261 = t110 * t121;
t223 = pkin(1) * t261;
t186 = qJD(2) * t223;
t271 = t106 * t90;
t276 = (t101 * t245 + t211) * pkin(2);
t52 = t271 * t85 + 0.2e1 * t186 + t276;
t133 = -t168 * t80 + t313 * t52;
t40 = t133 * t116;
t76 = t80 ^ 2;
t66 = t74 * t76 + 0.1e1;
t62 = 0.1e1 / t66;
t73 = 0.1e1 / t79;
t146 = t62 * (t36 * t283 - t40 * t73);
t294 = pkin(4) * t105;
t240 = 0.2e1 * t294;
t339 = -0.4e1 * t79;
t219 = pkin(2) * t246;
t193 = pkin(1) * t219;
t338 = -0.2e1 * t193;
t243 = qJD(1) * qJD(2);
t206 = t111 * t243;
t337 = pkin(2) * t206;
t336 = t85 * t86;
t231 = t62 * t294;
t195 = t74 * t231;
t335 = t116 * t80 * t195;
t333 = -0.2e1 * t74;
t332 = 0.4e1 * t80;
t331 = 0.2e1 * t81;
t330 = t116 ^ 2;
t329 = -0.2e1 * t102;
t327 = qJD(1) ^ 2;
t114 = qJD(2) ^ 2;
t197 = pkin(2) * t227;
t314 = -t102 / 0.2e1;
t279 = t90 * t86;
t55 = t94 + (-t255 + (t239 - t279) * t111) * pkin(2);
t45 = (t197 * t79 + t314 * t55) * t116;
t59 = t86 * t271 + 0.2e1 * t223 + (t101 * t112 + t257) * pkin(2);
t48 = (-t197 * t80 + t313 * t59) * t116;
t145 = t62 * (-t283 * t45 - t48 * t73);
t63 = 0.1e1 / t66 ^ 2;
t207 = 0.4e1 * t63 * t283;
t286 = t63 * t73;
t75 = t73 * t74;
t282 = t75 * t76;
t29 = -t282 * t51 + t283 * t52;
t230 = t29 * t286;
t236 = t75 * t332;
t238 = 0.2e1 * t62 * t74;
t267 = 0.4e1 * t90 / t278;
t212 = -t267 / 0.4e1;
t166 = t212 * t336;
t312 = -t111 / 0.2e1;
t215 = t121 * t260;
t135 = t112 * t191 - 0.4e1 * t215;
t83 = t135 * qJD(2);
t132 = (t111 * t166 + 0.2e1 * (-t112 * t85 + t83 * t312) * t90) * t102;
t104 = t102 * t103;
t182 = t104 * t215;
t159 = t182 * t332;
t162 = qJD(2) * t182;
t169 = t116 * t73 * t231;
t258 = t111 * t121;
t214 = t112 * t258;
t172 = 0.6e1 * pkin(1) * t214;
t179 = t106 * t267 / 0.4e1;
t265 = t102 * t112;
t216 = t106 * t265;
t221 = t102 * t280;
t301 = pkin(1) * t103;
t277 = -0.2e1 * ((0.4e1 * t186 + t276) * t314 + t162 * t339 + (-t132 / 0.2e1 + (t51 * t262 + (-t216 + (t111 * t55 + t112 * t79) * t103) * qJD(2)) * pkin(1)) * pkin(2)) * t335 - 0.2e1 * ((qJD(2) * t172 + t179 * t336 + t271 * t83) * t313 + qJD(2) * t159 + ((t221 / 0.2e1 - t52 * t301) * t111 + ((t255 + (-t101 + t279) * t111) * t313 + (-t111 * t59 - t112 * t80) * t301) * qJD(2)) * pkin(2)) * t169;
t3 = (0.4e1 * t145 * t193 + ((t238 * t51 + 0.4e1 * t230) * t48 + (t29 * t207 + (t236 * t51 + t333 * t52) * t62) * t45) * t105) * pkin(4) + t277;
t326 = t3 / 0.2e1;
t264 = t102 * t116;
t319 = t80 / 0.2e1;
t320 = -t79 / 0.2e1;
t141 = (Ifges(5,1) * t319 + Ifges(5,4) * t320) * t264;
t20 = t146 * t240;
t17 = Ifges(5,5) * t20 + qJD(1) * t141;
t324 = -t17 / 0.2e1;
t321 = t78 / 0.2e1;
t311 = t111 / 0.2e1;
t310 = -t112 / 0.2e1;
t307 = Ifges(5,5) * t80;
t304 = Ifges(5,6) * t79;
t303 = Ifges(5,3) * t20;
t300 = pkin(1) * t118;
t299 = pkin(2) * t102;
t298 = pkin(2) * t103;
t295 = pkin(3) * t105;
t155 = 0.8e1 * t162;
t233 = 0.2e1 * t103;
t54 = t95 + (-t255 + (t241 - t279) * t111) * pkin(1);
t266 = t54 * t111;
t269 = t112 * t78;
t293 = (((0.4e1 * t185 + t275) * t102 + t78 * t155) * t118 + (t132 + (-0.2e1 * t50 * t262 + (t265 * t328 + (-t266 - t269) * t233) * qJD(2)) * pkin(2)) * t300) * t81;
t281 = t81 * t53;
t72 = t70 * t71;
t284 = t72 * t77;
t65 = 0.1e1 / t67 ^ 2;
t291 = (t281 * t71 - t284 * t50) * t65;
t263 = t102 * t118;
t208 = -t263 / 0.2e1;
t180 = t78 * t208;
t250 = qJD(1) * t111;
t202 = t250 / 0.2e1;
t249 = qJD(1) * t112;
t57 = t81 * t202 * t263 + t180 * t249;
t274 = Ifges(3,4) * t111;
t273 = Ifges(3,4) * t112;
t272 = pkin(1) * qJD(1);
t268 = t112 * t81;
t259 = t111 * t112;
t256 = t112 * t114;
t254 = t114 * t123;
t251 = qJD(1) * t102;
t248 = qJD(1) * t118;
t247 = qJD(2) * t103;
t244 = qJD(2) * t118;
t242 = -0.2e1 * t298;
t237 = t72 * t331;
t235 = pkin(1) * t298;
t234 = pkin(1) * t297;
t232 = t70 * t295;
t226 = pkin(1) * t258;
t225 = pkin(2) * t263;
t220 = pkin(2) * t249;
t218 = qJD(2) * t280;
t217 = -t299 / 0.2e1;
t213 = t112 * t254;
t210 = -t264 / 0.4e1;
t209 = t264 / 0.4e1;
t205 = t112 * t243;
t204 = -t251 / 0.4e1;
t203 = t251 / 0.4e1;
t194 = 0.2e1 * t218;
t192 = -0.2e1 * t214;
t190 = t103 * t226;
t189 = t114 * t226;
t187 = qJD(2) * t225;
t184 = pkin(2) * t208;
t84 = t85 ^ 2;
t183 = t84 * t212;
t181 = qJD(2) * t217;
t178 = -0.2e1 * t197;
t171 = 0.6e1 * pkin(2) * t122 * t259;
t170 = t118 * t64 * t232;
t167 = t103 * t189;
t165 = Ifges(3,5) * t112 - Ifges(3,6) * t111;
t164 = t118 * t181;
t163 = t187 / 0.2e1;
t47 = (t102 * t54 + t178 * t78) * t118;
t58 = -t86 * t270 + 0.2e1 * t222 + (t100 * t112 + t257) * pkin(1);
t49 = (t102 * t58 + t178 * t81) * t118;
t157 = t285 * t47 - t49 * t70;
t156 = t116 * t168;
t153 = t118 * t121 * t188;
t152 = t78 * t163;
t150 = t268 / 0.2e1 + t78 * t311;
t149 = -t101 * t114 + t194;
t148 = t111 * (Ifges(3,1) * t112 - t274);
t147 = (t112 * Ifges(3,2) + t274) * qJD(1);
t143 = t157 * t64;
t142 = (t81 * t311 - t269 / 0.2e1) * t102;
t140 = (Ifges(5,4) * t319 + Ifges(5,2) * t320) * t264;
t139 = (t307 / 0.2e1 - t304 / 0.2e1) * t264;
t60 = t150 * t263;
t82 = t135 * t114;
t138 = -t90 * t82 + t183 + t254;
t137 = (t110 * t78 + t259 * t81) * t235;
t136 = (-t110 * t81 + t259 * t78) * t235;
t25 = ((t310 * t53 + t312 * t50) * t102 + (t142 + t137) * qJD(2)) * t118;
t26 = ((t310 * t50 + t311 * t53) * t102 + (t102 * t150 + t136) * qJD(2)) * t118;
t109 = Ifges(3,4) * t249;
t97 = Ifges(3,1) * t250 + Ifges(3,5) * qJD(2) + t109;
t96 = Ifges(3,6) * qJD(2) + t147;
t56 = qJD(1) * t60;
t46 = qJD(1) * t48;
t44 = (t190 * t81 + t217 * t58) * t244;
t43 = qJD(1) * t45;
t42 = (t190 * t78 + t217 * t54) * t244;
t39 = (t167 * t81 + t181 * t53) * t118;
t37 = qJD(1) * t40;
t35 = (t167 * t78 + t181 * t50) * t118;
t34 = qJD(1) * t36;
t33 = -mrSges(4,1) * t57 - mrSges(4,2) * t56;
t31 = -t282 * t55 + t283 * t59;
t28 = (t136 + ((t81 / 0.2e1 - t54 / 0.2e1) * t112 + (t58 / 0.2e1 + t321) * t111) * t102) * t248;
t27 = ((-t266 / 0.2e1 + t58 * t310) * t102 + t137) * t248 + t57;
t24 = qJD(1) * t26;
t23 = qJD(1) * t25;
t19 = t144 * t295 + qJD(2);
t16 = Ifges(5,6) * t20 + qJD(1) * t140;
t12 = ((qJD(2) * t171 + t107 * t166 - t270 * t83) * t102 + t81 * t155 + ((t242 * t53 + t221) * t111 + ((t255 + (-t100 + t279) * t111) * t102 + (-t58 * t111 - t268) * pkin(2) * t233) * qJD(2)) * pkin(1)) * t170;
t9 = -t56 * Ifges(4,1) + t57 * Ifges(4,4) + t19 * Ifges(4,5);
t8 = -t56 * Ifges(4,4) + t57 * Ifges(4,2) + t19 * Ifges(4,6);
t6 = ((t107 * t183 + t114 * t171 - t82 * t270) * t102 + 0.8e1 * t81 * t114 * t182 + ((t213 + (-t100 * t114 + t194) * t111) * t102 + (-0.4e1 * t246 * t53 - 0.2e1 * t256 * t81) * t298) * pkin(1)) * t170 - 0.2e1 * t41 * t232 * t291 + pkin(3) * t193 * t342 + ((t50 * t64 * t72 + t291 * t71) * t38 * t331 + (-t41 * t50 - ((0.8e1 * t104 * t121 * t78 + 0.4e1 * t299) * t118 * t114 * t260 + ((t218 * t329 + (t102 * t198 + t242 * t78) * t114) * t112 + (-0.4e1 * pkin(2) * t247 * t50 + t102 * t138) * t111) * t300) * t81 - t38 * t53) * t64 * t71) * t295;
t5 = -0.2e1 * ((t114 * t172 + t84 * t179 + t82 * t271) * t313 + t114 * t159 + ((t111 * t149 + t213) * t313 + (-0.2e1 * t246 * t52 - t256 * t80) * t301) * pkin(2)) * t169 - 0.2e1 * ((t104 * t122 * t339 + pkin(1) * t329) * t114 * t261 + ((t111 * t138 - t112 * t149) * t314 + (-t114 * t216 + (0.2e1 * t246 * t51 + t256 * t79) * t103) * pkin(1)) * pkin(2)) * t335 + 0.2e1 * (t36 * t52 + t40 * t51) * t195 - 0.2e1 * pkin(4) * t338 * t146 + 0.2e1 * (t40 * t230 + (-t29 * t63 * t74 - t51 * t62 * t75) * t36 * t80) * t240;
t4 = t12 + (t143 * t338 + (0.2e1 * t157 * t291 + (t47 * t50 * t237 + (-t47 * t53 - t49 * t50 - t293) * t71) * t64) * t105) * pkin(3);
t2 = t12 + (t234 * t342 + (-0.2e1 * t334 * t65 * (-t284 * t54 + t285 * t58) + (t38 * t54 * t237 + (-t38 * t58 - t41 * t54 - t293) * t71) * t64) * t105) * pkin(3);
t1 = (0.4e1 * t146 * t234 + ((t238 * t55 + 0.4e1 * t286 * t31) * t40 - (t31 * t207 + (t236 * t55 + t333 * t59) * t62) * t36) * t105) * pkin(4) + t277;
t7 = [-(-mrSges(4,1) * t26 + mrSges(4,2) * t25) * t220 - (-mrSges(4,1) * t24 + mrSges(4,2) * t23) * t296 + (m(4) * t192 + t148) * t243 + (-(mrSges(5,1) * t134 + mrSges(5,2) * t133) * t272 + t20 * (Ifges(5,5) * t133 - Ifges(5,6) * t134) / 0.2e1) * t116 + t114 * t165 / 0.2e1 + t33 * t219 + (-Ifges(3,2) * t111 + t273) * t205 + (t164 * t26 * t81 + t152 * t25) * mrSges(4,3) - t34 * t140 / 0.2e1 + t37 * t141 / 0.2e1 + t5 * t139 / 0.2e1 + t52 * t17 * t209 + t51 * t16 * t210 + (t330 * (Ifges(5,1) * t133 - Ifges(5,4) * t134) * t203 + (Ifges(5,1) * t37 - Ifges(5,4) * t34 + Ifges(5,5) * t5) * t209 + t156 * t324) * t80 + (t330 * (Ifges(5,4) * t133 - Ifges(5,2) * t134) * t204 + (Ifges(5,4) * t37 - Ifges(5,2) * t34 + Ifges(5,6) * t5) * t210 + t16 * t156 / 0.2e1) * t79 - pkin(1) * (mrSges(5,1) * t34 + mrSges(5,2) * t37) + (t97 + qJD(1) * (t111 * Ifges(3,1) + t273)) * t245 / 0.2e1 - (t147 + t96) * t246 / 0.2e1 - (mrSges(4,2) * t337 - mrSges(4,3) * t35 + Ifges(4,1) * t23 + Ifges(4,4) * t24 + Ifges(4,5) * t6) * t60 + (-mrSges(4,1) * t337 + mrSges(4,3) * t39 + Ifges(4,4) * t23 + Ifges(4,2) * t24 + Ifges(4,6) * t6) * t118 * t142 + t25 * t9 / 0.2e1 + t26 * t8 / 0.2e1 + t19 * (Ifges(4,5) * t25 + Ifges(4,6) * t26) / 0.2e1 - t56 * (Ifges(4,1) * t25 + Ifges(4,4) * t26) / 0.2e1 + t57 * (Ifges(4,4) * t25 + Ifges(4,2) * t26) / 0.2e1; -t165 * t243 / 0.2e1 + m(4) * ((-t35 * t78 - t39 * t81) * t225 + ((t281 / 0.2e1 + t50 * t321) * t121 * t247 + (-t128 - t77) * pkin(2) * t104 * t189) / t340) / 0.2e1 - m(4) * (t327 * t192 + (-t42 * t78 - t44 * t81) * t187) / 0.2e1 - t327 * t148 / 0.2e1 + pkin(2) * (mrSges(4,1) * t6 - mrSges(4,3) * t23) * t180 + t81 * (-mrSges(4,2) * t6 + mrSges(4,3) * t24) * t184 + (-mrSges(4,1) * t28 + mrSges(4,2) * t27) * t220 + t96 * t202 - pkin(2) * t33 * t250 + Ifges(3,5) * t205 + ((t79 * (Ifges(5,4) * t46 + Ifges(5,2) * t43 + Ifges(5,6) * t1) + t3 * t307) * t203 + (t80 * (Ifges(5,1) * t46 + Ifges(5,4) * t43 + Ifges(5,5) * t1) + t3 * t304) * t204) * t116 + (-t2 / 0.2e1 + t4) * (-t56 * Ifges(4,5) + t57 * Ifges(4,6) + t19 * Ifges(4,3)) + (t153 * t78 + t184 * t50 - t42) * (mrSges(4,1) * t19 + mrSges(4,3) * t56) + (t153 * t81 + t184 * t53 - t44) * (-mrSges(4,2) * t19 + mrSges(4,3) * t57) + (-t1 / 0.2e1 + t326) * (qJD(1) * t139 + t303) - (-Ifges(3,2) * t250 + t109 + t97) * t249 / 0.2e1 + (Ifges(5,5) * t37 - Ifges(5,6) * t34 + Ifges(5,3) * t5) * t145 * t240 + (mrSges(4,3) * t28 + (-t2 + t4) * mrSges(4,2)) * t81 * t163 - t27 * t9 / 0.2e1 - t28 * t8 / 0.2e1 - t19 * (Ifges(4,5) * t27 + Ifges(4,6) * t28 + Ifges(4,3) * t2) / 0.2e1 - t43 * t16 / 0.2e1 - t20 * (Ifges(5,5) * t46 + Ifges(5,6) * t43 + Ifges(5,3) * t1) / 0.2e1 + t56 * (Ifges(4,1) * t27 + Ifges(4,4) * t28 + Ifges(4,5) * t2) / 0.2e1 + (-mrSges(5,1) * t43 + mrSges(5,2) * t46) * t272 + t46 * t324 + t303 * t326 - t57 * (Ifges(4,4) * t27 + Ifges(4,2) * t28 + Ifges(4,6) * t2) / 0.2e1 + (mrSges(4,1) * t35 - mrSges(4,2) * t39 + Ifges(4,5) * t23 + Ifges(4,6) * t24 + Ifges(4,3) * t6) * (-t143 * t295 + 0.1e1) + (mrSges(4,1) * t2 - mrSges(4,3) * t27) * t152 + t78 * t4 * mrSges(4,1) * t164 - Ifges(3,6) * t206;];
tauc = t7(:);
