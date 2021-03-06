% Calculate inertial parameters regressor of coriolis joint torque vector for
% fourbar1turnDE1
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
% 
% Output:
% tauc_reg [2x(2*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = fourbar1turnDE1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:27:15
% EndTime: 2020-04-12 19:27:37
% DurationCPUTime: 8.94s
% Computational Cost: add. (110229->273), mult. (158035->688), div. (5778->22), fcn. (42461->8), ass. (0->259)
t357 = pkin(1) * pkin(2);
t136 = pkin(1) ^ 2;
t124 = cos(qJ(2));
t319 = pkin(2) * t124;
t266 = -0.2e1 * pkin(1) * t319 + t136;
t334 = -pkin(3) - pkin(4);
t110 = (pkin(2) - t334) * (pkin(2) + t334) + t266;
t333 = pkin(4) - pkin(3);
t111 = (pkin(2) - t333) * (pkin(2) + t333) + t266;
t281 = t110 * t111;
t139 = sqrt(-t281);
t104 = 0.1e1 / t139;
t123 = sin(qJ(2));
t207 = (-t110 - t111) * t357;
t100 = t123 * t207;
t99 = qJD(2) * t100;
t356 = 0.4e1 * t104 / t281 * t99 ^ 2;
t135 = pkin(2) ^ 2;
t117 = t135 + t266;
t348 = pkin(3) ^ 2;
t349 = pkin(4) ^ 2;
t264 = t348 - t349;
t113 = t117 - t264;
t118 = pkin(1) - t319;
t275 = t123 * t139;
t93 = -pkin(2) * t275 + t113 * t118;
t88 = 0.1e1 / t93 ^ 2;
t320 = pkin(2) * t123;
t108 = t113 * t320;
t94 = t118 * t139 + t108;
t298 = t88 * t94;
t129 = 0.1e1 / pkin(4);
t115 = 0.1e1 / t117 ^ 2;
t260 = qJD(2) * t123;
t237 = pkin(2) * t260;
t210 = t115 * t237;
t191 = pkin(1) * t210;
t114 = 0.1e1 / t117;
t328 = -t114 / 0.2e1;
t289 = t104 * t99;
t236 = t123 * t289;
t254 = 0.2e1 * t118 * pkin(1);
t274 = t124 * t139;
t59 = (-t236 + (-t274 + (t113 + t254) * t123) * qJD(2)) * pkin(2);
t41 = (t93 * t191 + t59 * t328) * t129;
t327 = t114 / 0.2e1;
t224 = t139 * t260;
t259 = qJD(2) * t124;
t284 = t104 * t118;
t121 = t123 ^ 2;
t322 = pkin(1) * t135;
t354 = 0.2e1 * t121 * t322;
t60 = qJD(2) * t354 + t99 * t284 + (t113 * t259 + t224) * pkin(2);
t43 = (-t94 * t191 + t60 * t327) * t129;
t90 = t94 ^ 2;
t80 = t88 * t90 + 0.1e1;
t76 = 0.1e1 / t80;
t87 = 0.1e1 / t93;
t159 = t76 * (t41 * t298 + t43 * t87);
t126 = qJD(2) ^ 2;
t355 = t126 * t123;
t249 = 0.2e1 * t115;
t112 = t117 + t264;
t119 = pkin(1) * t124 - pkin(2);
t92 = -pkin(1) * t275 - t112 * t119;
t82 = t92 ^ 2;
t83 = 0.1e1 / t92;
t85 = t83 / t82;
t326 = pkin(1) * t112;
t109 = t123 * t326;
t95 = -t119 * t139 + t109;
t91 = t95 ^ 2;
t299 = t85 * t91;
t84 = 0.1e1 / t92 ^ 2;
t300 = t84 * t95;
t256 = -0.2e1 * pkin(2) * t119;
t220 = t112 + t256;
t58 = (-t236 + (t123 * t220 - t274) * qJD(2)) * pkin(1);
t283 = t104 * t119;
t277 = t121 * t136;
t345 = 0.2e1 * pkin(2) * t277;
t61 = pkin(1) * t224 + qJD(2) * t345 + t259 * t326 - t99 * t283;
t81 = t84 * t91 + 0.1e1;
t79 = 0.1e1 / t81 ^ 2;
t314 = (-t58 * t299 + t61 * t300) * t79;
t353 = -0.2e1 * t314;
t273 = t126 * t124;
t116 = t114 * t115;
t227 = t135 * t277;
t352 = t126 * t116 * t227;
t351 = t322 * t355;
t346 = 0.6e1 * t124;
t317 = pkin(4) * t117;
t255 = 0.2e1 * t317;
t130 = 0.1e1 / t349;
t86 = t93 ^ 2;
t292 = t86 + t90;
t74 = t292 * t130 * t115;
t66 = t74 ^ (-0.1e1 / 0.2e1);
t344 = t115 * t66;
t18 = t159 * t255;
t250 = pkin(1) * t320;
t219 = t115 * t250;
t285 = t104 * t100;
t63 = t108 + (-t274 + (t254 - t285) * t123) * pkin(2);
t45 = (t93 * t219 + t63 * t328) * t129;
t65 = t100 * t284 + t354 + (t124 * t113 + t275) * pkin(2);
t47 = (-t94 * t219 + t65 * t327) * t129;
t158 = t76 * (-t45 * t298 - t47 * t87);
t19 = t158 * t255;
t206 = 0.2e1 * t250;
t343 = (qJD(2) * t19 + t18) * t206 * t344;
t132 = 0.1e1 / pkin(3);
t179 = -0.2e1 * t191;
t42 = (t114 * t58 + t179 * t92) * t132;
t44 = (t114 * t61 + t179 * t95) * t132;
t342 = -t42 * t300 + t44 * t83;
t70 = 0.1e1 / t74;
t133 = 0.1e1 / t348;
t293 = t82 + t91;
t75 = t293 * t133 * t115;
t72 = 0.1e1 / t75;
t341 = 0.2e1 * t59;
t340 = -0.2e1 * t88;
t339 = 0.4e1 * t94;
t338 = -0.2e1 * t95;
t337 = 0.2e1 * t95;
t68 = t75 ^ (-0.1e1 / 0.2e1);
t205 = 0.4e1 * t250;
t156 = t292 * t116 * t205;
t29 = ((t63 * t93 + t65 * t94) * t249 - t156) * t130;
t335 = -t29 / 0.2e1;
t77 = 0.1e1 / t80 ^ 2;
t222 = 0.4e1 * t77 * t298;
t301 = t77 * t87;
t89 = t87 / t86;
t297 = t89 * t90;
t31 = -t59 * t297 + t60 * t298;
t246 = t31 * t301;
t251 = t89 * t339;
t253 = 0.2e1 * t76 * t88;
t33 = -t63 * t297 + t65 * t298;
t332 = ((-qJD(2) * t158 - t159) * t205 + ((t253 * t63 + 0.4e1 * t33 * t301) * t43 + (t33 * t222 + (t251 * t63 + t65 * t340) * t76) * t41 - (t253 * t59 + 0.4e1 * t246) * t47 - (t31 * t222 + (t251 * t59 + t60 * t340) * t76) * t45) * t117) * pkin(4);
t202 = -0.2e1 * t219;
t62 = t109 + (-t274 + (t256 - t285) * t123) * pkin(1);
t46 = (t114 * t62 + t202 * t92) * t132;
t64 = -t100 * t283 + t345 + (t112 * t124 + t275) * pkin(1);
t48 = (t114 * t64 + t202 * t95) * t132;
t175 = t46 * t300 - t48 * t83;
t295 = t95 * t64;
t78 = 0.1e1 / t81;
t331 = ((t175 * qJD(2) + t342) * t78 * t206 + (-0.2e1 * t342 * t79 * (t84 * t295 - t62 * t299) + t175 * t353 + ((t62 * t42 - t58 * t46) * t85 * t337 + (-t64 * t42 - t44 * t62 + t61 * t46 + t48 * t58) * t84) * t78) * t117) * pkin(3);
t203 = -t356 / 0.4e1;
t269 = t139 * t126;
t96 = (t124 * t207 - 0.4e1 * t227) * t126;
t153 = -t104 * t96 + t203 + t269;
t229 = qJD(2) * t289;
t212 = 0.2e1 * t229;
t163 = -t113 * t126 + t212;
t213 = pkin(1) * t237;
t247 = t76 * t317;
t215 = t88 * t247;
t225 = t124 * t269;
t272 = t126 * t135;
t325 = pkin(1) * t114;
t5 = 0.2e1 * (-t41 * t60 + t43 * t59) * t215 - 0.4e1 * pkin(4) * t213 * t159 + 0.2e1 * (t43 * t246 + (t88 * t77 * t31 + t89 * t76 * t59) * t41 * t94) * t255 + 0.2e1 * (-((t118 * t356 / 0.4e1 + t346 * t351 + t96 * t284) * t327 + t339 * t352 + ((t123 * t163 + t225) * t327 + (-0.2e1 * t260 * t60 - t94 * t273) * t115 * pkin(1)) * pkin(2)) * t87 * t247 - ((-0.4e1 * t116 * t136 * t93 - 0.2e1 * t325) * t121 * t272 + ((t123 * t153 - t124 * t163) * t328 + (-t114 * t118 * t273 + (t260 * t341 + t93 * t273) * t115) * pkin(1)) * pkin(2)) * t94 * t215) * t129;
t330 = t5 * t66;
t318 = pkin(3) * t117;
t216 = t78 * t84 * t318;
t196 = t95 * t216;
t248 = t83 * t318;
t217 = t78 * t248;
t240 = t136 * t320;
t290 = pkin(2) * qJD(2);
t321 = pkin(2) * t115;
t6 = t44 * t248 * t353 + (-t42 * t61 - t44 * t58) * t216 + (t85 * t78 * t58 + t84 * t314) * t42 * t318 * t337 + 0.2e1 * t342 * pkin(3) * t78 * t213 + (((t126 * t240 * t346 + t119 * t203 - t96 * t283) * t114 + 0.8e1 * t95 * t352 + ((t225 + (-t112 * t126 + t212) * t123) * t114 + (-0.4e1 * t260 * t61 + t273 * t338) * t321) * pkin(1)) * t217 - ((0.8e1 * t116 * t135 * t92 + 0.4e1 * pkin(2) * t114) * t126 * t277 + ((-0.2e1 * t114 * t229 + (t114 * t220 - 0.2e1 * t321 * t92) * t126) * t124 + (-0.4e1 * t58 * t115 * t290 + t114 * t153) * t123) * pkin(1)) * t196) * t132;
t329 = t68 * t6;
t27 = ((t59 * t93 + t60 * t94) * t249 - qJD(2) * t156) * t130;
t71 = 0.1e1 / t74 ^ 2;
t315 = t27 * t71;
t183 = t123 * t92 + t124 * t95;
t279 = t114 * t132;
t234 = t68 * t279;
t39 = t183 * t234;
t313 = t39 * t92;
t287 = t124 * t92;
t288 = t123 * t95;
t40 = (-t287 + t288) * t234;
t312 = t40 * t95;
t309 = t59 * t66;
t308 = t60 * t66;
t67 = t66 * t70;
t307 = t67 * t93;
t306 = t67 * t94;
t69 = t68 * t72;
t305 = t69 * t92;
t304 = t69 * t95;
t303 = t70 * t93;
t302 = t70 * t94;
t296 = t93 * t94;
t291 = t92 + t64;
t276 = t123 * t124;
t127 = qJD(1) ^ 2;
t271 = t127 * t129;
t270 = t127 * t130;
t265 = -t124 ^ 2 + t121;
t263 = qJD(1) * t129;
t262 = qJD(1) * t130;
t261 = qJD(1) * t132;
t258 = qJD(1) * qJD(2);
t243 = t71 * t296;
t242 = t27 * t307;
t233 = t29 * t71 / 0.2e1;
t232 = -t306 / 0.2e1;
t231 = -t305 / 0.2e1;
t230 = -t304 / 0.2e1;
t226 = t127 * t276;
t218 = t116 * t250;
t173 = t293 * t218;
t157 = -0.4e1 * t173;
t160 = 0.2e1 * t58 * t92 + 0.2e1 * t61 * t95;
t28 = (qJD(2) * t157 + t115 * t160) * t133;
t211 = t28 * t230;
t198 = t258 * t276;
t197 = t70 * t218;
t190 = t92 * t62 + t295;
t189 = t66 * t210;
t188 = t234 * t288;
t187 = t68 * t249 * t357;
t185 = -0.2e1 * t198;
t184 = 0.2e1 * t197;
t178 = -0.2e1 * t240 * t344;
t174 = t115 * t68 * t351;
t170 = t19 * t27 / 0.2e1 - t18 * t335;
t17 = -t196 * t42 + t217 * t44 + qJD(2);
t20 = -t196 * t46 + t217 * t48 + 0.1e1;
t30 = (t190 * t249 + t157) * t133;
t169 = t20 * t28 / 0.2e1 - t17 * t30 / 0.2e1;
t56 = qJD(1) * t188;
t166 = t197 * t296;
t165 = -0.4e1 * qJD(2) * t197;
t164 = t231 * t28 + t58 * t68;
t162 = pkin(1) * t18 * t189;
t161 = 0.4e1 * t136 * t189;
t155 = t69 * (-t288 / 0.2e1 + t287 / 0.2e1);
t154 = t132 * t20 * t174;
t151 = (-t121 * t92 - t95 * t276) * t187;
t150 = (-t121 * t95 + t92 * t276) * t187;
t148 = (qJD(2) * t151 + (t164 * t123 + (t211 + (qJD(2) * t92 + t61) * t68) * t124) * t114) * t132;
t10 = (qJD(2) * t150 + (t28 * t155 + (qJD(2) * t183 + t123 * t61 - t124 * t58) * t68) * t114) * t132;
t73 = 0.1e1 / t75 ^ 2;
t38 = -qJD(1) * t234 * t287 + t56;
t37 = qJD(1) * t39;
t12 = (t150 + (t30 * t155 + ((-t62 + t95) * t124 + t291 * t123) * t68) * t114) * t261;
t11 = -t56 + (t151 + ((t231 * t30 + t62 * t68) * t123 + (t30 * t230 + t291 * t68) * t124) * t114) * t261;
t9 = -qJD(2) * t188 + t148;
t8 = qJD(1) * t10;
t7 = qJD(1) * t148 - qJD(2) * t56;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t198, -0.2e1 * t265 * t258, t273, t185, -t355, 0, 0, 0, 0, 0, t37 * t9 + t39 * t7, -t10 * t37 - t38 * t9 - t39 * t8 - t40 * t7, -t17 * t9 - t39 * t6, t10 * t38 + t40 * t8, t10 * t17 + t40 * t6, 0, (-t38 * t260 + t124 * t8 + (t10 * t124 - t260 * t40) * qJD(1)) * pkin(2), (-t37 * t260 + t124 * t7 + (t124 * t9 - t260 * t39) * qJD(1)) * pkin(2), (0.2e1 * (t312 + t313) * t174 + ((t312 / 0.2e1 + t313 / 0.2e1) * t69 * t28 + (-t10 * t95 - t39 * t58 - t40 * t61 - t9 * t92) * t68) * t114 * t290) * t132, t135 * t185, (t90 * t165 + (0.2e1 * t302 * t60 - t90 * t315) * t115) * t262, (0.8e1 * qJD(2) * t166 + (t27 * t243 + (-t59 * t94 - t60 * t93) * t70) * t249) * t262, (0.2e1 * t94 * t162 + (t94 * t330 - (t232 * t27 + t308) * t18) * t114) * t129, (t86 * t165 + (t303 * t341 - t86 * t315) * t115) * t262, (-0.2e1 * t93 * t162 + (-t93 * t330 - (-t309 + t242 / 0.2e1) * t18) * t114) * t129, 0, (t93 * t161 + (t242 - 0.2e1 * t309) * t325) * t263, (t94 * t161 + (t27 * t306 - 0.2e1 * t308) * t325) * t263, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t226, t265 * t127, 0, t226, 0, 0, 0, 0, 0, 0, -t37 * t11, t11 * t38 + t12 * t37, t11 * t17 - t20 * t7 + t331 * t37, -t38 * t12, -t12 * t17 + t20 * t8 - t331 * t38, -t331 * t17 + t20 * t6, 0.2e1 * t92 * t154 + ((-t124 * t12 + t123 * t38) * qJD(1) + (-t92 * t329 - t164 * t17 + (t169 * t305 + (t62 * t17 - t58 * t20 + t331 * t92) * t68) * qJD(2)) * t279) * pkin(2), t154 * t338 + ((-t11 * t124 + t123 * t37) * qJD(1) + (t95 * t329 + (t61 * t68 + t211) * t17 + (-t169 * t304 + (-t64 * t17 + t61 * t20 - t331 * t95) * t68) * qJD(2)) * t279) * pkin(2), ((-t58 * t37 - t61 * t38 - t92 * t7 - t95 * t8 + (t11 * t92 + t12 * t95 + t37 * t62 + t38 * t64) * qJD(2)) * t68 + (-qJD(2) * t30 + t28) * t69 * (t95 * t38 / 0.2e1 + t92 * t37 / 0.2e1)) * pkin(2) * t279, t135 * t226 + (-0.2e1 * t72 * t173 * t272 + ((-t190 * t72 + (t91 / 0.2e1 + t82 / 0.2e1) * t73 * t30) * t126 + (-t28 * t293 * t73 + t160 * t72) * qJD(2)) * t135 * t115) * t133, (t90 * t184 + (t233 * t90 - t302 * t65) * t115) * t270, (-0.4e1 * t166 + (-t29 * t243 + (t63 * t94 + t65 * t93) * t70) * t115) * t270, (-t94 * t343 + (-t170 * t306 + (t65 * t18 + t60 * t19 - t332 * t94) * t66) * t114) * t263, (t86 * t184 + (t233 * t86 - t303 * t63) * t115) * t270, (t93 * t343 + (t170 * t307 + (-t63 * t18 - t59 * t19 + t332 * t93) * t66) * t114) * t263, t332 * t18 + t19 * t5, (t93 * t178 + (t307 * t335 + t63 * t66) * t325) * t271, (t94 * t178 + (t232 * t29 + t65 * t66) * t325) * t271, 0, 0;];
tauc_reg = t1;
