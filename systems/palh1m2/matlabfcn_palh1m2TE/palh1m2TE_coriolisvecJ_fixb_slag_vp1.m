% Calculate vector of centrifugal and Coriolis load on the joints for
% palh1m2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = palh1m2TE_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2TE_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_coriolisvecJ_fixb_slag_vp1: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2TE_coriolisvecJ_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2TE_coriolisvecJ_fixb_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m2TE_coriolisvecJ_fixb_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:36:37
% EndTime: 2020-05-01 20:36:53
% DurationCPUTime: 7.45s
% Computational Cost: add. (1414->363), mult. (2346->518), div. (0->0), fcn. (406->88), ass. (0->247)
t354 = rSges(6,2) * m(6);
t141 = rSges(6,1) * t354 - Icges(6,4);
t208 = qJD(1) ^ 2;
t101 = t141 * t208;
t199 = rSges(6,1) * pkin(9);
t189 = sin(qJ(4));
t317 = t189 * t208;
t107 = m(6) * t317 * t199;
t347 = m(6) * t208;
t295 = rSges(6,2) * t347;
t151 = pkin(9) * t295;
t193 = cos(qJ(4));
t85 = (rSges(6,1) ^ 2 - rSges(6,2) ^ 2) * m(6) - Icges(6,1) + Icges(6,2);
t283 = t85 * t317;
t248 = -t107 + t101 + (-t151 - t283) * t193;
t278 = t193 ^ 2 * t101;
t46 = -0.2e1 * t278;
t6 = t248 + t46;
t183 = qJ(3) + qJ(2);
t157 = sin(t183);
t202 = m(5) + m(6);
t367 = pkin(5) * t202 + m(11) * rSges(11,1) + m(4) * rSges(4,1);
t384 = t157 * t367;
t196 = cos(pkin(18));
t181 = t196 ^ 2;
t257 = -t248 + 0.2e1 * t278;
t383 = t257 * t181;
t197 = -rSges(6,3) - pkin(11);
t104 = m(5) * rSges(5,2) + t197 * m(6);
t172 = pkin(22) + pkin(21);
t265 = pkin(18) - t172;
t231 = -pkin(20) + t265;
t95 = qJ(2) + t231;
t74 = qJ(3) + t95;
t54 = cos(t74);
t96 = -qJ(2) + t231;
t75 = -qJ(3) + t96;
t55 = cos(t75);
t382 = (t54 / 0.2e1 - t55 / 0.2e1) * t104;
t210 = 0.2e1 * qJ(2);
t184 = t210 + qJ(3);
t381 = t367 * cos(t184);
t194 = cos(qJ(3));
t380 = t367 * t194;
t177 = qJ(3) + pkin(19);
t152 = qJ(2) + t177;
t108 = 0.2e1 * t152;
t138 = 0.2e1 * t183;
t341 = m(11) * rSges(11,2);
t356 = rSges(4,2) * m(4);
t359 = m(9) * rSges(9,2);
t379 = (rSges(11,1) * t341 + rSges(4,1) * t356 - Icges(11,4) - Icges(4,4)) * cos(t138) + (rSges(9,1) * t359 - Icges(9,4)) * cos(t108);
t240 = qJ(2) + t265;
t102 = qJ(3) + t240;
t241 = -qJ(2) + t265;
t103 = -qJ(3) + t241;
t378 = (cos(t102) + cos(t103)) * rSges(11,2);
t332 = (pkin(5) ^ 2 * t202 - Icges(11,1) - Icges(4,1) + Icges(11,2) + Icges(4,2) + (rSges(11,1) ^ 2 - rSges(11,2) ^ 2) * m(11) + (rSges(4,1) ^ 2 - rSges(4,2) ^ 2) * m(4)) * sin(t138);
t377 = t332 / 0.2e1;
t168 = t197 * rSges(6,1);
t112 = m(6) * t168 + Icges(6,5);
t320 = t112 * t189;
t336 = rSges(6,2) * t197;
t321 = (m(6) * t336 + Icges(6,6)) * t193;
t376 = (t321 / 0.2e1 + t320 / 0.2e1) * t208;
t375 = (t320 + t321) * t208;
t115 = t341 + t356;
t190 = sin(qJ(3));
t316 = t190 * t115;
t374 = -t316 + t380;
t167 = pkin(17) + qJ(2) - pkin(18);
t114 = 0.2e1 * t167;
t353 = rSges(10,2) * m(10);
t357 = rSges(3,2) * m(3);
t373 = (rSges(3,1) * t357 + rSges(10,1) * t353 - Icges(3,4) - Icges(10,4)) * cos(t210) + (m(7) * rSges(7,1) * rSges(7,2) - Icges(7,4)) * cos(t114);
t319 = t115 * sin(t184);
t372 = t319 - t381;
t334 = rSges(10,2) * sin(t177);
t337 = rSges(10,1) * cos(t177);
t371 = t334 + t337;
t155 = t210 + t177;
t335 = rSges(10,2) * sin(t155);
t338 = rSges(10,1) * cos(t155);
t370 = t335 - t338;
t176 = qJD(2) + qJD(3);
t343 = qJD(3) / 0.2e1;
t369 = (2 * qJD(2)) + 0.2e1 * t343;
t368 = m(6) * (t189 * rSges(6,1) + t193 * rSges(6,2));
t144 = m(8) + m(11) + m(4) + t202;
t366 = (pkin(1) ^ 2 * t144 - Icges(3,1) - Icges(10,1) + Icges(3,2) + Icges(10,2) + (rSges(10,1) ^ 2 - rSges(10,2) ^ 2) * m(10) + (rSges(3,1) ^ 2 - rSges(3,2) ^ 2) * m(3)) * sin(t210);
t363 = pkin(2) * m(10);
t362 = pkin(4) * m(11);
t361 = m(5) * rSges(5,3);
t360 = m(7) * rSges(7,3);
t358 = m(9) * rSges(9,3);
t200 = rSges(3,1) * m(3);
t355 = rSges(6,2) * pkin(9);
t352 = rSges(4,3) * m(4);
t351 = rSges(10,3) * m(10);
t350 = pkin(14) * m(7);
t349 = pkin(5) * t190;
t348 = pkin(5) * t194;
t346 = (m(10) * pkin(2) ^ 2 + (rSges(9,1) ^ 2 - rSges(9,2) ^ 2) * m(9) - Icges(9,1) + Icges(9,2)) * sin(t108);
t345 = ((rSges(7,1) ^ 2 - rSges(7,2) ^ 2) * m(7) - Icges(7,1) + Icges(7,2)) * sin(t114);
t344 = t181 * t375;
t340 = m(11) * rSges(11,3);
t339 = m(6) * qJD(1);
t333 = pkin(15) * t208;
t171 = t176 ^ 2;
t195 = cos(qJ(2));
t207 = qJD(2) ^ 2;
t327 = (pkin(1) * t207 + t171 * t349) * t195;
t322 = pkin(15) * qJD(1);
t186 = sin(pkin(20));
t187 = cos(pkin(20));
t318 = t186 * t187;
t192 = sin(pkin(18));
t315 = t192 * t196;
t314 = qJ(4) - qJ(2);
t311 = rSges(10,1) * m(10) + t200;
t309 = qJD(4) * Icges(6,5);
t308 = qJD(4) * Icges(6,6);
t307 = qJD(4) * t197;
t306 = qJD(1) * qJD(4);
t304 = pkin(2) * t351;
t303 = pkin(18) - pkin(22);
t175 = qJD(2) - qJD(4);
t174 = qJD(2) + qJD(4);
t150 = qJD(3) + t175;
t302 = pkin(5) * m(6) * t150;
t301 = rSges(9,1) * t358;
t300 = rSges(9,2) * t358;
t299 = pkin(15) * t359;
t298 = rSges(4,2) * t352;
t296 = t208 * t362;
t294 = 0.2e1 * t322;
t191 = sin(qJ(2));
t292 = t191 * t348;
t291 = qJD(1) * t350;
t290 = t6 * t315;
t289 = t315 * t376;
t288 = (t290 - t344 + t376) * t318;
t284 = rSges(11,2) * t340;
t277 = -t339 / 0.2e1;
t276 = t339 / 0.2e1;
t166 = qJ(4) + t183;
t131 = cos(t166);
t274 = t131 * t306;
t273 = -t307 / 0.2e1;
t272 = t307 / 0.2e1;
t271 = pkin(1) * m(6) * t174 * t175;
t270 = -t306 / 0.2e1;
t268 = 0.2e1 * t6;
t261 = t171 * t292;
t260 = pkin(1) * t276;
t259 = -t296 / 0.2e1;
t258 = pkin(5) * t277;
t252 = 0.2e1 * t231;
t243 = t271 / 0.2e1;
t149 = qJD(3) + t174;
t242 = -t149 * t302 / 0.2e1;
t237 = t174 * t260;
t236 = t175 * t260;
t234 = pkin(5) * t149 * t276;
t233 = t150 * t258;
t165 = qJ(3) - t314;
t127 = sin(t165);
t128 = sin(t166);
t232 = (-t127 - t128) * t306;
t230 = t268 * t181 + t257 + 0.4e1 * t289;
t94 = -qJ(4) + t231;
t93 = qJ(4) + t231;
t217 = (t176 * t292 + t195 * (pkin(1) * qJD(2) + t176 * t349)) * t339;
t118 = sin(t152);
t122 = cos(t152);
t130 = cos(t165);
t133 = pkin(9) * m(6) + m(5) * rSges(5,1);
t134 = m(9) * rSges(9,1) + t363;
t161 = cos(t183);
t63 = qJ(4) + t74;
t31 = sin(t63);
t64 = -qJ(4) + t74;
t32 = sin(t64);
t35 = cos(t63);
t36 = cos(t64);
t48 = sin(t74);
t49 = sin(t75);
t80 = sin(t102);
t81 = sin(t103);
t214 = (t115 * t161 + t134 * t118) * t333 + t259 * t378 + (t81 * t296 / 0.2e1 + t80 * t259) * rSges(11,1) + (t122 * t299 + t346 / 0.2e1 + t379) * t208 + (-(rSges(6,2) * t35 + (t31 + t32) * rSges(6,1)) * t347 / 0.4e1 + (-t48 / 0.2e1 + t49 / 0.2e1) * t133 * t208 + t130 * t306 * t354 + t36 * t295 / 0.4e1) * pkin(5);
t209 = 0.2e1 * qJ(4);
t206 = qJD(4) ^ 2;
t185 = qJ(2) + qJ(4);
t173 = t187 ^ 2;
t163 = cos(t185);
t160 = cos(t314);
t159 = sin(t185);
t156 = -sin(t314);
t154 = -qJ(2) + t303;
t153 = qJ(2) + t303;
t147 = 0.2e1 * t206 + t208;
t143 = cos(t172);
t142 = sin(t172);
t132 = cos(t167);
t129 = sin(t167);
t126 = t353 + t357;
t124 = cos(t154);
t123 = cos(t153);
t120 = sin(t154);
t119 = sin(t153);
t116 = 0.2e1 / 0.3e1 * t208 + 0.4e1 / 0.3e1 * t206;
t113 = t171 + t208;
t106 = t322 - qJD(4) * pkin(9) / 0.2e1;
t100 = cos(t241);
t99 = cos(t240);
t98 = -qJ(4) + t252;
t97 = qJ(4) + t252;
t79 = -qJ(2) + t94;
t78 = -qJ(2) + t93;
t77 = qJ(2) + t94;
t76 = qJ(2) + t93;
t73 = rSges(6,2) * t106;
t72 = t106 * rSges(6,1);
t71 = cos(t96);
t70 = cos(t95);
t69 = sin(t96);
t68 = sin(t95);
t66 = -qJ(4) + t75;
t65 = qJ(4) + t75;
t62 = 0.2e1 * t94;
t61 = 0.2e1 * t93;
t59 = cos(t79);
t58 = cos(t78);
t57 = cos(t77);
t56 = cos(t76);
t53 = sin(t79);
t52 = sin(t78);
t51 = sin(t77);
t50 = sin(t76);
t38 = cos(t66);
t37 = cos(t65);
t34 = sin(t66);
t33 = sin(t65);
t4 = (-t113 * t195 * t348 + (t113 * t349 + pkin(1) * (t207 + t208)) * t191 - t333) * t368;
t2 = t192 * t4 + t196 * (t261 + t327) * t368;
t1 = t4 * t196 - 0.2e1 * t192 * (t261 / 0.2e1 + t327 / 0.2e1) * t368;
t3 = [(((t199 + t336) * m(6) + Icges(6,6)) * sin(t98) + t85 * sin(t209) + ((-t168 - t355) * m(6) - Icges(6,5)) * cos(t97)) * t306 / 0.2e1 - 0.2e1 * (((-rSges(10,2) * t351 / 0.2e1 - rSges(3,3) * t357 / 0.2e1 + Icges(3,6) / 0.2e1 + Icges(10,6) / 0.2e1) * qJD(2) + (t144 * pkin(1) + t311) * t322) * t195 + (rSges(7,2) * t291 + (-rSges(7,1) * t360 / 0.2e1 + Icges(7,5) / 0.2e1) * qJD(2)) * t129) * qJD(2) + (((t80 - t81) * rSges(11,1) + t378) * t362 - 0.2e1 * t379 - t332 - t346 + ((t48 - t49) * t133 + (t54 - t55) * t104) * pkin(5)) * qJD(1) * t176 + (t374 * qJD(3) + t372 * t369) * pkin(1) * qJD(1) + (t37 + t36) * rSges(6,2) * t233 + (t38 + t35) * rSges(6,2) * t234 + (t58 + t57) * rSges(6,1) * t236 + (t59 + t56) * rSges(6,1) * t237 + (t131 + t130) * rSges(6,2) * t242 + (t159 + t156) * rSges(6,2) * t243 + (t174 * t50 + t175 * t52) * pkin(1) * rSges(6,2) * t277 + (cos(t209) * t306 + (cos(t61) + cos(t62)) * t270) * t141 - 0.2e1 * (((-t284 / 0.2e1 - t298 / 0.2e1 + Icges(11,6) / 0.2e1 + Icges(4,6) / 0.2e1) * qJD(2) + (Icges(11,6) + Icges(4,6) - t284 - t298) * t343 + t367 * t322) * t157 + ((t301 / 0.2e1 + t304 / 0.2e1 - Icges(9,5) / 0.2e1) * qJD(2) + qJD(1) * t299 + (-Icges(9,5) + t301 + t304) * t343) * t122 + ((-t300 / 0.2e1 + Icges(9,6) / 0.2e1) * qJD(2) + t134 * t322 - qJD(3) * (-Icges(9,6) + t300) / 0.2e1) * t118) * t176 - (t115 * t294 + t176 * (pkin(5) * t361 + rSges(11,1) * t340 + rSges(4,1) * t352 - Icges(11,5) - Icges(4,5))) * t176 * t161 + (t149 * t34 * t258 + t33 * t233 + t31 * t234 + t128 * t242 + t160 * t243 - t163 * t271 / 0.2e1) * rSges(6,1) + (((rSges(6,2) * t272 + t72) * m(6) + t308 / 0.2e1) * sin(t93) - ((rSges(6,2) * t273 + t72) * m(6) - t308 / 0.2e1) * sin(t94) + ((rSges(6,1) * t272 + t73) * m(6) + t309 / 0.2e1) * cos(t94) + ((rSges(6,1) * t273 + t73) * m(6) - t309 / 0.2e1) * cos(t93) + (-pkin(9) * qJD(1) + qJD(4) * pkin(15)) * t368) * qJD(4) + (((-t168 + t355) * m(6) - Icges(6,5)) * cos(t98) + ((t199 - t336) * m(6) - Icges(6,6)) * sin(t97)) * t270 + (t371 * qJD(3) + t370 * t369) * qJD(1) * t363 + 0.2e1 * qJD(2) * (rSges(7,1) * t291 + (rSges(7,2) * t360 / 0.2e1 - Icges(7,6) / 0.2e1) * qJD(2)) * t132 + qJD(2) * (t126 * t294 + ((rSges(8,3) * m(8) + t340 + t352 + t361) * pkin(1) + rSges(10,1) * t351 + rSges(3,3) * t200 - Icges(3,5) - Icges(10,5)) * qJD(2)) * t191 + (qJD(1) * t32 + t149 * t127) * rSges(6,1) * t302 / 0.2e1 + (t51 * t236 + t53 * t237) * rSges(6,2) + (sin(t62) / 0.4e1 - sin(t61) / 0.4e1) * t85 * t306 + (((t100 + t99) * t362 + (t70 + t71) * t133 + (-t68 - t69) * t104 + ((-t119 - t120) * rSges(8,2) + (t123 + t124) * rSges(8,1)) * m(8)) * pkin(1) + t345 + t366 + 0.2e1 * t373) * qJD(2) * qJD(1); (((-t274 + (0.3e1 / 0.8e1 * t37 - 0.3e1 / 0.8e1 * t38) * t116) * rSges(6,2) + ((0.3e1 / 0.8e1 * t33 + 0.3e1 / 0.8e1 * t34) * t116 + t232) * rSges(6,1)) * pkin(5) + (((0.3e1 / 0.8e1 * t52 - 0.3e1 / 0.8e1 * t53) * t116 + (-t156 + t159) * t306) * rSges(6,2) + ((-0.3e1 / 0.8e1 * t58 - 0.3e1 / 0.8e1 * t59) * t116 + (-t160 - t163) * t306) * rSges(6,1)) * pkin(1)) * m(6) + (pkin(1) * t374 + t371 * t363) * qJD(3) * t369 + t214 + (-t345 / 0.2e1 + t377 - t366 / 0.2e1 + (-rSges(7,1) * t132 + rSges(7,2) * t129) * t350 - pkin(5) * t382 - t370 * t363 + (-t126 * t191 + t311 * t195 + t384) * pkin(15) + (pkin(15) * t144 * t195 + (-t70 / 0.2e1 - t71 / 0.2e1) * t133 - (-t68 / 0.2e1 - t69 / 0.2e1) * t104 + (-t99 / 0.2e1 - t100 / 0.2e1) * t362 + ((t119 / 0.2e1 + t120 / 0.2e1) * rSges(8,2) + (-t123 / 0.2e1 - t124 / 0.2e1) * rSges(8,1)) * m(8) + ((t50 / 0.4e1 - t51 / 0.4e1) * rSges(6,2) + (-t56 / 0.4e1 - t57 / 0.4e1) * rSges(6,1)) * m(6) - t372) * pkin(1) - t373) * t208; (-t208 * t382 + ((-t274 + (t37 / 0.4e1 - t38 / 0.4e1) * t147) * rSges(6,2) + ((t33 / 0.4e1 + t34 / 0.4e1) * t147 + t232) * rSges(6,1)) * m(6)) * pkin(5) + ((-t334 / 0.2e1 - t337 / 0.2e1) * t363 + (t316 / 0.2e1 - t380 / 0.2e1) * pkin(1)) * ((2 * t207) + t208) + t214 + (t377 + pkin(15) * t384 + (-t335 / 0.2e1 + t338 / 0.2e1) * t363 + (-t319 / 0.2e1 + t381 / 0.2e1) * pkin(1)) * t208; (-t186 * t1 + t2 * t187) * t142 + t230 * t173 + 0.2e1 * t288 + t383 - 0.2e1 * t289 + t46 + (0.2e1 * rSges(6,1) * t217 - t283) * t193 - 0.2e1 * rSges(6,2) * t189 * t217 + t101 + (((-0.4e1 * t290 + 0.4e1 * t344 - 0.2e1 * t375) * t173 + 0.4e1 * (t6 * t181 + 0.2e1 * t289 + t278 + (t283 / 0.2e1 + t151 / 0.2e1) * t193 + t107 / 0.2e1 - t101 / 0.2e1) * t318 - 0.2e1 * t344 + 0.2e1 * t290 + t208 * t321 + t112 * t317) * t142 + t1 * t187 + t186 * t2 + ((t268 - 0.8e1 * t289 + 0.4e1 * t383) * t173 - 0.4e1 * t288 + t230) * t143) * t143;];
tauc = t3(:);
