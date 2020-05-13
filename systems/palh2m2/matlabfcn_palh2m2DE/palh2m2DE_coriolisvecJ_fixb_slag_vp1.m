% Calculate vector of centrifugal and Coriolis load on the joints for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = palh2m2DE_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_coriolisvecJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2DE_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'palh2m2DE_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'palh2m2DE_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:28
% EndTime: 2020-05-03 01:06:41
% DurationCPUTime: 10.90s
% Computational Cost: add. (3016->485), mult. (7556->724), div. (0->0), fcn. (5581->8), ass. (0->302)
t209 = qJD(1) + qJD(4);
t213 = sin(qJ(1));
t217 = cos(qJ(1));
t170 = rSges(7,1) * t217 - rSges(7,2) * t213;
t172 = rSges(7,1) * t213 + rSges(7,2) * t217;
t210 = sin(qJ(4));
t214 = cos(qJ(4));
t287 = -t170 * t214 + t172 * t210;
t434 = t209 * t287;
t216 = cos(qJ(2));
t219 = qJD(2) ^ 2;
t353 = t216 * t219;
t317 = pkin(4) * t353;
t432 = 2 * qJD(2);
t327 = qJD(2) * t216;
t192 = pkin(4) * t327;
t215 = cos(qJ(3));
t321 = qJD(3) * t215;
t148 = pkin(5) * t321 + t192;
t211 = sin(qJ(3));
t212 = sin(qJ(2));
t412 = pkin(4) * t212;
t318 = t219 * t412;
t409 = pkin(5) * qJD(3) ^ 2;
t149 = -t211 * t409 - t318;
t316 = qJD(2) * t412;
t411 = pkin(4) * t216;
t431 = t148 * t316 + t149 * t411;
t323 = qJD(3) * t211;
t315 = pkin(5) * t323;
t147 = t315 + t316;
t116 = t147 * t217;
t293 = pkin(1) + t411;
t178 = pkin(2) + t293;
t410 = pkin(5) * t215;
t245 = t178 + t410;
t131 = rSges(6,1) + t245;
t90 = t217 * rSges(6,3) - t131 * t213;
t430 = qJD(1) * t90 - t116;
t176 = t213 * t315;
t177 = t213 * t316;
t339 = t177 + t176;
t89 = rSges(6,3) * t213 + t131 * t217;
t376 = qJD(1) * t89;
t69 = -t339 + t376;
t164 = -t210 * rSges(7,1) - t214 * rSges(7,2);
t173 = rSges(7,1) * t214 - rSges(7,2) * t210;
t264 = t164 * t217 - t173 * t213;
t179 = rSges(4,1) + t293;
t114 = t217 * rSges(4,3) - t179 * t213;
t326 = qJD(2) * t217;
t310 = t212 * t326;
t429 = -pkin(4) * t310 + qJD(1) * t114;
t428 = 0.2e1 * qJD(3);
t425 = t148 * t315 + t149 * t410;
t199 = Icges(5,4) * t215;
t267 = -Icges(5,2) * t211 + t199;
t158 = Icges(5,1) * t211 + t199;
t200 = Icges(3,4) * t216;
t268 = -Icges(3,2) * t212 + t200;
t160 = Icges(3,1) * t212 + t200;
t166 = rSges(3,1) * t212 + rSges(3,2) * t216;
t356 = t213 * t216;
t314 = rSges(3,1) * t356;
t359 = t212 * t213;
t337 = rSges(3,2) * t359 + t217 * rSges(3,3);
t55 = -t166 * t326 + (-pkin(1) * t213 - t314 + t337) * qJD(1);
t203 = t213 * rSges(3,3);
t354 = t216 * t217;
t399 = rSges(3,2) * t212;
t244 = rSges(3,1) * t354 + t203 + (pkin(1) - t399) * t217;
t328 = qJD(2) * t213;
t56 = qJD(1) * t244 - t166 * t328;
t423 = t213 * t55 - t217 * t56;
t153 = Icges(3,5) * t216 - Icges(3,6) * t212;
t152 = Icges(3,5) * t212 + Icges(3,6) * t216;
t249 = qJD(2) * t152;
t386 = Icges(3,4) * t212;
t161 = Icges(3,1) * t216 - t386;
t103 = Icges(3,5) * t213 + t161 * t217;
t99 = Icges(3,6) * t213 + t217 * t268;
t394 = t212 * t99;
t271 = t103 * t216 - t394;
t378 = Icges(3,3) * t217;
t422 = t217 * t249 + (t153 * t213 + t271 - t378) * qJD(1);
t185 = Icges(3,4) * t359;
t384 = Icges(3,5) * t217;
t102 = Icges(3,1) * t356 - t185 - t384;
t380 = Icges(3,6) * t217;
t98 = Icges(3,4) * t356 - Icges(3,2) * t359 - t380;
t395 = t212 * t98;
t272 = t102 * t216 - t395;
t95 = Icges(3,3) * t213 + t153 * t217;
t373 = qJD(1) * t95;
t421 = qJD(1) * t272 + t213 * t249 - t373;
t151 = Icges(5,5) * t215 - Icges(5,6) * t211;
t150 = Icges(5,5) * t211 + Icges(5,6) * t215;
t246 = qJD(3) * t150;
t385 = Icges(5,4) * t211;
t159 = Icges(5,1) * t215 - t385;
t101 = Icges(5,5) * t213 + t159 * t217;
t97 = Icges(5,6) * t213 + t217 * t267;
t396 = t211 * t97;
t273 = t101 * t215 - t396;
t377 = Icges(5,3) * t217;
t420 = t217 * t246 + (t151 * t213 + t273 - t377) * qJD(1);
t361 = t211 * t213;
t184 = Icges(5,4) * t361;
t357 = t213 * t215;
t383 = Icges(5,5) * t217;
t100 = Icges(5,1) * t357 - t184 - t383;
t379 = Icges(5,6) * t217;
t96 = Icges(5,4) * t357 - Icges(5,2) * t361 - t379;
t397 = t211 * t96;
t274 = t100 * t215 - t397;
t93 = Icges(5,3) * t213 + t151 * t217;
t374 = qJD(1) * t93;
t419 = qJD(1) * t274 + t213 * t246 - t374;
t94 = Icges(3,5) * t356 - Icges(3,6) * t359 - t378;
t31 = t213 * t272 - t217 * t94;
t92 = Icges(5,5) * t357 - Icges(5,6) * t361 - t377;
t29 = t213 * t274 - t217 * t92;
t154 = Icges(5,2) * t215 + t385;
t266 = t154 * t211 - t158 * t215;
t418 = t266 * qJD(1) + t151 * qJD(3);
t156 = Icges(3,2) * t216 + t386;
t265 = t156 * t212 - t160 * t216;
t417 = t265 * qJD(1) + t153 * qJD(2);
t416 = t213 * (-t156 * t217 + t103) - t217 * (-Icges(3,2) * t356 + t102 - t185);
t415 = t213 * (-t154 * t217 + t101) - t217 * (-Icges(5,2) * t357 + t100 - t184);
t414 = t213 / 0.2e1;
t413 = -t217 / 0.2e1;
t407 = -qJD(1) / 0.2e1;
t406 = qJD(1) / 0.2e1;
t355 = t215 * t217;
t405 = -t100 * t355 - t213 * t92;
t404 = t101 * t355 + t213 * t93;
t403 = -t102 * t354 - t213 * t94;
t402 = t103 * t354 + t213 * t95;
t401 = rSges(3,1) * t216;
t400 = rSges(5,1) * t215;
t398 = rSges(5,2) * t211;
t202 = t213 * rSges(5,3);
t392 = t217 * t55;
t130 = pkin(3) + t245;
t372 = t130 * t213;
t371 = t130 * t217;
t369 = t147 * t213;
t368 = t150 * t213;
t367 = t150 * t217;
t366 = t152 * t213;
t365 = t152 * t217;
t360 = t211 * t217;
t358 = t212 * t217;
t51 = -t213 * t266 - t367;
t352 = t51 * qJD(1);
t52 = -t213 * t265 - t365;
t351 = t52 * qJD(1);
t346 = -t154 + t159;
t345 = t158 + t267;
t344 = -t156 + t161;
t343 = t160 + t268;
t163 = qJD(1) * t177;
t342 = qJD(1) * t176 + t163;
t331 = qJD(1) * t213;
t312 = t211 * t331;
t330 = qJD(1) * t217;
t341 = rSges(5,2) * t312 + rSges(5,3) * t330;
t311 = t212 * t331;
t340 = rSges(3,2) * t311 + rSges(3,3) * t330;
t338 = rSges(5,2) * t361 + t217 * rSges(5,3);
t336 = rSges(5,1) * t355 + t202;
t220 = qJD(1) ^ 2;
t335 = pkin(1) * t220;
t333 = qJD(1) * t151;
t332 = qJD(1) * t153;
t329 = qJD(2) * t148;
t325 = qJD(3) * t148;
t165 = rSges(5,1) * t211 + rSges(5,2) * t215;
t324 = qJD(3) * t165;
t168 = -t398 + t400;
t142 = qJD(3) * t168;
t322 = qJD(3) * t213;
t320 = qJD(3) * t217;
t313 = rSges(5,1) * t357;
t309 = t331 / 0.2e1;
t308 = t330 / 0.2e1;
t307 = -t328 / 0.2e1;
t304 = t326 / 0.2e1;
t303 = -t322 / 0.2e1;
t300 = t320 / 0.2e1;
t299 = -t92 + t396;
t298 = -t94 + t394;
t79 = t101 * t357;
t297 = t217 * t93 - t79;
t80 = t103 * t356;
t296 = t217 * t95 - t80;
t292 = -t178 - t400;
t291 = t178 - t398;
t290 = -pkin(1) - t401;
t289 = t264 * qJD(4);
t78 = t164 * t213 + t173 * t217;
t113 = rSges(4,3) * t213 + t179 * t217;
t75 = qJD(1) * t113 - t177;
t169 = -t399 + t401;
t77 = t170 * t210 + t172 * t214;
t37 = -t130 * t331 - t209 * t77 - t116;
t263 = -t130 * t330 + t339;
t38 = t209 * t78 - t263;
t277 = -t213 * t38 - t217 * t37;
t41 = (-t316 - t324) * t217 + (-t178 * t213 - t313 + t338) * qJD(1);
t42 = -t165 * t322 - t177 + (t291 * t217 + t336) * qJD(1);
t276 = -t213 * t42 - t217 * t41;
t275 = -t213 * t69 - t217 * t430;
t47 = t100 * t211 + t215 * t96;
t48 = t101 * t211 + t215 * t97;
t49 = t102 * t212 + t216 * t98;
t50 = t103 * t212 + t216 * t99;
t262 = 0.2e1 * t166;
t261 = 0.2e1 * t165;
t260 = qJD(1) * t166;
t259 = qJD(1) * t165;
t30 = -t361 * t97 - t297;
t257 = (t213 * t30 - t217 * t29) * qJD(3);
t32 = -t359 * t99 - t296;
t256 = (t213 * t32 - t217 * t31) * qJD(2);
t33 = -t360 * t96 - t405;
t34 = -t360 * t97 + t404;
t255 = (t213 * t34 - t217 * t33) * qJD(3);
t35 = -t358 * t98 - t403;
t36 = -t358 * t99 + t402;
t254 = (t213 * t36 - t217 * t35) * qJD(2);
t253 = -t215 * t409 - t317;
t251 = qJD(2) * t160;
t250 = qJD(2) * t156;
t248 = qJD(3) * t158;
t247 = qJD(3) * t154;
t243 = -t213 * t97 + t217 * t96;
t242 = -t213 * t99 + t217 * t98;
t240 = -qJD(3) * t142 - t178 * t220 - t317;
t239 = -t130 * t220 + t253;
t238 = (-t211 * t345 + t215 * t346) * qJD(1);
t237 = (-t212 * t343 + t216 * t344) * qJD(1);
t23 = t240 * t213 + (-qJD(1) * t313 + (-qJD(3) * t261 - 0.2e1 * t316) * t217 + t341) * qJD(1);
t24 = 0.2e1 * t163 + t240 * t217 + ((-t168 * t217 - t202) * qJD(1) + t261 * t322) * qJD(1);
t236 = -t213 * t23 - t217 * t24 + (t213 * t41 - t217 * t42) * qJD(1);
t39 = qJD(1) * t264 + t289;
t21 = -0.2e1 * t147 * t330 + t209 * t39 + t213 * t239;
t22 = t147 * t331 + t209 * t434 + t217 * t239 + t342;
t235 = (t213 * t37 - t217 * t38) * qJD(1) - t21 * t213 - t217 * t22;
t71 = -t376 + t369;
t27 = qJD(1) * t71 + t217 * t253 + t342;
t28 = t253 * t213 + (t430 - t116) * qJD(1);
t234 = (t213 * t430 - t217 * t69) * qJD(1) - t213 * t28 - t217 * t27;
t62 = qJD(1) * t97 - t213 * t247;
t66 = qJD(1) * t101 - t213 * t248;
t229 = qJD(1) * t92 - qJD(3) * t47 - t211 * t62 + t215 * t66;
t61 = -t217 * t247 + (-t213 * t267 + t379) * qJD(1);
t65 = -t217 * t248 + (-t159 * t213 + t383) * qJD(1);
t228 = -qJD(3) * t48 - t211 * t61 + t215 * t65 + t374;
t64 = qJD(1) * t99 - t213 * t250;
t68 = qJD(1) * t103 - t213 * t251;
t227 = qJD(1) * t94 - qJD(2) * t49 - t212 * t64 + t216 * t68;
t63 = -t217 * t250 + (-t213 * t268 + t380) * qJD(1);
t67 = -t217 * t251 + (-t161 * t213 + t384) * qJD(1);
t226 = -qJD(2) * t50 - t212 * t63 + t216 * t67 + t373;
t136 = t267 * qJD(3);
t138 = t159 * qJD(3);
t225 = qJD(1) * t150 - t136 * t211 + t138 * t215 + (-t154 * t215 - t158 * t211) * qJD(3);
t137 = t268 * qJD(2);
t139 = t161 * qJD(2);
t224 = qJD(1) * t152 - t137 * t212 + t139 * t216 + (-t156 * t216 - t160 * t212) * qJD(2);
t223 = -t211 * t415 + t243 * t215;
t222 = -t212 * t416 + t242 * t216;
t143 = qJD(2) * t169;
t110 = (-t212 * t330 - t213 * t327) * pkin(4);
t109 = (-t216 * t326 + t311) * pkin(4);
t108 = (-t211 * t330 - t213 * t321) * pkin(5);
t107 = (-t215 * t320 + t312) * pkin(5);
t106 = t192 + t142;
t91 = -qJD(3) * t324 - t318;
t54 = -t217 * t265 + t366;
t53 = -t217 * t266 + t368;
t46 = t54 * qJD(1);
t45 = t53 * qJD(1);
t44 = -qJD(1) * t75 - t217 * t317 + t163;
t43 = qJD(1) * t429 + (-qJD(1) * t310 - t213 * t353) * pkin(4);
t26 = (-qJD(2) * t143 - t335) * t217 + ((-t169 * t217 - t203) * qJD(1) + t262 * t328) * qJD(1);
t25 = -t213 * t335 + qJD(1) * (-qJD(1) * t314 + t340) + (-t213 * t143 - t262 * t330) * qJD(2);
t20 = t224 * t213 - t217 * t417;
t19 = t225 * t213 - t217 * t418;
t18 = t213 * t417 + t224 * t217;
t17 = t213 * t418 + t225 * t217;
t16 = qJD(2) * t271 + t212 * t67 + t216 * t63;
t15 = t272 * qJD(2) + t212 * t68 + t216 * t64;
t14 = qJD(3) * t273 + t211 * t65 + t215 * t61;
t13 = t274 * qJD(3) + t211 * t66 + t215 * t62;
t12 = t46 + t254;
t11 = t45 + t255;
t10 = t256 + t351;
t9 = t257 + t352;
t1 = [m(3) * (t26 * (t213 * t290 + t337) + t25 * t244 + t56 * t340 + ((-pkin(1) - t169) * t392 + (-t55 * rSges(3,3) + t290 * t56) * t213) * qJD(1) + t423 * qJD(2) * t166) + (t45 + ((t30 - t79 + (t93 + t397) * t217 + t405) * t217 + t404 * t213) * qJD(3)) * t300 + (t46 + ((t32 - t80 + (t95 + t395) * t217 + t403) * t217 + t402 * t213) * qJD(2)) * t304 + (-t265 * qJD(2) + t137 * t216 + t139 * t212) * qJD(1) + m(6) * (t27 * t90 + t28 * t89 + (t69 + t71) * t430) + (-t266 * qJD(3) + t136 * t215 + t138 * t211) * qJD(1) + m(4) * (t113 * t43 + t114 * t44) + m(5) * (t24 * t338 + t41 * t177 + t23 * t336 + t42 * t341 + (t24 * t292 + t41 * t324 + (-t41 * rSges(5,3) + t292 * t42) * qJD(1)) * t213 + (t23 * t291 + t42 * (-rSges(5,1) * t323 - rSges(5,2) * t321 - t316) + t41 * (-t168 - t178) * qJD(1)) * t217) + (-t351 + ((t217 * t298 + t36 - t402) * t217 + (t213 * t298 + t296 + t35) * t213) * qJD(2) + t10) * t307 + (t16 + t18) * t328 / 0.2e1 + (-t352 + ((t217 * t299 + t34 - t404) * t217 + (t213 * t299 + t297 + t33) * t213) * qJD(3) + t9) * t303 + (t14 + t17) * t322 / 0.2e1 + (t22 * (-t77 - t372) + t21 * (t78 + t371) + (-(-qJD(1) * t130 - t173 * t209) * t213 - (t164 * t209 - t147) * t217 - t116 + t289 + (t264 - t372) * qJD(1)) * t38 + (-t371 * qJD(1) - t263 + t369) * t37) * m(7) - (t15 + t20 + t12) * t326 / 0.2e1 - (t13 + t19 + t11) * t320 / 0.2e1 + (((t49 + t52) * t213 + (t50 + t54) * t217) * qJD(2) + ((t47 + t51) * t213 + (t48 + t53) * t217) * qJD(3)) * t406; (-t15 * t217 + t16 * t213 + (t213 * t49 + t217 * t50) * qJD(1)) * t406 + ((-t328 * t365 + t332) * t213 + (t237 + (t213 * t366 + t222) * qJD(2)) * t217) * t307 + ((-t326 * t366 - t332) * t217 + (t237 + (t217 * t365 + t222) * qJD(2)) * t213) * t304 + ((t212 * t344 + t216 * t343) * qJD(1) + (t242 * t212 + t216 * t416) * qJD(2)) * t407 + (qJD(1) * t18 + (t213 * (-t213 * t422 + t226 * t217) - t217 * (-t213 * t421 + t227 * t217) + (t35 * t213 + t36 * t217) * qJD(1)) * t432) * t414 + (qJD(1) * t20 + (t213 * (t226 * t213 + t217 * t422) - t217 * (t227 * t213 + t217 * t421) + (t31 * t213 + t32 * t217) * qJD(1)) * t432) * t413 + (t10 + t256) * t309 + (t12 + t254) * t308 + ((t277 * t327 + (t235 - t329) * t212) * pkin(4) - t109 * t37 - t110 * t38 + t431) * m(7) + ((t275 * t327 + (t234 - t329) * t212) * pkin(4) - t109 * t430 - t110 * t69 + t431) * m(6) + (((qJD(2) * t276 + t91) * t216 + (-qJD(2) * t106 + t236) * t212) * pkin(4) + t106 * t316 - t109 * t41 - t110 * t42) * m(5) + ((-t213 * t56 - t392) * t143 - t55 * (-t169 * t326 + t213 * t260) - t56 * (-t169 * t328 - t217 * t260) + (qJD(1) * t423 - t169 * t219 - t25 * t213 - t26 * t217) * t166) * m(3) + (-t109 * t429 - t110 * t75 + ((-t213 * t75 - t217 * t429) * t327 + (-t213 * t43 - t217 * t44 + (t213 * t429 - t217 * t75) * qJD(1) - t317) * t212) * pkin(4)) * m(4); (-t13 * t217 + t14 * t213 + (t47 * t213 + t217 * t48) * qJD(1)) * t406 + ((-t322 * t367 + t333) * t213 + (t238 + (t213 * t368 + t223) * qJD(3)) * t217) * t303 + ((-t320 * t368 - t333) * t217 + (t238 + (t217 * t367 + t223) * qJD(3)) * t213) * t300 + ((t211 * t346 + t215 * t345) * qJD(1) + (t243 * t211 + t215 * t415) * qJD(3)) * t407 + (qJD(1) * t17 + (-(-t213 * t419 + t229 * t217) * t217 + (-t213 * t420 + t228 * t217) * t213 + (t33 * t213 + t34 * t217) * qJD(1)) * t428) * t414 + (qJD(1) * t19 + (t213 * (t228 * t213 + t217 * t420) - t217 * (t229 * t213 + t217 * t419) + (t29 * t213 + t30 * t217) * qJD(1)) * t428) * t413 + (t9 + t257) * t309 + (t11 + t255) * t308 + ((t277 * t321 + (t235 - t325) * t211) * pkin(5) - t107 * t37 - t108 * t38 + t425) * m(7) + ((t275 * t321 + (t234 - t325) * t211) * pkin(5) - t107 * t430 - t108 * t69 + t425) * m(6) + (t276 * t142 + t236 * t165 + t168 * t91 - t41 * (-t168 * t320 + t213 * t259) - t42 * (-t168 * t322 - t217 * t259)) * m(5); (t21 * t78 - t22 * t77 + t37 * t434 + t38 * t39 - (t38 * t264 + t37 * t287) * t209) * m(7);];
tauc = t1(:);
