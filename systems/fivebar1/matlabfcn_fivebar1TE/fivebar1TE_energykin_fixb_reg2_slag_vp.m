% Calculate inertial parameters regressor of fixed base kinetic energy for
% fivebar1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% 
% Output:
% T_reg [1x(2*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 10:28
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = fivebar1TE_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fivebar1TE_energykin_fixb_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fivebar1TE_energykin_fixb_reg2_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1TE_energykin_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 09:40:38
% EndTime: 2020-04-27 09:40:50
% DurationCPUTime: 13.05s
% Computational Cost: add. (117198->592), mult. (354982->950), div. (802->14), fcn. (44577->6), ass. (0->418)
t259 = (pkin(2) ^ 2);
t494 = -6 * t259;
t215 = sin(qJ(2));
t218 = cos(qJ(1));
t403 = qJD(1) * t218;
t157 = pkin(2) * t403;
t336 = pkin(3) * t157;
t216 = sin(qJ(1));
t217 = cos(qJ(2));
t399 = qJD(2) * t217;
t343 = t216 * t399;
t468 = pkin(2) * pkin(3);
t492 = -t215 * t336 - t343 * t468;
t255 = pkin(3) ^ 2;
t169 = 0.10e2 / 0.3e1 * t255;
t267 = t259 ^ 2;
t253 = t255 ^ 2;
t261 = (pkin(1) ^ 2);
t260 = t261 ^ 2;
t411 = t253 + t260;
t242 = 2 * t261;
t247 = (pkin(5) ^ 2);
t414 = t242 - t247;
t424 = t261 * t247;
t246 = t247 ^ 2;
t250 = pkin(4) ^ 2;
t249 = t250 ^ 2;
t480 = -t246 / 0.6e1 + t249 / 0.6e1;
t99 = t255 * t414 + t411 - t424 - t480;
t288 = t267 + t99;
t491 = (t169 + t414) * t494 - 0.6e1 * t288;
t493 = -2 * t261;
t410 = t255 - t261;
t154 = t410 * t259;
t171 = pkin(3) * t217;
t151 = pkin(1) + t171;
t404 = qJD(1) * t216;
t348 = t151 * t404;
t401 = qJD(2) * t215;
t373 = pkin(3) * t401;
t458 = pkin(2) * t218;
t482 = pkin(2) * t348 + t373 * t458;
t397 = 0.2e1 * pkin(1);
t394 = 0.4e1 * pkin(3);
t196 = -t247 / 0.3e1;
t484 = t196 - t250 / 0.3e1;
t356 = t261 + t484;
t325 = t255 + t356;
t204 = t250 / 0.3e1;
t358 = t247 / 0.3e1 + t204 + t242;
t490 = 0.10e2 / 0.3e1 * t267 - 0.2e1 * (-t255 + t358) * t259 + t325 * t493;
t489 = 0.2e1 * t157;
t172 = pkin(2) * t216;
t488 = 0.2e1 * t172;
t487 = 0.2e1 * t255;
t347 = t215 * t403;
t289 = t343 + t347;
t486 = -0.4e1 * t289;
t431 = t216 * t217;
t433 = t215 * t218;
t113 = -t431 + t433;
t374 = pkin(2) * t404;
t148 = pkin(1) * t374;
t150 = pkin(1) - t458;
t434 = t215 * t216;
t378 = pkin(3) * t434;
t340 = pkin(2) * t378;
t141 = -0.2e1 * t340;
t392 = pkin(1) * t458;
t153 = -0.2e1 * t392;
t408 = t259 + t261;
t351 = t255 + t408;
t446 = t150 * t217;
t379 = pkin(3) * t446;
t72 = t141 + t153 + t351 + 0.2e1 * t379;
t70 = 0.1e1 / t72 ^ 2;
t453 = (t148 + (-t150 * t401 + (-qJD(1) * t113 - t343) * pkin(2)) * pkin(3)) * t70;
t485 = -0.2e1 * t453;
t457 = pkin(2) * t255;
t182 = t217 ^ 2;
t455 = pkin(3) * t182;
t165 = -t255 / 0.3e1 + t261;
t483 = t492 * t165;
t481 = qJD(1) - qJD(2);
t114 = t217 * t218 + t434;
t126 = t150 + t171;
t406 = t260 + t267;
t137 = ((t261 * t494) + t406) * t253;
t405 = t261 - t247;
t350 = t255 + t405;
t146 = t250 + t350;
t170 = -30 * t247 + 60 * t261;
t180 = t182 ^ 2;
t185 = t218 ^ 2;
t183 = t185 ^ 2;
t226 = -2 * t247;
t231 = 0.6e1 * t255;
t240 = 4 * t261;
t266 = pkin(2) * t259;
t256 = t266 ^ 2;
t263 = pkin(1) * t261;
t270 = pkin(3) * t255;
t195 = -t247 / 0.4e1;
t211 = t255 / 0.2e1;
t102 = -0.2e1 / 0.3e1 * t340 + t261 + t211 + t195;
t162 = (t240 + t247) * t255;
t167 = t261 - 0.2e1 / 0.3e1 * t259;
t197 = -t247 / 0.2e1;
t143 = t197 + t351;
t329 = -0.4e1 * t340;
t301 = t143 * t329;
t181 = t217 * t182;
t439 = t181 * t270;
t377 = 0.16e2 * t439;
t426 = t253 * t180;
t386 = 0.8e1 * t426;
t425 = t255 * t182;
t429 = (pkin(1) + pkin(3)) * (pkin(1) - pkin(3));
t132 = t261 + t259 / 0.4e1 + t255 / 0.4e1 - t247 / 0.8e1;
t416 = 0.4e1 / 0.7e1 * t261 - t247 / 0.7e1;
t46 = -0.32e2 / 0.21e2 * t132 * t340 + t267 / 0.7e1 + (0.16e2 / 0.21e2 * t255 + t416) * t259 + t253 / 0.7e1 + t416 * t255 + t260 - 0.3e1 / 0.7e1 * t424 + t246 / 0.42e2 - t249 / 0.42e2;
t210 = t255 / 0.3e1;
t213 = t259 / 0.2e1;
t417 = t213 + t261;
t133 = t195 + t210 + t417;
t198 = -0.2e1 / 0.3e1 * t247;
t209 = 0.4e1 / 0.3e1 * t255;
t462 = 0.4e1 / 0.3e1 * t259;
t50 = -0.8e1 / 0.3e1 * t133 * t340 + t267 / 0.3e1 + (t209 + t196) * t259 + t260 - t253 / 0.3e1 + (t462 + 0.2e1 / 0.3e1 * t255 + t198) * t261 + t246 / 0.18e2 - t249 / 0.18e2;
t31 = t167 * t386 + 0.14e2 * t46 * t425 + t165 * t301 - t410 * t267 + (t162 - 0.10e2 / 0.3e1 * t253 + 0.2e1 * t260 - t424) * t259 + t99 * t429 + (t102 * t377 + 0.6e1 * t171 * t50) * pkin(1);
t413 = t246 - t249;
t327 = -(6 * t424) + 0.6e1 * t260 + t413;
t236 = -3 * t259;
t178 = t236 + t261;
t389 = 0.8e1 * t439;
t341 = pkin(1) * t389;
t129 = t178 * t341;
t391 = pkin(1) * t171;
t361 = 0.6e1 * t391;
t375 = 0.12e2 * t425;
t223 = 0.15e2 * t253;
t224 = 0.15e2 * t255;
t238 = 0.3e1 * t260;
t421 = t246 / 0.2e1 - t249 / 0.2e1;
t328 = -(3 * t424) + t238 + t421;
t199 = -0.3e1 / 0.2e1 * t247;
t241 = 3 * t261;
t418 = t199 + t241;
t175 = t255 + t261;
t450 = t256 + t175 * ((t199 + t242) * t255 - 0.3e1 / 0.2e1 * t424 + t411 + t421);
t39 = t340 * t491 + (t224 + t418) * t267 + t450 + (t223 + (-9 * t247 + 18 * t261) * t255 + t328) * t259;
t56 = t301 + (t231 + t414) * t259 + t288;
t166 = t261 - t259 / 0.3e1;
t110 = t166 * t141;
t428 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
t77 = t143 * t428 + t110;
t33 = t361 * t56 + t375 * t77 + t129 + t39;
t184 = t218 * t185;
t436 = t184 * t266;
t364 = t151 * t436;
t292 = -t340 + t417;
t100 = t195 - t255 + t292;
t177 = -0.3e1 * t255 + t261;
t393 = 0.4e1 * t171;
t239 = 8 * t261;
t380 = t270 * t172;
t474 = 0.4e1 * t253;
t475 = -0.4e1 * t215;
t94 = t380 * t475 + t474 + (4 * t259 + t226 + t239) * t255;
t42 = t141 * t429 + t143 * t177 + t182 * t94 + (t100 * t393 + t389) * pkin(1);
t407 = t260 - t253;
t44 = t166 * t301 - t256 + (-t169 - t405) * t267 + (t162 + t407 + t480) * t259 + t99 * t261;
t168 = t175 ^ 2;
t442 = t168 * (-t250 + t350);
t445 = t151 * t218;
t227 = -5 * t247;
t228 = 0.7e1 * t253;
t230 = 0.7e1 * t255;
t49 = (t230 + t418) * t267 + (t228 + (t227 + 10 * t261) * t255 + t328) * t259 + t450;
t149 = -0.12e2 * pkin(1) * t270 + t263 * t394;
t423 = t261 * t255;
t164 = -0.8e1 * t253 + 0.12e2 * t423;
t338 = pkin(1) * t377;
t63 = t149 * t217 + t164 * t182 + t338 + t386 + t411 - 0.6e1 * t423;
t80 = t141 * t428 + t143 * t178;
t17 = -0.32e2 * t42 * t364 + 0.16e2 * t137 * t180 + 0.24e2 * t44 * t425 + (t226 + t240 + 0.28e2 * t255) * t256 + t146 * t442 + (t170 * t253 + 0.24e2 * t31 * t185 + t231 * t327 + t242 * t413 - 0.6e1 * t260 * t247 + 0.4e1 * t263 ^ 2 + 0.28e2 * t270 ^ 2) * t259 + 0.8e1 * (-t33 * t445 - t378 * t49) * pkin(2) + (0.8e1 * t171 * t39 + 0.32e2 * t439 * t80) * pkin(1) + (t170 * t255 + 0.16e2 * t183 * t63 + 0.70e2 * t253 + t267 + t327) * t267;
t152 = 0.2e1 * t391;
t103 = 0.4e1 / 0.3e1 * t425 + t152 + t165;
t437 = t183 * t267;
t332 = -0.24e2 * t103 * t437;
t334 = 0.8e1 * t364;
t456 = pkin(3) * t215;
t369 = -t456 / 0.2e1;
t384 = 0.4e1 * t425;
t396 = 0.4e1 * pkin(1);
t440 = t181 * t253;
t194 = -t247 / 0.6e1;
t420 = t194 - t250 / 0.6e1;
t357 = t261 + t420;
t117 = t462 + t211 + t357;
t326 = t213 + t357;
t118 = t209 + t326;
t64 = -t117 * t456 + t118 * t172;
t342 = 0.20e2 / 0.3e1 * t255;
t203 = 0.2e1 / 0.3e1 * t250;
t359 = 0.2e1 / 0.3e1 * t247 + t203 + t240;
t360 = 0.4e1 / 0.3e1 * t247 + 0.4e1 / 0.3e1 * t250 + t493;
t75 = -t267 + (-t342 + t359) * t259 - 0.3e1 * t253 + t360 * t255 + t260;
t121 = t255 + t326;
t472 = 2 * t259;
t160 = t472 - t410;
t78 = t121 * t172 + t160 * t369;
t81 = -t154 - 0.5e1 / 0.3e1 * t253 + t358 * t255 + t261 * t356;
t34 = t64 * t384 + t75 * t369 + (-0.8e1 / 0.3e1 * t426 + t81) * t172 + (t171 * t78 - t215 * t440) * t396;
t229 = 0.5e1 * t253;
t225 = 0.10e2 * t255;
t415 = t225 + t242;
t208 = -0.2e1 / 0.3e1 * t250;
t419 = t198 + t208;
t173 = -t247 - t250;
t158 = t241 + t173;
t444 = t158 * t255;
t65 = t267 + (t415 + t419) * t259 + t229 + 0.2e1 * t444 + t261 * (t261 + t419);
t234 = 5 * t267;
t353 = t198 + t175;
t79 = t234 + (0.5e1 * t255 + t158) * t472 + (t208 + t353) * t175;
t43 = t172 * t65 - t456 * t79;
t235 = 3 * t259;
t123 = t235 + t325;
t232 = 0.3e1 * t255;
t161 = t232 + t408;
t124 = t161 + t484;
t68 = -t123 * t456 + t124 * t172;
t131 = t211 + t259 + t420;
t362 = t215 * t428;
t83 = pkin(3) * t362 + t131 * t488;
t35 = -0.4e1 * t83 * t425 + (-0.8e1 * t181 * t380 + t393 * t68) * pkin(1) + t43;
t435 = t185 * t259;
t376 = -0.12e2 * t435;
t385 = -0.6e1 * t425;
t76 = -(3 * t267) + (-t342 + t360) * t259 + t359 * t255 + t407;
t45 = t172 * t76 + t456 * t490;
t459 = pkin(1) * t217;
t136 = 0.10e2 * t444;
t441 = t173 * t261;
t47 = t256 + (0.21e2 * t255 + t158) * t267 + (t136 + t238 + 0.35e2 * t253 + 0.2e1 * t441) * t259 + (t228 + (t227 + t239 - 0.5e1 * t250) * t255 + t261 * (-t250 + t405)) * t175;
t470 = 0.24e2 * t166;
t355 = t197 - t250 / 0.2e1 + t261;
t120 = 0.3e1 / 0.2e1 * t259 + t232 + t355;
t135 = t172 + 0.2e1 * t456;
t395 = 0.2e1 * pkin(3);
t381 = t216 * t457;
t127 = -t215 * t270 + t381;
t447 = t127 * t182;
t48 = t177 * t172 + 0.4e1 * t447 + (t120 * t215 + t135 * t459) * t395;
t51 = 0.7e1 * t256 + (t230 + t158) * t234 + (t136 + 0.21e2 * t253 + 0.9e1 * t260 + 0.6e1 * t441) * t259 + t442;
t122 = t235 + 0.3e1 / 0.2e1 * t255 + t355;
t370 = pkin(3) * t178 / 0.2e1;
t82 = t122 * t172 + t215 * t370;
t23 = t82 * t338 + t48 * t334 + t45 * t385 + t34 * t376 + (-0.6e1 * t43 * t459 + (t51 + t332) * t215) * pkin(3) + (0.6e1 * t35 * t445 + (t426 * t470 - t47) * t216) * pkin(2);
t349 = t259 + t405;
t324 = t255 + t349;
t139 = -t250 + t324;
t106 = t153 + t139;
t387 = 0.2e1 * t435;
t409 = -t259 + t261;
t109 = t153 + t387 + t409;
t477 = -0.2e1 * t182;
t310 = t109 * t477 - t405;
t382 = pkin(2) * t434;
t463 = pkin(5) - pkin(4);
t464 = -pkin(4) - pkin(5);
t86 = t141 + t106;
t262 = sqrt(0.4e1 * t154 * t185 + 0.4e1 * t139 * t392 - t253 - (t261 + (pkin(2) - t464) * (pkin(2) + t464)) * (t261 + (pkin(2) - t463) * (pkin(2) + t463)) + (t236 + t250 + t310) * t487 + (t106 * t382 - t446 * t86) * t394);
t16 = t126 * t17 + t23 * t262;
t105 = t153 + t250 + t324;
t57 = t105 * t150 + 0.2e1 * t109 * t171;
t476 = 0.4e1 * t182;
t61 = t105 * t217 + (t476 - 0.2e1) * t150 * pkin(3);
t140 = -pkin(3) + t382;
t96 = -t140 + t446;
t30 = t172 * t61 + t215 * t57 + t262 * t96;
t138 = t232 + t250 + t349;
t87 = t138 + t153 + t329;
t98 = pkin(2) * t431 + t150 * t215;
t27 = -t87 * t446 + t98 * t262 + (-0.2e1 * pkin(1) * t140 * t218 + t138 * t434) * pkin(2) + (-t235 - t250 - t255 + t310 + t387) * pkin(3);
t108 = t152 + t384 + t177;
t335 = -0.8e1 * t364;
t412 = -t247 + t250;
t352 = t241 + t412;
t125 = t196 + t204 + t351;
t89 = -0.4e1 / 0.9e1 * t340 + t261 + t259 / 0.3e1 + t210 + t250 / 0.9e1 - t247 / 0.9e1;
t97 = t250 / 0.6e1 + t194 + t292;
t40 = t165 * t141 + t125 * t429 + 0.6e1 * t89 * t425 + (t171 * t97 + t439) * t396;
t155 = t409 * t255;
t354 = t198 + t203 + t242;
t422 = (t203 + t353) * t175 + t267;
t53 = t125 * t329 + (t231 + t354) * t259 + t422;
t95 = t141 + t125;
t41 = t155 * t476 + 0.4e1 * t391 * t95 + t53;
t71 = t125 * t428 + t110;
t73 = (t169 + t354) * t259 + t422;
t25 = t108 * t335 + t129 + t71 * t375 + t53 * t361 + t256 + (t224 + t352) * t267 + t168 * t146 + (0.12e2 * t185 * t40 + t231 * t352 + t242 * t412 + t223 + t238) * t259 + 0.6e1 * (-t378 * t73 - t41 * t445) * pkin(2);
t107 = t415 * t259 + t229 + t406 + 0.6e1 * t423;
t115 = t234 + (t225 + (6 * t261)) * t259 + t168;
t388 = -0.4e1 * t435;
t130 = t161 * t172;
t134 = t172 - t456;
t159 = t235 + t175;
t443 = t159 * t215;
t54 = t381 * t477 + t130 + (0.2e1 * t134 * t459 - t443) * pkin(3);
t174 = t487 + t259;
t55 = t174 * t456 + 0.2e1 * t447 + (-t410 + t152) * t172;
t88 = -pkin(3) * t443 + t130;
t156 = pkin(2) * t474 + 0.8e1 * t255 * t266;
t471 = 0.4e1 * t270;
t91 = t156 * t216 + t362 * t471;
t32 = t55 * t388 + t182 * t91 + (-0.4e1 * t88 * t459 + (t115 + t334) * t215) * pkin(3) + (0.4e1 * t54 * t445 + (-t107 + t341) * t216) * pkin(2);
t22 = t126 * t25 + t262 * t32;
t467 = 0.1e1 / t22 / 0.4e1;
t367 = t27 * t467;
t461 = -t262 / 0.4e1;
t290 = t16 * t367 + t30 * t461;
t252 = 0.1e1 / pkin(4) ^ 2;
t427 = 0.1e1 / pkin(5) * t252;
t69 = 0.1e1 / t72;
t371 = t69 * t427;
t13 = t290 * t371;
t366 = t30 * t467;
t460 = t262 / 0.4e1;
t291 = t16 * t366 + t27 * t460;
t14 = t291 * t371;
t6 = t113 * t13 + t114 * t14;
t479 = t6 ^ 2;
t478 = 0.8e1 * pkin(1);
t473 = -0.8e1 * t255;
t469 = pkin(1) * pkin(2);
t466 = t27 / 0.2e1;
t465 = t69 / 0.2e1;
t286 = t217 * t289;
t285 = pkin(3) * t286;
t284 = t150 * t285;
t344 = t215 * t399;
t321 = t109 * t344;
t402 = qJD(1) * t259;
t145 = 0.2e1 * t148;
t322 = t216 * t218 * t402;
t309 = -0.4e1 * t322;
t101 = t145 + t309;
t449 = t101 * t182;
t454 = ((0.8e1 * t321 - 0.4e1 * t449) * t255 + (-0.4e1 * t139 * t469 - 0.8e1 * t154 * t218) * t404 + (t215 * t216 ^ 2 * t402 * t478 + (-0.8e1 * t148 * t217 + 0.4e1 * t401 * t86) * t150 + 0.4e1 * (-t217 * t404 * t86 + t106 * t289 + 0.2e1 * t284) * pkin(2)) * pkin(3)) / t262;
t287 = t289 * pkin(3);
t90 = pkin(2) * t287;
t452 = t217 * t90;
t282 = t90 * t425;
t281 = t166 * t282;
t398 = qJD(2) * t270;
t316 = t182 * t215 * t398;
t302 = -0.24e2 * t316;
t294 = pkin(1) * t302;
t451 = t178 * t294 - 0.24e2 * t281;
t448 = (t157 * t255 - t217 * t398) * t182;
t438 = t182 * t270;
t432 = t215 * t262;
t430 = t217 * t262;
t400 = qJD(2) * t216;
t116 = -t373 + t374;
t337 = pkin(1) * t373;
t144 = -0.2e1 * t337;
t179 = t215 ^ 2;
t283 = t255 * t286 * t469;
t300 = t185 * t266 * t348;
t293 = -0.24e2 * t300;
t307 = pkin(1) * t316;
t295 = -0.48e2 * t307;
t372 = pkin(3) * t399;
t299 = t157 * t65 - t372 * t79;
t317 = t255 * t344;
t303 = -0.24e2 * t317;
t304 = t157 * t439;
t305 = t157 * t426;
t345 = qJD(2) * t436;
t330 = pkin(3) * t345;
t306 = t215 * t330;
t308 = -0.6e1 * t317;
t315 = -0.6e1 * t337;
t318 = t127 * t344;
t319 = t399 * t428;
t346 = t181 * t401;
t320 = t253 * t346;
t323 = t184 * t267 * t404;
t333 = -0.2e1 * t344;
t363 = t215 * t437;
t365 = t454 / 0.2e1;
t383 = pkin(2) * t445;
t10 = (t332 * t372 - 0.24e2 * pkin(3) * (-0.8e1 / 0.3e1 * t317 + t144) * t363 + 0.96e2 * t103 * t323 * t456 + (t177 * t157 + 0.2e1 * t120 * t372 + (t217 * t489 + (-0.2e1 * t135 * t215 + 0.4e1 * t455) * qJD(2)) * pkin(3) * pkin(1) - 0.8e1 * t318 + 0.4e1 * t448) * t334 + (-0.8e1 / 0.3e1 * t305 + (-t117 * t372 + t118 * t157) * t384 + t81 * t157 - t75 * t372 / 0.2e1 + (0.32e2 / 0.3e1 * t440 * t172 + t217 * t64 * t473) * t401 + (0.4e1 * t121 * t217 * t336 + ((0.12e2 * t179 * t182 - 0.4e1 * t180) * t253 + (-0.2e1 * t160 * t455 + t475 * t78) * pkin(3)) * qJD(2)) * pkin(1)) * t376 + 0.24e2 * t34 * t322 + 0.6e1 * ((-0.4e1 * (pkin(3) * t319 + t131 * t489) * t182 + 0.8e1 * t83 * t344) * t255 + (-0.8e1 * t304 + (-t123 * t372 + t124 * t157) * t393 + (-0.4e1 * pkin(3) * t68 + 0.24e2 * t182 * t380) * t401) * pkin(1) + t299) * t383 + t305 * t470 - 0.96e2 * t166 * t320 * t172 + (t122 * t157 + t370 * t399) * t338 + t82 * t295 + (t157 * t76 + t372 * t490) * t385 + 0.12e2 * t45 * t317 - 0.6e1 * t299 * t391 + 0.6e1 * t43 * t337 - t47 * t157 + t51 * t372 + (-0.8e1 * t306 + t293) * t48 - 0.6e1 * t482 * t35) * t262 + t23 * t365 + 0.8e1 * (0.2e1 * (-0.48e2 * pkin(1) * t438 - 0.2e1 * t164 * t217 - t149 - 0.32e2 * t440) * qJD(2) * t363 - 0.8e1 * t63 * t323 - 0.4e1 * (t94 * t333 + (t302 + (-t100 * t401 - t452) * t394) * pkin(1) + (-0.2e1 * t287 * t429 + t438 * t486) * pkin(2)) * t364 + 0.3e1 * (-0.32e2 * t167 * t320 + (-0.2e1 / 0.3e1 * t343 - 0.2e1 / 0.3e1 * t347) * t338 * t468 + t102 * t295 - 0.64e2 / 0.3e1 * t132 * t282 - 0.28e2 * t46 * t317 + 0.6e1 * (-0.8e1 / 0.3e1 * t343 - 0.8e1 / 0.3e1 * t347) * t133 * t459 * t457 + t50 * t315 + 0.4e1 * t483 * t143) * t435 - 0.6e1 * t31 * t322 - (t77 * t303 + (-0.6e1 * pkin(1) * t56 * t401 + (-0.24e2 * pkin(1) * t143 * t285 + t289 * t491) * pkin(2)) * pkin(3) + t451) * t383 - 0.8e1 * t137 * t346 - t90 * t341 * t428 - 0.12e2 * t80 * t307 - 0.12e2 * t143 * t281 + t44 * t308 + t283 * t491 - t39 * t337 + t492 * t49 + (0.4e1 * t306 + 0.12e2 * t300) * t42 + t482 * t33) * t126 + t17 * t116;
t18 = t96 * t365 + (-t150 * t432 + t217 * t57) * qJD(2) + (t101 * t217 - t109 * t401) * t215 * t395 + ((-t430 + (-t105 - 0.8e1 * t379) * t215) * t400 + ((t61 - t432) * t218 + (t430 + (t150 * t397 + t105) * t215 + (-pkin(3) + t459 + 0.2e1 * t455) * t488) * t216) * qJD(1)) * pkin(2);
t19 = t98 * t365 + (-t145 * t217 + (t215 * t87 + t430) * qJD(2)) * t150 + (t309 + 0.4e1 * t321 - 0.2e1 * t449) * pkin(3) + ((t138 * t217 - t432) * t400 + (t114 * t262 + t138 * t433 - t431 * t87) * qJD(1) + 0.4e1 * t284 + (t140 * t404 - t289 * t458) * t397) * pkin(2);
t298 = t157 * t161 - t159 * t372;
t368 = -((0.8e1 * t217 * t151 * t330 + t179 * t345 * t473 + t293 * t456 + (t174 * t372 - 0.4e1 * t318 + 0.2e1 * t448 + (-t410 * t403 + (-t215 * t400 + t217 * t403) * pkin(1) * t395) * pkin(2)) * t388 + 0.8e1 * t55 * t322 + 0.4e1 * ((0.4e1 * t215 * t343 + t403 * t477) * t457 + ((t157 - t372) * t171 - t134 * t373) * t397 + t298) * t383 + t304 * t478 + t294 * t172 + (t156 * t403 + t319 * t471) * t182 + t91 * t333 - 0.4e1 * t298 * t391 + 0.4e1 * t88 * t337 - t107 * t157 + t115 * t372 - 0.4e1 * t482 * t54) * t262 + t32 * t365 + ((t144 - 0.8e1 * t317) * t335 + 0.24e2 * (-0.6e1 * t307 + 0.3e1 * (-0.4e1 / 0.9e1 * t343 - 0.4e1 / 0.9e1 * t347) * t425 * t468 + t89 * t308 - t90 * t152 + t97 * t144 + t483) * t435 - 0.24e2 * t40 * t322 - 0.6e1 * (-0.8e1 * t155 * t344 + (t125 * pkin(2) * t486 + (-0.4e1 * t401 * t95 - 0.8e1 * t452) * pkin(1)) * pkin(3)) * t383 + t71 * t303 - 0.24e2 * t125 * t283 + t53 * t315 + t451 + (0.8e1 * t306 + 0.24e2 * t300) * t108) * t126 + t25 * t116 + 0.6e1 * (t482 * t41 + t492 * t73) * t126) / t22 ^ 2 / 0.4e1;
t2 = (t290 * t485 + (t10 * t367 + t18 * t461 - t30 * t454 / 0.8e1 + (t19 * t467 + t27 * t368) * t16) * t69) * t427;
t3 = (t291 * t485 + (t19 * t460 + t27 * t454 / 0.8e1 + t10 * t366 + (t18 * t467 + t30 * t368) * t16) * t69) * t427;
t8 = -t113 * t14 + t114 * t13;
t280 = t8 ^ 2;
t5 = 0.1e1 / t280;
t66 = t481 * t113;
t67 = t481 * t114;
t1 = qJD(1) + ((t113 * t2 + t114 * t3 - t13 * t67 + t14 * t66) / t8 - (-t113 * t3 + t114 * t2 + t13 * t66 + t14 * t67) * t6 * t5) / (t479 * t5 + 0.1e1);
t390 = pkin(2) * qJD(1) * t1;
t251 = 0.1e1 / pkin(4);
t279 = t27 ^ 2;
t26 = 0.1e1 / t279;
t28 = t30 ^ 2;
t11 = qJD(2) + ((t18 * t465 - t30 * t453) / t466 - 0.2e1 * (t19 * t465 - t27 * t453) * t30 * t26) * pkin(4) / (t26 * t28 + 0.1e1) * t251 * t72;
t311 = pkin(3) * qJD(2) * t11 * t251 * t69;
t245 = qJD(1) ^ 2;
t244 = qJD(2) ^ 2;
t4 = [0, 0, 0, 0, 0, t245 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 ^ 2 / 0.2e1, -t8 * t390, t6 * t390, 0, (t479 / 0.2e1 + t280 / 0.2e1) * t259 * t245, 0, 0, 0, 0, 0, t244 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 ^ 2 / 0.2e1, t311 * t466, -t30 * t311 / 0.2e1, 0, (t28 / 0.8e1 + t279 / 0.8e1) * t70 * t255 * t252 * t244;];
T_reg = t4;
