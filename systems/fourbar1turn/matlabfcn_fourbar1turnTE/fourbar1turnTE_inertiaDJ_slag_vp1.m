% Calculate time derivative of joint inertia matrix for
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
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [2x2]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1turnTE_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_inertiaDJ_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnTE_inertiaDJ_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_inertiaDJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnTE_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnTE_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnTE_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:18:37
% EndTime: 2020-04-12 19:19:13
% DurationCPUTime: 17.14s
% Computational Cost: add. (136704->541), mult. (193902->1013), div. (6920->14), fcn. (54247->6), ass. (0->367)
t209 = pkin(2) ^ 2;
t210 = pkin(1) ^ 2;
t202 = cos(qJ(2));
t411 = pkin(2) * t202;
t350 = -0.2e1 * pkin(1) * t411 + t210;
t184 = t209 + t350;
t348 = pkin(3) ^ 2 - pkin(4) ^ 2;
t176 = t184 - t348;
t189 = pkin(1) - t411;
t200 = sin(qJ(2));
t430 = -pkin(3) - pkin(4);
t169 = (pkin(2) - t430) * (pkin(2) + t430) + t350;
t429 = pkin(4) - pkin(3);
t170 = (pkin(2) - t429) * (pkin(2) + t429) + t350;
t365 = t169 * t170;
t211 = sqrt(-t365);
t359 = t200 * t211;
t138 = -pkin(2) * t359 + t176 * t189;
t133 = 0.1e1 / t138 ^ 2;
t412 = pkin(2) * t200;
t167 = t176 * t412;
t139 = t189 * t211 + t167;
t372 = t133 * t139;
t132 = 0.1e1 / t138;
t302 = pkin(1) * pkin(2) * (-t169 - t170);
t149 = t200 * t302;
t197 = t200 ^ 2;
t414 = pkin(1) * t209;
t333 = t197 * t414;
t155 = 0.1e1 / t211;
t367 = t155 * t189;
t105 = t149 * t367 + 0.2e1 * t333 + (t176 * t202 + t359) * pkin(2);
t206 = 0.1e1 / pkin(4);
t181 = 0.1e1 / t184 ^ 2;
t413 = pkin(2) * t181;
t338 = pkin(1) * t413;
t310 = t200 * t338;
t180 = 0.1e1 / t184;
t423 = t180 / 0.2e1;
t78 = (t105 * t423 - t139 * t310) * t206;
t403 = t132 * t78;
t339 = 0.2e1 * t189 * pkin(1);
t354 = t202 * t211;
t368 = t155 * t149;
t103 = t167 + (-t354 + (t339 - t368) * t200) * pkin(2);
t424 = -t180 / 0.2e1;
t76 = (t103 * t424 + t138 * t310) * t206;
t459 = -0.2e1 * t76 * t372 - 0.2e1 * t403;
t175 = t184 + t348;
t190 = pkin(1) * t202 - pkin(2);
t137 = -pkin(1) * t359 - t175 * t190;
t129 = 0.1e1 / t137;
t130 = 0.1e1 / t137 ^ 2;
t417 = pkin(1) * t175;
t168 = t200 * t417;
t140 = -t190 * t211 + t168;
t373 = t130 * t140;
t442 = -0.2e1 * t190;
t341 = pkin(2) * t442;
t102 = t168 + (-t354 + (t341 - t368) * t200) * pkin(1);
t208 = 0.1e1 / pkin(3);
t300 = -0.2e1 * t310;
t77 = (t102 * t180 + t137 * t300) * t208;
t361 = t197 * t210;
t366 = t155 * t190;
t104 = -t149 * t366 + 0.2e1 * pkin(2) * t361 + (t175 * t202 + t359) * pkin(1);
t79 = (t104 * t180 + t140 * t300) * t208;
t458 = t129 * t79 - t77 * t373;
t396 = t139 * rSges(5,2);
t401 = t138 * rSges(5,1);
t457 = t396 + t401;
t344 = qJD(2) * t200;
t331 = pkin(2) * t344;
t309 = pkin(1) * t331;
t293 = t181 * t309;
t456 = t293 * t457;
t363 = t180 * t206;
t425 = t139 / 0.2e1;
t427 = -t138 / 0.2e1;
t109 = (Icges(5,5) * t425 + Icges(5,6) * t427) * t363;
t201 = sin(qJ(1));
t203 = cos(qJ(1));
t380 = Icges(5,4) * t139;
t110 = (t380 / 0.2e1 + Icges(5,2) * t427) * t363;
t370 = t139 * t110;
t381 = Icges(5,4) * t138;
t316 = -t381 / 0.2e1;
t111 = (Icges(5,1) * t425 + t316) * t363;
t371 = t138 * t111;
t259 = t371 / 0.2e1 + t370 / 0.2e1;
t238 = t259 * t363;
t384 = Icges(5,1) * t138;
t241 = (-t384 / 0.2e1 - t380 / 0.2e1) * t363;
t92 = -Icges(5,5) * t203 + t201 * t241;
t393 = t139 * t92;
t377 = Icges(5,2) * t139;
t240 = (t316 - t377 / 0.2e1) * t363;
t90 = -Icges(5,6) * t203 + t201 * t240;
t400 = t138 * t90;
t455 = -t203 * t109 - t201 * t238 + (t393 / 0.2e1 - t400 / 0.2e1) * t363;
t93 = Icges(5,5) * t201 + t203 * t241;
t392 = t139 * t93;
t91 = Icges(5,6) * t201 + t203 * t240;
t399 = t138 * t91;
t454 = t201 * t109 - t203 * t238 + (t392 / 0.2e1 - t399 / 0.2e1) * t363;
t382 = Icges(3,4) * t202;
t283 = -Icges(3,2) * t200 + t382;
t160 = Icges(3,6) * t201 + t203 * t283;
t383 = Icges(3,4) * t200;
t284 = Icges(3,1) * t202 - t383;
t162 = Icges(3,5) * t201 + t203 * t284;
t276 = t160 * t200 - t162 * t202;
t452 = t201 * t276;
t159 = -Icges(3,6) * t203 + t201 * t283;
t161 = -Icges(3,5) * t203 + t201 * t284;
t278 = t159 * t200 - t161 * t202;
t451 = t203 * t278;
t356 = t202 * t140;
t421 = t200 / 0.2e1;
t257 = t356 / 0.2e1 + t137 * t421;
t245 = t257 * t180;
t450 = t208 * t245;
t449 = (t137 * t197 + t200 * t356) * t338;
t448 = t201 * t456;
t343 = qJD(2) * t202;
t345 = qJD(1) * t203;
t447 = -t200 * t345 - t201 * t343;
t282 = Icges(3,5) * t202 - Icges(3,6) * t200;
t157 = -Icges(3,3) * t203 + t201 * t282;
t445 = qJD(1) * t157;
t142 = qJD(2) * t149;
t178 = Icges(3,2) * t202 + t383;
t179 = Icges(3,1) * t200 + t382;
t273 = t178 * t200 - t179 * t202;
t444 = qJD(1) * t273 + t282 * qJD(2);
t443 = 2 * m(3);
t441 = m(4) / 0.2e1;
t440 = m(5) / 0.2e1;
t369 = t142 * t155;
t324 = t200 * t369;
t83 = (-t324 + (-t354 + (t176 + t339) * t200) * qJD(2)) * pkin(2);
t435 = -t83 / 0.2e1;
t304 = qJD(2) * t333;
t318 = t211 * t344;
t352 = (t176 * t343 + t318) * pkin(2);
t86 = t142 * t367 + 0.2e1 * t304 + t352;
t434 = -t86 / 0.2e1;
t433 = t86 / 0.2e1;
t432 = -t90 / 0.2e1;
t431 = -t91 / 0.2e1;
t136 = t140 ^ 2;
t118 = t130 * t136 + 0.1e1;
t115 = 0.1e1 / t118;
t410 = pkin(3) * t184;
t336 = t115 * t410;
t428 = t336 * t458 + 0.1e1;
t426 = -t139 / 0.2e1;
t422 = -t200 / 0.2e1;
t420 = t201 / 0.2e1;
t419 = -t202 / 0.2e1;
t418 = -t203 / 0.2e1;
t416 = pkin(1) * t181;
t415 = pkin(1) * t203;
t409 = pkin(4) * t184;
t135 = t139 ^ 2;
t117 = t133 * t135 + 0.1e1;
t113 = 0.1e1 / t117;
t141 = (t202 * t302 - 0.4e1 * t209 * t361) * qJD(2);
t301 = 0.4e1 / t365 * t142 * t368;
t285 = -t301 / 0.4e1;
t220 = (t200 * t285 + 0.2e1 * (t141 * t422 - t142 * t202) * t155) * t180;
t320 = qJD(2) * t361;
t286 = t180 * t181 * t209 * t320;
t340 = 0.2e1 * t409;
t294 = t139 * t76 * t340;
t305 = 0.6e1 * t200 * t343;
t326 = t180 * t369;
t334 = t113 * t409;
t362 = t181 * t200;
t364 = t180 * t202;
t402 = t132 * t133 * t83;
t5 = 0.2e1 * (t133 * t294 + t340 * t403) * (-t135 * t402 + t372 * t86) / t117 ^ 2 + 0.2e1 * ((-t76 * t86 + t78 * t83) * t133 + (-((t189 * t301 / 0.4e1 + t141 * t367 + t305 * t414) * t423 + 0.4e1 * t139 * t286 + ((t326 / 0.2e1 - t86 * t416) * t200 + ((t354 + (-t176 + t368) * t200) * t423 + (-t105 * t200 - t139 * t202) * t416) * qJD(2)) * pkin(2)) * t132 - ((0.4e1 * t304 + t352) * t424 - 0.4e1 * t138 * t286 + (-t220 / 0.2e1 + (t83 * t362 + (-t189 * t364 + (t103 * t200 + t138 * t202) * t181) * qJD(2)) * pkin(1)) * pkin(2)) * t372) * t206) * t334 + 0.2e1 * (pkin(4) * t309 * t459 + t294 * t402) * t113;
t408 = t203 * t5;
t407 = t5 * t201;
t406 = rSges(3,3) * t201;
t82 = (-t324 + (-t354 + (t175 + t341) * t200) * qJD(2)) * pkin(1);
t404 = t129 * t130 * t82;
t398 = t138 * t92;
t397 = t138 * t93;
t395 = t139 * t90;
t394 = t139 * t91;
t30 = t334 * t459;
t391 = t180 * t30;
t194 = t201 * rSges(5,3);
t390 = t201 * t30;
t375 = Icges(5,6) * t139;
t378 = Icges(5,5) * t138;
t239 = (-t378 / 0.2e1 - t375 / 0.2e1) * t363;
t89 = Icges(5,3) * t201 + t203 * t239;
t389 = t201 * t89;
t388 = t203 * rSges(3,3);
t387 = t203 * rSges(5,3);
t386 = t203 * t30;
t88 = -Icges(5,3) * t203 + t201 * t239;
t385 = t203 * t88;
t360 = t200 * t201;
t358 = t201 * t202;
t357 = t202 * t137;
t355 = t202 * t203;
t353 = t203 * t206;
t351 = pkin(1) * t318 + t343 * t417;
t349 = t201 ^ 2 + t203 ^ 2;
t158 = Icges(3,3) * t201 + t203 * t282;
t347 = qJD(1) * t158;
t346 = qJD(1) * t201;
t342 = qJD(2) * t203;
t337 = 0.2e1 * t181;
t335 = t129 * t410;
t303 = pkin(2) * t320;
t87 = -t142 * t366 + 0.2e1 * t303 + t351;
t261 = t419 * t82 + t421 * t87;
t223 = qJD(2) * t257 + t261;
t236 = (-t140 * t197 + t200 * t357) * t338;
t231 = qJD(2) * t236;
t258 = t357 / 0.2e1 + t140 * t422;
t34 = (t203 * t231 + (t203 * t223 + t258 * t346) * t180) * t208;
t262 = t419 * t87 + t422 * t82;
t222 = qJD(2) * t258 - t262;
t232 = qJD(2) * t449;
t35 = (-t203 * t232 + (t203 * t222 - t257 * t346) * t180) * t208;
t332 = t34 * rSges(4,1) + t35 * rSges(4,2) + rSges(4,3) * t345;
t246 = t258 * t180;
t107 = t208 * t246;
t100 = t203 * t107;
t101 = t203 * t450;
t71 = -t100 * rSges(4,1) + t101 * rSges(4,2) + t201 * rSges(4,3);
t325 = t201 * t363;
t313 = -t353 / 0.2e1;
t297 = t180 * t313;
t287 = rSges(5,2) * t297;
t288 = rSges(5,1) * t297;
t97 = t138 * t288 + t139 * t287 + t194;
t322 = qJD(1) * t428;
t319 = t201 * t344;
t315 = -t363 / 0.2e1;
t314 = t363 / 0.2e1;
t312 = t346 / 0.2e1;
t311 = t345 / 0.2e1;
t299 = 0.4e1 * t349;
t295 = 0.2e1 * t140 * t77 * t410;
t292 = rSges(3,1) * t202 - rSges(3,2) * t200;
t183 = rSges(3,1) * t200 + rSges(3,2) * t202;
t290 = -t385 + t389;
t265 = t401 / 0.2e1 + t396 / 0.2e1;
t96 = -t265 * t325 - t387;
t289 = t201 * t96 + t203 * t97;
t177 = Icges(3,5) * t200 + Icges(3,6) * t202;
t279 = t159 * t202 + t161 * t200;
t277 = t160 * t202 + t162 * t200;
t163 = rSges(3,1) * t358 - rSges(3,2) * t360 - t388;
t164 = t203 * t292 + t406;
t275 = t163 * t203 - t164 * t201;
t274 = t163 * t201 + t164 * t203;
t269 = t138 * t293;
t268 = t139 * t293;
t267 = 0.8e1 * t286;
t264 = -t398 / 0.2e1 - t395 / 0.2e1;
t263 = t397 / 0.2e1 + t394 / 0.2e1;
t255 = t312 * t363 * t457 + rSges(5,3) * t345 + t86 * t287 + t83 * t288 + t353 * t456;
t254 = qJD(2) * t179;
t253 = qJD(2) * t178;
t252 = qJD(2) * t177;
t84 = t387 + (t265 * t363 - pkin(1)) * t201;
t85 = t97 + t415;
t251 = -0.2e1 * t201 * t85 - 0.2e1 * t203 * t84;
t98 = t201 * t107;
t99 = t201 * t450;
t70 = -rSges(4,1) * t98 + rSges(4,2) * t99 - rSges(4,3) * t203;
t250 = t138 * t312 + t418 * t83;
t249 = t139 * t312 + t418 * t86;
t248 = t138 * t311 + t420 * t83;
t247 = t139 * t311 + t420 * t86;
t36 = (t201 * t231 + (t201 * t223 - t258 * t345) * t180) * t208;
t37 = (-t201 * t232 + (t201 * t222 + t257 * t345) * t180) * t208;
t244 = rSges(4,1) * t36 + rSges(4,2) * t37 + rSges(4,3) * t346;
t243 = t264 * t363;
t242 = t263 * t363;
t235 = t424 * t83 + t269;
t234 = t423 * t86 - t268;
t233 = t203 * t243;
t230 = (t380 + t384) * t293;
t229 = (t377 + t381) * t293;
t228 = (t375 + t378) * t293;
t227 = (t370 + t371) * t293;
t43 = Icges(5,6) * t346 + (t201 * t229 + (-Icges(5,4) * t248 - Icges(5,2) * t247) * t180) * t206;
t45 = Icges(5,5) * t346 + (t201 * t230 + (-Icges(5,1) * t248 - Icges(5,4) * t247) * t180) * t206;
t226 = t426 * t43 + t427 * t45 + t432 * t86 + t435 * t92;
t42 = Icges(5,6) * t345 + (t203 * t229 + (Icges(5,4) * t250 + Icges(5,2) * t249) * t180) * t206;
t44 = Icges(5,5) * t345 + (t203 * t230 + (Icges(5,1) * t250 + Icges(5,4) * t249) * t180) * t206;
t225 = t42 * t426 + t427 * t44 + t431 * t86 + t435 * t93;
t57 = (Icges(5,4) * t234 + Icges(5,2) * t235) * t206;
t58 = (Icges(5,1) * t234 + Icges(5,4) * t235) * t206;
t224 = t110 * t434 + t111 * t435 + t426 * t57 + t427 * t58;
t172 = t283 * qJD(2);
t173 = t284 * qJD(2);
t221 = qJD(1) * t177 - t172 * t200 + t173 * t202 + (-t178 * t202 - t179 * t200) * qJD(2);
t219 = ((t394 + t397) * t201 - (t395 + t398) * t203) * t30 * t293;
t174 = t292 * qJD(2);
t151 = -rSges(3,1) * t319 + (rSges(3,1) * t355 + t406) * qJD(1) + t447 * rSges(3,2);
t150 = -t183 * t342 + (-t201 * t292 + t388) * qJD(1);
t144 = -t201 * t252 + t347;
t143 = -t203 * t252 - t445;
t124 = t158 * t201 - t203 * t276;
t123 = t157 * t201 - t451;
t122 = -t158 * t203 - t452;
t121 = -t157 * t203 - t201 * t278;
t112 = (rSges(5,1) * t425 + rSges(5,2) * t427) * t363;
t75 = -rSges(4,1) * t450 - rSges(4,2) * t107;
t74 = -Icges(4,1) * t450 - Icges(4,4) * t107;
t73 = -Icges(4,4) * t450 - Icges(4,2) * t107;
t72 = -Icges(4,5) * t450 - Icges(4,6) * t107;
t69 = -Icges(4,1) * t100 + Icges(4,4) * t101 + Icges(4,5) * t201;
t68 = -Icges(4,1) * t98 + Icges(4,4) * t99 - Icges(4,5) * t203;
t67 = -Icges(4,4) * t100 + Icges(4,2) * t101 + Icges(4,6) * t201;
t66 = -Icges(4,4) * t98 + Icges(4,2) * t99 - Icges(4,6) * t203;
t65 = -Icges(4,5) * t100 + Icges(4,6) * t101 + Icges(4,3) * t201;
t64 = -Icges(4,5) * t98 + Icges(4,6) * t99 - Icges(4,3) * t203;
t63 = pkin(2) * t355 + t71;
t62 = -pkin(2) * t358 - t70;
t59 = (rSges(5,1) * t234 + rSges(5,2) * t235) * t206;
t56 = (Icges(5,5) * t234 + Icges(5,6) * t235) * t206;
t51 = -t203 * t242 + t389;
t50 = t201 * t88 + t233;
t49 = -t201 * t242 - t203 * t89;
t48 = t201 * t243 - t385;
t47 = (t261 * t180 + (t245 + t236) * qJD(2)) * t208;
t46 = (t262 * t180 + (-t246 + t449) * qJD(2)) * t208;
t41 = Icges(5,3) * t346 + (t201 * t228 + (-Icges(5,5) * t248 - Icges(5,6) * t247) * t180) * t206;
t40 = Icges(5,3) * t345 + (t203 * t228 + (Icges(5,5) * t250 + Icges(5,6) * t249) * t180) * t206;
t39 = (-t194 - t415) * qJD(1) + (-t448 + (rSges(5,1) * t248 + rSges(5,2) * t247) * t180) * t206;
t38 = -pkin(1) * t346 + t255;
t28 = t428 * t203;
t27 = t428 * t201;
t26 = -t203 * t412 - t28 * t75;
t25 = -pkin(2) * t360 - t27 * t75;
t24 = -t100 * t69 + t101 * t67 + t201 * t65;
t23 = -t100 * t68 + t101 * t66 + t201 * t64;
t22 = -t203 * t65 + t67 * t99 - t69 * t98;
t21 = -t203 * t64 + t66 * t99 - t68 * t98;
t20 = rSges(4,1) * t46 + rSges(4,2) * t47;
t19 = Icges(4,1) * t46 + Icges(4,4) * t47;
t18 = Icges(4,4) * t46 + Icges(4,2) * t47;
t17 = Icges(4,5) * t46 + Icges(4,6) * t47;
t16 = Icges(4,1) * t36 + Icges(4,4) * t37 + Icges(4,5) * t346;
t15 = Icges(4,1) * t34 + Icges(4,4) * t35 + Icges(4,5) * t345;
t14 = Icges(4,4) * t36 + Icges(4,2) * t37 + Icges(4,6) * t346;
t13 = Icges(4,4) * t34 + Icges(4,2) * t35 + Icges(4,6) * t345;
t12 = Icges(4,5) * t36 + Icges(4,6) * t37 + Icges(4,3) * t346;
t11 = Icges(4,5) * t34 + Icges(4,6) * t35 + Icges(4,3) * t345;
t10 = (-t202 * t345 + t319) * pkin(2) - t244;
t9 = (-t200 * t342 - t202 * t346) * pkin(2) + t332;
t6 = (t130 * t295 - 0.2e1 * t335 * t79) / t118 ^ 2 * (-t136 * t404 + t373 * t87) + (-t79 * t82 - ((0.4e1 * t303 + t351) * t180 + t137 * t267 + (t220 + (-0.2e1 * t82 * t362 + (t364 * t442 + (-t102 * t200 - t357) * t337) * qJD(2)) * pkin(2)) * pkin(1)) * t208 * t140 - t77 * t87) * t130 * t336 + (0.2e1 * t458 * pkin(3) * t309 + ((pkin(2) * t210 * t305 - t141 * t366 + t190 * t285) * t180 + t140 * t267 + ((-0.2e1 * t413 * t87 + t326) * t200 + ((t354 + (-t175 + t368) * t200) * t180 + (-t104 * t200 - t356) * pkin(2) * t337) * qJD(2)) * pkin(1)) * t208 * t335 + t295 * t404) * t115;
t4 = t201 * t6 + t203 * t322;
t3 = t201 * t322 - t203 * t6;
t2 = pkin(2) * t447 - t20 * t27 - t4 * t75;
t1 = -t20 * t28 + t3 * t75 + (t200 * t346 - t202 * t342) * pkin(2);
t7 = [t46 * t74 - t450 * t19 + t47 * t73 - t107 * t18 + t139 * t58 * t314 + t138 * t57 * t315 + t179 * t343 + t200 * t173 - t178 * t344 + t202 * t172 + (t150 * t164 + t151 * t163) * t443 + 0.2e1 * m(4) * (t10 * t62 + t63 * t9) + 0.2e1 * m(5) * (t38 * t85 + t39 * t84) + (-t206 * t268 + t314 * t86) * t111 + (t206 * t269 + t315 * t83) * t110; -(t109 * t346 - t203 * t56 + (t201 * t227 + (t201 * t224 - t259 * t345) * t180) * t206) * t386 / 0.2e1 + t5 * t251 * t112 * t440 + 0.2e1 * (t1 * t62 + t10 * t26 + t2 * t63 + t25 * t9) * t441 + m(3) * (t275 * t174 + (-qJD(1) * t274 - t150 * t201 + t151 * t203) * t183) + (-t107 * t66 - t203 * t72 - t450 * t68 + t73 * t99 - t74 * t98) * t3 / 0.2e1 + (-t100 * t74 + t101 * t73 - t107 * t67 + t201 * t72 - t450 * t69) * t4 / 0.2e1 + (-t100 * t19 + t101 * t18 - t107 * t13 - t15 * t450 + t17 * t201 + t34 * t74 + t345 * t72 + t35 * t73 + t46 * t69 + t47 * t67) * t27 / 0.2e1 - (-t107 * t14 - t16 * t450 - t17 * t203 + t18 * t99 - t19 * t98 + t346 * t72 + t36 * t74 + t37 * t73 + t46 * t68 + t47 * t66) * t28 / 0.2e1 + (-qJD(2) * t276 + (-qJD(1) * t159 - t203 * t253) * t202 + (-qJD(1) * t161 - t203 * t254) * t200 + t201 * t444 + t221 * t203) * t420 + (-qJD(2) * t278 + (qJD(1) * t160 - t201 * t253) * t202 + (qJD(1) * t162 - t201 * t254) * t200 + t201 * t221 - t203 * t444) * t418 - t455 * t408 / 0.2e1 + t454 * t407 / 0.2e1 + (t109 * t345 + t201 * t56 + ((-t392 + t399) * t293 + t203 * t227 + (t203 * t224 + t259 * t346 + t42 * t427 + t425 * t44 + t431 * t83 + t433 * t93) * t180) * t206) * t390 / 0.2e1 + (-t177 * t203 - t201 * t273 + t279) * t312 + (t201 * t177 - t203 * t273 + t277) * t311 + (((-t393 + t400) * t293 + (t425 * t45 + t427 * t43 + t432 * t83 + t433 * t92) * t180) * t313 + (t59 * t251 + 0.2e1 * (-t201 * t38 - t203 * t39 + (t201 * t84 - t203 * t85) * qJD(1)) * t112) * t440 + t455 * t312 + t454 * t311) * t30; ((t51 * t5 + (t201 * t40 + (t263 * t325 + t50) * qJD(1)) * t30) * t201 + (-t50 * t5 + (-t201 * t41 + (t290 + t51) * qJD(1)) * t30 + (t219 + (-t226 * t203 + (qJD(1) * t264 + t225) * t201) * t391) * t206) * t203) * t390 - ((-t48 * t5 + (t203 * t41 + (t49 - t233) * qJD(1)) * t30) * t203 + (t49 * t5 + (-t203 * t40 + (t290 + t48) * qJD(1)) * t30 + (t219 + (t225 * t201 + (-qJD(1) * t263 - t226) * t203) * t391) * t206) * t201) * t386 + (-t123 * t203 + t124 * t201) * t345 + t201 * ((t201 * t143 + (t123 + t452) * qJD(1)) * t201 + (t124 * qJD(1) + (t159 * t343 + t161 * t344 - t445) * t203 + (-t144 - t277 * qJD(2) + (t158 - t278) * qJD(1)) * t201) * t203) + t4 * (-t23 * t28 + t24 * t27) + t27 * ((-t100 * t15 + t101 * t13 + t11 * t201 + t34 * t69 + t345 * t65 + t35 * t67) * t27 + t24 * t4 - (-t100 * t16 + t101 * t14 + t12 * t201 + t34 * t68 + t345 * t64 + t35 * t66) * t28 + t23 * t3) + t3 * (-t21 * t28 + t22 * t27) - t28 * ((-t11 * t203 + t13 * t99 - t15 * t98 + t346 * t65 + t36 * t69 + t37 * t67) * t27 + t22 * t4 - (-t12 * t203 + t14 * t99 - t16 * t98 + t346 * t64 + t36 * t68 + t37 * t66) * t28 + t21 * t3) + (-t121 * t203 + t122 * t201) * t346 - t203 * ((t203 * t144 + (t122 + t451) * qJD(1)) * t203 + (t121 * qJD(1) + (-t160 * t343 - t162 * t344 + t347) * t201 + (-t143 + t279 * qJD(2) + (-t157 - t276) * qJD(1)) * t203) * t201) + ((t112 ^ 2 * t299 + 0.4e1 * t289 ^ 2) * t5 + (t59 * t112 * t299 + 0.4e1 * t289 * ((qJD(1) * t96 + t255) * t203 + ((-t97 + t194) * qJD(1) + (t448 + ((rSges(5,1) * t435 + rSges(5,2) * t434) * t201 - t265 * t345) * t180) * t206) * t201)) * t30) * t30 * t440 + (t274 * (qJD(1) * t275 + t203 * t150 + t201 * t151) + t349 * t183 * t174) * t443 + 0.4e1 * (t26 * t1 + t25 * t2 + (t27 * t70 + t28 * t71 + t349 * t411) * (t244 * t27 + t28 * t332 - t3 * t71 - t331 * t349 + t4 * t70)) * t441 + (t30 * t345 + t407) * (t201 * t51 - t203 * t50) * t30 + (t30 * t346 - t408) * (t201 * t49 - t203 * t48) * t30;];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t7(1), t7(2); t7(2), t7(3);];
Mq = res;
