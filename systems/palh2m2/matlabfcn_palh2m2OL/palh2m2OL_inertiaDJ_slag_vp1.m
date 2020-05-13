% Calculate time derivative of joint inertia matrix for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m2OL_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'palh2m2OL_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_inertiaDJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2OL_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'palh2m2OL_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'palh2m2OL_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:58:20
% EndTime: 2020-05-03 01:58:42
% DurationCPUTime: 4.99s
% Computational Cost: add. (1900->380), mult. (2149->507), div. (0->0), fcn. (405->96), ass. (0->271)
t180 = qJD(4) + qJD(5);
t150 = qJD(3) + t180;
t116 = (qJD(2) + t150);
t376 = m(7) * rSges(7,2);
t133 = rSges(7,1) * t376 - Icges(7,4);
t378 = m(6) * rSges(6,3);
t192 = qJ(4) + qJ(5);
t171 = qJ(3) + t192;
t145 = qJ(2) + t171;
t86 = sin(t145);
t339 = t86 * (rSges(6,2) * t378 - Icges(6,6));
t90 = cos(t145);
t344 = (rSges(6,1) * t378 - Icges(6,5)) * t90;
t392 = -2 * qJD(6);
t81 = t392 + t116;
t96 = -(2 * qJ(6)) + t145;
t359 = t81 * cos(t96);
t78 = 2 * qJD(6) + t116;
t219 = 2 * qJ(6);
t93 = t219 + t145;
t361 = t78 * cos(t93);
t337 = rSges(7,2) * pkin(3);
t367 = rSges(7,1) * rSges(7,3);
t265 = m(7) * (t337 - t367) + Icges(7,5);
t325 = -qJ(6) + qJ(3);
t271 = qJ(4) + t325;
t139 = qJ(2) + t271;
t104 = qJ(5) + t139;
t75 = cos(t104);
t244 = t75 * t265;
t179 = qJD(5) - qJD(6);
t149 = qJD(4) + t179;
t113 = qJD(3) + t149;
t80 = qJD(2) + t113;
t386 = t80 / 0.2e1;
t204 = rSges(7,2) * rSges(7,3);
t338 = pkin(3) * rSges(7,1);
t369 = (m(7) * (t204 + t338));
t50 = -Icges(7,6) + t369;
t71 = sin(t104);
t409 = -t50 * t80 * t71 / 0.2e1 + t244 * t386;
t225 = rSges(7,2) ^ 2;
t226 = rSges(7,1) ^ 2;
t62 = m(7) * (-t225 + t226) - Icges(7,1) + Icges(7,2);
t388 = -t62 / 0.4e1;
t7 = sin(t96) * t81 * t388;
t427 = (t339 - t344) * t116 + t409 + (t359 / 0.2e1 + t361 / 0.2e1) * t133 + t7;
t182 = qJD(2) + qJD(6);
t151 = qJD(3) + t182;
t188 = qJ(6) + qJ(2);
t173 = qJ(3) + t188;
t426 = t151 * (rSges(7,2) * cos(t173) + rSges(7,1) * sin(t173));
t207 = sin(qJ(4));
t375 = rSges(5,1) * m(5);
t385 = m(6) + m(7);
t76 = t385 * pkin(5) + t375;
t343 = t207 * t76;
t274 = qJD(4) * t343;
t381 = (m(5) * rSges(5,2));
t291 = qJD(4) * t381;
t212 = cos(qJ(4));
t319 = pkin(2) * t212;
t415 = -pkin(2) * t274 - t291 * t319;
t205 = sin(qJ(6));
t210 = cos(qJ(6));
t56 = rSges(7,1) * t205 + rSges(7,2) * t210;
t324 = t56 * m(7) * qJD(6);
t257 = pkin(3) * t324;
t266 = qJD(6) * t133 * cos(t219);
t275 = t62 * qJD(6) * sin(t219);
t422 = -t275 - 0.2e1 * t266 - 0.2e1 * t257;
t425 = 0.2e1 * t415 + t422;
t52 = m(7) * (t337 + t367) - Icges(7,5);
t138 = qJ(4) + t173;
t103 = qJ(5) + t138;
t74 = cos(t103);
t424 = t52 * t74 / 0.2e1;
t423 = (qJD(4) + t151) * (rSges(7,2) * cos(t138) + rSges(7,1) * sin(t138));
t198 = pkin(3) * m(7);
t111 = m(6) * rSges(6,1) + t198;
t379 = m(6) * rSges(6,2);
t118 = m(7) * rSges(7,3) + t379;
t122 = sin(t171);
t127 = cos(t171);
t181 = qJD(3) + qJD(4);
t193 = qJ(3) + qJ(4);
t235 = (-cos(t193) * t381 - t76 * sin(t193)) * t181;
t230 = t235 + (-t111 * t122 - t118 * t127) * t150;
t383 = m(4) * rSges(4,2);
t357 = m(5) + t385;
t47 = t357 * pkin(2) + m(4) * rSges(4,1);
t234 = (-cos(qJ(3)) * t383 - t47 * sin(qJ(3))) * qJD(3);
t419 = (t234 + t230) * pkin(4);
t178 = qJD(5) + qJD(6);
t148 = qJD(4) + t178;
t112 = qJD(3) + t148;
t79 = qJD(2) + t112;
t418 = t79 * t424;
t417 = t62 * sin(t93) * t78 / 0.4e1;
t380 = m(5) * rSges(5,3);
t119 = -t378 - t380;
t195 = qJ(2) + qJ(3);
t160 = sin(t195);
t167 = cos(t195);
t184 = qJD(2) + qJD(3);
t353 = pkin(2) * m(7);
t272 = t353 / 0.2e1;
t382 = m(4) * rSges(4,3);
t183 = qJD(2) - qJD(6);
t152 = qJD(3) + t183;
t174 = qJ(2) + t325;
t273 = -t353 / 0.2e1;
t254 = rSges(7,1) * t273;
t405 = (cos(t174) * rSges(7,2) * t272 + sin(t174) * t254) * t152;
t412 = t405 + t272 * t426 + ((pkin(2) * t119 - rSges(4,1) * t382 + Icges(4,5)) * t167 + t160 * (rSges(4,2) * t382 - Icges(4,6))) * t184;
t390 = pkin(4) * m(7);
t142 = qJ(5) + t271;
t345 = t113 * cos(t142);
t264 = pkin(4) * rSges(7,2) * t345;
t377 = m(7) * rSges(7,1);
t286 = -t377 / 0.2e1;
t346 = t113 * sin(t142);
t303 = pkin(4) * t346;
t391 = m(7) / 0.2e1;
t407 = t264 * t391 + t286 * t303;
t263 = pkin(5) * t286;
t389 = pkin(5) * m(7);
t299 = t389 / 0.2e1;
t406 = (cos(t139) * rSges(7,2) * t299 + sin(t139) * t263) * (qJD(4) + t152);
t190 = qJ(5) + qJ(6);
t156 = sin(t190);
t163 = cos(t190);
t403 = -(rSges(7,1) * t156 + rSges(7,2) * t163) * t178 * t389 / 0.2e1;
t237 = ((m(7) * t204 - Icges(7,6)) * t205 - (m(7) * t367 - Icges(7,5)) * t210) * qJD(6);
t211 = cos(qJ(5));
t317 = qJD(5) * t211;
t284 = pkin(5) * t317;
t255 = t118 * t284;
t206 = sin(qJ(5));
t318 = qJD(5) * t206;
t285 = pkin(5) * t318;
t256 = t111 * t285;
t191 = qJ(5) - qJ(6);
t46 = pkin(5) * t179 * cos(t191) * t376;
t401 = -0.2e1 * t256 - 0.2e1 * t255 + t46;
t153 = qJD(2) + t181;
t246 = t299 * t423 + t406 + t417;
t368 = m(7) * (-t204 + t338);
t51 = Icges(7,6) + t368;
t70 = sin(t103);
t364 = t51 * t70;
t281 = t364 / 0.2e1;
t175 = qJ(2) + t193;
t124 = sin(t175);
t335 = t124 * (rSges(5,2) * t380 - Icges(5,6));
t129 = cos(t175);
t341 = (pkin(5) * t378 + rSges(5,3) * t375 - Icges(5,5)) * t129;
t398 = (t335 - t341) * t153 + t418 + t79 * t281 + t246 + t427;
t10 = t118 * t180 + t324;
t1 = t10 * pkin(2) * t207;
t197 = pkin(5) * qJD(5);
t24 = t205 * t376 - t210 * t377 - t111;
t397 = (t1 + t24 * (t180 * t319 + t197)) * t206;
t396 = -0.2e1 * pkin(4);
t395 = -2 * pkin(1);
t394 = -0.2e1 * pkin(2);
t387 = -t79 / 0.2e1;
t13 = qJD(5) * t118 + t324;
t384 = pkin(5) * t13;
t374 = rSges(3,2) * m(3);
t372 = pkin(4) * t150;
t370 = pkin(5) * t206;
t358 = m(4) + t357;
t168 = qJ(4) + t190;
t141 = qJ(3) + t168;
t84 = sin(t141);
t88 = cos(t141);
t356 = -(t84 * rSges(7,1) + t88 * rSges(7,2)) * t112 * t390 / 0.2e1;
t354 = pkin(1) * m(7);
t169 = qJ(4) + t191;
t121 = sin(t169);
t336 = t121 * t149;
t334 = sin(t188) * t182;
t189 = -qJ(6) + qJ(2);
t333 = sin(t189) * t183;
t332 = sin(t191) * t179;
t330 = cos(t188) * t182;
t329 = cos(t189) * t183;
t328 = t163 * t178;
t327 = t180 * t207;
t326 = t180 * t212;
t223 = 2 * qJ(2);
t194 = qJ(3) + t223;
t323 = 2 * qJ(3) + t223;
t321 = pkin(2) * t149;
t320 = pkin(2) * t180;
t315 = qJD(6) * t206;
t314 = qJD(6) * t207;
t313 = qJD(6) * t212;
t177 = qJD(2) + qJD(3) / 0.2e1;
t202 = qJD(4) / 0.2e1;
t117 = t202 + t177;
t201 = qJD(5) / 0.2e1;
t67 = t201 + t117;
t311 = t67 * t396;
t91 = t201 + t153;
t310 = -0.2e1 * pkin(5) * t91;
t199 = qJD(6) / 0.2e1;
t309 = (t199 + t67) * t390;
t200 = -qJD(6) / 0.2e1;
t308 = (t200 + t67) * t390;
t307 = (t199 + t91) * t389;
t306 = (t200 + t91) * t389;
t305 = t117 * t396;
t146 = t202 + t184;
t82 = t201 + t146;
t300 = t82 * t394;
t297 = t79 * t354;
t296 = t80 * t354;
t295 = (t199 + t82) * t353;
t294 = (t200 + t82) * t353;
t293 = t122 * t372;
t292 = t127 * t372;
t290 = t116 * t395;
t289 = t153 * t395;
t287 = t146 * t394;
t283 = t24 * t327;
t277 = rSges(7,1) * t317;
t276 = rSges(7,2) * t317;
t172 = qJ(4) + t194;
t170 = qJ(4) + t323;
t270 = 2 * qJ(4) + t323;
t158 = sin(t192);
t269 = t158 * t320;
t165 = cos(t192);
t268 = t165 * t320;
t120 = sin(t168);
t125 = cos(t168);
t262 = t356 + (rSges(7,1) * t120 + rSges(7,2) * t125) * t148 * t273;
t251 = 2 * qJ(5) + t270;
t140 = qJ(5) + t270;
t144 = qJ(5) + t170;
t249 = t118 * t268;
t248 = t111 * t269;
t240 = -t156 * t178 - t332;
t35 = cos(t169) * t321 * t376;
t236 = t254 * t336 + t262 + t35 / 0.2e1 + t407;
t233 = (-t13 * t211 + t24 * t318) * pkin(5);
t232 = t46 / 0.2e1 + t263 * t332 + t236 + t403;
t229 = -0.2e1 * t248 - 0.2e1 * t249 + t35 + ((-pkin(2) * t125 * t148 - pkin(5) * t328) * rSges(7,2) + (t240 * pkin(5) + (-t120 * t148 - t336) * pkin(2)) * rSges(7,1)) * m(7) + t401 + t425;
t224 = -0.2e1 * Icges(7,6);
t215 = -Icges(7,1) / 0.2e1;
t214 = cos(qJ(2));
t209 = sin(qJ(2));
t196 = t226 / 0.2e1;
t143 = qJ(5) + t172;
t130 = 2 * t195;
t102 = t223 + t142;
t101 = t223 + t141;
t100 = -qJ(6) + t251;
t99 = qJ(6) + t251;
t98 = -qJ(6) + t140;
t97 = qJ(6) + t140;
t95 = -qJ(6) + t144;
t94 = qJ(6) + t144;
t92 = 2 * t175;
t77 = 0.2e1 * t145;
t58 = 2 * t104;
t57 = 2 * t103;
t9 = pkin(2) * t283;
t2 = [(-sin(t94) * t295 - sin(t95) * t294 - sin(t97) * t307 - sin(t98) * t306 - sin(t101) * t309 - sin(t102) * t308 - t71 * t296 - t70 * t297) * rSges(7,1) + (-t51 * sin(t99) - t52 * cos(t99)) * (t199 + t116) + (-t50 * sin(t100) + t265 * cos(t100)) * (t200 + t116) + (cos(t140) * t310 + cos(t143) * t311 + cos(t144) * t300 + t90 * t290 - t292) * t118 + (sin(t170) * t287 + sin(t172) * t305 + t124 * t289) * t76 + (cos(t194) * t383 + t47 * sin(t194)) * t177 * t396 + (-0.2e1 * (rSges(5,2) * t375 - Icges(5,4)) * cos(t92) - ((rSges(5,1) ^ 2 - rSges(5,2) ^ 2) * m(5) - Icges(5,1) + Icges(5,2) + t385 * pkin(5) ^ 2) * sin(t92)) * t153 + (t79 * sin(t57) + t80 * sin(t58)) * t388 + (cos(t170) * t287 + cos(t172) * t305 + t129 * t289) * t381 + (-cos(t94) * t295 + cos(t95) * t294 - cos(t97) * t307 + cos(t98) * t306 - cos(t101) * t309 + cos(t102) * t308 + t75 * t296 - t74 * t297) * rSges(7,2) - t255 - 0.2e1 * ((t214 * t374 + (t358 * pkin(4) + m(3) * rSges(3,1)) * t209) * pkin(1) + (rSges(3,1) * t374 - Icges(3,4)) * cos(t223)) * qJD(2) + (-((t196 + t225 / 0.2e1 - rSges(7,3) ^ 2 + pkin(3) ^ 2) * m(7) + (rSges(6,1) ^ 2 - rSges(6,2) ^ 2) * m(6) - Icges(6,1) + t215 + Icges(6,2) - Icges(7,2) / 0.2e1 + Icges(7,3)) * sin(t77) - 0.2e1 * (rSges(6,1) * t379 + rSges(7,3) * t198 - Icges(6,4)) * cos(t77)) * t116 + t232 + (((-rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * m(4) + Icges(4,1) - Icges(4,2) - t357 * pkin(2) ^ 2) * sin(t130) - 0.2e1 * (rSges(4,1) * t383 - Icges(4,4)) * cos(t130) + (t47 * t160 + t167 * t383) * t395) * t184 - t256 - t257 - t248 - t249 + t415 + t275 / 0.2e1 + (sin(t140) * t310 + sin(t143) * t311 + sin(t144) * t300 + t86 * t290 - t293) * t111 + t266 + pkin(4) * t234 + pkin(4) * t235 + (cos(t57) * t387 + cos(t58) * t386) * t133 - qJD(2) * ((rSges(3,1) ^ 2 - rSges(3,2) ^ 2) * m(3) - Icges(3,1) + Icges(3,2) + t358 * pkin(4) ^ 2) * sin(t223); 0.2e1 * (t361 / 0.4e1 + t359 / 0.4e1) * t133 + (t424 - (t224 - 0.2e1 * t368) * t70 / 0.4e1) * t79 + ((t330 / 0.2e1 + t329 / 0.2e1) * rSges(7,2) + (t334 / 0.2e1 - t333 / 0.2e1) * rSges(7,1)) * t390 + 0.4e1 * (t339 / 0.4e1 - t344 / 0.4e1) * t116 + t246 + 0.4e1 * (-t341 / 0.4e1 + t335 / 0.4e1) * t153 + (t244 / 0.2e1 - (t224 + (2 * t369)) * t71 / 0.4e1) * t80 + t7 + (Icges(3,5) * t214 - t209 * Icges(3,6) + (-rSges(3,1) * t214 + rSges(3,2) * t209) * rSges(3,3) * m(3) + (t119 - t382) * pkin(4) * t214) * qJD(2) + t412; ((-t112 * t88 + t345) * rSges(7,2) + (-t112 * t84 - t346) * rSges(7,1)) * t390 + t229 + 0.2e1 * t419; t398 + t412; t229 + t356 + t407 + t419; 0.2e1 * (-t10 * t319 - t384 + t9) * t211 + 0.2e1 * t397 + t425; t398; t230 * pkin(4) + (rSges(7,1) * t240 - rSges(7,2) * t328) * t389 + t236 + ((-t111 * t158 - t118 * t165) * t180 + (-t212 * t381 - t343) * qJD(4)) * pkin(2) + t422 + t401; (t9 - 0.2e1 * t384) * t211 + (0.2e1 * t197 * t24 + t1) * t206 + (-t274 + (t180 * t206 * t24 - t10 * t211 - t291) * t212) * pkin(2) + t422; 0.2e1 * t233 + t422; t417 + (t281 + t424) * t79 + t427; (-t269 - t285 - t293) * t111 + (-t268 - t284 - t292) * t118 + t232 + t422; (-t384 + (-t10 * t212 + t283) * pkin(2)) * t211 + t397 + t422; t233 + t422; -0.4e1 * qJD(6) * (Icges(7,4) / 0.2e1 + ((t215 + Icges(7,2) / 0.2e1) * t205 + t133 * t210) * t210 + (((t196 - t225 / 0.2e1) * t205 + t337 / 0.2e1) * t210 + (pkin(3) * t205 / 0.2e1 - rSges(7,2) / 0.2e1) * rSges(7,1)) * m(7)); t364 * t387 - t418 - ((t225 + t226) * m(7) + Icges(7,3)) * t116 * t86 + (pkin(1) * t56 * t392 - pkin(5) * t423 - pkin(2) * t426 + ((t329 - t330) * rSges(7,2) + (-t333 - t334) * rSges(7,1)) * pkin(4)) * t391 + t405 + t406 + t409; -t35 / 0.2e1 - t46 / 0.2e1 - t237 + (-t264 / 0.2e1 + (t303 / 0.2e1 + t121 * t321 / 0.2e1 + pkin(5) * t332 / 0.2e1) * rSges(7,1)) * m(7) + t262 + t403; -t237 + (((-rSges(7,1) * t315 - t276) * t210 - (-rSges(7,2) * t315 + t277) * t205) * pkin(5) + ((-(rSges(7,1) * t314 + rSges(7,2) * t326) * t211 - (rSges(7,1) * t313 - rSges(7,2) * t327) * t206) * t210 - ((rSges(7,1) * t326 - rSges(7,2) * t314) * t211 - (rSges(7,1) * t327 + rSges(7,2) * t313) * t206) * t205) * pkin(2)) * m(7); (-Icges(7,5) * t210 + Icges(7,6) * t205) * qJD(6) + ((-pkin(5) * t276 - qJD(6) * (rSges(7,1) * t370 - t367)) * t210 - (pkin(5) * t277 - qJD(6) * (rSges(7,2) * t370 - t204)) * t205) * m(7); -t237; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11), t2(16); t2(2), t2(3), t2(5), t2(8), t2(12), t2(17); t2(4), t2(5), t2(6), t2(9), t2(13), t2(18); t2(7), t2(8), t2(9), t2(10), t2(14), t2(19); t2(11), t2(12), t2(13), t2(14), t2(15), t2(20); t2(16), t2(17), t2(18), t2(19), t2(20), t2(21);];
Mq = res;
