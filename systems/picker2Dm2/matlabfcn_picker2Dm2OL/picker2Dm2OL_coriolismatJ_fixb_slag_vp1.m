% Calculate matrix of centrifugal and coriolis load on the joints for
% picker2Dm2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% qJD [12x1]
%   Generalized joint velocities
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
% Cq [12x12]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 23:20
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = picker2Dm2OL_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2OL_coriolismatJ_fixb_slag_vp1: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm2OL_coriolismatJ_fixb_slag_vp1: qJD has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2OL_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2OL_coriolismatJ_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm2OL_coriolismatJ_fixb_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'picker2Dm2OL_coriolismatJ_fixb_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 23:18:55
% EndTime: 2020-05-09 23:19:07
% DurationCPUTime: 2.73s
% Computational Cost: add. (21986->205), mult. (10688->267), div. (0->0), fcn. (7272->16), ass. (0->169)
t309 = qJ(1) + qJ(2);
t325 = qJ(4) + t309;
t320 = qJ(10) + t325;
t286 = sin(t320);
t316 = cos(t320);
t252 = -t286 * rSges(11,1) - rSges(11,2) * t316;
t319 = sin(t325);
t238 = pkin(4) * t319 + t252;
t301 = sin(t309);
t291 = pkin(3) * t301;
t220 = t291 + t238;
t306 = sin(qJ(1)) * pkin(1);
t214 = t306 + t220;
t253 = t316 * rSges(11,1) - t286 * rSges(11,2);
t296 = cos(t325);
t239 = pkin(4) * t296 - t253;
t303 = cos(t309);
t417 = pkin(3) * t303;
t221 = t239 + t417;
t385 = cos(qJ(1)) * pkin(1);
t215 = t221 + t385;
t103 = t214 * t239 - t238 * t215;
t105 = t220 * t239 - t238 * t221;
t261 = t319 * rSges(5,1) + t296 * rSges(5,2);
t246 = t291 + t261;
t240 = t306 + t246;
t264 = t296 * rSges(5,1) - t319 * rSges(5,2);
t247 = t264 + t417;
t242 = t247 + t385;
t121 = t240 * t264 - t261 * t242;
t129 = t246 * t264 - t261 * t247;
t433 = m(11) / 0.2e1;
t435 = m(5) / 0.2e1;
t383 = (-t105 + t103) * t433 + (-t129 + t121) * t435;
t384 = (t105 + t103) * t433 + (t129 + t121) * t435;
t2 = t384 - t383;
t458 = t2 * qJD(1);
t305 = qJ(3) + t309;
t321 = qJ(9) + t305;
t287 = cos(t321);
t317 = sin(t321);
t254 = t317 * rSges(10,1) + t287 * rSges(10,2);
t294 = sin(t305);
t244 = -pkin(6) * t294 + t254;
t292 = pkin(2) * t301;
t228 = t244 + t292;
t218 = t228 + t306;
t255 = t287 * rSges(10,1) - t317 * rSges(10,2);
t297 = cos(t305);
t245 = pkin(6) * t297 - t255;
t418 = pkin(2) * t303;
t229 = -t245 + t418;
t219 = t229 + t385;
t104 = -t245 * t218 - t244 * t219;
t262 = -t294 * rSges(4,1) - t297 * rSges(4,2);
t248 = t262 + t292;
t241 = t248 + t306;
t265 = t297 * rSges(4,1) - t294 * rSges(4,2);
t249 = -t265 + t418;
t243 = t249 + t385;
t122 = -t265 * t241 - t262 * t243;
t328 = t265 * t248 + t262 * t249;
t334 = t245 * t228 + t244 * t229;
t434 = m(10) / 0.2e1;
t436 = m(4) / 0.2e1;
t381 = (t334 + t104) * t434 + (t328 + t122) * t436;
t382 = (-t334 + t104) * t434 + (-t328 + t122) * t436;
t6 = t381 - t382;
t457 = t6 * qJD(1);
t331 = t253 * t238 + t252 * t239;
t456 = t331 * m(11) * qJD(4);
t333 = t253 * t220 + t252 * t221;
t455 = t333 * m(11) * qJD(2);
t447 = t253 * t214 + t252 * t215;
t454 = t447 * m(11) * qJD(1);
t304 = qJ(6) + t309;
t293 = sin(t304);
t295 = cos(t304);
t448 = m(7) * (-t385 * (-t293 * rSges(7,1) - t295 * rSges(7,2)) - (t295 * rSges(7,1) - t293 * rSges(7,2)) * t306);
t318 = t448 * qJD(1);
t405 = m(3) * (-t385 * (t301 * rSges(3,1) + t303 * rSges(3,2)) + t306 * (t303 * rSges(3,1) - t301 * rSges(3,2)));
t368 = m(10) * qJD(9);
t327 = m(11) * qJD(10);
t453 = t447 * t327;
t452 = t333 * t327;
t451 = t331 * t327;
t99 = t214 * t221 - t220 * t215;
t100 = t218 * t229 - t228 * t219;
t446 = -t218 * t255 + t254 * t219;
t115 = t240 * t247 - t246 * t242;
t116 = t241 * t249 - t248 * t243;
t308 = qJ(1) + qJ(8);
t300 = sin(t308);
t302 = cos(t308);
t445 = t385 * (-t300 * rSges(9,1) - t302 * rSges(9,2)) + (t302 * rSges(9,1) - t300 * rSges(9,2)) * t306;
t118 = t244 * t255 + t245 * t254;
t114 = t228 * t255 - t254 * t229;
t444 = (qJD(2) + qJD(6)) * t448;
t442 = 0.4e1 * qJD(1);
t440 = 0.4e1 * qJD(2);
t439 = 2 * qJD(3);
t438 = 0.2e1 * qJD(4);
t426 = m(10) * (t114 - t446);
t425 = m(10) * (-t114 - t446);
t423 = m(10) * (t118 - t446);
t422 = m(10) * (-t118 - t446);
t421 = m(10) * (t118 + t114);
t420 = m(10) * (-t118 + t114);
t414 = m(11) * (-t333 - t447);
t413 = m(11) * (t333 - t447);
t411 = m(11) * (-t331 - t447);
t410 = m(11) * (t331 - t447);
t409 = m(11) * (-t331 - t333);
t408 = m(11) * (t331 - t333);
t406 = m(11) * t99;
t403 = m(4) * t116;
t401 = m(4) * t122;
t400 = m(4) * t328;
t398 = m(5) * t115;
t396 = m(5) * t121;
t395 = m(5) * t129;
t391 = m(10) * t100;
t389 = m(10) * t104;
t388 = m(10) * t334;
t379 = m(11) * t103;
t378 = m(11) * t105;
t373 = m(9) * qJD(1);
t372 = m(9) * qJD(8);
t371 = m(10) * qJD(1);
t370 = m(10) * qJD(2);
t369 = m(10) * qJD(3);
t72 = t420 / 0.2e1;
t71 = t421 / 0.2e1;
t66 = t408 / 0.2e1;
t65 = t409 / 0.2e1;
t64 = -t388 - t400;
t62 = t422 / 0.2e1;
t61 = t423 / 0.2e1;
t58 = t378 + t395;
t55 = t410 / 0.2e1;
t54 = t411 / 0.2e1;
t47 = t425 / 0.2e1;
t46 = t426 / 0.2e1;
t43 = t389 + t401;
t42 = t413 / 0.2e1;
t41 = t414 / 0.2e1;
t40 = t379 + t396;
t27 = t72 - t421 / 0.2e1;
t26 = t71 - t420 / 0.2e1;
t25 = t71 + t72;
t22 = t62 - t423 / 0.2e1;
t21 = t61 - t422 / 0.2e1;
t20 = t61 + t62;
t19 = t66 - t409 / 0.2e1;
t18 = t65 - t408 / 0.2e1;
t17 = t65 + t66;
t16 = t55 - t411 / 0.2e1;
t15 = t54 - t410 / 0.2e1;
t14 = t54 + t55;
t13 = t391 + t448 + t398 + t403 + t405 + t406;
t12 = t47 - t426 / 0.2e1;
t11 = t46 - t425 / 0.2e1;
t10 = t46 + t47;
t9 = t42 - t414 / 0.2e1;
t8 = t41 - t413 / 0.2e1;
t7 = t41 + t42;
t5 = t381 + t382;
t1 = t383 + t384;
t3 = [t13 * qJD(2) + t43 * qJD(3) + t40 * qJD(4) + qJD(6) * t448 - t368 * t446 - t372 * t445 - t453, t13 * qJD(1) + t5 * qJD(3) + t1 * qJD(4) + t10 * qJD(9) + t7 * qJD(10) + 0.2e1 * (t405 / 0.2e1 + t100 * t434 + t115 * t435 + t116 * t436 + t99 * t433) * qJD(2) + t444, t43 * qJD(1) + t5 * qJD(2) + t20 * qJD(9) + (t104 * t434 + t122 * t436) * t439, t40 * qJD(1) + t1 * qJD(2) + t14 * qJD(10) + (t103 * t433 + t121 * t435) * t438, 0, t444 + t318, 0, -(t372 + t373) * t445, t10 * qJD(2) + t20 * qJD(3) + (-t368 - t371) * t446, t7 * qJD(2) + t14 * qJD(4) - t453 - t454, 0, 0; -t6 * qJD(3) + t2 * qJD(4) + t11 * qJD(9) + t8 * qJD(10) + (-t398 / 0.4e1 - t403 / 0.4e1 - t405 / 0.4e1 - t406 / 0.4e1 - t391 / 0.4e1) * t442 - t318, t64 * qJD(3) + t58 * qJD(4) + t114 * t368 - t452, -t457 + t64 * qJD(2) + t25 * qJD(9) + (-t328 * t436 - t334 * t434) * t439, t458 + t58 * qJD(2) + t17 * qJD(10) + (t105 * t433 + t129 * t435) * t438, 0, 0, 0, 0, t11 * qJD(1) + t25 * qJD(3) + (t368 + t370) * t114, t8 * qJD(1) + t17 * qJD(4) - t452 - t455, 0, 0; t6 * qJD(2) + t21 * qJD(9) + (-t401 / 0.4e1 - t389 / 0.4e1) * t442, t457 + t26 * qJD(9) + (t388 / 0.4e1 + t400 / 0.4e1) * t440, t118 * t368, 0, 0, 0, 0, 0, t21 * qJD(1) + t26 * qJD(2) + (t368 + t369) * t118, 0, 0, 0; t15 * qJD(10) - t2 * qJD(2) + (-t396 / 0.4e1 - t379 / 0.4e1) * t442, t18 * qJD(10) - t458 + (-t378 / 0.4e1 - t395 / 0.4e1) * t440, 0, -t451, 0, 0, 0, 0, 0, t15 * qJD(1) + t18 * qJD(2) - t451 - t456, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t318, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t445 * t373, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; t12 * qJD(2) + t22 * qJD(3) + t446 * t371, t12 * qJD(1) + t27 * qJD(3) - t114 * t370, t22 * qJD(1) + t27 * qJD(2) - t118 * t369, 0, 0, 0, 0, 0, 0, 0, 0, 0; t9 * qJD(2) + t16 * qJD(4) + t454, t9 * qJD(1) + t19 * qJD(4) + t455, 0, t16 * qJD(1) + t19 * qJD(2) + t456, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
Cq = t3;
