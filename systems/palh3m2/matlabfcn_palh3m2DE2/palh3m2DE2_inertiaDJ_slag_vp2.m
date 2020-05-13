% Calculate time derivative of joint inertia matrix for
% palh3m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [9x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh3m2DE2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE2_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_inertiaDJ_slag_vp2: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE2_inertiaDJ_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2DE2_inertiaDJ_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2DE2_inertiaDJ_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:22:15
% EndTime: 2020-05-07 04:22:16
% DurationCPUTime: 1.03s
% Computational Cost: add. (620->196), mult. (790->242), div. (0->0), fcn. (392->80), ass. (0->108)
t355 = pkin(17) + pkin(18);
t390 = (pkin(15) + t355);
t386 = (pkin(16) + t390);
t354 = qJD(2) + qJD(3);
t375 = m(5) + m(6);
t406 = pkin(8) * mrSges(6,2);
t351 = qJD(4) + qJD(2);
t405 = -t351 / 0.2e1;
t380 = 2 * qJ(2);
t404 = sin(t380);
t352 = qJD(4) - qJD(2);
t340 = qJD(3) - t352;
t403 = pkin(4) * t340;
t368 = sin(pkin(15));
t402 = pkin(8) * t368;
t372 = cos(pkin(15));
t401 = pkin(8) * t372;
t345 = pkin(10) * mrSges(6,2) - Ifges(6,6);
t366 = sin(qJ(3));
t367 = sin(qJ(2));
t399 = t367 * t366;
t371 = cos(qJ(2));
t398 = t371 * t366;
t358 = qJ(2) + qJ(4);
t397 = qJ(4) - qJ(2);
t329 = t375 * pkin(4) + mrSges(4,1) + mrSges(9,1);
t364 = mrSges(4,2) + mrSges(9,2);
t370 = cos(qJ(3));
t290 = (t329 * t366 + t364 * t370) * pkin(1) * qJD(3);
t333 = -pkin(4) * t370 + pkin(1);
t293 = pkin(4) * t399 + t333 * t371;
t396 = qJD(2) * t367;
t395 = qJD(2) * t371;
t393 = pkin(15) + pkin(18);
t392 = 2 * pkin(12);
t339 = qJD(3) + t351;
t391 = pkin(4) * t339 / 0.2e1;
t346 = pkin(10) * mrSges(6,1) - Ifges(6,5);
t389 = 2 * t386;
t388 = -qJ(2) + t390;
t387 = qJ(2) + t390;
t365 = sin(qJ(4));
t369 = cos(qJ(4));
t385 = mrSges(6,1) * t369 - t365 * mrSges(6,2);
t311 = mrSges(6,1) * t365 + mrSges(6,2) * t369;
t301 = t370 * t367 + t398;
t384 = -t371 * t370 + t399;
t323 = -qJ(2) + t386;
t322 = qJ(2) + t386;
t321 = -qJ(4) + t386;
t320 = qJ(4) + t386;
t383 = mrSges(6,1) * t402 + t346 * t372;
t316 = -qJ(2) + t321;
t315 = qJ(2) + t321;
t314 = -qJ(2) + t320;
t313 = qJ(2) + t320;
t359 = qJ(2) + qJ(3);
t341 = sin(t359);
t342 = cos(t359);
t347 = qJ(3) + t358;
t348 = qJ(3) - t397;
t382 = (t341 * (Ifges(4,6) + Ifges(9,6)) + (pkin(4) * mrSges(5,3) - Ifges(4,5) - Ifges(9,5)) * t342) * t354 + (cos(t347) * t391 + cos(t348) * t403 / 0.2e1) * mrSges(6,2) + (sin(t347) * t391 - sin(t348) * t403 / 0.2e1) * mrSges(6,1);
t379 = 2 * qJ(4);
t374 = pkin(8) * mrSges(6,1);
t361 = cos(pkin(16));
t360 = sin(pkin(16));
t357 = qJ(3) + t380;
t350 = pkin(8) * m(6) + mrSges(5,1);
t349 = -qJ(2) - pkin(15) + pkin(14);
t344 = -qJ(2) + t393;
t343 = qJ(2) + t393;
t338 = cos(t355);
t337 = sin(t355);
t335 = -pkin(10) * m(6) + mrSges(5,2) - mrSges(6,3);
t334 = 0.2e1 * t359;
t332 = cos(t349);
t331 = sin(t349);
t328 = 0.2e1 * t349;
t327 = qJ(3) + t387;
t326 = -qJ(3) + t388;
t325 = -qJ(4) + t389;
t324 = qJ(4) + t389;
t318 = qJ(3) + t322;
t317 = -qJ(3) + t323;
t310 = -qJ(3) + t316;
t309 = -qJ(3) + t314;
t308 = qJ(3) + t315;
t307 = qJ(3) + t313;
t306 = 2 * t321;
t305 = 2 * t320;
t304 = t385 * qJD(4);
t300 = mrSges(6,2) * t402 + t345 * t372;
t299 = mrSges(6,1) * t401 - t368 * t346;
t298 = mrSges(6,2) * t401 - t368 * t345;
t292 = -pkin(4) * t398 + t333 * t367;
t289 = t354 * t301;
t288 = t354 * t384;
t287 = -t368 * t301 - t372 * t384;
t286 = t301 * t372 - t368 * t384;
t285 = -t333 * t396 + (t301 * qJD(3) + t366 * t395) * pkin(4);
t284 = t333 * t395 + (t384 * qJD(3) + t366 * t396) * pkin(4);
t283 = -t368 * t292 + t293 * t372;
t282 = t292 * t372 + t368 * t293;
t281 = t368 * t288 - t289 * t372;
t280 = -t288 * t372 - t368 * t289;
t279 = -t368 * t284 + t285 * t372;
t278 = t284 * t372 + t368 * t285;
t1 = [(((cos(t308) / 0.2e1 + cos(t309) / 0.2e1) * t340 + (-cos(t307) / 0.2e1 - cos(t310) / 0.2e1) * t339) * mrSges(6,2) + ((-sin(t308) / 0.2e1 + sin(t309) / 0.2e1) * t340 + (-sin(t307) / 0.2e1 + sin(t310) / 0.2e1) * t339) * mrSges(6,1)) * pkin(4) + ((t329 * sin(t357) + t364 * cos(t357)) * (qJD(3) + 0.2e1 * qJD(2)) + ((cos(t315) / 0.2e1 + cos(t314) / 0.2e1) * t352 + (cos(t313) / 0.2e1 + cos(t316) / 0.2e1) * t351) * mrSges(6,2) + ((-sin(t315) / 0.2e1 + sin(t314) / 0.2e1) * t352 + (sin(t313) / 0.2e1 - sin(t316) / 0.2e1) * t351) * mrSges(6,1)) * pkin(1) + ((t329 * t341 + t364 * t342) * t392 + ((-sin(t318) + sin(t317)) * t350 + (cos(t318) - cos(t317)) * t335) * pkin(4) - (pkin(4) ^ 2 * t375 - Ifges(4,1) - Ifges(9,1) + Ifges(4,2) + Ifges(9,2)) * sin(t334) + 0.2e1 * (Ifges(4,4) + Ifges(9,4)) * cos(t334) + ((-cos(t327) - cos(t326)) * mrSges(9,2) + (-sin(t327) + sin(t326)) * mrSges(9,1)) * pkin(3)) * t354 + (-t311 * pkin(8) + (-sin(t305) / 0.4e1 + sin(t306) / 0.4e1 + sin(t379) / 0.2e1) * (Ifges(6,2) - Ifges(6,1)) + (cos(t305) / 0.2e1 + cos(t306) / 0.2e1 - cos(t379)) * Ifges(6,4) + ((cos(t320) + cos(t321)) * mrSges(6,2) + (sin(t320) - sin(t321)) * mrSges(6,1)) * pkin(12) - (t374 - t345) * sin(t324) / 0.2e1 + (t374 + t345) * sin(t325) / 0.2e1 + (-t346 - t406) * cos(t324) / 0.2e1 - (-t346 + t406) * cos(t325) / 0.2e1) * qJD(4) + (0.2e1 * Ifges(7,4) * cos(t328) + 0.2e1 * Ifges(3,4) * cos(t380) - (-Ifges(3,1) + Ifges(3,2)) * t404 + (Ifges(7,2) - Ifges(7,1)) * sin(t328) + (-mrSges(3,1) * t367 - mrSges(3,2) * t371) * t392 + 0.2e1 * (-mrSges(7,1) * t331 + mrSges(7,2) * t332) * pkin(6) + ((sin(t322) - sin(t323)) * t350 + (-t404 * pkin(1) - 0.2e1 * pkin(12) * t367) * (m(8) + m(9) + m(4) + t375) + (-cos(t322) + cos(t323)) * t335 + (-cos(t343) + cos(t344)) * mrSges(8,2) + (sin(t343) - sin(t344)) * mrSges(8,1) + (sin(t387) - sin(t388)) * m(9) * pkin(3)) * pkin(1)) * qJD(2) + t290; (Ifges(3,5) * t371 + Ifges(7,5) * t332 - Ifges(3,6) * t367 + Ifges(7,6) * t331) * qJD(2) + ((-mrSges(4,3) - mrSges(5,3) - mrSges(8,3) - mrSges(9,3)) * t395 + (t352 * cos(-t397) / 0.2e1 + cos(t358) * t405) * mrSges(6,2) + (-t352 * sin(-t397) / 0.2e1 + sin(t358) * t405) * mrSges(6,1)) * pkin(1) + t382; 0.2e1 * t290; t382; t290; 0; -t285 * t385 + ((-(t299 * t361 - t383 * t360) * t365 + (-t298 * t361 + t300 * t360) * t369) * t338 + (-(-t360 * t299 - t383 * t361) * t365 + (t298 * t360 + t300 * t361) * t369) * t337 + (pkin(12) + t293) * t311) * qJD(4); -((t278 * t361 + t279 * t360) * t338 + t337 * (-t360 * t278 + t279 * t361)) * t311 - ((t282 * t361 + t283 * t360) * t338 + t337 * (-t360 * t282 + t283 * t361)) * t304; (((t280 * t361 + t281 * t360) * t338 + t337 * (-t280 * t360 + t281 * t361)) * t311 + ((t286 * t361 + t287 * t360) * t338 + t337 * (-t286 * t360 + t287 * t361)) * t304) * pkin(4); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
