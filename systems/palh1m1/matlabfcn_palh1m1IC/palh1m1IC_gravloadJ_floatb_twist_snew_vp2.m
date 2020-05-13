% Calculate Gravitation load with newton euler on the joints for
% palh1m1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 20:03
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh1m1IC_gravloadJ_floatb_twist_snew_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1IC_gravloadJ_floatb_twist_snew_vp2: qJ has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1IC_gravloadJ_floatb_twist_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1IC_gravloadJ_floatb_twist_snew_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1IC_gravloadJ_floatb_twist_snew_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1IC_gravloadJ_floatb_twist_snew_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 20:03:10
% EndTime: 2020-04-15 20:03:11
% DurationCPUTime: 0.78s
% Computational Cost: add. (2698->161), mult. (5034->233), div. (20->8), fcn. (4738->34), ass. (0->109)
t448 = m(4) + m(11) + m(8);
t447 = pkin(2) * m(10);
t445 = cos(qJ(1));
t405 = sin(pkin(19));
t444 = pkin(4) * t405;
t406 = cos(pkin(19));
t443 = pkin(4) * t406;
t416 = sin(qJ(1));
t386 = -t445 * g(1) - t416 * g(2);
t415 = sin(qJ(2));
t424 = cos(qJ(2));
t371 = -t424 * g(3) - t415 * t386;
t373 = t415 * g(3) - t424 * t386;
t414 = sin(qJ(3));
t423 = cos(qJ(3));
t363 = t414 * t371 - t423 * t373;
t366 = t423 * t371 + t414 * t373;
t413 = sin(qJ(4));
t422 = cos(qJ(4));
t354 = t422 * t363 + t413 * t366;
t385 = t416 * g(1) - t445 * g(2);
t412 = sin(qJ(5));
t421 = cos(qJ(5));
t345 = -t412 * t354 - t421 * t385;
t352 = t413 * t363 - t422 * t366;
t328 = mrSges(6,2) * t352 - mrSges(6,3) * t345;
t346 = t421 * t354 - t412 * t385;
t329 = -mrSges(6,1) * t352 + mrSges(6,3) * t346;
t439 = m(6) * t421;
t440 = m(6) * t412;
t430 = -t345 * t440 + t346 * t439;
t311 = pkin(11) * t430 - mrSges(5,2) * t354 + t412 * t328 + t421 * t329 + (-m(6) * pkin(9) - mrSges(5,1)) * t352;
t442 = pkin(8) * t311;
t410 = sin(qJ(7));
t419 = cos(qJ(7));
t362 = -t410 * t371 + t419 * t373;
t365 = t419 * t371 + t410 * t373;
t407 = cos(qJ(10));
t438 = sin(qJ(10));
t376 = t405 * t407 - t406 * t438;
t377 = -t405 * t438 - t406 * t407;
t340 = t377 * t362 - t376 * t365;
t341 = t376 * t362 + t377 * t365;
t323 = mrSges(11,1) * t340 - mrSges(11,2) * t341;
t441 = pkin(10) * t323;
t437 = m(11) * t376;
t436 = m(11) * t377;
t318 = m(5) * t354 + t430;
t332 = (-m(5) - m(6)) * t352;
t435 = t413 * t318 + t422 * t332;
t434 = t340 * t436 + t341 * t437;
t433 = -t340 * t437 + t341 * t436;
t325 = t345 * t439 + t346 * t440;
t402 = -qJ(7) + pkin(19);
t432 = -m(5) * t385 + t325;
t409 = sin(qJ(8));
t418 = cos(qJ(8));
t361 = -t409 * t371 + t418 * t373;
t364 = t418 * t371 + t409 * t373;
t408 = sin(qJ(9));
t417 = cos(qJ(9));
t351 = -t417 * t361 + t408 * t364;
t353 = -t408 * t361 - t417 * t364;
t431 = mrSges(10,1) * t351 - mrSges(10,2) * t353;
t429 = mrSges(6,1) * t345 - mrSges(6,2) * t346;
t427 = -mrSges(9,2) * t364 + (-t351 * t417 - t353 * t408) * t447 + mrSges(9,1) * t361 + t431;
t426 = mrSges(8,1) * t362 - mrSges(8,2) * t365 + t433 * t444 + t434 * t443 + t323;
t425 = pkin(5) * t435 + mrSges(4,1) * t366 - mrSges(4,2) * t363 + t311;
t420 = cos(qJ(6));
t411 = sin(qJ(6));
t404 = qJ(8) + qJ(9);
t403 = qJ(3) + pkin(17);
t400 = cos(t404);
t399 = sin(t404);
t398 = pkin(20) + qJ(7) + qJ(2);
t397 = qJ(3) + qJ(4) + pkin(18);
t396 = -qJ(10) + t402;
t395 = cos(t403);
t394 = sin(t403);
t393 = cos(t398);
t392 = cos(t397);
t391 = sin(t398);
t390 = sin(t397);
t389 = cos(t396);
t388 = sin(t396);
t383 = -t418 * pkin(2) + t400 * pkin(12);
t382 = t409 * pkin(2) - t399 * pkin(12);
t381 = t424 * pkin(1) + pkin(3) * t393;
t380 = t415 * pkin(1) + pkin(3) * t391;
t379 = t423 * pkin(5) + pkin(10) * t392;
t378 = t414 * pkin(5) + pkin(10) * t390;
t375 = t389 * pkin(8) - pkin(4) * cos(t402);
t374 = t388 * pkin(8) - pkin(4) * sin(t402);
t372 = -t411 * g(3) + t420 * t386;
t370 = -t420 * g(3) - t411 * t386;
t368 = 0.1e1 / (-t388 * t390 + t389 * t392) / pkin(10) / pkin(8);
t348 = mrSges(10,1) * t385 + mrSges(10,3) * t353;
t347 = -mrSges(10,2) * t385 - mrSges(10,3) * t351;
t338 = mrSges(11,1) * t385 + mrSges(11,3) * t341;
t337 = -mrSges(11,2) * t385 - mrSges(11,3) * t340;
t320 = -mrSges(9,2) * t385 - mrSges(9,3) * t361 - t417 * t347 + t408 * t348;
t319 = mrSges(9,3) * t364 - t408 * t347 - t417 * t348 + (mrSges(9,1) + t447) * t385;
t316 = mrSges(8,3) * t365 + t376 * t337 + t377 * t338 + (m(11) * t443 + mrSges(8,1)) * t385;
t315 = -mrSges(8,3) * t362 + t377 * t337 - t376 * t338 + (m(11) * t444 - mrSges(8,2)) * t385;
t313 = -pkin(9) * t325 + mrSges(5,1) * t385 + mrSges(5,3) * t354 - t429;
t312 = -pkin(11) * t325 - mrSges(5,2) * t385 + mrSges(5,3) * t352 + t421 * t328 - t412 * t329;
t310 = -mrSges(4,2) * t385 - mrSges(4,3) * t366 + t422 * t312 - t413 * t313;
t309 = -pkin(5) * t432 + mrSges(4,1) * t385 + mrSges(4,3) * t363 + t413 * t312 + t422 * t313;
t1 = [-mrSges(2,2) * t386 + t424 * (-mrSges(3,3) * t373 + t423 * t309 + t414 * t310 + t419 * t315 - t410 * t316 - t409 * t319 + t418 * t320) - t415 * (-pkin(1) * t432 + mrSges(3,3) * t371 + t414 * t309 - t423 * t310 + t410 * t315 + t419 * t316 + t418 * t319 + t409 * t320) - pkin(15) * t432 + (-t411 * t370 + t420 * t372) * mrSges(7,3) + (mrSges(2,1) - t424 * mrSges(3,2) - t415 * (t448 * pkin(1) + mrSges(3,1)) + pkin(15) * (m(3) + m(9) + m(10) + t448) - t411 * mrSges(7,2) + t420 * mrSges(7,1) - m(7) * pkin(14)) * t385; (t410 * (m(8) * t365 + t433) + t419 * (m(8) * t362 + t434) - t423 * (m(4) * t363 + t422 * t318 - t413 * t332) + t414 * (m(4) * t366 + t435)) * pkin(1) + t425 + ((t380 * t393 - t381 * t391) * (mrSges(7,1) * t370 - mrSges(7,2) * t372) * pkin(3) + (t426 + (-(-t374 * t390 + t375 * t392) * t441 - (t374 * t389 - t375 * t388) * t442) * t368) * (t380 * t411 + t381 * t420) * pkin(7)) / (-t391 * t411 - t393 * t420) / pkin(7) / pkin(3) - mrSges(3,2) * t371 + mrSges(3,1) * t373 + t426 + t427; ((-t378 * t392 + t379 * t390) * t441 + (t378 * t388 - t379 * t389) * t442) * t368 + ((-t382 * t394 + t383 * t395) * t431 + (-t394 * t399 - t395 * t400) * t427 * pkin(12)) / (t382 * t400 + t383 * t399) / pkin(12) * pkin(6) + t425; t429;];
taug = t1(:);
