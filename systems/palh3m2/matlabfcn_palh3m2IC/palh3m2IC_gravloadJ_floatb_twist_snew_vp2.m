% Calculate Gravitation load with newton euler on the joints for
% palh3m2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 05:00
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh3m2IC_gravloadJ_floatb_twist_snew_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(3,1),zeros(16,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2IC_gravloadJ_floatb_twist_snew_vp2: qJ has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2IC_gravloadJ_floatb_twist_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2IC_gravloadJ_floatb_twist_snew_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2IC_gravloadJ_floatb_twist_snew_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2IC_gravloadJ_floatb_twist_snew_vp2: mrSges has to be [9x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 05:00:17
% EndTime: 2020-05-07 05:00:17
% DurationCPUTime: 0.74s
% Computational Cost: add. (2288->151), mult. (4282->178), div. (14->8), fcn. (4085->39), ass. (0->85)
t373 = m(4) + m(8) + m(9);
t350 = qJ(4) + pkin(14);
t348 = (qJ(3) + t350);
t345 = (pkin(15) + t348);
t357 = (qJ(2) - qJ(6));
t342 = (t345 - t357);
t338 = -2 * qJ(7) - pkin(16) + t342;
t343 = (t345 + t357);
t339 = pkin(16) + t343;
t372 = -cos((qJ(8) - t338)) + cos((qJ(8) - t339));
t368 = cos(qJ(1));
t367 = sin(qJ(8));
t316 = sin(pkin(15));
t366 = pkin(3) * t316;
t317 = cos(pkin(15));
t365 = pkin(3) * t317;
t320 = sin(qJ(5));
t364 = m(6) * t320;
t328 = cos(qJ(5));
t363 = m(6) * t328;
t325 = cos(qJ(8));
t304 = t316 * t325 - t317 * t367;
t362 = m(9) * t304;
t305 = -t316 * t367 - t317 * t325;
t361 = m(9) * t305;
t324 = sin(qJ(1));
t309 = -t368 * g(1) - t324 * g(2);
t323 = sin(qJ(2));
t331 = cos(qJ(2));
t300 = -g(3) * t331 - t309 * t323;
t302 = -g(3) * t323 + t309 * t331;
t322 = sin(qJ(3));
t330 = cos(qJ(3));
t295 = -t300 * t330 + t302 * t322;
t297 = -t300 * t322 - t302 * t330;
t321 = sin(qJ(4));
t329 = cos(qJ(4));
t290 = t295 * t321 + t297 * t329;
t308 = g(1) * t324 - t368 * g(2);
t286 = -t290 * t320 - t308 * t328;
t289 = -t329 * t295 + t297 * t321;
t270 = mrSges(6,2) * t289 - mrSges(6,3) * t286;
t287 = t290 * t328 - t308 * t320;
t271 = -mrSges(6,1) * t289 + mrSges(6,3) * t287;
t346 = -t286 * t364 + t287 * t363;
t255 = pkin(10) * t346 - mrSges(5,2) * t290 + t320 * t270 + t328 * t271 + (-m(6) * pkin(8) - mrSges(5,1)) * t289;
t360 = t255 / pkin(9);
t318 = sin(qJ(7));
t326 = cos(qJ(7));
t294 = t300 * t326 - t302 * t318;
t296 = t300 * t318 + t302 * t326;
t281 = t294 * t305 - t296 * t304;
t282 = t294 * t304 + t296 * t305;
t265 = mrSges(9,1) * t281 - mrSges(9,2) * t282;
t359 = t265 / pkin(7);
t356 = qJ(7) + qJ(8);
t262 = m(5) * t290 + t346;
t273 = (-m(5) - m(6)) * t289;
t355 = t321 * t262 + t329 * t273;
t354 = t281 * t361 + t282 * t362;
t353 = -t281 * t362 + t282 * t361;
t267 = t286 * t363 + t287 * t364;
t351 = qJ(7) + pkin(16);
t349 = -m(5) * t308 + t267;
t347 = t351 + t357;
t344 = mrSges(6,1) * t286 - mrSges(6,2) * t287;
t341 = -qJ(7) + t342;
t340 = -qJ(7) + t343;
t337 = mrSges(8,1) * t294 - mrSges(8,2) * t296 + t353 * t366 + t354 * t365 + t265;
t336 = pkin(4) * t355 + mrSges(4,1) * t295 - mrSges(4,2) * t297 + t255;
t335 = 0.1e1 / pkin(2);
t327 = cos(qJ(6));
t319 = sin(qJ(6));
t313 = sin(t347);
t301 = -g(3) * t319 + t309 * t327;
t299 = -g(3) * t327 - t309 * t319;
t279 = mrSges(9,1) * t308 + mrSges(9,3) * t282;
t278 = -mrSges(9,2) * t308 - mrSges(9,3) * t281;
t260 = mrSges(8,3) * t296 + t278 * t304 + t279 * t305 + (m(9) * t365 + mrSges(8,1)) * t308;
t259 = -mrSges(8,3) * t294 + t278 * t305 - t279 * t304 + (m(9) * t366 - mrSges(8,2)) * t308;
t257 = -pkin(8) * t267 + mrSges(5,1) * t308 + mrSges(5,3) * t290 - t344;
t256 = -pkin(10) * t267 - mrSges(5,2) * t308 + mrSges(5,3) * t289 + t270 * t328 - t271 * t320;
t254 = -mrSges(4,2) * t308 - mrSges(4,3) * t295 + t256 * t329 - t257 * t321;
t253 = -pkin(4) * t349 + mrSges(4,1) * t308 + mrSges(4,3) * t297 + t321 * t256 + t329 * t257;
t1 = [-mrSges(2,2) * t309 + t323 * (-mrSges(3,3) * t300 + t253 * t322 - t254 * t330 + t259 * t326 - t260 * t318) + t331 * (-pkin(1) * t349 + mrSges(3,3) * t302 - t330 * t253 - t322 * t254 + t318 * t259 + t326 * t260) - pkin(12) * t349 + (-t319 * t299 + t327 * t301) * mrSges(7,3) + (mrSges(2,1) - t323 * mrSges(3,2) + t331 * (t373 * pkin(1) + mrSges(3,1)) + pkin(12) * (m(3) + t373) - t319 * mrSges(7,2) + t327 * mrSges(7,1) - m(7) * pkin(6)) * t308; mrSges(3,1) * t300 - mrSges(3,2) * t302 + t336 + t337 + (t318 * (m(8) * t296 + t353) + t326 * (m(8) * t294 + t354) - t322 * (m(4) * t297 + t262 * t329 - t273 * t321) - t330 * (m(4) * t295 + t355)) * pkin(1) + (pkin(3) * (pkin(1) * (cos((-qJ(8) + t357)) - cos((qJ(8) + t357))) + (cos((-qJ(8) + t347)) - cos((qJ(8) + t347))) * pkin(2)) * t360 + (((cos(t341) - cos(t340)) * pkin(1) + (cos(t338) - cos(t339)) * pkin(2)) * pkin(3) + ((-cos((-qJ(8) + t341)) + cos((-qJ(8) + t340))) * pkin(1) + t372 * pkin(2)) * pkin(7)) * t359) / t372 * t335 + ((-pkin(2) * t313 - pkin(1) * sin(t357)) * t335 * t337 + pkin(1) * sin(t351) / pkin(5) * (mrSges(7,1) * t299 - mrSges(7,2) * t301)) / t313; -pkin(9) * t360 + t336 + (-(cos(0.2e1 * qJ(3) + t350) - cos((2 * qJ(7)) - 0.2e1 * pkin(15) + (2 * qJ(8)) + t350)) / (cos((2 * t348)) - cos(0.2e1 * pkin(15) - (2 * t356))) * t360 + sin(t350) / sin((t345 - t356)) * t359) * pkin(4); t344;];
taug = t1(:);
