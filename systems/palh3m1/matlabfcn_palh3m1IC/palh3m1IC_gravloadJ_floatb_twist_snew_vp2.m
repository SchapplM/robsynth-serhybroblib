% Calculate Gravitation load with newton euler on the joints for
% palh3m1IC
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
% Datum: 2020-04-20 17:32
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh3m1IC_gravloadJ_floatb_twist_snew_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(3,1),zeros(16,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1IC_gravloadJ_floatb_twist_snew_vp2: qJ has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1IC_gravloadJ_floatb_twist_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1IC_gravloadJ_floatb_twist_snew_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1IC_gravloadJ_floatb_twist_snew_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1IC_gravloadJ_floatb_twist_snew_vp2: mrSges has to be [9x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:32:09
% EndTime: 2020-04-20 17:32:10
% DurationCPUTime: 0.61s
% Computational Cost: add. (2362->127), mult. (4416->187), div. (16->6), fcn. (4180->26), ass. (0->86)
t357 = m(4) + m(8) + m(9);
t355 = cos(qJ(1));
t354 = sin(qJ(8));
t321 = sin(pkin(15));
t353 = pkin(3) * t321;
t322 = cos(pkin(15));
t352 = pkin(3) * t322;
t329 = sin(qJ(1));
t308 = -t355 * g(1) - t329 * g(2);
t328 = sin(qJ(2));
t336 = cos(qJ(2));
t295 = -g(3) * t336 - t328 * t308;
t297 = -t328 * g(3) + t308 * t336;
t327 = sin(qJ(3));
t335 = cos(qJ(3));
t289 = -t295 * t335 + t297 * t327;
t291 = -t295 * t327 - t297 * t335;
t326 = sin(qJ(4));
t334 = cos(qJ(4));
t282 = t289 * t326 + t291 * t334;
t307 = g(1) * t329 - t355 * g(2);
t325 = sin(qJ(5));
t333 = cos(qJ(5));
t278 = -t282 * t325 - t307 * t333;
t281 = -t334 * t289 + t291 * t326;
t262 = mrSges(6,2) * t281 - mrSges(6,3) * t278;
t279 = t282 * t333 - t307 * t325;
t263 = -mrSges(6,1) * t281 + mrSges(6,3) * t279;
t348 = m(6) * t333;
t349 = m(6) * t325;
t341 = -t278 * t349 + t279 * t348;
t247 = pkin(10) * t341 - mrSges(5,2) * t282 + t325 * t262 + t333 * t263 + (-m(6) * pkin(8) - mrSges(5,1)) * t281;
t351 = pkin(7) * t247;
t323 = sin(qJ(7));
t331 = cos(qJ(7));
t288 = t295 * t331 - t297 * t323;
t290 = t295 * t323 + t297 * t331;
t330 = cos(qJ(8));
t300 = t321 * t330 - t322 * t354;
t301 = -t321 * t354 - t322 * t330;
t273 = t288 * t301 - t290 * t300;
t274 = t288 * t300 + t290 * t301;
t257 = mrSges(9,1) * t273 - mrSges(9,2) * t274;
t350 = pkin(9) * t257;
t347 = m(9) * t300;
t346 = m(9) * t301;
t254 = m(5) * t282 + t341;
t265 = (-m(5) - m(6)) * t281;
t345 = t326 * t254 + t334 * t265;
t344 = t273 * t346 + t274 * t347;
t343 = -t273 * t347 + t274 * t346;
t259 = t278 * t348 + t279 * t349;
t320 = -qJ(7) + pkin(15);
t342 = -m(5) * t307 + t259;
t340 = mrSges(6,1) * t278 - mrSges(6,2) * t279;
t338 = mrSges(8,1) * t288 - mrSges(8,2) * t290 + t343 * t353 + t344 * t352 + t257;
t337 = pkin(4) * t345 + mrSges(4,1) * t289 - mrSges(4,2) * t291 + t247;
t332 = cos(qJ(6));
t324 = sin(qJ(6));
t318 = pkin(16) + qJ(7) + qJ(2);
t317 = -qJ(8) + t320;
t316 = qJ(3) + qJ(4) + pkin(14);
t315 = cos(t318);
t314 = cos(t317);
t313 = sin(t318);
t312 = sin(t317);
t310 = cos(t316);
t309 = sin(t316);
t305 = t336 * pkin(1) + pkin(2) * t315;
t304 = -pkin(1) * t328 - pkin(2) * t313;
t303 = pkin(4) * t335 + pkin(9) * t310;
t302 = -pkin(4) * t327 - pkin(9) * t309;
t299 = -t314 * pkin(7) + pkin(3) * cos(t320);
t298 = -t312 * pkin(7) + pkin(3) * sin(t320);
t296 = -g(3) * t324 + t308 * t332;
t294 = -g(3) * t332 - t308 * t324;
t292 = 0.1e1 / (t309 * t314 + t310 * t312) / pkin(9) / pkin(7);
t271 = mrSges(9,1) * t307 + mrSges(9,3) * t274;
t270 = -mrSges(9,2) * t307 - mrSges(9,3) * t273;
t252 = mrSges(8,3) * t290 + t270 * t300 + t271 * t301 + (m(9) * t352 + mrSges(8,1)) * t307;
t251 = -mrSges(8,3) * t288 + t270 * t301 - t271 * t300 + (m(9) * t353 - mrSges(8,2)) * t307;
t249 = -pkin(8) * t259 + mrSges(5,1) * t307 + mrSges(5,3) * t282 - t340;
t248 = -pkin(10) * t259 - mrSges(5,2) * t307 + mrSges(5,3) * t281 + t262 * t333 - t263 * t325;
t246 = -mrSges(4,2) * t307 - mrSges(4,3) * t289 + t248 * t334 - t249 * t326;
t245 = -pkin(4) * t342 + mrSges(4,1) * t307 + mrSges(4,3) * t291 + t326 * t248 + t334 * t249;
t1 = [-mrSges(2,2) * t308 + t328 * (-mrSges(3,3) * t295 + t245 * t327 - t246 * t335 + t251 * t331 - t252 * t323) + t336 * (-pkin(1) * t342 + mrSges(3,3) * t297 - t335 * t245 - t327 * t246 + t323 * t251 + t331 * t252) - pkin(12) * t342 + (-t324 * t294 + t332 * t296) * mrSges(7,3) + (mrSges(2,1) - t328 * mrSges(3,2) + t336 * (t357 * pkin(1) + mrSges(3,1)) + pkin(12) * (m(3) + t357) - t324 * mrSges(7,2) + t332 * mrSges(7,1) - m(7) * pkin(6)) * t307; (-t327 * (m(4) * t291 + t254 * t334 - t326 * t265) - t335 * (m(4) * t289 + t345) + t323 * (m(8) * t290 + t343) + t331 * (m(8) * t288 + t344)) * pkin(1) + t337 + mrSges(3,1) * t295 - mrSges(3,2) * t297 + ((-t304 * t315 - t305 * t313) * (mrSges(7,1) * t294 - mrSges(7,2) * t296) * pkin(2) + (-t338 + ((-t298 * t314 + t299 * t312) * t351 + (-t298 * t310 - t299 * t309) * t350) * t292) * (t304 * t332 + t305 * t324) * pkin(5)) / (-t313 * t332 + t315 * t324) / pkin(5) / pkin(2) + t338; ((t302 * t310 + t303 * t309) * t350 + (t302 * t314 - t303 * t312) * t351) * t292 + t337; t340;];
taug = t1(:);
