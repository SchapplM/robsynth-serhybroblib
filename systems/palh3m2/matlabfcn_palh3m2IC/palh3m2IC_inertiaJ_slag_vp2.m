% Calculate joint inertia matrix for
% palh3m2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
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
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 05:00
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh3m2IC_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2IC_inertiaJ_slag_vp2: qJ has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2IC_inertiaJ_slag_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2IC_inertiaJ_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2IC_inertiaJ_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2IC_inertiaJ_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 05:00:11
% EndTime: 2020-05-07 05:00:13
% DurationCPUTime: 1.33s
% Computational Cost: add. (3157->275), mult. (4235->388), div. (109->12), fcn. (4294->37), ass. (0->159)
t341 = qJ(4) + pkin(14);
t331 = (qJ(3) + t341);
t324 = (pkin(15) + t331);
t348 = (qJ(2) - qJ(6));
t318 = (t324 - t348);
t313 = -2 * qJ(7) - pkin(16) + t318;
t319 = (t324 + t348);
t314 = pkin(16) + t319;
t382 = -cos((qJ(8) - t313)) + cos((qJ(8) - t314));
t286 = sin(pkin(15));
t287 = cos(pkin(15));
t294 = cos(qJ(8));
t374 = sin(qJ(8));
t245 = t286 * t294 - t287 * t374;
t246 = -t286 * t374 - t287 * t294;
t288 = sin(qJ(7));
t373 = pkin(1) * t288;
t261 = pkin(3) * t286 + t373;
t295 = cos(qJ(7));
t371 = pkin(1) * t295;
t262 = pkin(3) * t287 + t371;
t208 = -t245 * t261 + t246 * t262;
t204 = t208 * mrSges(9,1);
t381 = 0.2e1 * mrSges(8,1) * t371 - 0.2e1 * mrSges(8,2) * t373 + 0.2e1 * t204;
t290 = sin(qJ(5));
t284 = t290 ^ 2;
t364 = mrSges(6,3) * t284;
t275 = pkin(10) * t364;
t297 = cos(qJ(5));
t285 = t297 ^ 2;
t363 = mrSges(6,3) * t285;
t278 = pkin(10) * t363;
t322 = mrSges(6,1) * t297 - mrSges(6,2) * t290;
t367 = pkin(8) * t322;
t380 = t275 + t278 + t367;
t298 = cos(qJ(4));
t368 = pkin(4) * t298;
t273 = -pkin(8) - t368;
t238 = t273 * t322;
t291 = sin(qJ(4));
t272 = t291 * pkin(4) + pkin(10);
t255 = t272 * t364;
t256 = t272 * t363;
t279 = mrSges(5,1) * t368;
t379 = -t238 + t255 + t256 + t279;
t292 = sin(qJ(3));
t293 = sin(qJ(2));
t299 = cos(qJ(3));
t300 = cos(qJ(2));
t251 = t292 * t293 - t299 * t300;
t271 = -pkin(1) * t300 - pkin(12);
t228 = -pkin(4) * t251 + t271;
t378 = 0.2e1 * t228;
t372 = pkin(1) * t292;
t370 = pkin(1) * t299;
t342 = qJ(7) + pkin(16);
t325 = t342 + t348;
t369 = pkin(3) * (pkin(1) * (cos((qJ(8) - t348)) - cos((qJ(8) + t348))) + (cos((qJ(8) - t325)) - cos((qJ(8) + t325))) * pkin(2));
t366 = Ifges(8,3) + Ifges(9,3);
t365 = mrSges(5,2) * t291;
t362 = Ifges(6,4) * t290;
t361 = Ifges(6,4) * t297;
t360 = Ifges(9,3) / pkin(7) ^ 2;
t209 = t245 * t262 + t246 * t261;
t359 = t209 * mrSges(9,2);
t216 = (t245 * t287 + t246 * t286) * pkin(3);
t358 = t216 * mrSges(9,2);
t274 = pkin(4) - t370;
t236 = t291 * t274 - t298 * t372;
t357 = t236 * mrSges(5,2);
t215 = (-t245 * t286 + t246 * t287) * pkin(3);
t213 = t215 * mrSges(9,1);
t356 = Ifges(9,3) + t213;
t315 = -qJ(7) + t319;
t316 = -qJ(7) + t318;
t183 = ((cos(t316) - cos(t315)) * pkin(1) + (cos(t313) - cos(t314)) * pkin(2)) * pkin(3) + ((-cos((qJ(8) - t316)) + cos((qJ(8) - t315))) * pkin(1) + t382 * pkin(2)) * pkin(7);
t304 = 0.1e1 / pkin(7);
t355 = t183 * t304;
t347 = qJ(7) + qJ(8);
t345 = cos((2 * t331)) - cos(0.2e1 * pkin(15) - (2 * t347));
t354 = (t345 * pkin(9) + (cos(0.2e1 * qJ(3) + t341) - cos(-(2 * qJ(7)) + 0.2e1 * pkin(15) - (2 * qJ(8)) - t341)) * pkin(4)) / t345;
t253 = -t292 * t300 - t293 * t299;
t223 = t251 * t291 + t253 * t298;
t353 = t223 * t290;
t352 = t223 * t297;
t224 = 0.1e1 / t382;
t307 = 0.1e1 / pkin(2);
t351 = t224 * t307;
t267 = sin(t325);
t266 = 0.1e1 / t267;
t350 = (-pkin(2) * t267 - pkin(1) * sin(t348)) * t266;
t257 = sin((t324 - t347));
t280 = sin(t341);
t349 = 0.1e1 / t257 * t280;
t250 = -t288 * t293 + t295 * t300;
t252 = t288 * t300 + t293 * t295;
t202 = -t245 * t252 + t246 * t250;
t203 = t245 * t250 + t246 * t252;
t186 = Ifges(9,5) * t203 + Ifges(9,6) * t202;
t222 = -t298 * t251 + t253 * t291;
t346 = Ifges(6,5) * t352 + Ifges(6,3) * t222;
t344 = Ifges(6,5) * t290 + Ifges(6,6) * t297;
t343 = t284 + t285;
t339 = pkin(4) * t365;
t338 = mrSges(4,1) * t370;
t259 = Ifges(6,2) * t297 + t362;
t260 = Ifges(6,1) * t290 + t361;
t337 = t297 * t259 + t290 * t260 + Ifges(5,3);
t336 = t224 * t355;
t335 = t183 * t351;
t303 = 0.1e1 / pkin(9);
t334 = t303 * t354;
t333 = t307 * t350;
t332 = t304 * t349;
t330 = t343 * t272;
t329 = t351 * t369;
t328 = Ifges(4,3) + t337;
t326 = Ifges(8,5) * t252 + Ifges(8,6) * t250 + t186;
t190 = Ifges(9,3) + t204 - t359;
t323 = t303 * t329;
t321 = mrSges(6,1) * t290 + mrSges(6,2) * t297;
t196 = -mrSges(6,2) * t222 - mrSges(6,3) * t353;
t197 = mrSges(6,1) * t222 - mrSges(6,3) * t352;
t320 = t297 * t196 - t290 * t197;
t235 = t274 * t298 + t291 * t372;
t187 = Ifges(6,6) * t222 + (-Ifges(6,2) * t290 + t361) * t223;
t188 = Ifges(6,5) * t222 + (Ifges(6,1) * t297 - t362) * t223;
t317 = t290 * t188 / 0.2e1 + t297 * t187 / 0.2e1 - t259 * t353 / 0.2e1 + t260 * t352 / 0.2e1 + Ifges(5,5) * t223 + (t344 / 0.2e1 - Ifges(5,6)) * t222;
t312 = Ifges(4,5) * t253 + Ifges(4,6) * t251 + t317;
t233 = -pkin(8) - t235;
t225 = t233 * t322;
t234 = pkin(10) + t236;
t226 = t234 * t364;
t227 = t234 * t363;
t231 = t235 * mrSges(5,1);
t311 = -t225 + t226 + t227 + t231 + t337 - t357;
t179 = m(6) * (pkin(10) * t234 * t343 - pkin(8) * t233) + t311 + t380;
t182 = m(6) * (-pkin(8) * t273 + pkin(10) * t330) + t337 - t339 + t379 + t380;
t276 = mrSges(4,2) * t372;
t309 = Ifges(4,3) + m(6) * (t233 * t273 + t234 * t330) + m(5) * (t235 * t298 + t236 * t291) * pkin(4) + t276 + t311 - t338 + t379;
t296 = cos(qJ(6));
t289 = sin(qJ(6));
t281 = sin(t342);
t229 = -pkin(10) * t321 + t344;
t211 = (-t250 * t287 - t252 * t286) * pkin(3) + t271;
t199 = 0.2e1 * t278 + 0.2e1 * t275 + 0.2e1 * t367 + m(6) * (pkin(10) ^ 2 * t343 + pkin(8) ^ 2) + t337;
t195 = t321 * t223;
t194 = t356 - t358;
t192 = pkin(8) * t222 - pkin(10) * t223 + t228;
t189 = -t229 * t334 - t272 * t321 + t344;
t181 = t229 * t323 - t234 * t321 + t344;
t178 = -t199 * t334 + t182;
t177 = -Ifges(6,6) * t353 + t192 * t322 + t346;
t176 = (Ifges(9,3) * t336 + t194 * t350) * t307 + t190;
t175 = t199 * t323 + t179;
t174 = -pkin(8) * t195 + pkin(10) * t320 + t317;
t173 = -t174 * t334 + t195 * t273 + t320 * t272 + (t186 * t332 + (-t222 * t291 - t223 * t298) * mrSges(5,3)) * pkin(4) + t312;
t172 = Ifges(3,5) * t293 + Ifges(3,6) * t300 + t312 + t320 * t234 + t233 * t195 + (t281 / pkin(5) * t266 * (Ifges(7,5) * t289 + Ifges(7,6) * t296) + (t288 * t250 - t295 * t252) * mrSges(8,3) + (-t292 * t251 + t299 * t253) * mrSges(4,3)) * pkin(1) + (((t202 * t216 - t203 * t215) * mrSges(9,3) + t326) * t350 + (t174 * t303 * t369 + t186 * t355) * t224) * t307 + (t209 * t202 - t208 * t203) * mrSges(9,3) + (-t236 * t222 - t235 * t223) * mrSges(5,3) + t326;
t1 = [(m(4) + m(8)) * t271 ^ 2 + t296 * (Ifges(7,4) * t289 + Ifges(7,2) * t296) + t289 * (Ifges(7,1) * t289 + Ifges(7,4) * t296) + 0.2e1 * pkin(6) * (-mrSges(7,1) * t296 + mrSges(7,2) * t289) - 0.2e1 * pkin(12) * (-mrSges(3,1) * t300 + mrSges(3,2) * t293) + t293 * (Ifges(3,1) * t293 + Ifges(3,4) * t300) + t300 * (Ifges(3,4) * t293 + Ifges(3,2) * t300) + Ifges(8,1) * t252 ^ 2 + Ifges(4,1) * t253 ^ 2 + m(5) * t228 ^ 2 + 0.2e1 * t211 * (-mrSges(9,1) * t202 + mrSges(9,2) * t203) + m(9) * t211 ^ 2 + t203 * (Ifges(9,1) * t203 + Ifges(9,4) * t202) + t202 * (Ifges(9,4) * t203 + Ifges(9,2) * t202) + (mrSges(5,2) * t378 + Ifges(5,1) * t223 - t290 * t187 + t297 * t188 + (-Ifges(6,6) * t290 - (2 * Ifges(5,4))) * t222) * t223 + (mrSges(5,1) * t378 + Ifges(5,2) * t222 + t346) * t222 + Ifges(2,3) + m(7) * pkin(6) ^ 2 + m(3) * pkin(12) ^ 2 + (0.2e1 * Ifges(4,4) * t253 + Ifges(4,2) * t251) * t251 + (0.2e1 * Ifges(8,4) * t252 + Ifges(8,2) * t250) * t250 + 0.2e1 * (-mrSges(4,1) * t251 - mrSges(8,1) * t250 + mrSges(4,2) * t253 + mrSges(8,2) * t252) * t271 + (m(6) * t343 * t192 + 0.2e1 * t290 * t196 + 0.2e1 * t297 * t197) * t192, t172, t173, t177; t172, t381 + (t190 + t176) * t304 * t335 + (t175 + t179) * t323 - 0.2e1 * t338 + t328 + (t281 ^ 2 / pkin(5) ^ 2 / t267 ^ 2 * Ifges(7,3) + m(8) * (t288 ^ 2 + t295 ^ 2) + m(4) * (t292 ^ 2 + t299 ^ 2)) * pkin(1) ^ 2 + (((m(9) * (t215 ^ 2 + t216 ^ 2) - 0.2e1 * t358 + 0.2e1 * t213 + t366) * t350 + t194 * t336) * t307 + (2 * Ifges(8,3)) + 0.2e1 * m(9) * (t208 * t215 + t209 * t216) + 0.2e1 * (-t209 - t216) * mrSges(9,2) + 0.2e1 * t356 + t381) * t333 + m(5) * (t235 ^ 2 + t236 ^ 2) + m(9) * (t208 ^ 2 + t209 ^ 2) + m(6) * (t234 ^ 2 * t343 + t233 ^ 2) + 0.2e1 * t231 - 0.2e1 * t357 - 0.2e1 * t359 + t366 + 0.2e1 * t276 - 0.2e1 * t225 + 0.2e1 * t226 + 0.2e1 * t227 + Ifges(3,3), (t176 * t332 - t365) * pkin(4) + (-t175 * t354 + t182 * t329) * t303 + t309, t181; t173, t309 + (t178 * t329 - t179 * t354) * t303 + (-t365 + (t335 * t360 + (t194 * t333 + t190) * t304) * t349) * pkin(4), 0.2e1 * t256 + 0.2e1 * t255 - 0.2e1 * t339 + 0.2e1 * t279 - 0.2e1 * t238 + m(6) * (t272 ^ 2 * t343 + t273 ^ 2) - (t178 + t182) * t334 + (m(5) * (t291 ^ 2 + t298 ^ 2) + t280 ^ 2 / t257 ^ 2 * t360) * pkin(4) ^ 2 + t328, t189; t177, t181, t189, Ifges(6,3);];
Mq = t1;
