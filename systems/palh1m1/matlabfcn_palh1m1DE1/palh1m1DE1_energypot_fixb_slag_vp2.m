% Calculate potential energy for
% palh1m1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-14 19:47
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m1DE1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(23,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE1_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1DE1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE1_energypot_fixb_slag_vp2: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1DE1_energypot_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1DE1_energypot_fixb_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-13 14:47:41
% EndTime: 2020-04-13 14:47:51
% DurationCPUTime: 7.91s
% Computational Cost: add. (159809->201), mult. (240238->270), div. (11448->9), fcn. (152251->46), ass. (0->138)
t371 = -m(6) - m(5);
t370 = pkin(4) * m(11);
t369 = -m(11) - m(4) - m(8);
t319 = pkin(1) ^ 2;
t293 = sin(qJ(2));
t295 = sin(pkin(19));
t299 = cos(qJ(2));
t301 = cos(pkin(19));
t270 = t293 * t301 - t295 * t299;
t347 = pkin(7) * t270;
t339 = -0.2e1 * pkin(1) * t347 + t319;
t332 = pkin(7) ^ 2 + t339;
t336 = pkin(3) ^ 2 - pkin(8) ^ 2;
t252 = t332 - t336;
t265 = pkin(1) * t270 - pkin(7);
t362 = -pkin(8) + pkin(3);
t363 = -pkin(3) - pkin(8);
t250 = sqrt(-((pkin(7) - t362) * (pkin(7) + t362) + t339) * ((pkin(7) - t363) * (pkin(7) + t363) + t339));
t273 = t293 * t295 + t299 * t301;
t344 = t250 * t273;
t242 = -pkin(1) * t344 - t252 * t265;
t245 = pkin(1) * t252 * t273 - t250 * t265;
t296 = sin(pkin(18));
t255 = 0.1e1 / t332;
t343 = t255 / pkin(8);
t353 = cos(pkin(18)) / 0.2e1;
t236 = atan2((t245 * t353 + t242 * t296 / 0.2e1) * t343, (t242 * t353 - t296 * t245 / 0.2e1) * t343);
t233 = sin(t236);
t234 = cos(t236);
t368 = -m(3) * pkin(16) + m(7) * pkin(15) + mrSges(3,1) * t293 - mrSges(7,1) * t234 + mrSges(3,2) * t299 + mrSges(7,2) * t233 - mrSges(2,1);
t367 = -m(6) * pkin(12) + mrSges(5,2) - mrSges(6,3);
t291 = sin(qJ(4));
t297 = cos(qJ(4));
t366 = -m(6) * pkin(10) - t297 * mrSges(6,1) + t291 * mrSges(6,2) - mrSges(5,1);
t365 = t291 * mrSges(6,1) + t297 * mrSges(6,2) - mrSges(2,2) + mrSges(11,3) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(7,3) + mrSges(8,3) + mrSges(9,3) + mrSges(10,3);
t364 = -pkin(2) - pkin(13);
t361 = -pkin(9) - pkin(11);
t360 = -pkin(9) + pkin(11);
t359 = pkin(13) - pkin(2);
t253 = t332 + t336;
t264 = pkin(1) - t347;
t243 = -pkin(7) * t344 + t253 * t264;
t358 = -t243 / 0.2e1;
t244 = pkin(7) * t253 * t273 + t250 * t264;
t357 = t244 / 0.2e1;
t356 = sin(pkin(23)) / 0.2e1;
t355 = sin(pkin(21)) / 0.2e1;
t298 = cos(qJ(3));
t354 = -t298 / 0.2e1;
t352 = 0.1e1 / pkin(2) / 0.2e1;
t351 = pkin(1) * t293;
t292 = sin(qJ(3));
t342 = t255 / pkin(3);
t240 = (t244 * t354 + t292 * t358) * t342;
t241 = (t243 * t354 + t292 * t357) * t342;
t282 = pkin(23) + pkin(22);
t277 = sin(t282);
t278 = cos(t282);
t229 = t240 * t278 + t241 * t277;
t350 = pkin(5) * t229;
t286 = sin(pkin(20));
t290 = cos(pkin(20));
t349 = pkin(6) * (-t286 * t298 - t290 * t292);
t348 = pkin(6) * (t286 * t292 - t290 * t298);
t294 = sin(qJ(1));
t279 = t294 * pkin(16);
t300 = cos(qJ(1));
t280 = t300 * pkin(16);
t276 = t299 * pkin(1) + pkin(14);
t341 = -0.2e1 * pkin(4) * t350 + pkin(5) ^ 2;
t210 = sqrt(-((pkin(4) - t360) * (pkin(4) + t360) + t341) * ((pkin(4) - t361) * (pkin(4) + t361) + t341));
t230 = -t240 * t277 + t241 * t278;
t346 = t210 * t230;
t333 = pkin(4) ^ 2 + t341;
t225 = 0.1e1 / t333;
t345 = t225 / pkin(11);
t335 = pkin(1) * t349;
t263 = -0.2e1 * t335;
t312 = pkin(6) ^ 2;
t340 = t263 + t312;
t338 = pkin(9) ^ 2 - pkin(11) ^ 2;
t337 = t312 + t319;
t331 = -pkin(13) ^ 2 + t337;
t330 = t225 / pkin(9) / 0.2e1;
t329 = 0.1e1 / (t263 + t337) * t352;
t328 = 0.1e1 / pkin(13) * t352;
t274 = -t294 * t351 + t279;
t275 = -t300 * t351 + t280;
t287 = cos(pkin(23));
t235 = atan2((t243 * t356 + t287 * t357) * t342, (t244 * t356 + t287 * t358) * t342);
t231 = sin(t235);
t232 = cos(t235);
t215 = -t231 * t299 - t232 * t293;
t327 = t231 * t293 - t232 * t299;
t249 = sqrt(-((pkin(1) - t359) * (pkin(1) + t359) + t340) * ((pkin(1) - t364) * (pkin(1) + t364) + t340));
t317 = pkin(2) ^ 2;
t251 = t263 + t317 + t331;
t262 = -pkin(1) + t349;
t239 = atan2((-t249 * t262 + t251 * t348) * t329, (-t249 * t348 - t251 * t262) * t329);
t237 = sin(t239);
t238 = cos(t239);
t223 = -t237 * t299 - t238 * t293;
t326 = t237 * t293 - t238 * t299;
t271 = t292 * t299 + t293 * t298;
t272 = -t292 * t293 + t298 * t299;
t221 = t333 + t338;
t226 = -pkin(4) + t350;
t321 = atan2((pkin(5) * t221 * t230 - t210 * t226) * t330, (-pkin(5) * t346 - t221 * t226) * t330);
t320 = sin(t321);
t289 = cos(pkin(21));
t288 = cos(pkin(22));
t284 = sin(pkin(22));
t261 = t271 * t300;
t260 = t272 * t300;
t259 = t271 * t294;
t258 = t272 * t294;
t248 = atan2(t249 * t328, (t317 - t331 + 0.2e1 * t335) * t328);
t247 = cos(t248);
t246 = sin(t248);
t227 = -pkin(4) * t229 + pkin(5);
t222 = t333 - t338;
t220 = t223 * t300;
t219 = t326 * t300;
t218 = t223 * t294;
t217 = t326 * t294;
t214 = t215 * t300;
t213 = t327 * t300;
t212 = t215 * t294;
t211 = t327 * t294;
t209 = pkin(4) * t222 * t230 + t210 * t227;
t208 = -pkin(4) * t346 + t222 * t227;
t207 = cos(t321);
t205 = atan2((t208 * t355 + t209 * t289 / 0.2e1) * t345, (-t208 * t289 / 0.2e1 + t209 * t355) * t345);
t204 = cos(t205);
t203 = sin(t205);
t202 = -t288 * t207 - t284 * t320;
t201 = t207 * t284 - t288 * t320;
t1 = (m(7) * pkin(17) - (t201 * t215 - t202 * t327) * mrSges(11,1) - (t201 * t327 + t202 * t215) * mrSges(11,2) - t233 * mrSges(7,1) - t234 * mrSges(7,2) + mrSges(8,1) * t327 - mrSges(8,2) * t215 - mrSges(4,1) * t271 - mrSges(4,2) * t272 - mrSges(3,1) * t299 + mrSges(3,2) * t293 - (-t215 * t284 - t288 * t327) * t370 - mrSges(1,3) - mrSges(2,3) + t371 * (t271 * pkin(5) + t276) + t369 * t276 + (t246 * mrSges(10,1) + t247 * mrSges(10,2) - mrSges(9,2)) * t223 - (-m(10) * pkin(2) + t247 * mrSges(10,1) - t246 * mrSges(10,2) - mrSges(9,1)) * t326 + t366 * (t203 * t272 + t204 * t271) + t367 * (t203 * t271 - t272 * t204) + (-m(2) - m(9) - m(3) - m(7) - m(10)) * pkin(14)) * g(3) + (-m(10) * (pkin(2) * t218 + t279) - mrSges(9,1) * t218 - mrSges(9,2) * t217 - (t201 * t211 + t202 * t212) * mrSges(11,1) - (-t201 * t212 + t202 * t211) * mrSges(11,2) - mrSges(8,1) * t212 - mrSges(8,2) * t211 - (-t217 * t246 - t218 * t247) * mrSges(10,1) - (-t217 * t247 + t218 * t246) * mrSges(10,2) - (-t211 * t284 + t212 * t288) * t370 - m(9) * t279 - mrSges(4,1) * t258 + mrSges(4,2) * t259 - mrSges(1,2) + t371 * (t258 * pkin(5) + t274) + t369 * t274 + t367 * (t203 * t258 + t259 * t204) + t368 * t294 + t366 * (-t203 * t259 + t204 * t258) + t365 * t300) * g(2) + (-m(10) * (pkin(2) * t220 + t280) - mrSges(8,1) * t214 - mrSges(8,2) * t213 - mrSges(9,1) * t220 - mrSges(9,2) * t219 - (-t219 * t246 - t220 * t247) * mrSges(10,1) - (-t219 * t247 + t220 * t246) * mrSges(10,2) - (-t213 * t284 + t214 * t288) * t370 - m(9) * t280 - (t201 * t213 + t202 * t214) * mrSges(11,1) - (-t201 * t214 + t202 * t213) * mrSges(11,2) - mrSges(4,1) * t260 + mrSges(4,2) * t261 - mrSges(1,1) + t371 * (t260 * pkin(5) + t275) + t369 * t275 + t367 * (t203 * t260 + t261 * t204) + t368 * t300 + t366 * (-t203 * t261 + t204 * t260) - t365 * t294) * g(1);
U = t1;
