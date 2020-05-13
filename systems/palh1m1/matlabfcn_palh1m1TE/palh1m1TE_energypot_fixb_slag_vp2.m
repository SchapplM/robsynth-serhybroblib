% Calculate potential energy for
% palh1m1TE
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
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m1TE_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(23,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1TE_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1TE_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_energypot_fixb_slag_vp2: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1TE_energypot_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1TE_energypot_fixb_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 20:10:25
% EndTime: 2020-04-12 20:10:32
% DurationCPUTime: 6.31s
% Computational Cost: add. (80051->201), mult. (120442->280), div. (5724->9), fcn. (76231->28), ass. (0->135)
t355 = -m(5) - m(6);
t354 = pkin(4) * m(11);
t251 = sin(pkin(20));
t253 = cos(pkin(20));
t255 = sin(qJ(3));
t260 = cos(qJ(3));
t335 = pkin(6) * (-t251 * t260 - t253 * t255);
t319 = pkin(1) * t335;
t234 = -0.2e1 * t319;
t268 = pkin(6) ^ 2;
t324 = t234 + t268;
t342 = pkin(13) - pkin(2);
t347 = -pkin(2) - pkin(13);
t220 = sqrt(-((pkin(1) - t342) * (pkin(1) + t342) + t324) * ((pkin(1) - t347) * (pkin(1) + t347) + t324));
t269 = pkin(2) ^ 2;
t271 = pkin(1) ^ 2;
t322 = t268 + t271;
t312 = -pkin(13) ^ 2 + t322;
t223 = t234 + t269 + t312;
t233 = -pkin(1) + t335;
t334 = pkin(6) * (t251 * t255 - t253 * t260);
t218 = -t220 * t334 - t223 * t233;
t219 = -t220 * t233 + t223 * t334;
t261 = cos(qJ(2));
t270 = 0.1e1 / pkin(2);
t328 = 0.1e1 / (t234 + t322) * t270;
t256 = sin(qJ(2));
t338 = -t256 / 0.2e1;
t300 = (t261 * t218 / 0.2e1 + t219 * t338) * t328;
t353 = -m(11) - m(8) - m(4);
t258 = sin(pkin(19));
t263 = cos(pkin(19));
t240 = t256 * t263 - t258 * t261;
t333 = pkin(7) * t240;
t323 = -0.2e1 * pkin(1) * t333 + t271;
t313 = pkin(7) ^ 2 + t323;
t321 = -pkin(3) ^ 2 + pkin(8) ^ 2;
t224 = t313 + t321;
t235 = pkin(1) * t240 - pkin(7);
t345 = -pkin(8) + pkin(3);
t346 = -pkin(8) - pkin(3);
t221 = sqrt(-((pkin(7) - t345) * (pkin(7) + t345) + t323) * ((pkin(7) - t346) * (pkin(7) + t346) + t323));
t243 = t256 * t258 + t261 * t263;
t326 = t243 * t221;
t304 = -pkin(1) * t326 - t224 * t235;
t305 = pkin(1) * t224 * t243 - t221 * t235;
t317 = cos(pkin(18)) / 0.2e1;
t226 = 0.1e1 / t313;
t327 = t226 / pkin(8);
t337 = sin(pkin(18));
t212 = (t304 * t317 - t337 * t305 / 0.2e1) * t327;
t213 = (t305 * t317 + t304 * t337 / 0.2e1) * t327;
t352 = -m(3) * pkin(16) + m(7) * pkin(15) + mrSges(3,1) * t256 - mrSges(7,1) * t212 + mrSges(3,2) * t261 + mrSges(7,2) * t213 - mrSges(2,1);
t351 = m(6) * pkin(12) - mrSges(5,2) + mrSges(6,3);
t254 = sin(qJ(4));
t259 = cos(qJ(4));
t350 = -m(6) * pkin(10) - t259 * mrSges(6,1) + t254 * mrSges(6,2) - mrSges(5,1);
t349 = t254 * mrSges(6,1) + t259 * mrSges(6,2) - mrSges(2,2) + mrSges(11,3) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(7,3) + mrSges(8,3) + mrSges(9,3) + mrSges(10,3);
t348 = 0.1e1 / pkin(3);
t344 = -pkin(9) - pkin(11);
t343 = pkin(11) - pkin(9);
t341 = -t220 / 0.2e1;
t340 = t220 / 0.2e1;
t339 = -t269 / 0.2e1 + t312 / 0.2e1 - t319;
t336 = pkin(1) * t256;
t257 = sin(qJ(1));
t247 = t257 * pkin(16);
t262 = cos(qJ(1));
t248 = t262 * pkin(16);
t246 = t261 * pkin(1) + pkin(14);
t332 = cos(pkin(21));
t331 = cos(pkin(23));
t330 = sin(pkin(21));
t329 = sin(pkin(23));
t325 = 0.1e1 / pkin(13) * t270;
t320 = -pkin(9) ^ 2 + pkin(11) ^ 2;
t318 = pkin(23) + pkin(22);
t316 = mrSges(10,1) * t325;
t315 = mrSges(10,2) * t325;
t311 = pkin(1) - t333;
t310 = cos(t318);
t309 = sin(t318);
t244 = -t257 * t336 + t247;
t245 = -t262 * t336 + t248;
t308 = t313 - t321;
t298 = t348 * (-pkin(7) * t326 + t308 * t311);
t296 = -t298 / 0.2e1;
t299 = t348 * (pkin(7) * t243 * t308 + t221 * t311);
t297 = t299 / 0.2e1;
t210 = (t296 * t331 + t297 * t329) * t226;
t211 = (t331 * t297 + t329 * t298 / 0.2e1) * t226;
t202 = t210 * t261 - t211 * t256;
t203 = -t210 * t256 - t211 * t261;
t241 = t255 * t261 + t256 * t260;
t242 = -t255 * t256 + t260 * t261;
t209 = (-t261 * t219 / 0.2e1 + t218 * t338) * t328;
t295 = t255 * t296 - t260 * t299 / 0.2e1;
t294 = t255 * t297 + t260 * t296;
t293 = t226 * (t294 * t309 + t295 * t310);
t292 = t226 * (t294 * t310 - t295 * t309);
t291 = pkin(5) * t293;
t290 = pkin(4) - t291;
t289 = -pkin(4) * t293 + pkin(5);
t288 = -0.2e1 * pkin(4) * t291 + pkin(5) ^ 2;
t287 = pkin(4) ^ 2 + t288;
t286 = 0.1e1 / t287;
t285 = 0.1e1 / pkin(9) * t286;
t284 = 0.1e1 / pkin(11) * t286;
t283 = t287 + t320;
t282 = t287 - t320;
t281 = sqrt(-((pkin(4) - t343) * (pkin(4) + t343) + t288) * ((pkin(4) - t344) * (pkin(4) + t344) + t288));
t280 = t281 * t292;
t279 = (-pkin(5) * t280 + t282 * t290) * t285;
t278 = (-pkin(4) * t280 + t283 * t289) * t284;
t277 = -(pkin(5) * t282 * t292 + t281 * t290) * t285 / 0.2e1;
t276 = (pkin(4) * t283 * t292 + t281 * t289) * t284 / 0.2e1;
t252 = cos(pkin(22));
t250 = sin(pkin(22));
t232 = t241 * t262;
t231 = t242 * t262;
t230 = t241 * t257;
t229 = t242 * t257;
t207 = t262 * t300;
t206 = t262 * t209;
t205 = t257 * t300;
t204 = t257 * t209;
t201 = t202 * t262;
t200 = t203 * t262;
t199 = t202 * t257;
t198 = t203 * t257;
t193 = -t332 * t278 / 0.2e1 + t330 * t276;
t192 = t330 * t278 / 0.2e1 + t332 * t276;
t191 = -t252 * t279 / 0.2e1 + t250 * t277;
t190 = t250 * t279 / 0.2e1 + t252 * t277;
t1 = (m(7) * pkin(17) - mrSges(3,1) * t261 + mrSges(3,2) * t256 - t212 * mrSges(7,2) - t213 * mrSges(7,1) - (t202 * t252 - t203 * t250) * t354 - (t190 * t203 + t191 * t202) * mrSges(11,1) - (-t190 * t202 + t191 * t203) * mrSges(11,2) - mrSges(4,2) * t242 - mrSges(8,1) * t202 - mrSges(8,2) * t203 - mrSges(4,1) * t241 - mrSges(1,3) - mrSges(2,3) + t355 * (t241 * pkin(5) + t246) + t353 * t246 + (-t339 * t315 - t341 * t316 - mrSges(9,2)) * t209 + (-m(10) * pkin(2) - t340 * t315 - t339 * t316 - mrSges(9,1)) * t300 + t350 * (t192 * t242 + t193 * t241) + t351 * (-t192 * t241 + t193 * t242) + (-m(9) - m(7) - m(3) - m(2) - m(10)) * pkin(14)) * g(3) + (-m(10) * (pkin(2) * t204 + t247) - (t204 * t340 - t205 * t339) * t315 - (t204 * t339 - t205 * t341) * t316 - (-t190 * t199 + t191 * t198) * mrSges(11,1) - (-t190 * t198 - t191 * t199) * mrSges(11,2) - (t198 * t252 + t199 * t250) * t354 + mrSges(4,2) * t230 - mrSges(9,1) * t204 + mrSges(9,2) * t205 - mrSges(8,1) * t198 + mrSges(8,2) * t199 - mrSges(4,1) * t229 - mrSges(1,2) - m(9) * t247 + t355 * (t229 * pkin(5) + t244) + t353 * t244 + t351 * (-t192 * t229 - t193 * t230) + t352 * t257 + t350 * (-t192 * t230 + t193 * t229) + t349 * t262) * g(2) + (-m(10) * (pkin(2) * t206 + t248) - (t206 * t340 - t207 * t339) * t315 - (t206 * t339 - t207 * t341) * t316 - (-t190 * t201 + t191 * t200) * mrSges(11,1) - (-t190 * t200 - t191 * t201) * mrSges(11,2) - (t200 * t252 + t201 * t250) * t354 + mrSges(4,2) * t232 - mrSges(9,1) * t206 + mrSges(9,2) * t207 - mrSges(8,1) * t200 + mrSges(8,2) * t201 - mrSges(4,1) * t231 - mrSges(1,1) - m(9) * t248 + t355 * (t231 * pkin(5) + t245) + t353 * t245 + t351 * (-t192 * t231 - t193 * t232) + t352 * t262 + t350 * (-t192 * t232 + t193 * t231) - t349 * t257) * g(1);
U = t1;
