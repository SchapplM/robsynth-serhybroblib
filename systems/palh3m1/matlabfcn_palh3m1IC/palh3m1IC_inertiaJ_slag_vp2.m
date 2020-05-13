% Calculate joint inertia matrix for
% palh3m1IC
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
% Datum: 2020-04-20 17:32
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh3m1IC_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1IC_inertiaJ_slag_vp2: qJ has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1IC_inertiaJ_slag_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1IC_inertiaJ_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1IC_inertiaJ_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m1IC_inertiaJ_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:32:04
% EndTime: 2020-04-20 17:32:05
% DurationCPUTime: 1.38s
% Computational Cost: add. (3723->248), mult. (5276->375), div. (128->6), fcn. (5054->24), ass. (0->145)
t287 = pkin(16) + qJ(7) + qJ(2);
t273 = sin(t287);
t275 = cos(t287);
t296 = sin(qJ(6));
t303 = cos(qJ(6));
t228 = 0.1e1 / (-t273 * t303 + t275 * t296) / pkin(5) / pkin(2);
t300 = sin(qJ(2));
t257 = -pkin(1) * t300 - pkin(2) * t273;
t307 = cos(qJ(2));
t347 = pkin(1) * t307;
t258 = pkin(2) * t275 + t347;
t195 = (t257 * t303 + t258 * t296) * pkin(5) * t228;
t293 = sin(pkin(15));
t294 = cos(pkin(15));
t301 = cos(qJ(8));
t352 = sin(qJ(8));
t247 = t293 * t301 - t294 * t352;
t248 = -t293 * t352 - t294 * t301;
t295 = sin(qJ(7));
t351 = pkin(1) * t295;
t264 = pkin(3) * t293 + t351;
t302 = cos(qJ(7));
t349 = pkin(1) * t302;
t265 = pkin(3) * t294 + t349;
t210 = -t247 * t264 + t248 * t265;
t206 = t210 * mrSges(9,1);
t218 = (-t247 * t293 + t248 * t294) * pkin(3);
t216 = t218 * mrSges(9,1);
t335 = Ifges(9,3) + t216;
t211 = t247 * t265 + t248 * t264;
t337 = t211 * mrSges(9,2);
t219 = (t247 * t294 + t248 * t293) * pkin(3);
t360 = t219 * mrSges(9,2);
t362 = Ifges(9,3) + t206 - t337 - t195 * (t335 - t360);
t290 = -qJ(7) + pkin(15);
t286 = -qJ(8) + t290;
t272 = sin(t286);
t241 = -t272 * pkin(7) + pkin(3) * sin(t290);
t274 = cos(t286);
t242 = -t274 * pkin(7) + pkin(3) * cos(t290);
t285 = qJ(3) + qJ(4) + pkin(14);
t269 = sin(t285);
t270 = cos(t285);
t226 = 0.1e1 / (t269 * t274 + t270 * t272) / pkin(9) / pkin(7);
t314 = t226 * t195;
t177 = (-t241 * t270 - t242 * t269) * pkin(9) * t314;
t173 = Ifges(9,3) * t177 + t362;
t299 = sin(qJ(3));
t255 = -pkin(4) * t299 - pkin(9) * t269;
t306 = cos(qJ(3));
t256 = pkin(4) * t306 + pkin(9) * t270;
t188 = (t255 * t270 + t256 * t269) * t226 * pkin(9);
t348 = pkin(1) * t306;
t279 = pkin(4) - t348;
t298 = sin(qJ(4));
t305 = cos(qJ(4));
t350 = pkin(1) * t299;
t237 = t279 * t305 + t298 * t350;
t235 = -pkin(8) - t237;
t238 = t298 * t279 - t305 * t350;
t236 = pkin(10) + t238;
t345 = pkin(4) * t305;
t278 = -pkin(8) - t345;
t281 = mrSges(4,2) * t350;
t304 = cos(qJ(5));
t297 = sin(qJ(5));
t339 = Ifges(6,4) * t297;
t262 = Ifges(6,2) * t304 + t339;
t338 = Ifges(6,4) * t304;
t263 = Ifges(6,1) * t297 + t338;
t325 = t304 * t262 + t297 * t263 + Ifges(5,3);
t323 = Ifges(4,3) + t325;
t346 = pkin(4) * t298;
t277 = pkin(10) + t346;
t291 = t297 ^ 2;
t292 = t304 ^ 2;
t329 = t291 + t292;
t324 = t329 * t277;
t326 = mrSges(4,1) * t348;
t321 = mrSges(6,1) * t304 - mrSges(6,2) * t297;
t240 = t278 * t321;
t341 = mrSges(6,3) * t291;
t259 = t277 * t341;
t340 = mrSges(6,3) * t292;
t260 = t277 * t340;
t284 = mrSges(5,1) * t345;
t355 = -t240 + t259 + t260 + t284;
t227 = t235 * t321;
t229 = t236 * t341;
t230 = t236 * t340;
t233 = t237 * mrSges(5,1);
t356 = -t227 + t229 + t230 + t233;
t361 = m(6) * (t235 * t278 + t236 * t324) + m(5) * (t237 * t305 + t238 * t298) * pkin(4) + t281 + (-t238 - t346) * mrSges(5,2) - t326 + t323 + t355 + t356 + t173 * t188;
t357 = 0.2e1 * mrSges(8,1) * t349 - 0.2e1 * mrSges(8,2) * t351 + 0.2e1 * t206;
t252 = t299 * t300 - t306 * t307;
t276 = -pkin(12) - t347;
t231 = -pkin(4) * t252 + t276;
t354 = 0.2e1 * t231;
t343 = pkin(8) * t321;
t342 = Ifges(8,3) + Ifges(9,3);
t336 = t238 * mrSges(5,2);
t254 = -t299 * t307 - t300 * t306;
t225 = t252 * t298 + t254 * t305;
t334 = t225 * t297;
t333 = t225 * t304;
t251 = -t295 * t300 + t302 * t307;
t253 = t295 * t307 + t300 * t302;
t204 = -t247 * t253 + t248 * t251;
t205 = t247 * t251 + t248 * t253;
t185 = Ifges(9,5) * t205 + Ifges(9,6) * t204;
t224 = -t305 * t252 + t254 * t298;
t331 = Ifges(6,5) * t333 + Ifges(6,3) * t224;
t330 = Ifges(6,5) * t297 + Ifges(6,6) * t304;
t327 = mrSges(5,2) * t346;
t322 = Ifges(8,5) * t253 + Ifges(8,6) * t251 + t185;
t320 = mrSges(6,1) * t297 + mrSges(6,2) * t304;
t198 = -mrSges(6,2) * t224 - mrSges(6,3) * t334;
t199 = mrSges(6,1) * t224 - mrSges(6,3) * t333;
t318 = t198 * t304 - t199 * t297;
t186 = Ifges(6,6) * t224 + (-Ifges(6,2) * t297 + t338) * t225;
t187 = Ifges(6,5) * t224 + (Ifges(6,1) * t304 - t339) * t225;
t316 = t297 * t187 / 0.2e1 + t304 * t186 / 0.2e1 - t262 * t334 / 0.2e1 + t263 * t333 / 0.2e1 + Ifges(5,5) * t225 + (t330 / 0.2e1 - Ifges(5,6)) * t224;
t280 = pkin(10) * t341;
t283 = pkin(10) * t340;
t315 = t280 + t283 + t325 + t343;
t312 = Ifges(4,5) * t254 + Ifges(4,6) * t252 + t316;
t179 = m(6) * (t329 * t236 * pkin(10) - pkin(8) * t235) + t315 - t336 + t356;
t182 = m(6) * (-pkin(8) * t278 + pkin(10) * t324) + t315 - t327 + t355;
t232 = -t320 * pkin(10) + t330;
t215 = (-t251 * t294 - t253 * t293) * pkin(3) + t276;
t201 = m(6) * (t329 * pkin(10) ^ 2 + pkin(8) ^ 2) + 0.2e1 * t343 + 0.2e1 * t283 + 0.2e1 * t280 + t325;
t197 = t320 * t225;
t193 = (-t257 * t275 - t258 * t273) * t228 * pkin(2);
t192 = pkin(8) * t224 - pkin(10) * t225 + t231;
t189 = (t255 * t274 - t256 * t272) * t226 * pkin(7);
t181 = t189 * t232 - t320 * t277 + t330;
t178 = (-t241 * t274 + t242 * t272) * pkin(7) * t314;
t176 = t189 * t201 + t182;
t175 = -Ifges(6,6) * t334 + t321 * t192 + t331;
t174 = t178 * t232 - t320 * t236 + t330;
t172 = t178 * t201 + t179;
t171 = -pkin(8) * t197 + t318 * pkin(10) + t316;
t170 = t171 * t189 + t185 * t188 + t197 * t278 + t318 * t277 + (-t224 * t298 - t225 * t305) * mrSges(5,3) * pkin(4) + t312;
t169 = t312 + ((t251 * t295 - t253 * t302) * mrSges(8,3) + (-t252 * t299 + t254 * t306) * mrSges(4,3)) * pkin(1) + Ifges(3,6) * t307 + t318 * t236 + Ifges(3,5) * t300 + t193 * (Ifges(7,5) * t296 + Ifges(7,6) * t303) + t235 * t197 + t177 * t185 + t178 * t171 - t195 * t322 + (-t195 * (t204 * t219 - t205 * t218) - t210 * t205 + t211 * t204) * mrSges(9,3) + (-t224 * t238 - t225 * t237) * mrSges(5,3) + t322;
t1 = [(mrSges(5,2) * t354 + Ifges(5,1) * t225 - t297 * t186 + t304 * t187 + (-Ifges(6,6) * t297 - (2 * Ifges(5,4))) * t224) * t225 + (mrSges(5,1) * t354 + Ifges(5,2) * t224 + t331) * t224 + (m(8) + m(4)) * t276 ^ 2 - 0.2e1 * pkin(12) * (-mrSges(3,1) * t307 + mrSges(3,2) * t300) + t300 * (Ifges(3,1) * t300 + Ifges(3,4) * t307) + t307 * (Ifges(3,4) * t300 + Ifges(3,2) * t307) + t303 * (Ifges(7,4) * t296 + Ifges(7,2) * t303) + t296 * (Ifges(7,1) * t296 + Ifges(7,4) * t303) + 0.2e1 * pkin(6) * (-mrSges(7,1) * t303 + mrSges(7,2) * t296) + Ifges(8,1) * t253 ^ 2 + Ifges(4,1) * t254 ^ 2 + m(5) * t231 ^ 2 + t205 * (Ifges(9,1) * t205 + Ifges(9,4) * t204) + t204 * (Ifges(9,4) * t205 + Ifges(9,2) * t204) + 0.2e1 * t215 * (-mrSges(9,1) * t204 + mrSges(9,2) * t205) + m(9) * t215 ^ 2 + Ifges(2,3) + m(7) * pkin(6) ^ 2 + m(3) * pkin(12) ^ 2 + (0.2e1 * Ifges(4,4) * t254 + Ifges(4,2) * t252) * t252 + (0.2e1 * Ifges(8,4) * t253 + Ifges(8,2) * t251) * t251 + 0.2e1 * (-mrSges(4,1) * t252 - mrSges(8,1) * t251 + mrSges(4,2) * t254 + mrSges(8,2) * t253) * t276 + (m(6) * t329 * t192 + 0.2e1 * t297 * t198 + 0.2e1 * t304 * t199) * t192, t169, t170, t175; t169, t342 - 0.2e1 * t336 - 0.2e1 * t337 + m(6) * (t329 * t236 ^ 2 + t235 ^ 2) - 0.2e1 * t326 + t323 - 0.2e1 * t227 + 0.2e1 * t229 + m(5) * (t237 ^ 2 + t238 ^ 2) + m(9) * (t210 ^ 2 + t211 ^ 2) + (m(8) * (t295 ^ 2 + t302 ^ 2) + m(4) * (t299 ^ 2 + t306 ^ 2)) * pkin(1) ^ 2 - (0.2e1 * (-t211 - t219) * mrSges(9,2) + 0.2e1 * Ifges(8,3) + 0.2e1 * m(9) * (t210 * t218 + t211 * t219) + 0.2e1 * t335 - (0.2e1 * t216 + m(9) * (t218 ^ 2 + t219 ^ 2) + t342 - 0.2e1 * t360) * t195 + t357) * t195 + (t173 + t362) * t177 + (t179 + t172) * t178 + t193 ^ 2 * Ifges(7,3) + 0.2e1 * t230 + 0.2e1 * t233 + 0.2e1 * t281 + t357 + Ifges(3,3), t172 * t189 + t178 * t182 + t361, t174; t170, t176 * t178 + t189 * t179 + t361, -0.2e1 * t327 + t188 ^ 2 * Ifges(9,3) - 0.2e1 * t240 + 0.2e1 * t259 + 0.2e1 * t260 + 0.2e1 * t284 + (t182 + t176) * t189 + m(6) * (t329 * t277 ^ 2 + t278 ^ 2) + m(5) * (t298 ^ 2 + t305 ^ 2) * pkin(4) ^ 2 + t323, t181; t175, t174, t181, Ifges(6,3);];
Mq = t1;
