% Calculate joint inertia matrix for
% palh1m2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [11x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:49
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m2IC_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2IC_inertiaJ_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2IC_inertiaJ_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2IC_inertiaJ_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2IC_inertiaJ_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2IC_inertiaJ_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:47:30
% EndTime: 2020-05-02 23:47:32
% DurationCPUTime: 1.80s
% Computational Cost: add. (3363->312), mult. (4451->447), div. (144->13), fcn. (4673->40), ass. (0->180)
t338 = qJ(4) + pkin(18);
t339 = (pkin(19) + qJ(3));
t320 = (t338 + t339);
t346 = (qJ(6) - qJ(2));
t312 = (t320 + t346);
t307 = -2 * qJ(7) - pkin(20) + t312;
t313 = (t320 - t346);
t308 = pkin(20) + t313;
t385 = cos((qJ(10) - t307)) + cos((qJ(10) - t308));
t271 = sin(pkin(19));
t272 = cos(pkin(19));
t273 = cos(qJ(10));
t369 = sin(qJ(10));
t222 = t271 * t273 - t272 * t369;
t223 = -t271 * t369 - t272 * t273;
t276 = sin(qJ(7));
t376 = pkin(1) * t276;
t241 = t271 * pkin(4) + t376;
t284 = cos(qJ(7));
t375 = pkin(1) * t284;
t242 = t272 * pkin(4) + t375;
t171 = -t222 * t241 + t223 * t242;
t170 = t171 * mrSges(11,1);
t384 = 0.2e1 * mrSges(8,1) * t375 - 0.2e1 * mrSges(8,2) * t376 + 0.2e1 * t170;
t278 = sin(qJ(5));
t269 = t278 ^ 2;
t366 = mrSges(6,3) * t269;
t257 = pkin(11) * t366;
t286 = cos(qJ(5));
t270 = t286 ^ 2;
t365 = mrSges(6,3) * t270;
t259 = pkin(11) * t365;
t316 = t286 * mrSges(6,1) - t278 * mrSges(6,2);
t371 = pkin(9) * t316;
t383 = t257 + t259 + t371;
t287 = cos(qJ(4));
t372 = pkin(5) * t287;
t256 = -pkin(9) - t372;
t212 = t256 * t316;
t279 = sin(qJ(4));
t255 = pkin(5) * t279 + pkin(11);
t236 = t255 * t366;
t237 = t255 * t365;
t260 = mrSges(5,1) * t372;
t382 = -t212 + t236 + t237 + t260;
t288 = cos(qJ(3));
t374 = pkin(1) * t288;
t337 = qJ(7) + pkin(20);
t321 = t337 - t346;
t373 = pkin(4) * (pkin(1) * (sin((qJ(10) - t346)) + sin((qJ(10) + t346))) + (sin((qJ(10) - t321)) + sin((qJ(10) + t321))) * pkin(3));
t280 = sin(qJ(3));
t370 = t280 * pkin(1);
t282 = cos(qJ(9));
t368 = mrSges(10,1) * t282;
t274 = sin(qJ(9));
t367 = mrSges(10,2) * t274;
t364 = Ifges(6,4) * t278;
t363 = Ifges(6,4) * t286;
t254 = pkin(5) + t370;
t208 = t254 * t279 - t287 * t374;
t362 = t208 * mrSges(5,2);
t361 = t279 * mrSges(5,2);
t360 = Ifges(11,3) + Ifges(8,3);
t281 = sin(qJ(2));
t251 = t281 * pkin(1) - pkin(15);
t359 = Ifges(11,3) / pkin(8) ^ 2;
t309 = -qJ(7) + t312;
t310 = -qJ(7) + t313;
t149 = (-(cos(t310) + cos(t309)) * pkin(1) + (-cos(t307) - cos(t308)) * pkin(3)) * pkin(4) + ((cos((qJ(10) - t310)) + cos((qJ(10) - t309))) * pkin(1) + t385 * pkin(3)) * pkin(8);
t295 = 0.1e1 / pkin(8);
t358 = t149 * t295;
t172 = t222 * t242 + t223 * t241;
t357 = t172 * mrSges(11,2);
t182 = (t222 * t272 + t223 * t271) * pkin(4);
t356 = t182 * mrSges(11,2);
t289 = cos(qJ(2));
t231 = t289 * t280 + t281 * t288;
t348 = t289 * t288;
t232 = -t281 * t280 + t348;
t195 = t287 * t231 + t232 * t279;
t355 = t195 * t278;
t354 = t195 * t286;
t196 = 0.1e1 / t385;
t298 = 0.1e1 / pkin(3);
t353 = t196 * t298;
t340 = qJ(10) + qJ(7);
t235 = cos((t320 - t340));
t234 = 0.1e1 / t235;
t352 = (-pkin(10) * t235 - pkin(5) * cos((t339 - t340))) * t234;
t327 = -qJ(8) + pkin(17) + qJ(3);
t246 = cos(qJ(9) - t327);
t351 = (-pkin(12) * t246 + pkin(2) * cos(t327)) / pkin(12);
t247 = cos(t321);
t245 = 0.1e1 / t247;
t350 = (-pkin(3) * t247 - pkin(1) * cos(t346)) * t245;
t264 = sin(t338);
t349 = t234 * t264;
t181 = (-t222 * t271 + t223 * t272) * pkin(4);
t179 = t181 * mrSges(11,1);
t347 = Ifges(11,3) + t179;
t228 = -t276 * t281 + t284 * t289;
t229 = -t276 * t289 - t284 * t281;
t168 = -t222 * t228 + t223 * t229;
t169 = t222 * t229 + t223 * t228;
t150 = Ifges(11,5) * t169 + Ifges(11,6) * t168;
t233 = Ifges(10,3) + (t367 - t368) * pkin(2);
t317 = Ifges(10,3) + Ifges(9,3) + (0.2e1 * t367 - 0.2e1 * t368 + m(10) * (t274 ^ 2 + t282 ^ 2) * pkin(2)) * pkin(2);
t334 = pkin(6) / t274 / pkin(2);
t318 = t334 * t351;
t345 = t246 * t317 * t334 + t233 * t318;
t194 = -t279 * t231 + t287 * t232;
t344 = Ifges(6,5) * t354 - Ifges(6,3) * t194;
t275 = sin(qJ(8));
t283 = cos(qJ(8));
t227 = -t275 * t289 - t283 * t281;
t230 = -t275 * t281 + t283 * t289;
t192 = -t282 * t227 + t274 * t230;
t193 = -t274 * t227 - t282 * t230;
t343 = Ifges(10,5) * t193 + Ifges(10,6) * t192;
t342 = Ifges(6,5) * t278 + Ifges(6,6) * t286;
t341 = t269 + t270;
t335 = pkin(5) * t361;
t239 = Ifges(6,2) * t286 + t364;
t240 = Ifges(6,1) * t278 + t363;
t333 = t286 * t239 + t278 * t240 + Ifges(5,3);
t332 = t196 * t358;
t331 = t149 * t353;
t294 = 0.1e1 / pkin(10);
t330 = t294 * t352;
t329 = t298 * t350;
t328 = t295 * t349;
t326 = t341 * t255;
t325 = t353 * t373;
t324 = Ifges(4,3) + t333;
t322 = Ifges(8,5) * t228 + Ifges(8,6) * t229 + t150;
t155 = Ifges(11,3) + t170 - t357;
t319 = t294 * t325;
t315 = mrSges(6,1) * t278 + mrSges(6,2) * t286;
t162 = t194 * mrSges(6,2) - mrSges(6,3) * t355;
t163 = -t194 * mrSges(6,1) - mrSges(6,3) * t354;
t314 = t162 * t286 - t163 * t278;
t204 = -pkin(5) * t348 - pkin(15) + (t280 * pkin(5) + pkin(1)) * t281;
t209 = t254 * t287 + t279 * t374;
t153 = -Ifges(6,6) * t194 + (-Ifges(6,2) * t278 + t363) * t195;
t154 = -Ifges(6,5) * t194 + (Ifges(6,1) * t286 - t364) * t195;
t311 = t278 * t154 / 0.2e1 + t286 * t153 / 0.2e1 - t239 * t355 / 0.2e1 + t240 * t354 / 0.2e1 + Ifges(5,5) * t195 + (-t342 / 0.2e1 + Ifges(5,6)) * t194;
t306 = Ifges(4,5) * t231 + Ifges(4,6) * t232 + t311;
t305 = Ifges(9,5) * t230 + Ifges(9,6) * t227 + t343 + (-t192 * t274 + t193 * t282) * pkin(2) * mrSges(10,3);
t206 = -pkin(9) - t209;
t197 = t206 * t316;
t205 = pkin(11) + t208;
t198 = t205 * t366;
t199 = t205 * t365;
t202 = t209 * mrSges(5,1);
t304 = -t197 + t198 + t199 + t202 + t333 - t362;
t145 = m(6) * (t341 * t205 * pkin(11) - pkin(9) * t206) + t304 + t383;
t148 = m(6) * (-pkin(9) * t256 + pkin(11) * t326) + t333 - t335 + t382 + t383;
t258 = mrSges(4,1) * t370;
t262 = mrSges(4,2) * t374;
t302 = Ifges(4,3) + m(6) * (t205 * t326 + t256 * t206) + m(5) * (t208 * t279 + t209 * t287) * pkin(5) + t258 + t262 + t304 + t345 + t382;
t285 = cos(qJ(6));
t277 = sin(qJ(6));
t263 = sin(t337);
t210 = -t227 * pkin(2) - pkin(15);
t201 = -t315 * pkin(11) + t342;
t177 = (-t228 * t271 - t229 * t272) * pkin(4) + t251;
t165 = 0.2e1 * t259 + 0.2e1 * t257 + 0.2e1 * t371 + m(6) * (t341 * pkin(11) ^ 2 + pkin(9) ^ 2) + t333;
t161 = t315 * t195;
t160 = t347 - t356;
t158 = t201 * t330 - t315 * t255 + t342;
t157 = -pkin(9) * t194 - pkin(11) * t195 + t204;
t147 = -t201 * t319 - t315 * t205 + t342;
t144 = t165 * t330 + t148;
t143 = -Ifges(6,6) * t355 + t316 * t157 + t344;
t142 = -t165 * t319 + t145;
t141 = (Ifges(11,3) * t332 + t160 * t350) * t298 + t155;
t140 = -pkin(9) * t161 + t314 * pkin(11) + t311;
t139 = t140 * t330 + t256 * t161 + t314 * t255 + (t246 * t305 + t343 * t351) * t334 + (t150 * t328 + (t194 * t279 - t195 * t287) * mrSges(5,3)) * pkin(5) + t306;
t138 = t314 * t205 + (t194 * t208 - t195 * t209) * mrSges(5,3) + Ifges(3,5) * t289 - Ifges(3,6) * t281 + t206 * t161 + t306 + t305 + (t263 / pkin(7) * t245 * (Ifges(7,5) * t277 + Ifges(7,6) * t285) + (-t228 * t284 + t229 * t276) * mrSges(8,3) + (-t231 * t280 - t232 * t288) * mrSges(4,3)) * pkin(1) + (((t168 * t182 - t169 * t181) * mrSges(11,3) + t322) * t350 + (-t140 * t294 * t373 + t150 * t358) * t196) * t298 + (t172 * t168 - t171 * t169) * mrSges(11,3) + t322;
t1 = [(Ifges(4,1) * t231 + 0.2e1 * Ifges(4,4) * t232) * t231 + (Ifges(8,1) * t228 + 0.2e1 * Ifges(8,4) * t229) * t228 + t277 * (Ifges(7,1) * t277 + Ifges(7,4) * t285) + t285 * (Ifges(7,4) * t277 + Ifges(7,2) * t285) + 0.2e1 * pkin(14) * (-mrSges(7,1) * t285 + mrSges(7,2) * t277) + (-0.2e1 * mrSges(5,1) * t204 + Ifges(5,2) * t194 - t344) * t194 + m(5) * t204 ^ 2 + 0.2e1 * t210 * (-mrSges(10,1) * t192 + mrSges(10,2) * t193) + m(10) * t210 ^ 2 + t193 * (Ifges(10,1) * t193 + Ifges(10,4) * t192) + t192 * (Ifges(10,4) * t193 + Ifges(10,2) * t192) + 0.2e1 * t177 * (-mrSges(11,1) * t168 + mrSges(11,2) * t169) + m(11) * t177 ^ 2 + t168 * (Ifges(11,4) * t169 + Ifges(11,2) * t168) + t169 * (Ifges(11,1) * t169 + Ifges(11,4) * t168) + (0.2e1 * t204 * mrSges(5,2) + Ifges(5,1) * t195 - t278 * t153 + t286 * t154 + (Ifges(6,6) * t278 + (2 * Ifges(5,4))) * t194) * t195 + (m(6) * t341 * t157 + 0.2e1 * t162 * t278 + 0.2e1 * t163 * t286) * t157 + 0.2e1 * (-mrSges(4,1) * t232 - mrSges(8,1) * t229 + mrSges(4,2) * t231 + mrSges(8,2) * t228) * t251 - 0.2e1 * (mrSges(3,1) * t281 - mrSges(9,1) * t227 + mrSges(3,2) * t289 + mrSges(9,2) * t230) * pkin(15) + Ifges(4,2) * t232 ^ 2 + Ifges(3,1) * t289 ^ 2 + (-0.2e1 * Ifges(3,4) * t289 + Ifges(3,2) * t281) * t281 + (0.2e1 * Ifges(9,4) * t230 + Ifges(9,2) * t227) * t227 + Ifges(9,1) * t230 ^ 2 + Ifges(8,2) * t229 ^ 2 + (m(8) + m(4)) * t251 ^ 2 + (m(9) + m(3)) * pkin(15) ^ 2 + Ifges(2,3) + m(7) * pkin(14) ^ 2, t138, t139, t143; t138, m(5) * (t208 ^ 2 + t209 ^ 2) + m(11) * (t171 ^ 2 + t172 ^ 2) + 0.2e1 * t258 + 0.2e1 * t262 + m(6) * (t341 * t205 ^ 2 + t206 ^ 2) + 0.2e1 * t202 + 0.2e1 * t199 - 0.2e1 * t197 + 0.2e1 * t198 + t324 + t317 + (-t142 - t145) * t319 + (t155 + t141) * t295 * t331 + (m(8) * (t276 ^ 2 + t284 ^ 2) + m(4) * (t280 ^ 2 + t288 ^ 2) + t263 ^ 2 / pkin(7) ^ 2 / t247 ^ 2 * Ifges(7,3)) * pkin(1) ^ 2 + (((m(11) * (t181 ^ 2 + t182 ^ 2) - 0.2e1 * t356 + 0.2e1 * t179 + t360) * t350 + t160 * t332) * t298 + (2 * Ifges(8,3)) + 0.2e1 * m(11) * (t171 * t181 + t172 * t182) + 0.2e1 * (-t172 - t182) * mrSges(11,2) + 0.2e1 * t347 + t384) * t329 + t360 - 0.2e1 * t362 + t384 - 0.2e1 * t357 + Ifges(3,3), (t142 * t352 - t148 * t325) * t294 + (t141 * t328 - t361) * pkin(5) + t302, t147; t139, (-t361 + (t331 * t359 + (t160 * t329 + t155) * t295) * t349) * pkin(5) + (-t144 * t325 + t145 * t352) * t294 + t302, 0.2e1 * t237 + 0.2e1 * t236 + 0.2e1 * t260 - 0.2e1 * t335 - 0.2e1 * t212 + m(6) * (t341 * t255 ^ 2 + t256 ^ 2) + (t144 + t148) * t330 + (m(5) * (t279 ^ 2 + t287 ^ 2) + t264 ^ 2 / t235 ^ 2 * t359) * pkin(5) ^ 2 + (t345 * t246 + (Ifges(10,3) * t351 + t233 * t246) * t318) * t334 + t324, t158; t143, t147, t158, Ifges(6,3);];
Mq = t1;
