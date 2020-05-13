% Calculate Gravitation load with newton euler on the joints for
% picker2Dm1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:55
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = picker2Dm1IC_gravloadJ_floatb_twist_snew_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(3,1),zeros(8,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1IC_gravloadJ_floatb_twist_snew_vp2: qJ has to be [12x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1IC_gravloadJ_floatb_twist_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1IC_gravloadJ_floatb_twist_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1IC_gravloadJ_floatb_twist_snew_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm1IC_gravloadJ_floatb_twist_snew_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 05:55:12
% EndTime: 2020-05-11 05:55:13
% DurationCPUTime: 0.88s
% Computational Cost: add. (2328->159), mult. (3917->182), div. (37->10), fcn. (3746->64), ass. (0->106)
t287 = 0.1e1 / pkin(5);
t277 = sin(qJ(1));
t286 = cos(qJ(1));
t239 = -g(1) * t277 + g(2) * t286;
t240 = g(1) * t286 + g(2) * t277;
t270 = sin(qJ(8));
t279 = cos(qJ(8));
t223 = -t239 * t279 + t240 * t270;
t225 = -t239 * t270 - t240 * t279;
t304 = mrSges(9,1) * t223 - mrSges(9,2) * t225;
t320 = t287 * t304;
t319 = t286 * pkin(1);
t263 = qJ(1) + qJ(8);
t308 = pkin(8) + qJ(5);
t315 = qJ(3) - qJ(6);
t300 = t308 + t315;
t292 = t300 - t263;
t301 = t308 - t315;
t293 = t301 - t263;
t318 = 0.1e1 / (pkin(6) * (cos(t293) - cos(t292)) + (-cos(qJ(9) - t293) + cos(qJ(9) + t292)) * pkin(2)) * t287;
t306 = -qJ(9) - t315;
t220 = 0.1e1 / ((-cos(qJ(4) - t315) + cos(qJ(4) + t315)) * pkin(6) + (cos(qJ(4) + t306) - cos(qJ(4) - t306)) * pkin(2));
t289 = 0.1e1 / pkin(3);
t317 = t220 * t289;
t316 = qJ(2) + qJ(4);
t276 = sin(qJ(2));
t285 = cos(qJ(2));
t224 = t239 * t285 - t240 * t276;
t226 = t239 * t276 + t240 * t285;
t274 = sin(qJ(4));
t283 = cos(qJ(4));
t214 = t224 * t283 - t226 * t274;
t217 = t224 * t274 + t226 * t283;
t265 = sin(qJ(10));
t267 = cos(qJ(10));
t205 = -t214 * t267 + t217 * t265;
t206 = -t214 * t265 - t217 * t267;
t290 = (-t205 * t267 - t206 * t265) * m(11);
t195 = m(5) * t214 + t290;
t196 = m(5) * t217 + (t205 * t265 - t206 * t267) * m(11);
t314 = t283 * t195 + t274 * t196;
t260 = qJ(1) + t316;
t254 = -qJ(9) + t260;
t241 = t254 - t315;
t253 = qJ(9) + t260;
t242 = t253 + t315;
t313 = sin(t242) - sin(t241);
t312 = cos(t242) - cos(t241);
t310 = sin(t254) - sin(t253);
t309 = cos(t254) - cos(t253);
t262 = 0.1e1 / t274;
t288 = 0.1e1 / pkin(4);
t307 = pkin(1) * t262 * t288;
t305 = -qJ(1) + t308;
t272 = sin(qJ(6));
t281 = cos(qJ(6));
t213 = -t224 * t281 + t226 * t272;
t216 = -t224 * t272 - t226 * t281;
t209 = mrSges(7,1) * t213 - mrSges(7,2) * t216;
t275 = sin(qJ(3));
t284 = cos(qJ(3));
t215 = -t224 * t284 + t226 * t275;
t218 = -t224 * t275 - t226 * t284;
t269 = sin(qJ(9));
t278 = cos(qJ(9));
t207 = -t215 * t278 + t218 * t269;
t208 = -t215 * t269 - t218 * t278;
t202 = mrSges(10,1) * t207 - mrSges(10,2) * t208;
t201 = mrSges(11,1) * t205 - mrSges(11,2) * t206;
t303 = -qJ(4) + t305;
t302 = qJ(4) + t305;
t299 = -qJ(8) + t303;
t298 = -qJ(8) + t302;
t291 = (-t207 * t278 - t208 * t269) * m(10);
t197 = m(4) * t215 + t291;
t198 = m(4) * t218 + (t207 * t269 - t208 * t278) * m(10);
t297 = -t197 * t284 - t198 * t275;
t251 = sin(t260);
t252 = cos(t260);
t264 = qJ(1) + qJ(2);
t256 = sin(t264);
t257 = cos(t264);
t296 = t251 * t257 - t252 * t256;
t271 = sin(qJ(7));
t280 = cos(qJ(7));
t295 = t251 * t271 + t252 * t280;
t192 = pkin(6) * t291 + mrSges(4,1) * t215 - mrSges(4,2) * t218 + t202;
t191 = pkin(4) * t290 + mrSges(5,1) * t214 - mrSges(5,2) * t217 + t201;
t188 = pkin(2) * t297 + pkin(3) * t314 + mrSges(3,1) * t224 - mrSges(3,2) * t226 + t191 + t192 + t209;
t282 = cos(qJ(5));
t273 = sin(qJ(5));
t268 = cos(pkin(8));
t266 = sin(pkin(8));
t261 = t277 * pkin(1);
t259 = -qJ(9) + t308;
t258 = qJ(9) + t308;
t255 = sin(t316);
t250 = qJ(9) + t300;
t249 = -qJ(9) + t301;
t234 = -pkin(5) * sin(t263) + t261;
t233 = t319 - pkin(5) * cos(t263);
t232 = -t266 * t273 + t268 * t282;
t231 = t266 * t282 + t268 * t273;
t230 = pkin(3) * t257 + pkin(4) * t252 + t319;
t229 = pkin(3) * t256 + pkin(4) * t251 + t261;
t1 = [-pkin(5) * t320 - t276 * t201 * t307 + mrSges(2,1) * t239 - mrSges(2,2) * t240 + t188 + t270 / sin(-qJ(8) + t305) * (mrSges(6,1) * (g(1) * t231 - g(2) * t232) - mrSges(6,2) * (-g(1) * t232 - g(2) * t231)) + t304 + (((t310 * t229 + t309 * t230) * t317 + (-(sin(t259) - sin(t258)) * t234 - (cos(t259) - cos(t258)) * t233) * t318) * t209 + ((-t313 * t229 - t312 * t230) * t317 + ((sin(t250) - sin(t249)) * t234 + (cos(t250) - cos(t249)) * t233) * t318) * t192) * pkin(2) + ((-cos(qJ(8) + 0.2e1 * qJ(1)) + cos(0.2e1 * pkin(8) + 0.2e1 * qJ(5) + qJ(8))) / (-cos(0.2e1 * t263) + cos(0.2e1 * t308)) * t320 + (-t223 * t279 - t225 * t270) * m(9) + t276 * (m(3) * t226 - t195 * t274 + t196 * t283 + t197 * t275 - t198 * t284) + t285 * (m(3) * t224 + t297 + t314) + (t276 * (t213 * t272 - t216 * t281) + t285 * (-t213 * t281 - t216 * t272)) * m(7)) * pkin(1) + ((-pkin(1) * t255 - pkin(3) * t274) * t262 * t188 + (pkin(3) * t276 + pkin(4) * t255) * t191 * t307 - pkin(1) * (pkin(3) * (cos(t303) - cos(t302)) + (cos(-qJ(2) + t299) - cos(qJ(2) + t298)) * pkin(5)) * t287 / (cos(t299) - cos(t298)) * t202) * t289; mrSges(8,1) * (-g(1) * t280 - g(2) * t271) - mrSges(8,2) * (-g(1) * t271 + g(2) * t280) + ((t313 * t192 - t310 * t209) * t280 + (-t312 * t192 + t309 * t209) * t271) * t220 * pkin(2) + (((t295 * t191 + t296 * t201) * pkin(4) + (t191 - t201) * pkin(3) * (t256 * t271 + t257 * t280)) * t288 - (t188 + t202) * t295) / t296;];
taug = t1(:);
