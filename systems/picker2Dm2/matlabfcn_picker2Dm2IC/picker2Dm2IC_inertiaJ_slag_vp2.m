% Calculate joint inertia matrix for
% picker2Dm2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
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
% Mq [2x2]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 09:21
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = picker2Dm2IC_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2IC_inertiaJ_slag_vp2: qJ has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2IC_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2IC_inertiaJ_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm2IC_inertiaJ_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'picker2Dm2IC_inertiaJ_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 09:21:02
% EndTime: 2020-05-11 09:21:03
% DurationCPUTime: 1.24s
% Computational Cost: add. (5151->286), mult. (4203->331), div. (313->10), fcn. (3408->60), ass. (0->185)
t253 = qJ(1) + qJ(2);
t243 = sin(t253);
t244 = cos(t253);
t258 = sin(qJ(7));
t265 = cos(qJ(7));
t348 = (t243 * t258 + t244 * t265) * pkin(3);
t269 = cos(qJ(2));
t249 = t269 * pkin(1);
t232 = t249 + pkin(3);
t260 = sin(qJ(4));
t267 = cos(qJ(4));
t262 = sin(qJ(2));
t344 = pkin(1) * t262;
t196 = t267 * t232 - t260 * t344;
t191 = pkin(4) + t196;
t199 = t260 * t232 + t267 * t344;
t254 = sin(qJ(10));
t255 = cos(qJ(10));
t169 = -t255 * t191 + t254 * t199;
t165 = t169 * mrSges(11,1);
t188 = t196 * mrSges(5,1);
t347 = Ifges(5,3) + t165 + t188;
t233 = t249 + pkin(2);
t261 = sin(qJ(3));
t268 = cos(qJ(3));
t197 = -t268 * t233 + t261 * t344;
t192 = pkin(6) + t197;
t200 = -t261 * t233 - t268 * t344;
t256 = sin(qJ(9));
t263 = cos(qJ(9));
t171 = -t263 * t192 + t256 * t200;
t166 = t171 * mrSges(10,1);
t189 = t197 * mrSges(4,1);
t332 = t200 * mrSges(4,2);
t346 = t166 + t189 - t332;
t257 = sin(qJ(8));
t345 = 0.2e1 * t257 * pkin(1) * mrSges(9,2);
t274 = 0.1e1 / pkin(4);
t343 = pkin(1) * t274;
t312 = qJ(3) - qJ(6);
t296 = -qJ(9) - t312;
t178 = 0.1e1 / ((-cos(qJ(4) - t312) + cos(qJ(4) + t312)) * pkin(6) + (cos(qJ(4) + t296) - cos(qJ(4) - t296)) * pkin(2));
t342 = pkin(2) * t178;
t341 = pkin(2) * t261;
t340 = t260 * pkin(3);
t339 = t267 * pkin(3);
t338 = t268 * pkin(2);
t337 = cos(qJ(1)) * pkin(1);
t336 = Ifges(4,3) + Ifges(10,3);
t172 = -t256 * t192 - t263 * t200;
t335 = t172 * mrSges(10,2);
t231 = pkin(6) - t338;
t198 = -t256 * t231 + t263 * t341;
t334 = t198 * mrSges(10,2);
t333 = t199 * mrSges(5,2);
t259 = sin(qJ(6));
t266 = cos(qJ(6));
t204 = (-t259 * t269 - t262 * t266) * pkin(1);
t331 = t204 * mrSges(7,2);
t330 = t262 * mrSges(3,2);
t329 = t263 * mrSges(10,1);
t328 = Ifges(11,3) + Ifges(5,3);
t195 = -t263 * t231 - t256 * t341;
t187 = t195 * mrSges(10,1);
t327 = Ifges(10,3) + t187;
t272 = 0.1e1 / pkin(5);
t252 = qJ(1) + qJ(8);
t306 = pkin(8) + qJ(5);
t289 = t306 + t312;
t281 = t289 - t252;
t290 = t306 - t312;
t282 = t290 - t252;
t326 = 0.1e1 / (pkin(6) * (cos(t282) - cos(t281)) + (-cos(qJ(9) - t282) + cos(qJ(9) + t281)) * pkin(2)) * t272;
t170 = -t254 * t191 - t255 * t199;
t325 = t170 * mrSges(11,2);
t295 = -qJ(1) + t306;
t291 = qJ(4) + t295;
t287 = -qJ(8) + t291;
t292 = -qJ(4) + t295;
t288 = -qJ(8) + t292;
t324 = (pkin(3) * (cos(t292) - cos(t291)) + (cos(qJ(2) - t288) - cos(qJ(2) + t287)) * pkin(5)) / (cos(t288) - cos(t287));
t276 = 0.1e1 / pkin(3);
t323 = t178 * t276;
t309 = -cos(0.2e1 * t252) + cos(0.2e1 * t306);
t322 = (-t309 * pkin(5) + (-cos(qJ(8) + 0.2e1 * qJ(1)) + cos((2 * pkin(8)) + (2 * qJ(5)) + qJ(8))) * pkin(1)) / t309;
t313 = qJ(2) + qJ(4);
t247 = qJ(1) + t313;
t228 = sin(t247);
t229 = cos(t247);
t286 = t228 * t244 - t243 * t229;
t182 = 0.1e1 / t286;
t285 = t228 * t258 + t229 * t265;
t321 = t182 * t285;
t230 = pkin(4) + t339;
t194 = -t254 * t230 - t255 * t340;
t320 = t194 * mrSges(11,2);
t242 = sin(t313);
t207 = -t242 * pkin(1) - t340;
t251 = 0.1e1 / t260;
t319 = t207 * t251;
t318 = t207 * t276;
t317 = (t262 * pkin(3) + t242 * pkin(4)) * t276;
t316 = t251 * t274;
t315 = t255 * mrSges(11,1);
t193 = -t255 * t230 + t254 * t340;
t184 = t193 * mrSges(11,1);
t314 = Ifges(11,3) + t184;
t235 = -qJ(9) + t247;
t216 = t235 - t312;
t234 = qJ(9) + t247;
t217 = t234 + t312;
t311 = sin(t217) - sin(t216);
t310 = cos(t217) - cos(t216);
t308 = sin(t235) - sin(t234);
t307 = cos(t235) - cos(t234);
t264 = cos(qJ(8));
t305 = -0.2e1 * mrSges(9,1) * t264;
t304 = mrSges(5,2) * t340;
t303 = pkin(6) * t329;
t302 = mrSges(4,1) * t338;
t161 = -t310 * t258 + t311 * t265;
t301 = t161 * t342;
t300 = pkin(4) * t315;
t299 = t276 * t324;
t298 = t272 * t322;
t297 = t251 * t318;
t294 = pkin(1) * t272 * t324;
t293 = t272 * t299;
t203 = (t259 * t262 - t266 * t269) * pkin(1);
t201 = t203 * mrSges(7,1);
t177 = Ifges(7,3) + t201 - t331;
t158 = Ifges(10,3) + t166 - t335;
t157 = Ifges(11,3) + t165 - t325;
t173 = t314 - t320;
t237 = t256 * pkin(6) * mrSges(10,2);
t206 = Ifges(10,3) + t237 - t303;
t283 = Ifges(3,3) + Ifges(7,3) + t328 + t336;
t236 = t254 * pkin(4) * mrSges(11,2);
t205 = Ifges(11,3) + t236 - t300;
t239 = mrSges(4,2) * t341;
t280 = Ifges(4,3) + t239 - t302 + t327;
t145 = Ifges(4,3) + m(10) * (-t171 * t263 - t172 * t256) * pkin(6) + t206 - t335 + t346;
t144 = m(11) * (-t169 * t255 - t170 * t254) * pkin(4) + t205 - t325 - t333 + t347;
t240 = mrSges(5,1) * t339;
t241 = mrSges(3,1) * t249;
t279 = t240 + t241 + t280 + m(10) * (t195 * t171 + t198 * t172) + m(11) * (t193 * t169 + t194 * t170) + m(4) * (-t197 * t268 - t200 * t261) * pkin(2) + t177 + (-t172 - t198) * mrSges(10,2) + (-t199 - t340) * mrSges(5,2) + (-t170 - t194) * mrSges(11,2) + m(5) * (t196 * t267 + t199 * t260) * pkin(3) + Ifges(3,3) + t314 + t346 + t347;
t250 = t257 ^ 2;
t248 = sin(qJ(1)) * pkin(1);
t246 = -qJ(9) + t306;
t245 = qJ(9) + t306;
t226 = qJ(9) + t289;
t225 = -qJ(9) + t290;
t219 = -sin(qJ(8) - t295);
t210 = -pkin(5) * sin(t252) + t248;
t209 = t337 - pkin(5) * cos(t252);
t186 = pkin(3) * t244 + pkin(4) * t229 + t337;
t185 = pkin(3) * t243 + pkin(4) * t228 + t248;
t181 = -0.2e1 * t303 + 0.2e1 * t237 + m(10) * (t256 ^ 2 + t263 ^ 2) * pkin(6) ^ 2 + t336;
t180 = -0.2e1 * t300 + 0.2e1 * t236 + m(11) * (t254 ^ 2 + t255 ^ 2) * pkin(4) ^ 2 + t328;
t176 = t285 * pkin(4) + t348;
t175 = t327 - t334;
t168 = t307 * t258 - t308 * t265;
t167 = t286 * pkin(4) - t348;
t160 = t168 * Ifges(7,3) * t342;
t159 = t175 * t321;
t154 = -t334 + t237 + (m(10) * (-t195 * t263 - t198 * t256) - t329) * pkin(6) + t280;
t153 = -t304 + Ifges(5,3) + t236 + t240 + (m(11) * (-t193 * t255 - t194 * t254) - t315) * pkin(4) + t173;
t152 = -Ifges(7,3) * t321 + t160;
t149 = (t173 * t318 + (-Ifges(11,3) * t262 + t205 * t317) * t343) * t251 + t157;
t148 = 0.2e1 * t184 + 0.2e1 * t239 + 0.2e1 * t240 + m(4) * (t261 ^ 2 + t268 ^ 2) * pkin(2) ^ 2 + m(5) * (t260 ^ 2 + t267 ^ 2) * pkin(3) ^ 2 - 0.2e1 * t320 - 0.2e1 * t334 + t283 - 0.2e1 * t304 - 0.2e1 * t302 + 0.2e1 * t187 + m(10) * (t195 ^ 2 + t198 ^ 2) + m(11) * (t193 ^ 2 + t194 ^ 2);
t147 = -Ifges(10,3) * t321 + t206 * t301 - t159;
t146 = (-t173 * t285 + (Ifges(11,3) * t167 + t176 * t205) * t274) * t182;
t143 = ((t308 * t185 + t307 * t186) * t323 + (-(sin(t246) - sin(t245)) * t210 - (cos(t246) - cos(t245)) * t209) * t326) * pkin(2);
t142 = t143 * Ifges(7,3);
t141 = ((-t311 * t185 - t310 * t186) * t323 + ((sin(t226) - sin(t225)) * t210 + (cos(t226) - cos(t225)) * t209) * t326) * pkin(2);
t140 = t181 * t301 - (t154 + t206) * t321;
t139 = (-t153 * t285 + (t167 * t205 + t176 * t180) * t274) * t182;
t138 = Ifges(7,3) * t297 + t142 + t177;
t137 = (t153 * t318 + (t180 * t317 - t205 * t262) * t343) * t251 + t144;
t136 = t141 * t206 + (-Ifges(10,3) * t294 + t175 * t319) * t276 + t158;
t135 = -pkin(1) * t330 + t279;
t134 = t141 * t181 + (t154 * t319 - t206 * t294) * t276 + t145;
t133 = t154 * t301 - t159 + t160 + (-t148 * t285 + (t153 * t176 + t167 * t173) * t274) * t182;
t132 = t142 + t279 + t141 * t154 + (-t175 * t293 - t330 + (t153 * t317 - t173 * t262) * t316) * pkin(1) + t148 * t297;
t1 = [(t345 + (0.2e1 + t298) * Ifges(9,3)) * t298 + (t132 + t135) * t297 + (t305 - 0.2e1 * t330 + (t305 * t322 + (-t136 - t158) * t299) * t272 + ((-t149 - t157) * t262 + (t137 + t144) * t317) * t316 + (m(3) * (t262 ^ 2 + t269 ^ 2) + m(9) * (t264 ^ 2 + t250)) * pkin(1)) * pkin(1) - 0.2e1 * t333 - 0.2e1 * t335 + 0.2e1 * t241 + m(4) * (t197 ^ 2 + t200 ^ 2) + m(7) * (t203 ^ 2 + t204 ^ 2) + m(11) * (t169 ^ 2 + t170 ^ 2) + m(10) * (t171 ^ 2 + t172 ^ 2) + m(5) * (t196 ^ 2 + t199 ^ 2) + 0.2e1 * t165 + 0.2e1 * t166 + t283 + t250 / t219 ^ 2 * Ifges(6,3) + (t138 + t177) * t143 + (t134 + t145) * t141 + 0.2e1 * t188 + 0.2e1 * t189 + 0.2e1 * t201 + t345 + Ifges(2,3) + Ifges(9,3) - 0.2e1 * t331 - 0.2e1 * t332 - 0.2e1 * t325, (t134 * t161 + t138 * t168) * t342 + ((t137 * t176 + t149 * t167) * t274 - (t132 + t136) * t285) * t182; t133 * t297 + t140 * t141 + t152 * t143 + (t145 * t161 + t168 * t177) * t342 + (-t147 * t293 + (t139 * t317 - t146 * t262) * t316) * pkin(1) + ((t144 * t176 + t157 * t167) * t274 - (t135 + t158) * t285) * t182, Ifges(8,3) + (t140 * t161 + t152 * t168) * t342 + ((t139 * t176 + t146 * t167) * t274 - (t133 + t147) * t285) * t182;];
Mq = t1;
