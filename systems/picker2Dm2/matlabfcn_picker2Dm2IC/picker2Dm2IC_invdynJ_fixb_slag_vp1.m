% Calculate vector of inverse dynamics joint torques with ic for
% picker2Dm2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% qJD [12x1]
%   Generalized joint velocities
% qJDD [12x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [11x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 09:21
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = picker2Dm2IC_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(12,1),zeros(3,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2IC_invdynJ_fixb_slag_vp1: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm2IC_invdynJ_fixb_slag_vp1: qJD has to be [12x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [12 1]), ...
  'picker2Dm2IC_invdynJ_fixb_slag_vp1: qJDD has to be [12x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2IC_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2IC_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2IC_invdynJ_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm2IC_invdynJ_fixb_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'picker2Dm2IC_invdynJ_fixb_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 09:20:39
% EndTime: 2020-05-11 09:20:44
% DurationCPUTime: 5.61s
% Computational Cost: add. (7398->348), mult. (4383->353), div. (37->10), fcn. (2146->59), ass. (0->226)
t231 = sin(qJ(1));
t218 = t231 * pkin(1);
t233 = cos(qJ(1));
t236 = qJD(1) ^ 2;
t302 = pkin(1) * qJDD(1);
t257 = t236 * t218 - t233 * t302;
t355 = t257 - g(2);
t323 = pkin(1) * t233;
t276 = t231 * t302 + t236 * t323;
t354 = t276 - g(1);
t221 = qJDD(1) + qJDD(2);
t199 = qJDD(6) + t221;
t223 = qJD(1) + qJD(2);
t202 = qJD(6) + t223;
t227 = qJ(1) + qJ(2);
t214 = qJ(6) + t227;
t181 = sin(t214);
t184 = cos(t214);
t102 = t184 * rSges(7,1) - rSges(7,2) * t181;
t59 = t202 * t102;
t99 = rSges(7,1) * t181 + rSges(7,2) * t184;
t353 = -t199 * t99 - t202 * t59 + t354;
t306 = t202 * t99;
t352 = t102 * t199 - t202 * t306 + t355;
t209 = sin(t227);
t211 = cos(t227);
t219 = t223 ^ 2;
t342 = t209 * t219 - t211 * t221;
t351 = pkin(3) * t342 + t355;
t350 = pkin(2) * t342 + t355;
t341 = t209 * t221 + t211 * t219;
t349 = pkin(2) * t341 + t354;
t348 = pkin(3) * t341 + t354;
t310 = pkin(1) * qJD(1);
t272 = t233 * t310;
t288 = t211 * t223;
t245 = pkin(3) * t288 + t272;
t287 = qJ(2) + qJ(4);
t215 = qJ(1) + t287;
t185 = cos(t215);
t203 = qJD(4) + t223;
t295 = t185 * t203;
t175 = qJD(10) + t203;
t187 = qJ(10) + t215;
t162 = sin(t187);
t163 = cos(t187);
t264 = t163 * rSges(11,1) - rSges(11,2) * t162;
t343 = t175 * t264;
t24 = pkin(4) * t295 + t245 - t343;
t196 = t231 * t310;
t291 = t209 * t223;
t285 = pkin(3) * t291 + t196;
t182 = sin(t215);
t297 = t182 * t203;
t86 = rSges(11,1) * t162 + rSges(11,2) * t163;
t54 = t175 * t86;
t23 = pkin(4) * t297 + t285 - t54;
t112 = -rSges(3,1) * t211 + t209 * rSges(3,2);
t331 = t223 * t112;
t69 = t272 - t331;
t200 = qJDD(4) + t221;
t165 = qJDD(10) + t200;
t197 = t203 ^ 2;
t340 = -t165 * t86 - t175 * t343 + (t182 * t200 + t185 * t197) * pkin(4) + t348;
t339 = t165 * t264 - t175 * t54 + (t182 * t197 - t185 * t200) * pkin(4) + t351;
t100 = rSges(5,1) * t182 + rSges(5,2) * t185;
t61 = rSges(5,1) * t295 - rSges(5,2) * t297;
t338 = t100 * t200 + t203 * t61 + t348;
t216 = qJ(3) + t227;
t183 = sin(t216);
t186 = cos(t216);
t101 = rSges(4,1) * t183 + rSges(4,2) * t186;
t201 = qJDD(3) + t221;
t204 = qJD(3) + t223;
t294 = t186 * t204;
t296 = t183 * t204;
t63 = rSges(4,1) * t294 - rSges(4,2) * t296;
t337 = -t101 * t201 - t204 * t63 + t349;
t110 = rSges(3,1) * t209 + rSges(3,2) * t211;
t347 = t110 * t221 - t223 * t331 + t354;
t267 = -rSges(5,1) * t185 + t182 * rSges(5,2);
t60 = rSges(5,1) * t297 + rSges(5,2) * t295;
t336 = t200 * t267 + t203 * t60 + t351;
t265 = t186 * rSges(4,1) - rSges(4,2) * t183;
t78 = t204 * t101;
t335 = t201 * t265 - t204 * t78 + t350;
t82 = rSges(3,1) * t291 + rSges(3,2) * t288;
t346 = t112 * t221 + t223 * t82 + t355;
t174 = qJDD(9) + t201;
t176 = qJD(9) + t204;
t198 = t204 ^ 2;
t271 = qJ(9) + t227;
t193 = qJ(3) + t271;
t171 = cos(t193);
t298 = t171 * t176;
t168 = sin(t193);
t299 = t168 * t176;
t45 = rSges(10,1) * t298 - rSges(10,2) * t299;
t90 = t168 * rSges(10,1) + t171 * rSges(10,2);
t334 = t174 * t90 + t176 * t45 + (-t183 * t201 - t186 * t198) * pkin(6) + t349;
t266 = -rSges(10,1) * t171 + t168 * rSges(10,2);
t44 = rSges(10,1) * t299 + rSges(10,2) * t298;
t333 = t174 * t266 + t176 * t44 + (-t183 * t198 + t186 * t201) * pkin(6) + t350;
t246 = pkin(2) * t288 + t272;
t79 = t204 * t265;
t31 = t246 - t79;
t344 = t31 * t78;
t284 = pkin(2) * t291 + t196;
t29 = t284 - t78;
t47 = -t59 + t272;
t226 = qJ(1) + qJ(8);
t208 = sin(t226);
t210 = cos(t226);
t111 = t210 * rSges(9,1) - rSges(9,2) * t208;
t222 = qJD(1) + qJD(8);
t81 = t222 * t111;
t68 = -t81 + t272;
t109 = rSges(9,1) * t208 + rSges(9,2) * t210;
t300 = t109 * t222;
t66 = t196 - t300;
t46 = t196 - t306;
t228 = sin(qJ(7));
t232 = cos(qJ(7));
t332 = pkin(3) * (t209 * t228 + t211 * t232);
t250 = -pkin(6) * t296 + t284;
t64 = t176 * t90;
t25 = t250 + t64;
t76 = t203 * t100;
t28 = t285 + t76;
t330 = -t29 * t63 + t344;
t77 = t203 * t267;
t30 = t245 - t77;
t329 = t28 * t61 - t30 * t60;
t244 = -pkin(6) * t294 + t246;
t65 = t176 * t266;
t26 = t244 - t65;
t328 = t25 * t45 - t26 * t44;
t326 = (t306 * t47 - t46 * t59 + (-t202 * t47 - t353) * t99 + (t202 * t46 + t352) * t102) * m(7);
t322 = pkin(2) * t211;
t321 = pkin(3) * t211;
t320 = pkin(4) * t185;
t237 = 0.1e1 / pkin(5);
t224 = pkin(8) + qJ(5);
t286 = qJ(3) - qJ(6);
t258 = t224 + t286;
t247 = t258 - t226;
t259 = t224 - t286;
t248 = t259 - t226;
t305 = t237 / (pkin(6) * (cos(t248) - cos(t247)) + (-cos(qJ(9) - t248) + cos(qJ(9) + t247)) * pkin(2));
t239 = 0.1e1 / pkin(3);
t270 = -qJ(9) - t286;
t37 = 0.1e1 / ((-cos(qJ(4) - t286) + cos(qJ(4) + t286)) * pkin(6) + (cos(qJ(4) + t270) - cos(qJ(4) - t270)) * pkin(2));
t303 = t239 * t37;
t192 = -qJ(9) + t215;
t151 = t192 - t286;
t191 = qJ(4) + t271;
t152 = t191 + t286;
t283 = sin(t152) - sin(t151);
t282 = cos(t152) - cos(t151);
t145 = Icges(11,3) * t165;
t281 = Icges(5,3) * t200 + t145;
t147 = Icges(10,3) * t174;
t280 = Icges(4,3) * t201 + t147;
t279 = -cos(0.2e1 * t226) + cos(0.2e1 * t224);
t278 = sin(t192) - sin(t191);
t277 = cos(t192) - cos(t191);
t229 = sin(qJ(4));
t225 = 0.1e1 / t229;
t238 = 0.1e1 / pkin(4);
t275 = pkin(1) * t225 * t238;
t268 = -qJ(1) + t224;
t67 = t110 * t223 + t196;
t261 = -qJ(4) + t268;
t260 = qJ(4) + t268;
t56 = -pkin(6) * t183 + t90;
t57 = pkin(6) * t186 + t266;
t256 = rSges(2,1) * t233 - rSges(2,2) * t231;
t142 = rSges(2,1) * t231 + rSges(2,2) * t233;
t205 = sin(t224);
t206 = cos(t224);
t108 = rSges(6,1) * t206 - rSges(6,2) * t205;
t107 = rSges(6,1) * t205 + rSges(6,2) * t206;
t143 = rSges(8,1) * t232 - t228 * rSges(8,2);
t141 = rSges(8,1) * t228 + rSges(8,2) * t232;
t255 = -qJ(8) + t261;
t254 = -qJ(8) + t260;
t253 = t182 * t211 - t185 * t209;
t252 = t182 * t228 + t185 * t232;
t180 = pkin(2) * t209;
t42 = t180 + t56;
t74 = t265 - t322;
t72 = t267 - t321;
t49 = t264 - t320;
t73 = -t101 + t180;
t179 = pkin(3) * t209;
t71 = t100 + t179;
t157 = Icges(7,3) * t199;
t249 = Icges(3,3) * t221 + t157 + t280 + t281;
t161 = pkin(4) * t182;
t48 = t161 - t86;
t43 = t57 - t322;
t38 = t179 + t48;
t39 = t49 - t321;
t230 = sin(qJ(2));
t220 = qJDD(1) + qJDD(8);
t213 = -qJ(9) + t224;
t212 = qJ(9) + t224;
t207 = sin(t287);
t189 = Icges(9,3) * t220;
t173 = qJ(9) + t258;
t172 = -qJ(9) + t259;
t125 = -pkin(5) * t208 + t218;
t124 = -pkin(5) * t210 + t323;
t93 = t111 - t323;
t91 = -t109 + t218;
t89 = t320 + t321 + t323;
t88 = t218 + t179 + t161;
t21 = t111 * t220 - t222 * t300 + t257;
t19 = -t109 * t220 - t222 * t81 + t276;
t6 = t157 + t326;
t5 = t147 + ((t176 * t26 + t334) * t90 + (t176 * t25 + t333) * t266 + t328) * m(10);
t4 = t145 + ((-t175 * t24 - t340) * t86 + (t175 * t23 + t339) * t264 - t23 * t343 + t24 * t54) * m(11);
t3 = t280 + (t333 * t57 + t334 * t56 + (-t44 + t64) * t26 + (t45 + t65) * t25) * m(10) + ((t204 * t29 + t335) * t265 + (-t204 * t31 - t337) * t101 + t330) * m(4);
t2 = t281 + ((t203 * t28 + t336) * t267 + (t203 * t30 + t338) * t100 + t329) * m(5) + (t339 * t49 + t340 * t48) * m(11);
t1 = t249 + (t25 * t65 + t26 * t64 + t333 * t43 + t334 * t42 + t328) * m(10) + (t28 * t77 + t30 * t76 + t336 * t72 + t338 * t71 + t329) * m(5) + (t29 * t79 + t335 * t74 + t337 * t73 + t330 - t344) * m(4) + (-t67 * t331 - t69 * t82 + (t223 * t67 + t346) * t112 + (t223 * t69 + t347) * t110) * m(3) + (t339 * t39 + t340 * t38) * m(11) + t326;
t7 = [sin(qJ(8)) / sin(-qJ(8) + t268) * (Icges(6,3) * qJDD(5) + ((t107 ^ 2 + t108 ^ 2) * qJDD(5) + g(1) * t107 - g(2) * t108) * m(6)) + t189 + (-t279 * pkin(5) + (-cos(qJ(8) + 0.2e1 * qJ(1)) + cos((2 * pkin(8)) + (2 * qJ(5)) + qJ(8))) * pkin(1)) * t237 / t279 * (t189 + (-t66 * t81 + t68 * t300 + (t66 * t222 - g(2) + t21) * t111 + (-t68 * t222 + g(1) - t19) * t109) * m(9)) - m(9) * (g(1) * t91 + g(2) * t93) + t249 + m(9) * (t19 * t91 + t21 * t93) - t230 * t4 * t275 + Icges(2,3) * qJDD(1) + (t25 * (t244 + t45) - t26 * (t250 + t44) + t333 * (t43 - t323) + t334 * (t218 + t42)) * m(10) + (t352 * (t102 - t323) + t353 * (t218 - t99)) * m(7) + (t28 * (t245 + t61) - t30 * (t285 + t60) + t336 * (t72 - t323) + t338 * (t218 + t71)) * m(5) + (t335 * (t74 - t323) + t337 * (t218 + t73) + (t246 - t63 - t31) * t29) * m(4) + (t346 * (t112 - t323) + t347 * (t110 + t218) + (t67 - t196 - t82) * t69) * m(3) + (t339 * (t39 - t323) + t340 * (t218 + t38)) * m(11) + ((qJDD(1) * t256 + g(2)) * t256 + (qJDD(1) * t142 - g(1)) * t142) * m(2) + (((t277 * t89 + t278 * t88) * t303 + (-(sin(t213) - sin(t212)) * t125 - (cos(t213) - cos(t212)) * t124) * t305) * t6 + ((-t282 * t89 - t283 * t88) * t303 + ((sin(t173) - sin(t172)) * t125 + (cos(t173) - cos(t172)) * t124) * t305) * t3) * pkin(2) + ((-pkin(1) * t207 - pkin(3) * t229) * t225 * t1 - (pkin(3) * (cos(t261) - cos(t260)) + (cos(-qJ(2) + t255) - cos(qJ(2) + t254)) * pkin(5)) * pkin(1) * t237 / (cos(t255) - cos(t254)) * t5 + (pkin(3) * t230 + pkin(4) * t207) * t2 * t275) * t239; Icges(8,3) * qJDD(7) + ((-t228 * t282 + t232 * t283) * t3 + (t228 * t277 - t232 * t278) * t6) * pkin(2) * t37 + ((t141 ^ 2 + t143 ^ 2) * qJDD(7) - g(1) * t143 - g(2) * t141) * m(8) + (-(t1 + t5) * t252 + ((pkin(4) * t252 + t332) * t2 + (pkin(4) * t253 - t332) * t4) * t238) / t253;];
tau = t7(:);
