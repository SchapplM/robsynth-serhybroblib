% Calculate minimal parameter regressor of Coriolis joint torque vector for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see palh2m1OL_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:46
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = palh2m1OL_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'palh2m1OL_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:45:20
% EndTime: 2020-06-30 17:45:28
% DurationCPUTime: 3.52s
% Computational Cost: add. (2943->288), mult. (7585->407), div. (0->0), fcn. (6176->8), ass. (0->132)
t305 = cos(qJ(5));
t358 = qJD(5) * t305;
t307 = cos(qJ(3));
t308 = cos(qJ(2));
t363 = qJD(1) * t308;
t344 = t307 * t363;
t303 = sin(qJ(3));
t304 = sin(qJ(2));
t364 = qJD(1) * t304;
t347 = t303 * t364;
t277 = -t344 + t347;
t322 = t303 * t308 + t304 * t307;
t278 = t322 * qJD(1);
t302 = sin(qJ(4));
t306 = cos(qJ(4));
t249 = t306 * t277 + t278 * t302;
t404 = t249 * t305;
t409 = t358 - t404;
t297 = qJD(2) + qJD(3);
t384 = qJD(2) * pkin(2);
t350 = t307 * t384;
t281 = pkin(3) * t297 + t350;
t351 = t303 * t384;
t335 = t302 * t351;
t265 = t281 * t306 - t335;
t295 = qJD(4) + t297;
t262 = -pkin(4) * t295 - t265;
t408 = t249 * t262;
t285 = qJD(1) * pkin(1) + pkin(2) * t363;
t264 = -pkin(3) * t277 + t285;
t326 = t277 * t302 - t306 * t278;
t361 = qJD(4) * t302;
t346 = t281 * t361;
t407 = -t326 * t264 - t346;
t301 = sin(qJ(5));
t359 = qJD(5) * t301;
t353 = qJD(1) * qJD(2);
t343 = t304 * t353;
t254 = qJD(3) * t347 - t297 * t344 + t303 * t343;
t318 = t322 * qJD(3);
t261 = qJD(2) * t322 + t318;
t312 = t261 * qJD(1);
t360 = qJD(4) * t306;
t220 = t306 * t254 + t277 * t360 + t278 * t361 + t302 * t312;
t369 = t305 * t220 + t295 * t358;
t210 = -t326 * t359 + t369;
t208 = t210 * t301;
t209 = t210 * t305;
t235 = t295 * t301 + t305 * t326;
t383 = t220 * t301;
t211 = qJD(5) * t235 + t383;
t221 = qJD(4) * t326 + t254 * t302 - t306 * t312;
t219 = t305 * t221;
t379 = t326 * t301;
t233 = -t305 * t295 + t379;
t354 = -qJD(5) + t249;
t217 = t301 * t221;
t370 = -t354 * t358 + t217;
t403 = t354 * t301;
t406 = (t295 * t326 - t221) * MDP(21) + (-t249 ^ 2 + t326 ^ 2) * MDP(19) + (-t249 * t295 + t220) * MDP(20) - t249 * t326 * MDP(18) + (t409 * t235 + t208) * MDP(25) + (-t235 * t326 + t354 * t404 + t370) * MDP(27) + (t233 * t326 - t354 * t403 + t219) * MDP(28) + (-t301 * t211 - t409 * t233 + t235 * t403 + t209) * MDP(26);
t405 = t249 * t264;
t373 = t303 * t306;
t324 = t302 * t307 + t373;
t345 = t303 * t360;
t388 = qJD(3) * t324 + t345;
t402 = -t384 * t388 + t407;
t401 = (t304 ^ 2 - t308 ^ 2) * MDP(5) + (t308 * MDP(10) + t304 * MDP(9)) * pkin(1);
t228 = pkin(4) * t326 - pkin(6) * t249;
t222 = -pkin(4) * t249 - pkin(6) * t326 + t264;
t266 = t281 * t302 + t306 * t351;
t263 = pkin(6) * t295 + t266;
t215 = t222 * t305 - t263 * t301;
t314 = t388 * pkin(2);
t239 = qJD(2) * t314 + t346;
t321 = -t215 * t326 - t239 * t305 + t262 * t359;
t400 = -t301 * t408 + t321;
t216 = t222 * t301 + t263 * t305;
t330 = t216 * t326 + t239 * t301 + t262 * t358;
t399 = -t262 * t404 + t330;
t357 = t326 * MDP(29);
t394 = t354 * t357;
t279 = t303 * t304 - t307 * t308;
t260 = t297 * t279;
t325 = t306 * t279 + t302 * t322;
t224 = qJD(4) * t325 + t260 * t306 + t261 * t302;
t259 = t279 * t302 - t306 * t322;
t293 = t308 * pkin(2) + pkin(1);
t269 = -pkin(3) * t279 + t293;
t227 = -pkin(4) * t325 - pkin(6) * t259 + t269;
t334 = qJD(3) * t350;
t367 = (qJD(3) + qJD(4)) * t335;
t238 = t306 * (qJD(4) * t281 + t334) - t367;
t389 = qJD(5) * t227 * t354 + t262 * t224 + t239 * t259 + (qJD(5) * t222 + t238) * t325;
t386 = pkin(3) * t278;
t385 = t304 * pkin(2);
t378 = t259 * t262;
t375 = t278 * t285;
t374 = t302 * t303;
t365 = MDP(11) * t278;
t362 = qJD(3) * t297;
t352 = pkin(2) * t364;
t349 = t304 * t384;
t226 = t228 - t386;
t290 = pkin(3) * t302 + pkin(6);
t336 = qJD(5) * t290 + t226;
t274 = t324 * t384;
t333 = pkin(3) * t361 - t274;
t323 = t306 * t307 - t374;
t275 = t323 * t384;
t332 = -pkin(3) * t360 + t275;
t331 = (-qJD(3) + t297) * t384;
t329 = t367 - t405;
t225 = qJD(4) * t259 + t260 * t302 - t306 * t261;
t251 = -pkin(3) * t261 - t349;
t328 = -(pkin(4) * t225 - pkin(6) * t224 + t251) * t354 + t227 * t221;
t327 = -t221 * t290 - t408;
t320 = t224 * t305 - t259 * t359;
t319 = t293 * t322;
t313 = (-t277 * t297 + t254) * MDP(13) + (-t278 * t297 + t312) * MDP(14) + (-t277 ^ 2 + t278 ^ 2) * MDP(12) + t406;
t240 = (-pkin(3) * t318 + (-pkin(3) * t322 - t385) * qJD(2)) * qJD(1);
t309 = qJD(2) ^ 2;
t292 = pkin(2) * t307 + pkin(3);
t291 = -pkin(3) * t306 - pkin(4);
t273 = pkin(2) * t373 + t292 * t302 + pkin(6);
t272 = pkin(2) * t374 - t292 * t306 - pkin(4);
t267 = -t352 - t386;
t253 = t292 * t361 + t314;
t252 = t292 * t360 + (qJD(3) * t323 - t303 * t361) * pkin(2);
t223 = t226 - t352;
t205 = t221 * pkin(4) - t220 * pkin(6) + t240;
t204 = t305 * t205;
t1 = [(t221 * t269 + t225 * t264 - t240 * t325 - t249 * t251) * MDP(23) + (t220 * t269 + t224 * t264 + t240 * t259 + t251 * t326) * MDP(24) + (-t221 * t325 - t225 * t354) * MDP(29) + (t220 * t325 - t221 * t259 + t224 * t249 - t225 * t326) * MDP(19) + (t220 * t259 + t224 * t326) * MDP(18) + (t277 * t349 - t285 * t261 + (-qJD(3) * t319 + (t279 * t385 - t319) * qJD(2)) * qJD(1)) * MDP(16) + (t209 * t259 + t235 * t320) * MDP(25) + ((-t233 * t305 - t235 * t301) * t224 + (-t208 - t211 * t305 + (t233 * t301 - t235 * t305) * qJD(5)) * t259) * MDP(26) + (-t259 * t217 + t211 * t325 - t225 * t233 - (-t224 * t301 - t259 * t358) * t354) * MDP(28) + (-t216 * t225 + t389 * t305 + ((-qJD(5) * t263 + t205) * t325 - qJD(5) * t378 - t328) * t301) * MDP(31) + (-t204 * t325 + t215 * t225 + ((t263 * t325 + t378) * qJD(5) + t328) * t305 + t389 * t301) * MDP(30) + (-t210 * t325 + t219 * t259 + t225 * t235 - t320 * t354) * MDP(27) + (t254 * t293 + t260 * t285 + 0.2e1 * t278 * t349) * MDP(17) + (t254 * t279 + t260 * t277 - t278 * t261 - t312 * t322) * MDP(12) + t309 * t304 * MDP(7) + (-t254 * t322 - t260 * t278) * MDP(11) + (0.2e1 * MDP(4) * t343 - MDP(6) * t309) * t308 + (MDP(13) * t260 + MDP(14) * t261) * t297 + (MDP(20) * t224 - MDP(21) * t225) * t295 - 0.2e1 * t401 * t353; (-MDP(17) * t285 + t365) * t277 + (-t304 * t308 * MDP(4) + t401) * qJD(1) ^ 2 + ((-t277 * t364 - t303 * t362) * MDP(16) + (-t278 * t364 - t307 * t362) * MDP(17) + (-MDP(23) * t345 + ((-MDP(23) * t306 - MDP(16)) * t303 + (-MDP(23) * t302 - MDP(24) * t306 - MDP(17)) * t307) * qJD(3)) * qJD(2)) * pkin(2) + (t249 * t267 - t253 * t295 + t407) * MDP(23) + (-t252 * t295 - t267 * t326 - t281 * t360 + t329) * MDP(24) + (t211 * t272 - t273 * t217 + t233 * t253 + t400) * MDP(30) + (t210 * t272 - t219 * t273 + t235 * t253 + t399) * MDP(31) - (-t357 + (-t223 * t305 - t252 * t301 - t273 * t358) * MDP(30) + (t223 * t301 - t252 * t305 + t273 * t359) * MDP(31)) * t354 + MDP(16) * t375 + t313; (t303 * t331 + t375) * MDP(16) + (-t277 * t285 + t307 * t331) * MDP(17) + (t326 * t386 + t275 * t295 + (-t334 + (-pkin(3) * t295 - t281) * qJD(4)) * t306 + t329) * MDP(24) + t394 + (t291 * t210 + t327 * t305 + t333 * t235 - (t301 * t336 + t305 * t332) * t354 + t330) * MDP(31) + (t274 * t295 + (-t249 * t278 - t295 * t361) * pkin(3) + t402) * MDP(23) + (t291 * t211 + t327 * t301 + t333 * t233 - (t301 * t332 - t305 * t336) * t354 + t321) * MDP(30) + t277 * t365 + t313; (t266 * t295 + t402) * MDP(23) + (t265 * t295 - t238 - t405) * MDP(24) + t394 + (-pkin(4) * t211 + (t228 * t305 - t265 * t301) * t354 - t266 * t233 - t370 * pkin(6) + t400) * MDP(30) + (-pkin(4) * t210 - (t228 * t301 + t265 * t305) * t354 - t266 * t235 + (-t354 * t359 - t219) * pkin(6) + t399) * MDP(31) + t406; t235 * t233 * MDP(25) + (-t233 ^ 2 + t235 ^ 2) * MDP(26) + (-t233 * t354 + t369) * MDP(27) + (-t235 * t354 - t383) * MDP(28) + t221 * MDP(29) + (-t216 * t354 - t235 * t262 - t238 * t301 + t204) * MDP(30) + (-t205 * t301 - t215 * t354 + t233 * t262 - t238 * t305) * MDP(31) + (-MDP(27) * t379 - MDP(28) * t235 - MDP(30) * t216 - MDP(31) * t215) * qJD(5);];
tauc = t1;
