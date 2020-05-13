% Calculate joint inertia matrix for
% palh3m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% rSges [9x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [9x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh3m2DE2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_inertiaJ_slag_vp1: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE2_inertiaJ_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2DE2_inertiaJ_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m2DE2_inertiaJ_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:11:13
% EndTime: 2020-05-07 04:11:20
% DurationCPUTime: 6.45s
% Computational Cost: add. (2957->447), mult. (1638->540), div. (0->0), fcn. (1159->137), ass. (0->264)
t283 = cos(qJ(2));
t443 = 0.2e1 * t283;
t430 = m(6) / 0.2e1;
t374 = -qJ(2) - pkin(15);
t241 = pkin(18) - t374;
t200 = pkin(17) + qJ(3) + t241;
t163 = atan2(-sin(t200), -cos(t200));
t170 = atan2(sin(t241), -cos(t241));
t362 = pkin(17) - t170;
t119 = -t163 + t362;
t117 = -qJ(2) + t119;
t112 = cos(t117);
t115 = cos(t119);
t297 = rSges(9,2) ^ 2;
t304 = rSges(9,1) ^ 2;
t250 = -t297 + t304;
t261 = t283 ^ 2;
t278 = sin(qJ(2));
t375 = t278 * t283;
t345 = rSges(9,2) * t375;
t370 = t304 / 0.2e1 - t297 / 0.2e1;
t442 = t250 * t112 ^ 2 + 0.2e1 * t115 ^ 2 * (0.2e1 * rSges(9,1) * t345 - t250 * t261 + t370);
t287 = rSges(6,1) * pkin(8);
t286 = pkin(10) + rSges(6,3);
t385 = rSges(6,2) * t286;
t441 = (t287 - t385) * t430 + Icges(6,6) / 0.2e1;
t440 = (t287 + t385) * t430 - Icges(6,6) / 0.2e1;
t361 = m(8) + m(4) + m(5);
t254 = m(9) + t361;
t439 = m(6) + t254;
t282 = cos(qJ(3));
t311 = pkin(4) ^ 2;
t387 = t311 * m(5);
t224 = t282 ^ 2 * t387;
t265 = qJ(2) + qJ(3);
t236 = cos(t265);
t302 = rSges(4,2) ^ 2;
t309 = rSges(4,1) ^ 2;
t253 = t302 + t309;
t289 = m(5) + m(6);
t392 = m(6) * t311;
t438 = t253 * m(4) + t311 * t289 * t236 ^ 2 + Icges(4,3) + Icges(9,3) + t224 + 0.2e1 * m(9) * t297 + t392 / 0.2e1;
t245 = qJ(1) + t265;
t210 = sin(t245);
t215 = cos(t245);
t246 = qJ(1) - t265;
t216 = cos(t246);
t211 = sin(t246);
t407 = -t211 / 0.2e1;
t436 = -t215 * t216 / 0.2e1 + t210 * t407;
t435 = cos(qJ(4)) * rSges(6,1) - sin(qJ(4)) * rSges(6,2);
t194 = pkin(4) * t289 + m(4) * rSges(4,1);
t188 = pkin(16) + t200;
t162 = atan2(-sin(t188), cos(t188));
t373 = t162 - qJ(4);
t154 = qJ(3) + t373;
t150 = qJ(2) + t154;
t140 = qJ(1) + t150;
t129 = cos(t140);
t158 = t162 + qJ(4);
t153 = qJ(3) + t158;
t149 = qJ(2) + t153;
t141 = qJ(1) - t149;
t130 = cos(t141);
t267 = qJ(1) - qJ(4);
t238 = cos(t267);
t322 = t129 / 0.4e1 + t130 / 0.4e1 - t238 / 0.2e1;
t138 = qJ(1) + t149;
t127 = cos(t138);
t139 = qJ(1) - t150;
t128 = cos(t139);
t266 = qJ(1) + qJ(4);
t237 = cos(t266);
t323 = -t127 / 0.4e1 - t128 / 0.4e1 - t237 / 0.2e1;
t125 = sin(t140);
t126 = sin(t141);
t231 = sin(t267);
t324 = -t125 / 0.4e1 - t126 / 0.4e1 + t231 / 0.2e1;
t123 = sin(t138);
t124 = sin(t139);
t230 = sin(t266);
t325 = t123 / 0.4e1 + t124 / 0.4e1 + t230 / 0.2e1;
t100 = (t322 - t323) * rSges(6,2) + (t324 + t325) * rSges(6,1);
t432 = 0.2e1 * t100;
t101 = (-t324 + t325) * rSges(6,2) + (t322 + t323) * rSges(6,1);
t431 = 0.2e1 * t101;
t428 = pkin(1) * m(6);
t427 = pkin(4) * m(5);
t426 = pkin(4) * m(6);
t425 = pkin(8) * m(6);
t424 = m(4) * rSges(4,2);
t423 = m(4) * rSges(4,3);
t422 = m(5) * rSges(5,2);
t421 = m(5) * rSges(5,3);
t420 = m(7) * rSges(7,3);
t419 = m(8) * rSges(8,2);
t417 = m(9) * rSges(9,1);
t416 = m(9) * rSges(9,2);
t415 = m(9) * rSges(9,3);
t413 = rSges(6,1) * m(6);
t412 = rSges(3,2) * m(3);
t411 = rSges(6,2) * pkin(8);
t410 = rSges(6,2) * m(6);
t269 = qJ(1) - qJ(2);
t240 = cos(t269);
t205 = pkin(1) * t240;
t408 = -t205 / 0.2e1;
t279 = sin(qJ(1));
t247 = t279 * rSges(6,2);
t406 = t247 / 0.2e1;
t405 = t261 / 0.2e1;
t268 = qJ(1) + qJ(2);
t232 = sin(t268);
t401 = pkin(1) * t232;
t233 = sin(t269);
t204 = pkin(1) * t233;
t239 = cos(t268);
t400 = pkin(1) * t239;
t398 = pkin(4) * t210;
t397 = pkin(4) * t215;
t396 = pkin(4) * t216;
t394 = pkin(4) * t236;
t229 = sin(t265);
t393 = m(4) * t229;
t391 = rSges(9,1) * rSges(9,2);
t390 = t250 * m(9);
t300 = rSges(6,2) ^ 2;
t307 = rSges(6,1) ^ 2;
t389 = (-t300 + t307) * m(6);
t388 = t286 * m(6);
t386 = rSges(9,1) * t112;
t111 = sin(t117);
t384 = rSges(9,2) * t111;
t382 = t279 * rSges(6,1);
t284 = cos(qJ(1));
t380 = t284 * rSges(6,1);
t249 = t286 * rSges(6,1);
t114 = sin(t119);
t379 = t114 * t115;
t378 = (rSges(9,1) * t375 + rSges(9,2) * t261 - rSges(9,2)) * t114;
t313 = pkin(1) ^ 2;
t376 = t254 * t313;
t160 = qJ(3) + t162;
t263 = -qJ(4) + qJ(2);
t262 = qJ(4) + qJ(2);
t190 = pkin(4) * t407;
t372 = t190 + t398 / 0.2e1;
t191 = -t397 / 0.2e1;
t371 = t396 / 0.2e1 + t191;
t292 = 0.2e1 * qJ(2);
t369 = 0.2e1 * qJ(3) + t292;
t264 = t292 + qJ(3);
t299 = rSges(7,2) ^ 2;
t306 = rSges(7,1) ^ 2;
t368 = (t299 + t306);
t251 = t300 + t307;
t303 = rSges(3,2) ^ 2;
t310 = rSges(3,1) ^ 2;
t367 = (t303 + t310);
t353 = t282 * t427;
t350 = -0.2e1 * t386;
t349 = m(9) * t384;
t348 = -t413 / 0.2e1;
t347 = -t410 / 0.2e1;
t346 = t410 / 0.2e1;
t344 = t283 * t384;
t343 = t250 * t375;
t157 = qJ(2) + t160;
t342 = 0.2e1 * t162 + t369;
t341 = 0.2e1 * t390 - 0.4e1 * t224 + 0.2e1 * t387;
t156 = t162 + t264;
t155 = t162 + t369;
t293 = -0.2e1 * qJ(2);
t338 = t293 + t362;
t335 = t251 * m(6) + Icges(6,3);
t331 = pkin(4) * t413 / 0.2e1;
t329 = pkin(4) * t346;
t277 = sin(qJ(3));
t328 = t282 * t277 * t387;
t327 = -t398 + t401;
t326 = -t397 + t400;
t274 = sin(pkin(16));
t275 = cos(pkin(16));
t280 = sin(pkin(15));
t285 = cos(pkin(15));
t171 = t274 * t285 + t275 * t280;
t172 = -t274 * t280 + t275 * t285;
t120 = t171 * t282 + t277 * t172;
t121 = -t277 * t171 + t172 * t282;
t321 = t120 * t283 + t278 * t121;
t175 = t278 * t277 - t283 * t282;
t176 = t283 * t277 + t278 * t282;
t320 = t175 * t285 + t280 * t176;
t257 = m(9) * t391;
t319 = -t257 + t328;
t316 = t442 * m(9) + t438;
t242 = qJ(3) + t262;
t207 = sin(t242);
t243 = qJ(3) + t263;
t208 = sin(t243);
t212 = cos(t242);
t213 = cos(t243);
t315 = t236 * (rSges(4,2) * t423 - Icges(4,6)) + (rSges(9,2) * t415 - Icges(9,6)) * t112 + (-rSges(9,1) * t415 + Icges(9,5)) * t111 + (pkin(4) * t421 + rSges(4,1) * t423 - Icges(4,5)) * t229 + pkin(4) * t212 * t348 + t213 * t331 + (t207 + t208) * t329;
t314 = (-(t194 * t283 + t278 * t424) * t236 + (-rSges(4,1) * t278 + rSges(4,2) * t283) * t393 - t353) * pkin(1);
t312 = pkin(3) ^ 2;
t308 = rSges(5,1) ^ 2;
t305 = (rSges(8,1) ^ 2);
t301 = rSges(5,2) ^ 2;
t298 = (rSges(8,2) ^ 2);
t295 = -4 * Icges(6,5);
t294 = -2 * Icges(6,2);
t290 = 0.2e1 * qJ(4);
t270 = 0.4e1 * t387;
t258 = pkin(17) + pkin(18);
t255 = t283 * pkin(1);
t248 = t284 * rSges(6,2);
t244 = pkin(14) + t374;
t235 = cos(t263);
t234 = cos(t262);
t228 = sin(t263);
t227 = sin(t262);
t226 = cos(t258);
t225 = sin(t258);
t221 = 0.2e1 * t265;
t220 = -t382 / 0.2e1;
t217 = m(5) * rSges(5,1) + t425;
t214 = cos(t244);
t209 = sin(t244);
t206 = pkin(1) * t277 * t427;
t203 = 0.4e1 * t249;
t202 = 0.4e1 * pkin(1) * t353;
t201 = 0.2e1 * t244;
t199 = t388 - t422;
t197 = 0.4e1 * t390;
t173 = -rSges(9,1) * t261 + rSges(9,1) + t345;
t169 = qJ(2) + t170;
t168 = t292 + t170;
t165 = -qJ(2) + t362;
t164 = 0.2e1 * t169;
t152 = qJ(1) - t157;
t151 = qJ(1) + t157;
t148 = -qJ(4) + t156;
t147 = -qJ(4) + t155;
t146 = qJ(4) + t155;
t145 = -qJ(4) + t342;
t144 = qJ(4) + t342;
t143 = qJ(4) + t156;
t142 = cos(t157);
t137 = 0.2e1 * t157;
t136 = cos(t150);
t135 = cos(t149);
t134 = sin(t150);
t133 = sin(t149);
t132 = 0.2e1 * t150;
t131 = 0.2e1 * t149;
t122 = t280 * t175 - t176 * t285;
t118 = t293 - t163 - 0.2e1 * t170 + 0.2e1 * pkin(17);
t116 = -t163 + t338;
t109 = 0.2e1 * t117;
t104 = t122 * t274 - t320 * t275;
t103 = t122 * t275 + t320 * t274;
t102 = -t278 * t120 + t121 * t283;
t99 = (t248 + t382) * t238 + (t247 - t380) * t231 + (t248 - t382) * t237 + t230 * (t247 + t380) + (t129 + t130) * (t220 - t248 / 0.2e1) + (t128 + t127) * (t220 + t248 / 0.2e1) + (-t125 - t126) * (t406 - t380 / 0.2e1) + (t124 + t123) * (t406 + t380 / 0.2e1) + ((cos(t152) - cos(t151)) * t284 + (sin(t152) - sin(t151)) * t279) * rSges(6,3);
t1 = [(-0.2e1 * rSges(7,1) * t214 - 0.2e1 * rSges(7,2) * t209 + pkin(6)) * m(7) * pkin(6) + (cos(t131) + cos(t132)) * (-Icges(6,1) / 0.8e1 + Icges(6,2) / 0.8e1 + t389 / 0.8e1) + (t312 * cos(0.2e1 * t165) + 0.2e1 * rSges(9,3) ^ 2 + t297 + t304 + t312) * m(9) / 0.2e1 + (sin(t290) / 0.2e1 + sin(t132) / 0.4e1 - sin(t131) / 0.4e1) * (rSges(6,1) * t410 - Icges(6,4)) + cos(t145) * t441 + cos(t144) * t440 + (0.2e1 * rSges(5,3) ^ 2 + t301 + t308) * m(5) / 0.2e1 + ((2 * rSges(8,3) ^ 2 + t298 + t305) * m(8)) / 0.2e1 + (m(6) * t313 - Icges(3,1) + Icges(3,2) + t376 + (-t303 + t310) * m(3)) * cos(t292) / 0.2e1 + ((-t299 + t306) * m(7) - Icges(7,1) + Icges(7,2)) * cos(t201) / 0.2e1 + ((-t298 + t305) * m(8) - Icges(8,1) + Icges(8,2)) * cos(t164) / 0.2e1 + t435 * t425 + (2 * rSges(3,3) ^ 2 + t367) * m(3) / 0.2e1 + (0.2e1 * rSges(4,3) ^ 2 + t253) * m(4) / 0.2e1 + ((sin(-t373) + sin(t146) + sin(t158)) * t347 + (sin(t155) + sin(t162)) * t199 + (cos(t155) + cos(t162)) * t217) * pkin(4) + (cos(-t373) + cos(t158) + cos(t147) + cos(t146)) * t331 + t376 / 0.2e1 + (0.2e1 * Icges(6,1) + t294 - 0.2e1 * t389) * cos(t290) / 0.8e1 + t387 / 0.2e1 + (2 * rSges(7,3) ^ 2 + t368) * m(7) / 0.2e1 + (pkin(8) * t388 - rSges(5,1) * t422 + Icges(5,4)) * sin(t137) + (-rSges(4,1) * t424 + Icges(4,4)) * sin(t221) + (-rSges(8,1) * t419 + Icges(8,4)) * sin(t164) + ((t203 + 0.4e1 * t411) * m(6) + t295) * sin(t145) / 0.8e1 + ((t203 - 0.4e1 * t411) * m(6) + t295) * sin(t144) / 0.8e1 + (-rSges(3,1) * t412 + Icges(3,4)) * sin(t292) + (0.2e1 * (-0.2e1 * (pkin(8) + t286) * (-pkin(8) + t286) + t251) * m(6) + 0.4e1 * (-t301 + t308) * m(5) - (4 * Icges(5,1)) - 0.2e1 * Icges(6,1) + (4 * Icges(5,2)) + t294 + 0.4e1 * Icges(6,3)) * cos(t137) / 0.8e1 + (0.4e1 * t392 + t270 + 0.4e1 * (-t302 + t309) * m(4) - (4 * Icges(4,1)) + (4 * Icges(4,2))) * cos(t221) / 0.8e1 + ((sin(t163) - sin(t118)) * t416 + (-cos(t163) - cos(t118)) * t417) * pkin(3) + (-Icges(9,4) + t257) * sin(t109) + (m(7) * rSges(7,1) * rSges(7,2) - Icges(7,4)) * sin(t201) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + sin(t147) * t329 + (t439 * pkin(12) * t443 + (cos(t168) + cos(t170)) * m(8) * rSges(8,1) + (cos(t362) + cos(t338)) * pkin(3) * m(9) + (-sin(t116) - t114) * t416 + (-cos(t116) - t115) * t417 + (-sin(t168) - sin(t170)) * t419 + (sin(t264) + t277) * t424 + (sin(t153) + sin(t143)) * t346 + (sin(t154) + sin(t148)) * t347 + (cos(t154) + cos(t153) + cos(t148) + cos(t143)) * t348 + (-sin(t156) - sin(t160)) * t199 + (-cos(t264) - t282) * t194 + (-cos(t156) - cos(t160)) * t217) * pkin(1) + (0.2e1 * (-rSges(8,2) * sin(t169) + rSges(8,1) * cos(t169)) * m(8) + (t133 - t134) * t410 + (-t135 - t136) * t413 + (0.2e1 * pkin(3) * cos(t165) + t350) * m(9) - 0.2e1 * t199 * sin(t157) - 0.2e1 * t217 * t142 - 0.2e1 * t194 * t236 - 0.2e1 * t278 * t412 - 0.2e1 * t349 + m(3) * rSges(3,1) * t443 + 0.2e1 * rSges(4,2) * t393 + (m(3) + t439) * pkin(12)) * pkin(12) + (0.4e1 * pkin(8) ^ 2 + 0.4e1 * t286 ^ 2 + 0.6e1 * t300 + 0.6e1 * t307 + 0.4e1 * t311 + 0.4e1 * t313) * m(6) / 0.8e1 + (t197 - (4 * Icges(9,1)) + (4 * Icges(9,2))) * cos(t109) / 0.8e1 + Icges(3,1) / 0.2e1 + Icges(4,1) / 0.2e1 + Icges(5,1) / 0.2e1 + Icges(6,1) / 0.4e1 + Icges(7,1) / 0.2e1 + Icges(8,1) / 0.2e1 + Icges(9,1) / 0.2e1 + Icges(3,2) / 0.2e1 + Icges(4,2) / 0.2e1 + Icges(5,2) / 0.2e1 + Icges(6,2) / 0.4e1 + Icges(7,2) / 0.2e1 + Icges(8,2) / 0.2e1 + Icges(9,2) / 0.2e1 + Icges(2,3) + Icges(6,3) / 0.2e1; Icges(3,5) * t278 + t283 * Icges(3,6) + (-rSges(3,1) * t278 - rSges(3,2) * t283) * rSges(3,3) * m(3) + (-rSges(7,2) * t420 + Icges(7,6)) * t214 + (rSges(7,1) * t420 - Icges(7,5)) * t209 + (-(m(8) * rSges(8,3)) - t415 - t421 - t423) * pkin(1) * t278 + ((-t228 / 0.2e1 - t227 / 0.2e1) * rSges(6,2) + (-t235 / 0.2e1 + t234 / 0.2e1) * rSges(6,1)) * t428 + t315; 0.2e1 * t314 + ((t255 - t384) * t350 + 0.2e1 * pkin(1) * t378 - 0.2e1 * pkin(1) * t344 + t313 + 0.2e1 * ((-t343 + (-0.2e1 * t261 + 0.1e1) * t391) * t114 - pkin(1) * t173) * t115) * m(9) + t361 * t313 + (t202 + t341) * t405 + ((t204 - t327) * t190 - t327 * t204 / 0.2e1 - (t205 - t326) * t396 / 0.2e1 + t326 * t408 - t398 * t401 / 0.2e1 + t191 * t400 + (0.1e1 / 0.2e1 + t261) * t313) * m(6) + t316 + t367 * m(3) + t368 * m(7) - 0.2e1 * (t206 - t319) * t375 + Icges(3,3) + Icges(7,3); t315; (t197 + t202 - 0.8e1 * t224 + t270) * t261 / 0.4e1 - (t206 + 0.2e1 * t257 - 0.2e1 * t328) * t375 + t314 + (-(t255 - 0.2e1 * t384) * t386 - (0.2e1 * t343 + (0.4e1 * t261 - 0.2e1) * t391) * t379 + (-t173 * t115 - t344 + t378) * pkin(1) + t442) * m(9) + (-t205 * t216 / 0.4e1 - t204 * t211 / 0.4e1 + t436 * pkin(4) + (t239 * t216 / 0.4e1 + t232 * t211 / 0.4e1 + (t240 / 0.4e1 - t239 / 0.4e1) * t215 + (t233 / 0.4e1 - t232 / 0.4e1) * t210) * pkin(1)) * t426 + t438; -0.4e1 * (t370 * t375 + (t261 - 0.1e1 / 0.2e1) * t391) * m(9) * t379 + 0.2e1 * t349 * t386 + t316 + 0.2e1 * t319 * t375 + t341 * t405 + t436 * t392; t136 * t441 + ((t249 + t411) * m(6) - Icges(6,5)) * t134 / 0.2e1 + t135 * t440 + ((t249 - t411) * m(6) - Icges(6,5)) * t133 / 0.2e1 + t335 * t142 - pkin(12) * t435 * m(6) - ((t207 - t208) * rSges(6,2) + (-t212 - t213) * rSges(6,1)) * t426 / 0.2e1 - ((-t227 + t228) * rSges(6,2) + (t234 + t235) * rSges(6,1)) * t428 / 0.2e1; ((t204 / 0.2e1 + t372) * t432 + (t408 + t371) * t431 + t99 * sin(atan2(t103 * t226 - t104 * t225, t103 * t225 + t104 * t226) + t265) * (t255 - t394) + (-t100 * t232 + t101 * t239) * pkin(1)) * t430; (t372 * t432 + t371 * t431 - t99 * sin(atan2(-t102 * t225 - t321 * t226, t102 * t226 - t225 * t321) + t265) * t394) * t430; t335;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
