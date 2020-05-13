% Calculate matrix of centrifugal and coriolis load on the joints for
% palh1m2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = palh1m2TE_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2TE_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_coriolismatJ_fixb_slag_vp1: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2TE_coriolismatJ_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2TE_coriolismatJ_fixb_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m2TE_coriolismatJ_fixb_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:32:32
% EndTime: 2020-05-01 20:32:55
% DurationCPUTime: 5.55s
% Computational Cost: add. (1909->380), mult. (2243->521), div. (0->0), fcn. (454->92), ass. (0->262)
t200 = qJ(3) + pkin(19);
t177 = qJ(2) + t200;
t145 = sin(t177);
t149 = cos(t177);
t380 = m(9) * rSges(9,2);
t386 = pkin(2) * m(10);
t425 = -(m(9) * rSges(9,1) + t386) * t145 - t149 * t380;
t333 = qJD(1) * pkin(15);
t137 = 0.2e1 * t177;
t205 = (qJ(3) + qJ(2));
t166 = 2 * t205;
t357 = m(11) * rSges(11,2);
t377 = rSges(4,2) * m(4);
t410 = (rSges(9,1) * t380 - Icges(9,4)) * cos(t137) + (rSges(11,1) * t357 + rSges(4,1) * t377 - Icges(11,4) - Icges(4,4)) * cos(t166);
t422 = -t410 * qJD(1) + t425 * t333;
t227 = m(5) + m(6);
t344 = (pkin(5) ^ 2 * t227 + (rSges(11,1) ^ 2 - rSges(11,2) ^ 2) * m(11) + (rSges(4,1) ^ 2 - rSges(4,2) ^ 2) * m(4) - Icges(4,1) + Icges(4,2) - Icges(11,1) + Icges(11,2)) * sin(t166);
t222 = -rSges(6,3) - pkin(11);
t135 = m(5) * rSges(5,2) + t222 * m(6);
t196 = pkin(22) + pkin(21);
t286 = pkin(18) - t196;
t252 = -pkin(20) + t286;
t126 = qJ(2) + t252;
t104 = qJ(3) + t126;
t82 = cos(t104);
t127 = -qJ(2) + t252;
t105 = -qJ(3) + t127;
t83 = cos(t105);
t414 = (t82 / 0.2e1 - t83 / 0.2e1) * t135;
t421 = -pkin(5) * t414 + t344 / 0.2e1;
t91 = qJ(4) + t104;
t61 = sin(t91);
t92 = -qJ(4) + t104;
t62 = sin(t92);
t419 = t61 + t62;
t93 = qJ(4) + t105;
t63 = sin(t93);
t94 = -qJ(4) + t105;
t64 = sin(t94);
t418 = t63 + t64;
t65 = cos(t91);
t68 = cos(t94);
t417 = t65 + t68;
t66 = cos(t92);
t67 = cos(t93);
t416 = t67 + t66;
t172 = m(8) + m(11) + m(4) + t227;
t220 = cos(qJ(2));
t415 = t172 * t220;
t199 = qJD(2) + qJD(3);
t413 = ((t67 + t68) * rSges(6,2) + (t63 - t64) * rSges(6,1)) * t199;
t232 = 0.2e1 * qJ(2);
t206 = t232 + qJ(3);
t411 = -pkin(5) * t227 - m(4) * rSges(4,1);
t399 = m(11) * rSges(11,1) - t411;
t340 = t399 * cos(t206);
t219 = cos(qJ(3));
t412 = t399 * t219;
t192 = pkin(17) + qJ(2) - pkin(18);
t156 = sin(t192);
t159 = cos(t192);
t142 = 0.2e1 * t192;
t375 = rSges(10,2) * m(10);
t378 = rSges(3,2) * m(3);
t238 = ((rSges(7,1) ^ 2 - rSges(7,2) ^ 2) * m(7) - Icges(7,1) + Icges(7,2)) * sin(t142) / 0.2e1 + (m(7) * rSges(7,1) * rSges(7,2) - Icges(7,4)) * cos(t142) + (rSges(3,1) * t378 + rSges(10,1) * t375 - Icges(3,4) - Icges(10,4)) * cos(t232);
t345 = (pkin(1) ^ 2 * t172 + (rSges(10,1) ^ 2 - rSges(10,2) ^ 2) * m(10) + (rSges(3,1) ^ 2 - rSges(3,2) ^ 2) * m(3) - Icges(3,1) - Icges(10,1) + Icges(3,2) + Icges(10,2)) * sin(t232);
t373 = pkin(14) * m(7);
t409 = t345 / 0.2e1 + (rSges(7,1) * t159 - rSges(7,2) * t156) * t373 + t238;
t100 = sin(t127);
t101 = cos(t126);
t102 = cos(t127);
t257 = qJ(2) + t286;
t130 = cos(t257);
t258 = -qJ(2) + t286;
t131 = cos(t258);
t320 = pkin(18) - pkin(22);
t178 = qJ(2) + t320;
t146 = sin(t178);
t179 = -qJ(2) + t320;
t147 = sin(t179);
t150 = cos(t178);
t151 = cos(t179);
t160 = pkin(9) * m(6) + m(5) * rSges(5,1);
t384 = pkin(4) * m(11);
t99 = sin(t126);
t408 = ((t150 / 0.2e1 + t151 / 0.2e1) * rSges(8,1) - (t147 / 0.2e1 + t146 / 0.2e1) * rSges(8,2)) * m(8) - (t99 / 0.2e1 + t100 / 0.2e1) * t135 + (t101 / 0.2e1 + t102 / 0.2e1) * t160 - pkin(15) * t415 + (t130 / 0.2e1 + t131 / 0.2e1) * t384;
t133 = qJ(3) + t257;
t112 = cos(t133);
t134 = -qJ(3) + t258;
t113 = cos(t134);
t407 = rSges(11,2) * (t112 + t113);
t110 = sin(t133);
t111 = sin(t134);
t190 = qJ(4) + t205;
t157 = cos(t190);
t191 = -qJ(4) + t205;
t158 = cos(t191);
t314 = qJD(1) * t384;
t282 = t314 / 0.2e1;
t255 = rSges(11,1) * t282;
t283 = -t314 / 0.2e1;
t256 = rSges(11,1) * t283;
t356 = m(6) * qJD(1);
t302 = t356 / 0.4e1;
t279 = rSges(6,2) * t302;
t259 = pkin(5) * t279;
t303 = -t356 / 0.4e1;
t280 = rSges(6,2) * t303;
t260 = pkin(5) * t280;
t281 = rSges(6,1) * t302;
t262 = pkin(5) * t281;
t263 = pkin(5) * rSges(6,1) * t303;
t367 = qJD(1) / 0.2e1;
t290 = t160 * t367;
t275 = pkin(5) * t290;
t358 = pkin(5) * qJD(1);
t305 = -t358 / 0.2e1;
t276 = t160 * t305;
t343 = ((rSges(9,1) ^ 2 - rSges(9,2) ^ 2) * m(9) + pkin(2) ^ 2 * m(10) - Icges(9,1) + Icges(9,2)) * sin(t137);
t301 = -t343 / 0.2e1;
t385 = pkin(5) * m(6);
t315 = t385 / 0.2e1;
t316 = -t385 / 0.2e1;
t155 = sin(t191);
t354 = rSges(6,1) * t155;
t154 = sin(t190);
t355 = rSges(6,1) * t154;
t368 = -qJD(1) / 0.2e1;
t379 = m(9) * rSges(9,3);
t76 = sin(t104);
t77 = sin(t105);
t406 = qJD(1) * t301 + t110 * t255 + t111 * t256 + t417 * t259 + t416 * t260 + t419 * t262 + t418 * t263 + t76 * t275 + t77 * t276 + t282 * t407 + t368 * t344 + ((rSges(9,2) * t379 - Icges(9,6)) * t145 + t315 * t354 + (-rSges(9,1) * t379 - rSges(10,3) * t386 + Icges(9,5)) * t149 + (t355 + (t157 + t158) * rSges(6,2)) * t316) * t199 + t422;
t143 = t357 + t377;
t215 = sin(qJ(3));
t327 = t215 * t143;
t249 = t327 - t412;
t404 = t249 * pkin(1);
t124 = qJ(4) + t252;
t108 = -qJ(2) + t124;
t86 = cos(t108);
t125 = -qJ(4) + t252;
t109 = -qJ(2) + t125;
t87 = cos(t109);
t403 = t86 + t87;
t180 = t232 + t200;
t148 = sin(t180);
t152 = cos(t180);
t400 = -rSges(10,1) * t152 + rSges(10,2) * t148;
t395 = (t67 - t68) * rSges(6,2) + t418 * rSges(6,1);
t173 = sin(t200);
t174 = cos(t200);
t393 = (t173 / 0.2e1 + t148 / 0.2e1) * rSges(10,2) + (t174 / 0.2e1 - t152 / 0.2e1) * rSges(10,1);
t391 = 0.8e1 * m(4);
t390 = 0.8e1 * m(11);
t389 = -0.8e1 * t199;
t388 = 0.8e1 * t199;
t387 = -pkin(9) / 0.2e1;
t383 = m(5) * rSges(5,3);
t214 = sin(qJ(4));
t218 = cos(qJ(4));
t98 = t214 * rSges(6,1) + t218 * rSges(6,2);
t382 = m(6) * t98;
t381 = m(7) * rSges(7,3);
t225 = rSges(3,1) * m(3);
t223 = rSges(10,1) * m(10);
t376 = rSges(6,2) * pkin(9);
t374 = rSges(4,3) * m(4);
t181 = sin(t205);
t372 = t181 / 0.8e1;
t185 = cos(t205);
t371 = t185 / 0.8e1;
t370 = pkin(5) * t219;
t369 = t215 * pkin(5);
t366 = qJD(3) / 0.2e1;
t360 = pkin(1) * qJD(1);
t359 = pkin(1) * qJD(2);
t207 = qJ(2) + qJ(4);
t187 = cos(t207);
t353 = rSges(6,1) * t187;
t351 = rSges(6,2) * t157;
t350 = rSges(6,2) * t222;
t348 = rSges(4,3) * t199;
t347 = pkin(15) * t181;
t193 = t222 * rSges(6,1);
t341 = t399 * t181;
t334 = rSges(11,3) * t199;
t337 = pkin(5) * t383 * t389 + (-rSges(11,1) * t334 - rSges(11,2) * t333) * t390;
t336 = rSges(6,1) * qJD(1);
t335 = rSges(6,2) * qJD(1);
t231 = 0.2e1 * qJ(4);
t201 = sin(t231);
t273 = (rSges(6,1) ^ 2 - rSges(6,2) ^ 2) * m(6) - Icges(6,1) + Icges(6,2);
t332 = t273 * t201;
t182 = sin(t206);
t331 = t143 * t182;
t169 = m(6) * rSges(6,1) * rSges(6,2) - Icges(6,4);
t203 = cos(t231);
t328 = t169 * t203;
t325 = t223 + t225;
t324 = qJD(1) * t143;
t323 = Icges(4,5) + Icges(11,5);
t322 = Icges(4,6) + Icges(11,6);
t216 = sin(qJ(2));
t38 = (pkin(1) + t369) * t216 - t220 * t370 - pkin(15);
t321 = (qJD(4) * pkin(9) + t38 * qJD(1)) * t382;
t198 = qJD(2) - qJD(4);
t197 = qJD(2) + qJD(4);
t319 = m(6) * t359;
t317 = qJD(1) * t386;
t313 = qJD(1) * t373;
t312 = rSges(4,2) * t333;
t106 = qJ(2) + t124;
t78 = sin(t106);
t81 = sin(t109);
t311 = -t78 / 0.4e1 + t81 / 0.4e1;
t107 = qJ(2) + t125;
t79 = sin(t107);
t80 = sin(t108);
t310 = t79 / 0.4e1 - t80 / 0.4e1;
t84 = cos(t106);
t308 = t84 / 0.4e1 + t87 / 0.4e1;
t85 = cos(t107);
t307 = t85 / 0.4e1 + t86 / 0.4e1;
t117 = t158 * t335;
t304 = pkin(9) * t367;
t292 = t135 * t367;
t274 = 0.2e1 * t252;
t183 = sin(t207);
t208 = qJ(2) - qJ(4);
t184 = sin(t208);
t188 = cos(t208);
t271 = (t355 * qJD(1) + t155 * t336 + t157 * t335) * pkin(5) + (t184 * t335 + t188 * t336 + (-rSges(6,2) * t183 + t353) * qJD(1)) * pkin(1);
t250 = t331 - t340;
t248 = (m(6) * t350 + Icges(6,6)) * t218 + (t193 * m(6) + Icges(6,5)) * t214;
t103 = 0.2e1 * t252;
t247 = t248 * sin(t103) * t367 + (-0.2e1 * pkin(9) * t382 - 0.2e1 * t328 - t332) * qJD(1) * cos(t103) / 0.4e1;
t242 = t117 + (-t351 - t355 - t354) * qJD(1);
t153 = t375 + t378;
t241 = t216 * t153 - t325 * t220;
t240 = t199 * t216 * t370 + (t199 * t369 + t359) * t220;
t239 = t331 / 0.2e1 - t340 / 0.2e1 - t249 / 0.2e1;
t236 = pkin(15) * t185 * t324 + t110 * t256 + t111 * t255 + t416 * t259 + t417 * t260 + t418 * t262 + t419 * t263 + t77 * t275 + t76 * t276 + t283 * t407 + t343 * t367 + (rSges(6,2) * t158 * t315 + (t351 + (t154 + t155) * rSges(6,1)) * t316) * qJD(4) - t422;
t224 = rSges(6,1) * pkin(9);
t221 = cos(pkin(18));
t217 = sin(pkin(18));
t212 = cos(pkin(19));
t211 = cos(pkin(20));
t210 = sin(pkin(19));
t209 = sin(pkin(20));
t195 = qJD(2) + t366;
t176 = qJD(3) + t198;
t175 = qJD(3) + t197;
t171 = cos(t196);
t170 = sin(t196);
t129 = -qJ(4) + t274;
t128 = qJ(4) + t274;
t90 = 0.2e1 * t125;
t89 = 0.2e1 * t124;
t50 = (-rSges(11,1) * t333 + rSges(11,2) * t334) * t390;
t9 = t248 * qJD(4);
t2 = t404 + ((-rSges(10,1) * t212 - rSges(10,2) * t210) * t219 + (rSges(10,1) * t210 - rSges(10,2) * t212) * t215) * t386;
t1 = [((t195 * t148 + t173 * t366) * rSges(10,2) + (-t195 * t152 + t174 * t366) * rSges(10,1)) * t386 + (t241 * pkin(15) + t409) * qJD(2) + (((-t66 / 0.4e1 - t67 / 0.4e1) * t176 + (t65 / 0.4e1 + t68 / 0.4e1) * t175) * rSges(6,2) + ((t62 / 0.4e1 - t63 / 0.4e1) * t176 + (t61 / 0.4e1 - t64 / 0.4e1) * t175) * rSges(6,1)) * t385 + (-((t224 - t350) * m(6) - Icges(6,6)) * sin(t128) / 0.4e1 + ((t224 + t350) * m(6) + Icges(6,6)) * sin(t129) / 0.4e1 + ((-t193 - t376) * m(6) - Icges(6,5)) * cos(t128) / 0.4e1 - ((-t193 + t376) * m(6) - Icges(6,5)) * cos(t129) / 0.4e1 + t382 * t387 + (-cos(t89) / 0.4e1 - cos(t90) / 0.4e1 + t203 / 0.2e1) * t169 + (-sin(t89) / 0.8e1 + sin(t90) / 0.8e1 + t201 / 0.4e1) * t273 + ((cos(t124) / 0.2e1 + cos(t125) / 0.2e1) * rSges(6,2) + (sin(t124) / 0.2e1 - sin(t125) / 0.2e1) * rSges(6,1)) * m(6) * pkin(15)) * qJD(4) + (t301 - t344 / 0.2e1 + (-t143 * t185 - t341 + t425) * pkin(15) + ((t76 / 0.2e1 - t77 / 0.2e1) * t160 + t414) * pkin(5) + ((t112 / 0.2e1 + t113 / 0.2e1) * rSges(11,2) + (t110 / 0.2e1 - t111 / 0.2e1) * rSges(11,1)) * t384 - t410) * t199 + (t250 * t195 + (-t327 / 0.2e1 + t412 / 0.2e1) * qJD(3) + ((t311 * t197 + t310 * t198) * rSges(6,2) + (t308 * t197 + t307 * t198) * rSges(6,1)) * m(6) + t408 * qJD(2)) * pkin(1), t400 * t317 + (rSges(7,1) * t313 + qJD(2) * (rSges(7,2) * t381 - Icges(7,6))) * t159 + (-rSges(7,2) * t313 + qJD(2) * (rSges(7,1) * t381 - Icges(7,5))) * t156 + t216 * (t153 * t333 + (rSges(3,3) * t225 + rSges(10,3) * t223 - Icges(3,5) - Icges(10,5)) * qJD(2)) + (-0.8e1 * m(4) * t312 + (rSges(4,1) * t374 - t323) * t389 + t337) * t371 + t345 * t367 - pkin(5) * t83 * t292 - t319 * t353 / 0.2e1 + (qJD(2) * (rSges(3,3) * t378 + rSges(10,3) * t375 - Icges(3,6) - Icges(10,6)) - t325 * t333) * t220 + t238 * qJD(1) - t135 * t82 * t305 + (t50 + (rSges(4,2) * t374 - t322) * t388 + 0.8e1 * t411 * t333) * t372 + ((t184 + t183) * rSges(6,2) + rSges(6,1) * t188) * t319 / 0.2e1 + (t216 * (rSges(11,3) * m(11) + rSges(8,3) * m(8) + t374 + t383) * qJD(2) + (t81 + t79) * t279 + (t80 + t78) * t280 + (t131 + t130) * t282 + (t85 + t84 + t403) * t281 - t333 * t415 + t182 * t324 + (t102 + t101) * t290 - (t100 + t99) * t292) * pkin(1) + t406 + (-t340 + (-(t147 + t146) * rSges(8,2) / 0.2e1 + (t151 + t150) * rSges(8,1) / 0.2e1) * m(8)) * t360, (t50 + (-rSges(4,1) * t333 + rSges(4,2) * t348) * t391 + t322 * t389) * t372 + ((-rSges(4,1) * t348 - t312) * t391 + t323 * t388 + t337) * t371 + t239 * t360 + (-t227 * t347 + t414) * t358 + t393 * t317 + t406, -(t38 * qJD(4) + t304) * t382 + (t209 * t170 - t211 * t171) * (-t217 * t9 + t221 * t321) + (t328 / 0.2e1 + t332 / 0.4e1) * qJD(1) + t247 + (-t211 * t170 - t209 * t171) * (t217 * t321 + t9 * t221); ((-t241 + t341) * pkin(15) - t400 * t386 + (((-t310 - t311) * rSges(6,2) + (-t307 - t308) * rSges(6,1)) * m(6) - t250 - t408) * pkin(1) - t409 + t421) * qJD(1) + t236 + ((t183 / 0.2e1 - t184 / 0.2e1) * rSges(6,2) + (-t187 / 0.2e1 - t188 / 0.2e1) * rSges(6,1)) * qJD(4) * m(6) * pkin(1), -qJD(3) * (t404 + (-rSges(10,1) * t174 - rSges(10,2) * t173) * t386), -t2 * t199, -m(6) * (-pkin(5) * t117 + (-t395 * pkin(5) + ((-t80 + t81) * rSges(6,2) + t403 * rSges(6,1)) * pkin(1)) * qJD(4) + t271) / 0.2e1; (-t239 * pkin(1) + t399 * t347 - t393 * t386 + t421) * qJD(1) + t236, t2 * qJD(2), 0, (t395 * qJD(4) + t242) * t315; t328 * t368 - qJD(1) * t332 / 0.4e1 - t247 + ((t209 * t217 + t211 * t221) * t171 + (-t209 * t221 + t211 * t217) * t170) * t38 * t98 * t356 + ((t240 * rSges(6,1) + rSges(6,2) * t304) * t218 - t214 * (t240 * rSges(6,2) + t336 * t387)) * m(6), m(6) * (((-t80 - t81) * rSges(6,2) + (t86 - t87) * rSges(6,1)) * t359 + (-t117 - t413) * pkin(5) + t271) / 0.2e1, (t242 + t413) * t316, 0;];
Cq = t1;
