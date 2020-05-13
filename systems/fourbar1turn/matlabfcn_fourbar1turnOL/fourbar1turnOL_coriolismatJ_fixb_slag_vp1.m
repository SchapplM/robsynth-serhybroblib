% Calculate matrix of centrifugal and coriolis load on the joints for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% m [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = fourbar1turnOL_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_coriolismatJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnOL_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnOL_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:49
% EndTime: 2020-04-12 19:40:55
% DurationCPUTime: 3.47s
% Computational Cost: add. (7626->311), mult. (11296->457), div. (0->0), fcn. (10540->8), ass. (0->194)
t226 = sin(qJ(1));
t330 = t226 / 0.2e1;
t229 = cos(qJ(1));
t327 = -t229 / 0.2e1;
t326 = t229 / 0.2e1;
t223 = qJ(2) + qJ(3);
t209 = sin(t223);
t210 = cos(t223);
t246 = rSges(4,1) * t209 + rSges(4,2) * t210;
t225 = sin(qJ(2));
t325 = pkin(2) * t225;
t264 = -t246 + t325;
t349 = t264 * t226;
t221 = t226 ^ 2;
t222 = t229 ^ 2;
t266 = t221 + t222;
t227 = cos(qJ(4));
t319 = rSges(5,1) * t227;
t265 = pkin(1) + t319;
t224 = sin(qJ(4));
t284 = t224 * t226;
t267 = rSges(5,2) * t284 + t229 * rSges(5,3);
t112 = -t226 * t265 + t267;
t283 = t224 * t229;
t260 = -rSges(5,2) * t283 + t226 * rSges(5,3);
t113 = t229 * t265 + t260;
t185 = rSges(5,1) * t224 + rSges(5,2) * t227;
t164 = t185 * t226;
t166 = t185 * t229;
t303 = Icges(5,4) * t224;
t177 = Icges(5,2) * t227 + t303;
t182 = Icges(5,1) * t227 - t303;
t348 = -(t182 / 0.2e1 - t177 / 0.2e1) * t224 - m(5) * (t112 * t164 - t113 * t166);
t285 = t210 * t229;
t289 = t209 * t229;
t151 = rSges(4,1) * t289 + rSges(4,2) * t285;
t281 = t225 * t229;
t101 = -pkin(2) * t281 + t151;
t228 = cos(qJ(2));
t279 = t226 * t228;
t282 = t225 * t226;
t142 = rSges(3,1) * t279 - rSges(3,2) * t282 - t229 * rSges(3,3);
t276 = t228 * t229;
t143 = rSges(3,1) * t276 - rSges(3,2) * t281 + t226 * rSges(3,3);
t186 = rSges(3,1) * t225 + rSges(3,2) * t228;
t165 = t186 * t226;
t167 = t186 * t229;
t306 = Icges(3,4) * t225;
t179 = Icges(3,2) * t228 + t306;
t184 = Icges(3,1) * t228 - t306;
t290 = t209 * t226;
t262 = -rSges(4,2) * t290 + rSges(4,3) * t229;
t320 = rSges(4,1) * t210;
t324 = pkin(2) * t228;
t88 = (t320 - t324) * t226 + t262;
t261 = -rSges(4,1) * t285 + t226 * rSges(4,3);
t318 = rSges(4,2) * t209;
t89 = (t318 + t324) * t229 + t261;
t347 = -(t184 / 0.2e1 - t179 / 0.2e1) * t225 - m(4) * (t101 * t89 + t349 * t88) - m(3) * (-t142 * t165 - t143 * t167);
t104 = t264 * t229;
t174 = t318 - t320;
t240 = Icges(4,5) * t209 + Icges(4,6) * t210;
t144 = t240 * t226;
t145 = t240 * t229;
t194 = Icges(4,4) * t289;
t111 = -Icges(4,1) * t285 + Icges(4,5) * t226 + t194;
t272 = Icges(4,2) * t285 + t111 + t194;
t193 = Icges(4,4) * t290;
t286 = t210 * t226;
t110 = Icges(4,1) * t286 + Icges(4,5) * t229 - t193;
t273 = -Icges(4,2) * t286 + t110 - t193;
t304 = Icges(4,4) * t210;
t170 = Icges(4,2) * t209 - t304;
t109 = Icges(4,6) * t226 + t170 * t229;
t244 = Icges(4,1) * t209 + t304;
t274 = -t244 * t229 + t109;
t108 = Icges(4,4) * t286 - Icges(4,2) * t290 + Icges(4,6) * t229;
t275 = t244 * t226 + t108;
t341 = (t272 * t226 + t229 * t273) * t209 + (t274 * t226 + t229 * t275) * t210;
t323 = (t221 * t145 + (-t226 * t144 + t341) * t229) * t330 + (t144 * t222 + (-t145 * t229 + t341) * t226) * t327;
t74 = -t226 * (rSges(4,1) * t286 + t262) + t229 * (rSges(4,2) * t289 + t261);
t67 = t266 * t324 + t74;
t150 = t246 * t226;
t82 = t226 * t150 + t229 * t151;
t7 = t323 + m(4) * (t67 * t82 + (t104 * t229 + t226 * t349) * t174);
t346 = t7 * qJD(3);
t217 = Icges(5,4) * t227;
t178 = -Icges(5,2) * t224 + t217;
t181 = Icges(5,1) * t224 + t217;
t218 = Icges(3,4) * t228;
t180 = -Icges(3,2) * t225 + t218;
t183 = Icges(3,1) * t225 + t218;
t305 = Icges(4,4) * t209;
t169 = -Icges(4,2) * t210 - t305;
t172 = -Icges(4,1) * t210 + t305;
t249 = (-t170 / 0.2e1 + t244 / 0.2e1) * t210 + (t169 / 0.2e1 - t172 / 0.2e1) * t209;
t106 = Icges(4,5) * t286 - Icges(4,6) * t290 + Icges(4,3) * t229;
t168 = -Icges(4,5) * t210 + Icges(4,6) * t209;
t294 = t168 * t229;
t107 = Icges(4,3) * t226 + t294;
t252 = t111 * t210 + t106;
t83 = t109 * t290;
t259 = t107 * t229 - t83;
t299 = t110 * t210;
t321 = t226 * t107 + t109 * t289;
t322 = t226 * t106 + t108 * t289;
t54 = -t111 * t286 - t259;
t55 = t110 * t285 - t322;
t56 = -t111 * t285 + t321;
t5 = ((t54 - t83 + (t107 - t299) * t229 + t322) * t229 + t321 * t226) * t327 + (t56 * t226 - t229 * t55) * t326 + ((t226 * t252 + t259 + t54 + t55) * t226 + (-t321 + t56 + t226 * (t108 * t209 - t299) + (t252 - t106) * t229) * t229) * t330;
t140 = Icges(3,5) * t226 + t184 * t229;
t268 = -t179 * t229 + t140;
t204 = Icges(3,4) * t282;
t139 = Icges(3,1) * t279 - Icges(3,5) * t229 - t204;
t269 = -Icges(3,2) * t279 + t139 - t204;
t136 = Icges(3,6) * t226 + t180 * t229;
t270 = -t183 * t229 - t136;
t135 = Icges(3,4) * t279 - Icges(3,2) * t282 - Icges(3,6) * t229;
t271 = t183 * t226 + t135;
t342 = (-t268 * t226 + t229 * t269) * t225 + (t270 * t226 + t229 * t271) * t228;
t233 = (-t226 * t89 - t229 * t88) * t174;
t337 = m(4) * (t104 * t150 - t151 * t349 + t233);
t336 = m(4) * (t233 - (-t101 * t226 - t229 * t349) * t246);
t333 = m(4) * (-t150 * t88 + t151 * t89);
t331 = -t226 / 0.2e1;
t280 = t226 * t227;
t129 = Icges(5,5) * t280 - Icges(5,6) * t284 - Icges(5,3) * t229;
t203 = Icges(5,4) * t284;
t137 = Icges(5,1) * t280 - Icges(5,5) * t229 - t203;
t277 = t227 * t229;
t312 = -t226 * t129 - t137 * t277;
t239 = Icges(5,5) * t227 - Icges(5,6) * t224;
t130 = Icges(5,3) * t226 + t229 * t239;
t138 = Icges(5,5) * t226 + t182 * t229;
t311 = t226 * t130 + t138 * t277;
t131 = Icges(3,5) * t279 - Icges(3,6) * t282 - Icges(3,3) * t229;
t310 = -t226 * t131 - t139 * t276;
t242 = Icges(3,5) * t228 - Icges(3,6) * t225;
t132 = Icges(3,3) * t226 + t229 * t242;
t309 = t226 * t132 + t140 * t276;
t133 = Icges(5,4) * t280 - Icges(5,2) * t284 - Icges(5,6) * t229;
t296 = t133 * t224;
t295 = t135 * t225;
t278 = t226 * t229;
t263 = -t174 - t324;
t94 = t138 * t280;
t258 = t130 * t229 - t94;
t95 = t140 * t279;
t257 = t132 * t229 - t95;
t134 = Icges(5,6) * t226 + t178 * t229;
t251 = t224 * t134 - t129;
t250 = t136 * t225 - t131;
t248 = t221 / 0.2e1 + t222 / 0.2e1;
t241 = -Icges(3,5) * t225 - Icges(3,6) * t228;
t238 = -Icges(5,5) * t224 - Icges(5,6) * t227;
t230 = (t169 - t172) * t210 + (t170 - t244) * t209;
t232 = -t5 + (t226 * t168 + t209 * t274 - t210 * t272 + t229 * t230) * t330 + (-t209 * t275 + t210 * t273 + t226 * t230 - t294) * t327;
t231 = -t249 + (t330 + t331) * (t108 * t210 + t110 * t209);
t189 = rSges(3,1) * t228 - rSges(3,2) * t225;
t188 = -rSges(5,2) * t224 + t319;
t155 = t241 * t229;
t154 = t241 * t226;
t153 = t238 * t229;
t152 = t238 * t226;
t105 = t263 * t229;
t103 = t263 * t226;
t73 = -t266 * t325 + t82;
t65 = -t136 * t281 + t309;
t64 = -t135 * t281 - t310;
t63 = -t134 * t283 + t311;
t62 = -t133 * t283 - t312;
t61 = -t136 * t282 - t257;
t59 = -t134 * t284 - t258;
t43 = (t181 / 0.2e1 + t178 / 0.2e1) * t227 - t348;
t42 = t65 * t226 - t229 * t64;
t41 = t63 * t226 - t229 * t62;
t40 = t61 * t226 - t229 * (-(-t139 * t228 + t295) * t226 - t131 * t229);
t39 = t59 * t226 - t229 * (-t226 * (-t137 * t227 + t296) - t129 * t229);
t35 = t249 + t333;
t34 = t336 / 0.2e1;
t33 = t337 / 0.2e1;
t20 = (t183 / 0.2e1 + t180 / 0.2e1) * t228 + t249 - t347;
t19 = (t61 - t95 + (t132 + t295) * t229 + t310) * t229 + t309 * t226;
t18 = (t59 - t94 + (t130 + t296) * t229 + t312) * t229 + t311 * t226;
t17 = (t250 * t229 - t309 + t65) * t229 + (t226 * t250 + t257 + t64) * t226;
t16 = (t229 * t251 - t311 + t63) * t229 + (t226 * t251 + t258 + t62) * t226;
t9 = m(4) * (-t266 * t246 * t174 + t74 * t82) + t323;
t8 = t9 * qJD(3);
t6 = (t41 / 0.2e1 - t18 / 0.2e1) * t229 + (t16 / 0.2e1 + t39 / 0.2e1) * t226;
t4 = t33 - t336 / 0.2e1 + t5;
t3 = t34 - t337 / 0.2e1 + t5;
t2 = t33 + t34 + t232;
t1 = (t42 / 0.2e1 - t19 / 0.2e1) * t229 + (t17 / 0.2e1 + t40 / 0.2e1) * t226 + t5;
t10 = [t20 * qJD(2) + t35 * qJD(3) + t43 * qJD(4), t20 * qJD(1) + t2 * qJD(3) + (t19 * t326 + (t225 * t270 + t228 * t268) * t330 + t232 + m(3) * ((t142 * t229 - t143 * t226) * t189 + (-t165 * t229 + t167 * t226) * t186) + m(4) * (t103 * t89 + t105 * t88 + (-t101 - t104) * t349) + t248 * t242 + (t17 + t40) * t331 + (-t225 * t271 + t228 * t269 + t42) * t327) * qJD(2), t35 * qJD(1) + t2 * qJD(2) + ((t233 - (t150 * t229 - t151 * t226) * t246) * m(4) + t232) * qJD(3), t43 * qJD(1) + (((-t177 * t229 + t138) * t227 + (-t181 * t229 - t134) * t224) * t330 + m(5) * ((-t112 * t229 - t113 * t226) * t188 + (-t164 * t229 + t166 * t226) * t185) + t18 * t326 + t248 * t239 + (t16 + t39) * t331 + ((-Icges(5,2) * t280 + t137 - t203) * t227 + (-t181 * t226 - t133) * t224 + t41) * t327) * qJD(4), 0; t1 * qJD(2) + t4 * qJD(3) + (t231 - (t183 + t180) * t228 / 0.2e1 + t347) * qJD(1), t1 * qJD(1) + (m(3) * ((t226 * t142 + t143 * t229) * (-t226 * t165 - t167 * t229) + t266 * t189 * t186) + (t221 * t155 + (-t226 * t154 + t342) * t229) * t330 + (t154 * t222 + (-t155 * t229 + t342) * t226) * t327 + m(4) * (-t103 * t349 - t104 * t105 + t67 * t73) + t323) * qJD(2) + t346, t4 * qJD(1) + t7 * qJD(2) + t346, 0, 0; (t231 - t333) * qJD(1) + t3 * qJD(2) + t5 * qJD(3), t3 * qJD(1) + ((t74 * t73 - (-t103 * t226 - t105 * t229) * t246) * m(4) + t323) * qJD(2) + t8, qJD(1) * t5 + qJD(2) * t9 + t8, 0, 0; (-(t181 + t178) * t227 / 0.2e1 + t348) * qJD(1) + t6 * qJD(4), 0, 0, t6 * qJD(1) + (m(5) * ((t226 * (rSges(5,1) * t280 - t267) + t229 * (rSges(5,1) * t277 + t260)) * (-t226 * t164 - t166 * t229) + t266 * t188 * t185) + (-t152 * t278 + t221 * t153) * t330 + (t222 * t152 - t153 * t278) * t327) * qJD(4), 0; 0, 0, 0, 0, 0;];
Cq = t10;
