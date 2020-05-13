% Calculate vector of inverse dynamics joint torques for
% picker2Dm2OL
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
% tau [12x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 23:20
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = picker2Dm2OL_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(12,1),zeros(3,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2OL_invdynJ_fixb_slag_vp1: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm2OL_invdynJ_fixb_slag_vp1: qJD has to be [12x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [12 1]), ...
  'picker2Dm2OL_invdynJ_fixb_slag_vp1: qJDD has to be [12x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2OL_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2OL_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2OL_invdynJ_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm2OL_invdynJ_fixb_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'picker2Dm2OL_invdynJ_fixb_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 23:18:54
% EndTime: 2020-05-09 23:19:09
% DurationCPUTime: 3.29s
% Computational Cost: add. (4248->257), mult. (2706->280), div. (0->0), fcn. (1280->20), ass. (0->170)
t195 = sin(qJ(1));
t185 = t195 * pkin(1);
t197 = cos(qJ(1));
t200 = qJD(1) ^ 2;
t240 = pkin(1) * qJDD(1);
t285 = t200 * t185 - t197 * t240 - g(2);
t257 = pkin(1) * t197;
t282 = t195 * t240 + t200 * t257 - g(1);
t188 = qJDD(1) + qJDD(2);
t169 = qJDD(6) + t188;
t190 = qJD(1) + qJD(2);
t172 = qJD(6) + t190;
t193 = qJ(1) + qJ(2);
t181 = qJ(6) + t193;
t153 = sin(t181);
t156 = cos(t181);
t90 = t156 * rSges(7,1) - rSges(7,2) * t153;
t51 = t172 * t90;
t87 = rSges(7,1) * t153 + rSges(7,2) * t156;
t289 = -t169 * t87 - t172 * t51 + t282;
t242 = t172 * t87;
t288 = t169 * t90 - t172 * t242 + t285;
t178 = sin(t193);
t180 = cos(t193);
t186 = t190 ^ 2;
t274 = t178 * t186 - t180 * t188;
t287 = t274 * pkin(3) + t285;
t286 = t274 * pkin(2) + t285;
t273 = t178 * t188 + t180 * t186;
t284 = t273 * pkin(2) + t282;
t283 = t273 * pkin(3) + t282;
t246 = pkin(1) * qJD(1);
t219 = t197 * t246;
t39 = -t51 + t219;
t189 = qJD(1) + qJD(8);
t192 = qJ(1) + qJ(8);
t177 = sin(t192);
t179 = cos(t192);
t99 = t179 * rSges(9,1) - rSges(9,2) * t177;
t72 = t189 * t99;
t60 = -t72 + t219;
t227 = t180 * t190;
t206 = pkin(3) * t227 + t219;
t182 = qJ(4) + t193;
t157 = cos(t182);
t173 = qJD(4) + t190;
t234 = t157 * t173;
t147 = qJD(10) + t173;
t159 = qJ(10) + t182;
t141 = sin(t159);
t142 = cos(t159);
t214 = t142 * rSges(11,1) - rSges(11,2) * t141;
t275 = t147 * t214;
t18 = pkin(4) * t234 + t206 - t275;
t166 = t195 * t246;
t230 = t178 * t190;
t226 = pkin(3) * t230 + t166;
t154 = sin(t182);
t236 = t154 * t173;
t77 = rSges(11,1) * t141 + rSges(11,2) * t142;
t46 = t147 * t77;
t17 = pkin(4) * t236 + t226 - t46;
t100 = -rSges(3,1) * t180 + t178 * rSges(3,2);
t264 = t190 * t100;
t61 = t219 - t264;
t170 = qJDD(4) + t188;
t53 = rSges(5,1) * t234 - rSges(5,2) * t236;
t88 = rSges(5,1) * t154 + rSges(5,2) * t157;
t272 = t170 * t88 + t173 * t53 + t283;
t171 = qJDD(3) + t188;
t174 = qJD(3) + t190;
t183 = qJ(3) + t193;
t158 = cos(t183);
t233 = t158 * t174;
t155 = sin(t183);
t235 = t155 * t174;
t55 = rSges(4,1) * t233 - rSges(4,2) * t235;
t89 = rSges(4,1) * t155 + rSges(4,2) * t158;
t271 = -t171 * t89 - t174 * t55 + t284;
t217 = -rSges(5,1) * t157 + t154 * rSges(5,2);
t52 = rSges(5,1) * t236 + rSges(5,2) * t234;
t270 = t170 * t217 + t173 * t52 + t287;
t143 = qJDD(10) + t170;
t167 = t173 ^ 2;
t269 = -t143 * t77 - t147 * t275 + (t154 * t170 + t157 * t167) * pkin(4) + t283;
t268 = t143 * t214 - t147 * t46 + (t154 * t167 - t157 * t170) * pkin(4) + t287;
t146 = qJDD(9) + t171;
t148 = qJD(9) + t174;
t168 = t174 ^ 2;
t163 = qJ(9) + t183;
t145 = cos(t163);
t237 = t145 * t148;
t144 = sin(t163);
t238 = t144 * t148;
t37 = rSges(10,1) * t237 - rSges(10,2) * t238;
t78 = t144 * rSges(10,1) + t145 * rSges(10,2);
t267 = t146 * t78 + t148 * t37 + (-t155 * t171 - t158 * t168) * pkin(6) + t284;
t216 = -rSges(10,1) * t145 + t144 * rSges(10,2);
t36 = rSges(10,1) * t238 + rSges(10,2) * t237;
t266 = t146 * t216 + t148 * t36 + (-t155 * t168 + t158 * t171) * pkin(6) + t286;
t187 = qJDD(1) + qJDD(8);
t97 = rSges(9,1) * t177 + rSges(9,2) * t179;
t281 = -t187 * t97 - t189 * t72 + t282;
t215 = t158 * rSges(4,1) - rSges(4,2) * t155;
t69 = t174 * t89;
t265 = t171 * t215 - t174 * t69 + t286;
t241 = t189 * t97;
t280 = t187 * t99 - t189 * t241 + t285;
t98 = rSges(3,1) * t178 + rSges(3,2) * t180;
t279 = t188 * t98 - t190 * t264 + t282;
t73 = rSges(3,1) * t230 + rSges(3,2) * t227;
t278 = t100 * t188 + t190 * t73 + t285;
t207 = pkin(2) * t227 + t219;
t70 = t174 * t215;
t24 = t207 - t70;
t276 = t24 * t69;
t225 = pkin(2) * t230 + t166;
t22 = t225 - t69;
t58 = t166 - t241;
t38 = t166 - t242;
t209 = pkin(6) * t235 - t225;
t56 = t148 * t78;
t19 = -t209 + t56;
t67 = t173 * t88;
t21 = t226 + t67;
t263 = -t22 * t55 + t276;
t68 = t173 * t217;
t23 = t206 - t68;
t262 = t21 * t53 - t23 * t52;
t205 = -pkin(6) * t233 + t207;
t57 = t148 * t216;
t20 = t205 - t57;
t261 = t19 * t37 - t20 * t36;
t259 = (t242 * t39 - t38 * t51 + (t38 * t172 + t288) * t90 + (-t39 * t172 - t289) * t87) * m(7);
t256 = pkin(2) * t180;
t255 = pkin(3) * t180;
t127 = Icges(11,3) * t143;
t224 = Icges(5,3) * t170 + t127;
t129 = Icges(10,3) * t146;
t223 = Icges(4,3) * t171 + t129;
t59 = t190 * t98 + t166;
t48 = -pkin(6) * t155 + t78;
t49 = pkin(6) * t158 + t216;
t210 = rSges(2,1) * t197 - rSges(2,2) * t195;
t124 = rSges(2,1) * t195 + rSges(2,2) * t197;
t191 = pkin(8) + qJ(5);
t175 = sin(t191);
t176 = cos(t191);
t96 = rSges(6,1) * t176 - rSges(6,2) * t175;
t95 = rSges(6,1) * t175 + rSges(6,2) * t176;
t194 = sin(qJ(7));
t196 = cos(qJ(7));
t125 = rSges(8,1) * t196 - t194 * rSges(8,2);
t123 = rSges(8,1) * t194 + rSges(8,2) * t196;
t152 = pkin(2) * t178;
t34 = t152 + t48;
t65 = t215 - t256;
t63 = t217 - t255;
t41 = -pkin(4) * t157 + t214;
t64 = t152 - t89;
t151 = pkin(3) * t178;
t62 = t151 + t88;
t136 = Icges(7,3) * t169;
t208 = Icges(3,3) * t188 + t136 + t223 + t224;
t40 = pkin(4) * t154 - t77;
t35 = t49 - t256;
t30 = t151 + t40;
t31 = t41 - t255;
t161 = Icges(9,3) * t187;
t1 = [Icges(2,3) * qJDD(1) + t161 + t208 + (t19 * (t205 + t37) - t20 * (-t209 + t36) + t266 * (t35 - t257) + t267 * (t185 + t34)) * m(10) + (t280 * (t99 - t257) + t281 * (t185 - t97)) * m(9) + (t288 * (t90 - t257) + t289 * (t185 - t87)) * m(7) + (t21 * (t206 + t53) - t23 * (t226 + t52) + t270 * (t63 - t257) + t272 * (t185 + t62)) * m(5) + (t265 * (t65 - t257) + t271 * (t185 + t64) + (t207 - t55 - t24) * t22) * m(4) + (t278 * (t100 - t257) + t279 * (t185 + t98) + (t59 - t166 - t73) * t61) * m(3) + (t268 * (t31 - t257) + t269 * (t185 + t30)) * m(11) + ((qJDD(1) * t210 + g(2)) * t210 + (qJDD(1) * t124 - g(1)) * t124) * m(2); t208 + (t19 * t57 + t20 * t56 + t266 * t35 + t267 * t34 + t261) * m(10) + (t21 * t68 + t23 * t67 + t270 * t63 + t272 * t62 + t262) * m(5) + (t22 * t70 + t265 * t65 + t271 * t64 + t263 - t276) * m(4) + (-t59 * t264 - t61 * t73 + (t190 * t61 + t279) * t98 + (t190 * t59 + t278) * t100) * m(3) + (t268 * t31 + t269 * t30) * m(11) + t259; t223 + (t266 * t49 + t267 * t48 + (-t36 + t56) * t20 + (t37 + t57) * t19) * m(10) + ((-t174 * t24 - t271) * t89 + (t174 * t22 + t265) * t215 + t263) * m(4); t224 + ((t173 * t23 + t272) * t88 + (t173 * t21 + t270) * t217 + t262) * m(5) + (t268 * t41 + t269 * t40) * m(11); Icges(6,3) * qJDD(5) + ((t95 ^ 2 + t96 ^ 2) * qJDD(5) + g(1) * t95 - g(2) * t96) * m(6); t136 + t259; Icges(8,3) * qJDD(7) + ((t123 ^ 2 + t125 ^ 2) * qJDD(7) - g(1) * t125 - g(2) * t123) * m(8); t161 + (-t58 * t72 + t60 * t241 + (t189 * t58 + t280) * t99 + (-t189 * t60 - t281) * t97) * m(9); t129 + ((t148 * t20 + t267) * t78 + (t148 * t19 + t266) * t216 + t261) * m(10); t127 + ((-t147 * t18 - t269) * t77 + (t147 * t17 + t268) * t214 - t17 * t275 + t18 * t46) * m(11); 0; 0;];
tau = t1;
