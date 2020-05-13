% Calculate vector of inverse dynamics joint torques for
% fourbarprisDE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% qJDD [1x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% m [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [1x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:45
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbarprisDE2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(1,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE2_invdynJ_fixb_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisDE2_invdynJ_fixb_slag_vp1: qJD has to be [1x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [1 1]), ...
  'fourbarprisDE2_invdynJ_fixb_slag_vp1: qJDD has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisDE2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE2_invdynJ_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisDE2_invdynJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbarprisDE2_invdynJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbarprisDE2_invdynJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:45:03
% EndTime: 2020-05-07 09:45:15
% DurationCPUTime: 5.57s
% Computational Cost: add. (4700->164), mult. (4844->325), div. (327->18), fcn. (240->2), ass. (0->136)
t87 = (qJ(1) + pkin(3));
t185 = -pkin(2) - t87;
t62 = pkin(1) - t185;
t51 = 1 / t62;
t184 = -pkin(2) + t87;
t64 = pkin(1) + t184;
t55 = 1 / t64;
t205 = (t51 * t55);
t65 = pkin(1) + t185;
t224 = 1 / t65;
t63 = pkin(1) - t184;
t225 = 1 / t63;
t231 = t224 * t225;
t168 = t205 * t231;
t108 = 0.1e1 / pkin(1) ^ 2;
t76 = (t87 ^ 2);
t196 = (t108 / t76 ^ 2);
t77 = 0.1e1 / t87;
t146 = (t77 * t168 * t196);
t165 = (t55 * t224 * t196);
t54 = 1 / (t63 ^ 2);
t56 = 1 / (t64 ^ 2);
t58 = 1 / (t65 ^ 2);
t241 = 4 * t146 + (-t54 * t165 + (t224 * t56 - t55 * t58) * t225 * t196) * t51;
t204 = t64 * t65;
t182 = t63 * t204;
t167 = t62 * t182;
t109 = sqrt(-t167);
t20 = 0.1e1 / t109;
t100 = (qJ(1) ^ 2);
t103 = (pkin(3) ^ 2);
t230 = (t103 + t100);
t155 = -2 * pkin(3) * qJ(1) - t230;
t105 = pkin(2) ^ 2;
t106 = (pkin(1) ^ 2);
t191 = (t106 - t105);
t33 = t155 + t191;
t208 = t20 * t33;
t116 = (t100 ^ 2);
t237 = 2 * qJ(1);
t102 = pkin(3) * t103;
t113 = t103 ^ 2;
t216 = 3 * t103;
t72 = (-t105 - t106);
t28 = 2 * t102 * rSges(3,3) + t72 * t216 + 3 * t113;
t29 = 3 * pkin(3) * t72 + rSges(3,3) * t216 + 7 * t102;
t90 = pkin(3) * rSges(3,3);
t35 = 8 * t103 + 2 * t90 + t72;
t93 = rSges(3,3) / 0.2e1;
t68 = t93 + 0.9e1 / 0.2e1 * pkin(3);
t88 = 3 * t100;
t99 = qJ(1) * t100;
t239 = (t29 * t237) + (t35 * t88) + 0.4e1 * t68 * t99 + (5 * t116) + t28;
t47 = (m(2) * rSges(2,2) - m(3) * rSges(3,1));
t48 = m(2) * rSges(2,1) + m(3) * rSges(3,3);
t24 = -t47 * g(1) + g(2) * t48;
t233 = t77 * t208;
t107 = 0.1e1 / pkin(1);
t232 = t20 * t107 ^ 2;
t61 = t90 + 2 * t103;
t69 = 0.5e1 / 0.2e1 * pkin(3) + t93;
t228 = t237 * t69 + t61 + t88;
t45 = (t87 * qJD(1));
t223 = -2 * t45;
t222 = 2 * t45;
t221 = -2 * t87;
t220 = 2 * t99;
t217 = 5 * t100;
t214 = m(3) * g(2);
t213 = t20 / 0.2e1;
t212 = -t33 / 0.2e1;
t211 = t33 / 0.2e1;
t13 = -t182 + (t204 + (t64 - t65) * t63) * t62;
t161 = 0.1e1 / t167 * t233 / 0.2e1;
t78 = 1 / t87 ^ 2;
t206 = t33 * t78;
t210 = (t13 * t161 + (t221 * t77 - t206) * t20) * qJD(1) * t77;
t209 = t113 / 0.2e1;
t121 = t87 * t76;
t79 = 1 / t121;
t96 = qJD(1) ^ 2;
t202 = t79 * t96;
t91 = pkin(1) + pkin(2);
t92 = pkin(1) - pkin(2);
t201 = t91 ^ 2 * t92 ^ 2;
t200 = t87 * t96;
t198 = rSges(3,1) * t109;
t197 = t106 * t24;
t195 = qJ(1) * t103;
t12 = t13 * qJD(1);
t194 = qJD(1) * t12;
t193 = t33 * qJD(1);
t192 = t45 * qJD(1);
t190 = 2 * m(4) * t121;
t181 = t20 * t193;
t180 = t78 * t193;
t179 = -rSges(2,1) * t20 / 0.2e1;
t178 = rSges(2,2) * t213;
t170 = t106 / 0.2e1 - t105 / 0.2e1;
t18 = t99 + t69 * t100 + (t61 * qJ(1)) + t102 / 0.2e1 + t170 * pkin(3) + (t103 / 0.2e1 - t170) * rSges(3,3);
t177 = t18 * t213;
t176 = -t208 / 0.2e1;
t175 = -t208 / 0.4e1;
t32 = -t155 + t191;
t174 = t32 * t211;
t173 = t211 - t32 / 0.2e1;
t171 = (rSges(4,1) * g(2) - rSges(4,2) * g(1)) * t190;
t162 = t12 * t177;
t160 = t12 * t175;
t159 = qJD(1) * t176;
t158 = t181 / 0.2e1;
t156 = t239 * qJD(1);
t154 = t225 * t165;
t11 = (t29 * t100) + (t72 * t102) + t68 * t116 + (t35 * t99) + ((t28 + t116) * qJ(1)) + (-t201 / 0.2e1 + t209) * rSges(3,3) + (t201 / 0.2e1 + t209) * pkin(3);
t5 = t212 * t198 + t11;
t152 = t173 * t222;
t151 = t78 * t159;
t150 = t78 * t158;
t141 = qJD(1) * t154;
t10 = rSges(3,1) * t174 + t18 * t109;
t31 = t33 ^ 2;
t52 = 0.1e1 / t62 ^ 2;
t138 = t31 * t52 * t54 * t56 * t58 * t194;
t137 = t51 * t141;
t1 = t161 * t194 + (-0.2e1 * t77 * t192 + (qJDD(1) * t77 - t78 * t96) * t33) * t20;
t129 = t1 * t77 / 0.2e1 + t176 * t202;
t110 = pkin(2) * t105;
t83 = t102 * t214;
t22 = t228 * qJD(1);
t17 = -rSges(2,1) * t109 - rSges(2,2) * t32;
t16 = -(rSges(2,1) * t32) + rSges(2,2) * t109;
t9 = t10 * qJD(1);
t7 = rSges(2,2) * t223 + t12 * t179;
t6 = (rSges(2,1) * t223) + t12 * t178;
t4 = t5 * qJD(1);
t2 = [((((rSges(4,1) * g(1) + rSges(4,2) * g(2)) * t190) + (-pkin(3) * t110 + (t102 + 4 * t195 + (t217 + t106) * pkin(3) + t220) * pkin(2)) * g(1) * m(3) + (t110 + (pkin(3) * t237 - t106 + t230) * pkin(2)) * (t48 * g(1) + t47 * g(2))) * t109 + 0.2e1 * (-t197 + t83 + (3 * t195 + (t88 + t106) * pkin(3) + t99) * t214) * t110 + ((t83 + (4 * qJ(1) * t214 + t24) * t103 + (t24 * t237 + (t217 - t106) * t214) * pkin(3) + t214 * t220 + t24 * t100 + t197) * pkin(2) + t171) * (pkin(1) + t87) * (pkin(1) - t87) + ((-pkin(3) * t214 + t24) * t110 + t171) * t105) * t107 * t78 / pkin(2) * t213 - m(2) * (-t17 * ((rSges(2,2) * t221 + t13 * t179) * t151 + (-t210 / 0.2e1 + t79 * t158) * t17) + (((rSges(2,1) * t221) + t13 * t178) * t150 + (t210 / 0.2e1 + t79 * t159) * t16) * t16) * t180 * t232 / 0.2e1 - m(3) * (-0.2e1 * (t9 * (0.2e1 * t173 * t87 * rSges(3,1) + t109 * t228 + t13 * t177) + t4 * ((t87 * t109 + t13 * t175) * rSges(3,1) + t239)) * t137 + (t4 ^ 2 + t9 ^ 2) * (t52 * t154 + t241)) / 0.2e1 + m(2) * ((-(-t129 * t17 + t151 * t7) * t17 + (t129 * t16 + t150 * t6) * t16) * t206 * t232 - (t6 * t16 + t7 * t17) * t31 * t137 / 0.2e1 + (t16 ^ 2 + t17 ^ 2) * (t33 * t45 * t137 + t31 * t96 * t146 - t138 * t196 / 0.4e1)) / 0.2e1 + m(3) * ((-((qJD(1) * t22 + qJDD(1) * t18) * t109 + qJD(1) * t162 + (qJD(1) * t152 + qJDD(1) * t174) * rSges(3,1)) * t10 - t9 * (rSges(3,1) * t152 + t22 * t109 + t162) - ((qJDD(1) * t212 + t192) * t198 + t11 * qJDD(1) + (t160 * rSges(3,1) + t156) * qJD(1)) * t5 - t4 * ((t45 * t109 + t160) * rSges(3,1) + t156)) * t51 * t154 + (t10 * t9 + t4 * t5) * (qJD(1) * t241 + t52 * t141)) + (Icges(3,2) + Icges(2,3)) * (t1 * t233 - t181 * t210 - t78 * t138 / 0.2e1 + (t180 * t222 + t31 * t202) * t168) + ((0.4e1 * t200 + (-0.32e2 * qJDD(1) * t76 - 0.64e2 * t200) * t106 * t105 * t108 / pkin(2) ^ 2 / 0.8e1) * t168 + (-(2 * (t224 * t54 + t225 * t58) * t205) + 0.2e1 * ((t51 * t56) + t52 * t55) * t231) * t76 * t96) * ((rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * m(4) + Icges(4,3));];
tau = t2;
