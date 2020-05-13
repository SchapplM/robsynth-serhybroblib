% Calculate vector of centrifugal and Coriolis load on the joints for
% picker2Dm1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% qJD [12x1]
%   Generalized joint velocities
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
% tauc [12x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:46
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = picker2Dm1OL_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1OL_coriolisvecJ_fixb_slag_vp1: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm1OL_coriolisvecJ_fixb_slag_vp1: qJD has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1OL_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1OL_coriolisvecJ_fixb_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm1OL_coriolisvecJ_fixb_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'picker2Dm1OL_coriolisvecJ_fixb_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 05:44:58
% EndTime: 2020-05-11 05:45:06
% DurationCPUTime: 2.58s
% Computational Cost: add. (3217->207), mult. (2129->221), div. (0->0), fcn. (948->16), ass. (0->151)
t138 = qJ(1) + qJ(2);
t129 = qJ(4) + t138;
t107 = sin(t129);
t135 = qJD(1) + qJD(2);
t120 = qJD(4) + t135;
t200 = t107 * t120;
t140 = sin(qJ(1));
t217 = pkin(1) * qJD(1);
t116 = t140 * t217;
t125 = sin(t138);
t195 = t125 * t135;
t86 = pkin(3) * t195;
t204 = t116 + t86;
t102 = qJD(10) + t120;
t112 = qJ(10) + t129;
t98 = sin(t112);
t99 = cos(t112);
t54 = rSges(11,1) * t98 + rSges(11,2) * t99;
t31 = t102 * t54;
t7 = pkin(4) * t200 + t204 - t31;
t142 = cos(qJ(1));
t185 = t142 * t217;
t127 = cos(t138);
t193 = t127 * t135;
t190 = pkin(3) * t193;
t110 = cos(t129);
t198 = t110 * t120;
t179 = rSges(11,1) * t99 - rSges(11,2) * t98;
t22 = t102 * t179;
t8 = pkin(4) * t198 + t185 + t190 - t22;
t121 = qJD(3) + t135;
t130 = qJ(3) + t138;
t108 = sin(t130);
t111 = cos(t130);
t61 = rSges(4,1) * t108 + rSges(4,2) * t111;
t48 = t121 * t61;
t87 = pkin(2) * t195;
t18 = t116 + t87 - t48;
t119 = qJD(6) + t135;
t128 = qJ(6) + t138;
t106 = sin(t128);
t109 = cos(t128);
t174 = rSges(7,1) * t109 - rSges(7,2) * t106;
t34 = t119 * t174;
t30 = -t34 + t185;
t199 = t108 * t121;
t188 = pkin(6) * t199;
t113 = qJ(9) + t130;
t101 = cos(t113);
t103 = qJD(9) + t121;
t201 = t101 * t103;
t100 = sin(t113);
t202 = t100 * t103;
t23 = rSges(10,1) * t202 + rSges(10,2) * t201;
t55 = t100 * rSges(10,1) + t101 * rSges(10,2);
t39 = t103 * t55;
t228 = t39 + t188 - t23;
t134 = qJD(1) + qJD(8);
t137 = qJ(1) + qJ(8);
t124 = sin(t137);
t126 = cos(t137);
t173 = rSges(9,1) * t126 - rSges(9,2) * t124;
t51 = t134 * t173;
t43 = -t51 + t185;
t68 = rSges(9,1) * t124 + rSges(9,2) * t126;
t209 = t134 * t68;
t41 = t116 - t209;
t59 = rSges(7,1) * t106 + rSges(7,2) * t109;
t210 = t119 * t59;
t29 = t116 - t210;
t227 = -t87 + t228;
t35 = rSges(5,1) * t200 + rSges(5,2) * t198;
t60 = t107 * rSges(5,1) + t110 * rSges(5,2);
t46 = t120 * t60;
t226 = t46 - t86 - t35;
t178 = rSges(3,1) * t127 - rSges(3,2) * t125;
t53 = t178 * t135;
t44 = t53 + t185;
t132 = t140 * pkin(1);
t145 = qJD(1) ^ 2;
t114 = t145 * t132;
t15 = -t119 * t210 + t114;
t224 = pkin(1) * t142;
t115 = t145 * t224;
t16 = -t119 * t34 + t115;
t225 = (t15 * t174 - t16 * t59 - t29 * t34 + t30 * t210 - (-t174 * t29 + t30 * t59) * t119) * m(7);
t223 = pkin(2) * t127;
t222 = pkin(3) * t127;
t221 = pkin(4) * t120 ^ 2;
t220 = pkin(6) * t121 ^ 2;
t172 = -rSges(10,1) * t101 + rSges(10,2) * t100;
t197 = t111 * t121;
t77 = pkin(6) * t197;
t219 = -t103 * t172 - t77;
t52 = rSges(3,1) * t195 + rSges(3,2) * t193;
t133 = t135 ^ 2;
t196 = t125 * t133;
t208 = pkin(3) * t196 + t114;
t207 = pkin(2) * t196 + t114;
t194 = t127 * t133;
t206 = pkin(3) * t194 + t115;
t205 = pkin(2) * t194 + t115;
t191 = pkin(2) * t193;
t104 = pkin(3) * t125;
t187 = t104 + t60;
t69 = rSges(3,1) * t125 + rSges(3,2) * t127;
t42 = t135 * t69 + t116;
t184 = -pkin(6) * t108 + t55;
t176 = -rSges(5,1) * t110 + rSges(5,2) * t107;
t183 = -t120 * t176 + t190;
t181 = t87 - t188;
t38 = rSges(4,1) * t197 - rSges(4,2) * t199;
t36 = rSges(5,1) * t198 - rSges(5,2) * t200;
t24 = rSges(10,1) * t201 - rSges(10,2) * t202;
t177 = rSges(4,1) * t111 - rSges(4,2) * t108;
t105 = pkin(2) * t125;
t170 = t105 + t184;
t169 = -t116 - t181;
t168 = t24 - t77;
t167 = pkin(4) * t107 - t54;
t166 = pkin(6) * t111 + t172;
t165 = t105 - t61;
t161 = t185 + t191;
t160 = t104 + t167;
t159 = t177 - t223;
t158 = t176 - t222;
t157 = -pkin(4) * t110 + t179;
t155 = -t38 + t191;
t154 = t36 + t190;
t150 = t166 - t223;
t149 = t168 + t191;
t147 = t157 - t222;
t49 = t121 * t177;
t28 = t135 * t53 + t115;
t27 = -t134 * t51 + t115;
t26 = t135 * t52 + t114;
t25 = -t134 * t209 + t114;
t20 = t161 - t49;
t19 = t183 + t185;
t17 = t204 + t46;
t14 = -t121 * t38 + t205;
t13 = t120 * t36 + t206;
t12 = -t121 * t48 + t207;
t11 = t120 * t35 + t208;
t10 = t161 + t219;
t9 = -t169 + t39;
t6 = t103 * t24 - t111 * t220 + t205;
t5 = t103 * t23 - t108 * t220 + t207;
t4 = -t102 * t22 + t110 * t221 + t206;
t3 = -t102 * t31 + t107 * t221 + t208;
t1 = [m(9) * (t27 * (t132 - t68) + t25 * (t173 - t224)) + m(7) * (t16 * (t132 - t59) + t15 * (t174 - t224)) + m(3) * (t28 * (t132 + t69) + t26 * (-t178 - t224) + (t42 - t116 - t52) * t44) + (t6 * (t132 + t170) + t9 * (t149 + t185) + t5 * (t150 - t224) + (-t169 - t9 - t116 + t227) * t10) * m(10) + (t13 * (t132 + t187) + t17 * (t154 + t185) + t11 * (t158 - t224) + (-t17 + t204 - t116 + t226) * t19) * m(5) + (t14 * (t132 + t165) + t12 * (t159 - t224) + (-t20 + t155 + t185) * t18) * m(4) + (t4 * (t132 + t160) + t3 * (t147 - t224)) * m(11); t225 + (t5 * t150 + t6 * t170 + (-t191 - t219 + t149) * t9 + (t181 + t227) * t10) * m(10) + (t11 * t158 + t13 * t187 + (t86 + t226) * t19 + (-t183 + t154) * t17) * m(5) + (t12 * t159 + t14 * t165 + (t49 - t191 + t155) * t18) * m(4) + (-(t178 * t42 - t44 * t69) * t135 - t26 * t178 + t28 * t69 + t42 * t53 - t44 * t52) * m(3) + (t3 * t147 + t4 * t160) * m(11); (t5 * t166 + t6 * t184 + (t168 - t219) * t9 + (-t188 + t228) * t10) * m(10) + (t12 * t177 - t14 * t61 - t18 * t38 + t20 * t48 - (-t177 * t18 + t20 * t61) * t121) * m(4); (t11 * t176 + t13 * t60 + t17 * t36 - t19 * t35 - (-t17 * t176 - t19 * t60) * t120) * m(5) + (t3 * t157 + t4 * t167) * m(11); 0; t225; 0; (t25 * t173 - t27 * t68 - t41 * t51 + t43 * t209 - (-t173 * t41 + t43 * t68) * t134) * m(9); (-t10 * t23 + t5 * t172 + t9 * t24 + t6 * t55 - (-t10 * t55 - t172 * t9) * t103) * m(10); (t3 * t179 + t8 * t31 - t7 * t22 - t4 * t54 - (-t179 * t7 + t8 * t54) * t102) * m(11); 0; 0;];
tauc = t1(:);
