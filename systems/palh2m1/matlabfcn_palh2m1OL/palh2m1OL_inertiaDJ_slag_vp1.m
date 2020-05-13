% Calculate time derivative of joint inertia matrix for
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
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m1OL_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1OL_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'palh2m1OL_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'palh2m1OL_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 00:10:43
% EndTime: 2020-05-03 00:10:53
% DurationCPUTime: 2.43s
% Computational Cost: add. (931->264), mult. (1205->379), div. (0->0), fcn. (231->68), ass. (0->190)
t132 = rSges(6,2) * pkin(4);
t130 = rSges(6,3) + pkin(6);
t98 = t130 * rSges(6,1);
t22 = (t98 - t132) * m(6) - Icges(6,5);
t197 = qJ(5) + qJ(3);
t91 = qJ(4) + t197;
t76 = qJ(2) + t91;
t46 = cos(t76);
t101 = qJD(4) + qJD(5);
t82 = qJD(3) + t101;
t55 = qJD(2) + t82;
t257 = t22 * t55 * t46;
t122 = sin(qJ(5));
t126 = cos(qJ(5));
t135 = pkin(4) * m(6);
t161 = qJD(5) * (t122 * rSges(6,1) + t126 * rSges(6,2)) * t135;
t136 = 0.2e1 * qJ(5);
t227 = rSges(6,2) * m(6);
t65 = rSges(6,1) * t227 - Icges(6,4);
t165 = qJD(5) * t65 * cos(t136);
t141 = rSges(6,2) ^ 2;
t142 = rSges(6,1) ^ 2;
t27 = (-t141 + t142) * m(6) - Icges(6,1) + Icges(6,2);
t252 = t27 * qJD(5) * sin(t136);
t256 = -t252 - 0.2e1 * t161 - 0.2e1 * t165;
t236 = -t27 / 0.4e1;
t103 = qJD(3) + qJD(4);
t86 = qJD(2) + t103;
t54 = 0.2e1 * qJD(5) + t86;
t57 = -0.2e1 * qJD(5) + t86;
t113 = qJ(3) + qJ(4);
t97 = qJ(2) + t113;
t69 = t136 + t97;
t78 = -0.2e1 * qJ(5) + t97;
t255 = sin(t78) * t57 * t236 + t27 * t54 * sin(t69) / 0.4e1;
t124 = sin(qJ(3));
t128 = cos(qJ(3));
t233 = m(4) * rSges(4,2);
t232 = rSges(4,1) * m(4);
t234 = m(6) + m(5);
t31 = t234 * pkin(3) + t232;
t146 = (-t31 * t124 - t128 * t233) * qJD(3);
t254 = pkin(2) * t146;
t253 = t257 / 0.2e1;
t63 = m(5) * rSges(5,1) + t135;
t10 = (rSges(6,1) * t126 - rSges(6,2) * t122) * m(6) + t63;
t123 = sin(qJ(4));
t219 = t130 * m(6);
t228 = rSges(5,2) * m(5);
t39 = -t219 + t228;
t201 = t39 * t103;
t189 = qJD(5) * t126;
t192 = qJD(5) * t122;
t231 = rSges(6,1) * m(6);
t206 = t189 * t227 + t192 * t231;
t6 = t201 + t206;
t2 = pkin(2) * t6 * t124;
t203 = qJD(4) * pkin(3);
t199 = t103 * t128;
t38 = pkin(2) * t199;
t251 = (t2 - t10 * (t38 + t203)) * t123;
t104 = qJD(2) + qJD(5);
t224 = pkin(3) * (qJD(3) + t104);
t96 = qJ(2) + t197;
t59 = sin(t96);
t61 = cos(t96);
t250 = (t59 * rSges(6,1) + t61 * rSges(6,2)) * m(6) * t224 / 0.2e1;
t105 = qJD(2) - qJD(5);
t170 = t227 / 0.2e1;
t156 = pkin(3) * t170;
t171 = -t231 / 0.2e1;
t157 = pkin(3) * t171;
t198 = -qJ(5) + qJ(3);
t94 = qJ(2) + t198;
t249 = (cos(t94) * t156 + sin(t94) * t157) * (qJD(3) + t105);
t205 = rSges(6,2) * t130;
t40 = m(6) * t205 - Icges(6,6);
t41 = m(6) * t98 - Icges(6,5);
t246 = (t40 * t122 - t41 * t126) * qJD(5);
t21 = (t98 + t132) * m(6) - Icges(6,5);
t102 = qJD(4) - qJD(5);
t83 = qJD(3) + t102;
t56 = qJD(2) + t83;
t235 = t56 / 0.2e1;
t230 = rSges(6,1) * pkin(4);
t222 = (-t205 + t230) * m(6);
t24 = Icges(6,6) + t222;
t92 = qJ(4) + t198;
t77 = qJ(2) + t92;
t43 = sin(t77);
t47 = cos(t77);
t245 = t21 * t47 * t235 - t24 * t56 * t43 / 0.2e1;
t244 = -0.2e1 * pkin(2);
t243 = -0.2e1 * pkin(3);
t242 = -0.2e1 * t86;
t106 = qJD(2) + qJD(3);
t241 = -0.2e1 * t106;
t240 = -0.2e1 * qJD(2);
t239 = m(6) * pkin(2);
t238 = pkin(1) * m(6);
t237 = pkin(3) * m(6);
t229 = rSges(3,2) * m(3);
t226 = rSges(4,3) * m(4);
t225 = rSges(5,3) * m(5);
t223 = (t205 + t230) * m(6);
t221 = pkin(2) * t104;
t127 = cos(qJ(4));
t200 = t103 * t124;
t148 = -t10 * t200 - t128 * t6;
t9 = qJD(4) * t39 + t206;
t220 = (t148 * pkin(2) - pkin(3) * t9) * t127;
t116 = qJ(2) + qJ(3);
t90 = cos(t116);
t218 = (rSges(4,1) * t226 + pkin(3) * t225 - Icges(4,5)) * t90;
t216 = t54 * cos(t69);
t213 = t57 * cos(t78);
t60 = sin(t97);
t212 = t60 * (rSges(5,2) * t225 - Icges(5,6));
t62 = cos(t97);
t211 = (rSges(5,1) * t225 - Icges(5,5)) * t62;
t88 = sin(t116);
t210 = t88 * (rSges(4,2) * t226 - Icges(4,6));
t209 = t9 * t127;
t208 = m(4) + t234;
t139 = 0.2e1 * qJ(2);
t196 = 0.2e1 * qJ(3) + t139;
t114 = t139 + qJ(3);
t195 = qJD(4) * t123;
t194 = qJD(4) * t127;
t191 = qJD(5) * t123;
t190 = qJD(5) * t124;
t188 = qJD(5) * t128;
t187 = pkin(1) * t242;
t100 = qJD(2) + qJD(3) / 0.2e1;
t120 = qJD(4) / 0.2e1;
t58 = t120 + t100;
t186 = t58 * t244;
t79 = t120 + t106;
t185 = t79 * t243;
t184 = t55 * t238;
t183 = t56 * t238;
t118 = qJD(5) / 0.2e1;
t182 = (t118 + t58) * t239;
t119 = -qJD(5) / 0.2e1;
t181 = (t119 + t58) * t239;
t180 = (t118 + t79) * t237;
t179 = (t119 + t79) * t237;
t176 = t123 * t237;
t23 = -Icges(6,6) + t223;
t42 = sin(t76);
t173 = t23 * t55 * t42;
t169 = -t224 / 0.2e1;
t168 = pkin(3) * t194;
t167 = -t221 / 0.2e1;
t163 = qJ(5) + t196;
t162 = -qJ(5) + t196;
t160 = pkin(2) * t171;
t159 = pkin(2) * t170;
t158 = -t82 * t239 / 0.2e1;
t110 = -qJ(5) + qJ(2);
t154 = t249 + (cos(t110) * t159 + sin(t110) * t160) * t105;
t153 = m(6) * t168;
t152 = -t101 * t237 / 0.2e1;
t147 = (-t10 * t195 - t209) * pkin(3);
t145 = t173 / 0.2e1 + (t213 / 0.2e1 + t216 / 0.2e1) * t65 + (t211 - t212) * t86 + t245 + t255;
t140 = 0.2e1 * Icges(6,6);
t137 = 0.2e1 * qJ(4);
t131 = -Icges(6,1) / 0.2e1;
t129 = cos(qJ(2));
t125 = sin(qJ(2));
t117 = t142 / 0.2e1;
t115 = qJ(2) + qJ(5);
t112 = qJ(4) - qJ(5);
t111 = qJ(4) + qJ(5);
t109 = 0.2e1 * t132;
t95 = qJ(4) + t114;
t93 = qJ(4) + t196;
t89 = cos(t115);
t87 = sin(t115);
t75 = t137 + t162;
t74 = t137 + t163;
t73 = qJ(4) + t162;
t72 = qJ(4) + t163;
t71 = t139 + t92;
t70 = t139 + t91;
t64 = 0.2e1 * t116;
t51 = 0.2e1 * t97;
t33 = 0.2e1 * t77;
t32 = 0.2e1 * t76;
t1 = [-((rSges(4,1) ^ 2 - rSges(4,2) ^ 2) * m(4) - Icges(4,1) + Icges(4,2) + t234 * pkin(3) ^ 2) * t106 * sin(t64) - ((rSges(3,1) ^ 2 - rSges(3,2) ^ 2) * m(3) - Icges(3,1) + Icges(3,2) + t208 * pkin(2) ^ 2) * qJD(2) * sin(t139) - t86 * ((t117 + t141 / 0.2e1 - (pkin(4) + t130) * (-pkin(4) + t130)) * m(6) + (rSges(5,1) ^ 2 - rSges(5,2) ^ 2) * m(5) - Icges(5,1) + t131 + Icges(5,2) - Icges(6,2) / 0.2e1 + Icges(6,3)) * sin(t51) - pkin(2) * cos(t113) * t201 + (rSges(3,1) * t229 - Icges(3,4)) * cos(t139) * t240 + (rSges(4,2) * t232 - Icges(4,4)) * cos(t64) * t241 + (rSges(5,1) * t228 - pkin(4) * t219 - Icges(5,4)) * cos(t51) * t242 + t252 / 0.2e1 - t161 + t165 + (sin(t92) * t160 + cos(t92) * t159) * t83 + (-t55 * cos(t32) / 0.2e1 + cos(t33) * t235) * t65 + (t21 * cos(t75) - t24 * sin(t75)) * (t119 + t86) + (t22 * cos(t74) - t23 * sin(t74)) * (t118 + t86) + (t55 * sin(t32) + t56 * sin(t33)) * t236 + t254 + (t31 * sin(t114) + cos(t114) * t233) * t100 * t244 + (sin(t112) * t157 + cos(t112) * t156) * t102 + (cos(t93) * t185 + cos(t95) * t186 + t62 * t187 - t168) * t39 + (-pkin(2) * t103 * sin(t113) + sin(t93) * t185 + sin(t95) * t186 - pkin(3) * t195 + t60 * t187) * t63 + (-cos(t70) * t182 + cos(t71) * t181 - cos(t72) * t180 + cos(t73) * t179 + cos(t91) * t158 + cos(t111) * t152 + t47 * t183 - t46 * t184) * rSges(6,2) + (-sin(t70) * t182 - sin(t71) * t181 - sin(t72) * t180 - sin(t73) * t179 + sin(t91) * t158 + sin(t111) * t152 - t43 * t183 - t42 * t184) * rSges(6,1) + ((t90 * t233 + t31 * t88) * t241 + ((rSges(3,1) * m(3) + t208 * pkin(2)) * t125 + t129 * t229) * t240) * pkin(1); 0.2e1 * (t213 / 0.4e1 + t216 / 0.4e1) * t65 + (-(t140 + 0.2e1 * t222) * t43 / 0.4e1 + ((t109 + 0.2e1 * t98) * m(6) - 0.2e1 * Icges(6,5)) * t47 / 0.4e1) * t56 + (-(t140 - 0.2e1 * t223) * t42 / 0.4e1 + ((t109 - 0.2e1 * t98) * m(6) + 0.2e1 * Icges(6,5)) * t46 / 0.4e1) * t55 + (rSges(6,1) * t87 / 0.2e1 + rSges(6,2) * t89 / 0.2e1) * m(6) * t221 + 0.4e1 * (-t212 / 0.4e1 + t211 / 0.4e1) * t86 + 0.4e1 * (t218 / 0.4e1 - t210 / 0.4e1) * t106 + t154 + (-Icges(3,5) * t129 + t125 * Icges(3,6) + (rSges(3,1) * t129 - rSges(3,2) * t125) * rSges(3,3) * m(3) + (t225 + t226) * pkin(2) * t129) * qJD(2) + t250 + t255; t256 + 0.2e1 * t220 + 0.2e1 * t251 + 0.2e1 * t254; -t253 - t106 * (t210 - t218) + t145 + t249 + t250; t209 * t243 + (t2 - (t38 + 0.2e1 * t203) * t10) * t123 + (t148 * t127 + t146) * pkin(2) + t256; 0.2e1 * t147 + t256; -t257 / 0.2e1 + t145; t220 + t251 + t256; t147 + t256; -0.4e1 * qJD(5) * (Icges(6,4) / 0.2e1 + ((t131 + Icges(6,2) / 0.2e1) * t122 + t65 * t126) * t126 + (((t117 - t141 / 0.2e1) * t122 + t132 / 0.2e1) * t126 + (t122 * pkin(4) / 0.2e1 - rSges(6,2) / 0.2e1) * rSges(6,1)) * m(6)); -t173 / 0.2e1 + t253 - ((t141 + t142) * m(6) + Icges(6,3)) * t86 * t60 + m(6) * ((-pkin(1) * t189 + t89 * t167 + t61 * t169) * rSges(6,2) + (-pkin(1) * t192 + t87 * t167 + t59 * t169) * rSges(6,1)) + t154 + t245; t246 + (((-rSges(6,1) * t191 - rSges(6,2) * t194) * t126 - (rSges(6,1) * t194 - rSges(6,2) * t191) * t122) * pkin(3) + ((-(rSges(6,1) * t190 + rSges(6,2) * t199) * t127 + (-rSges(6,1) * t188 + rSges(6,2) * t200) * t123) * t126 - ((rSges(6,1) * t199 - rSges(6,2) * t190) * t127 - (rSges(6,1) * t200 + rSges(6,2) * t188) * t123) * t122) * pkin(2)) * m(6); (-rSges(6,2) * t153 - qJD(5) * (rSges(6,1) * t176 + t41)) * t126 - t122 * (rSges(6,1) * t153 - qJD(5) * (rSges(6,2) * t176 + t40)); t246; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
