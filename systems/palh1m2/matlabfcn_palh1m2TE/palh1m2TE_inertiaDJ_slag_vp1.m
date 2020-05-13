% Calculate time derivative of joint inertia matrix for
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
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m2TE_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2TE_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_inertiaDJ_slag_vp1: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2TE_inertiaDJ_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2TE_inertiaDJ_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m2TE_inertiaDJ_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:18:14
% EndTime: 2020-05-01 20:18:31
% DurationCPUTime: 4.49s
% Computational Cost: add. (741->256), mult. (1040->336), div. (0->0), fcn. (408->88), ass. (0->147)
t100 = cos(qJ(3));
t200 = pkin(2) * m(10);
t196 = m(5) + m(6);
t28 = pkin(5) * t196 + m(11) * rSges(11,1) + m(4) * rSges(4,1);
t168 = m(11) * rSges(11,2);
t188 = rSges(4,2) * m(4);
t56 = t168 + t188;
t88 = qJ(3) + pkin(19);
t96 = sin(qJ(3));
t163 = ((-t28 * t100 + t96 * t56) * pkin(1) - (rSges(10,2) * sin(t88) + rSges(10,1) * cos(t88)) * t200) * qJD(3);
t202 = 2 * qJD(2);
t84 = pkin(22) + pkin(21);
t127 = pkin(18) - t84;
t114 = (-pkin(20) + t127);
t87 = qJD(2) + qJD(3);
t201 = pkin(1) * m(6);
t199 = pkin(5) * m(6);
t86 = qJD(2) - qJD(4);
t198 = t86 / 0.2e1;
t97 = sin(qJ(2));
t197 = -t97 / 0.2e1;
t194 = pkin(4) * m(11);
t192 = m(5) * rSges(5,3);
t191 = m(6) * rSges(6,2);
t190 = m(9) * rSges(9,2);
t189 = rSges(3,2) * m(3);
t187 = rSges(6,2) * pkin(9);
t186 = rSges(10,2) * m(10);
t185 = rSges(3,3) * m(3);
t184 = rSges(4,3) * m(4);
t183 = rSges(7,3) * m(7);
t182 = rSges(9,3) * m(9);
t181 = rSges(10,3) * m(10);
t180 = m(8) + m(11) + m(4) + t196;
t165 = rSges(11,3) * m(11);
t89 = qJ(3) + qJ(2);
t77 = cos(t89);
t177 = (-pkin(5) * t192 - rSges(11,1) * t165 - rSges(4,1) * t184 + Icges(11,5) + Icges(4,5)) * t77;
t72 = qJ(2) + t88;
t58 = cos(t72);
t176 = (pkin(2) * t181 + rSges(9,1) * t182 - Icges(9,5)) * t58;
t57 = sin(t72);
t175 = (rSges(9,2) * t182 - Icges(9,6)) * t57;
t76 = sin(t89);
t174 = t76 * (rSges(11,2) * t165 + rSges(4,2) * t184 - Icges(11,6) - Icges(4,6));
t167 = m(6) * qJD(4);
t103 = pkin(11) + rSges(6,3);
t166 = rSges(6,2) * t103;
t164 = t100 * t97;
t82 = t103 * rSges(6,1);
t162 = qJD(2) * t97;
t160 = qJD(4) * (rSges(6,1) * t191 - Icges(6,4));
t101 = cos(qJ(2));
t159 = t100 * t101;
t158 = qJD(2) * t101;
t155 = pkin(18) - pkin(22);
t85 = qJD(2) + qJD(4);
t149 = t201 / 0.2e1;
t148 = -t199 / 0.2e1;
t147 = t199 / 0.2e1;
t95 = sin(qJ(4));
t99 = cos(qJ(4));
t35 = t95 * rSges(6,1) + t99 * rSges(6,2);
t18 = pkin(9) * t35 * t167;
t135 = t85 * t149;
t134 = t86 * t149;
t70 = qJD(3) + t85;
t133 = t70 * t148;
t132 = t70 * t147;
t71 = qJD(3) + t86;
t131 = t71 * t148;
t120 = rSges(6,2) * t131;
t122 = rSges(6,1) * t71 * t147;
t123 = rSges(6,1) * t133;
t79 = qJ(4) + t89;
t80 = -qJ(4) + t89;
t130 = cos(t80) * t120 + sin(t80) * t122 + rSges(6,2) * cos(t79) * t133 + sin(t79) * t123;
t126 = 2 * t114;
t62 = t96 * pkin(5) + pkin(1);
t16 = -pkin(5) * t159 + t62 * t97;
t119 = -qJ(2) + t127;
t118 = qJ(2) + t127;
t116 = t99 * rSges(6,1) - t95 * rSges(6,2);
t102 = cos(pkin(18));
t15 = (-(-t82 * m(6) + Icges(6,5)) * t95 - (-t166 * m(6) + Icges(6,6)) * t99) * qJD(4);
t98 = sin(pkin(18));
t115 = t18 * t102 + t98 * t15;
t25 = t101 * t96 + t164;
t24 = -t97 * t96 + t159;
t48 = -qJ(2) + t114;
t47 = qJ(2) + t114;
t46 = -qJ(4) + t114;
t45 = qJ(4) + t114;
t113 = t24 * qJD(3);
t38 = -qJ(3) + t48;
t37 = qJ(3) + t47;
t109 = 0.2e1 * qJ(2);
t108 = 0.2e1 * qJ(4);
t104 = rSges(6,1) * pkin(9);
t94 = cos(pkin(20));
t93 = sin(pkin(20));
t92 = qJ(2) - qJ(4);
t91 = qJ(2) + qJ(4);
t90 = t109 + qJ(3);
t81 = pkin(17) + qJ(2) - pkin(18);
t75 = t109 + t88;
t74 = -qJ(2) + t155;
t73 = qJ(2) + t155;
t67 = cos(t84);
t66 = sin(t84);
t63 = 0.2e1 * t89;
t61 = pkin(9) * m(6) + m(5) * rSges(5,1);
t60 = cos(t81);
t59 = sin(t81);
t55 = 0.2e1 * t81;
t54 = 0.2e1 * t72;
t53 = -m(5) * rSges(5,2) + t103 * m(6);
t52 = -qJ(3) + t119;
t51 = qJ(3) + t118;
t50 = -qJ(4) + t126;
t49 = qJ(4) + t126;
t42 = -qJ(2) + t46;
t41 = -qJ(2) + t45;
t40 = qJ(2) + t46;
t39 = qJ(2) + t45;
t34 = -qJ(4) + t38;
t33 = qJ(4) + t38;
t32 = -qJ(4) + t37;
t31 = qJ(4) + t37;
t30 = 0.2e1 * t46;
t29 = 0.2e1 * t45;
t26 = t116 * qJD(4);
t17 = pkin(5) * t164 + t62 * t101;
t14 = t87 * t25;
t13 = qJD(2) * t24 + t113;
t12 = t25 * t102 - t98 * t24;
t11 = t24 * t102 + t98 * t25;
t10 = t62 * t158 + (qJD(3) * t25 + t100 * t162) * pkin(5);
t9 = -t62 * t162 + (t100 * t158 + t113) * pkin(5);
t8 = t17 * t102 + t16 * t98;
t7 = t16 * t102 - t98 * t17;
t5 = t15 * t102 - t18 * t98;
t4 = -t14 * t102 + t98 * t13;
t3 = t13 * t102 + t98 * t14;
t2 = t10 * t102 - t98 * t9;
t1 = t98 * t10 + t9 * t102;
t6 = [(((rSges(7,1) ^ 2 - rSges(7,2) ^ 2) * m(7) - Icges(7,1) + Icges(7,2)) * sin(t55) + ((rSges(10,1) ^ 2 - rSges(10,2) ^ 2) * m(10) + (rSges(3,1) ^ 2 - rSges(3,2) ^ 2) * m(3) - Icges(3,1) - Icges(10,1) + Icges(3,2) + Icges(10,2) + t180 * pkin(1) ^ 2) * sin(t109)) * qJD(2) + 0.2e1 * (-(pkin(1) * t180 + rSges(3,1) * m(3) + rSges(10,1) * m(10)) * t158 + (t186 + t189) * t162) * pkin(15) + 0.2e1 * ((-rSges(10,1) * cos(t75) + rSges(10,2) * sin(t75)) * t200 + (-t28 * cos(t90) + t56 * sin(t90)) * pkin(1)) * (qJD(2) + qJD(3) / 0.2e1) - t18 - (cos(t30) + cos(t29)) * t160 / 0.2e1 + (((t104 - t166) * m(6) + Icges(6,6)) * sin(t50) + ((t82 - t187) * m(6) - Icges(6,5)) * cos(t49)) * qJD(4) / 0.2e1 - (((t104 + t166) * m(6) - Icges(6,6)) * sin(t49) + ((t82 + t187) * m(6) - Icges(6,5)) * cos(t50)) * qJD(4) / 0.2e1 - (t86 * sin(t41) + t85 * sin(t39)) * pkin(1) * t191 / 0.2e1 + (cos(t41) + cos(t40)) * rSges(6,1) * t134 + (cos(t42) + cos(t39)) * rSges(6,1) * t135 + (cos(t34) + cos(t31)) * rSges(6,2) * t132 + (-sin(t29) / 0.4e1 + sin(t30) / 0.4e1 + sin(t108) / 0.2e1) * qJD(4) * ((rSges(6,1) ^ 2 - rSges(6,2) ^ 2) * m(6) - Icges(6,1) + Icges(6,2)) + (cos(t33) + cos(t32)) * t120 + ((cos(t45) + cos(t46)) * rSges(6,2) + (sin(t45) - sin(t46)) * rSges(6,1)) * pkin(15) * t167 + (sin(t40) * t134 + sin(t42) * t135) * rSges(6,2) - t163 + (sin(t31) * t132 + sin(t33) * t131) * rSges(6,1) + sin(t32) * t122 + sin(t34) * t123 + cos(t108) * t160 + (-0.2e1 * (rSges(11,1) * t168 + rSges(4,1) * t188 - Icges(11,4) - Icges(4,4)) * cos(t63) - 0.2e1 * (rSges(9,1) * t190 - Icges(9,4)) * cos(t54) - 0.2e1 * (t58 * t190 + t56 * t77 + t28 * t76 + (m(9) * rSges(9,1) + t200) * t57) * pkin(15) - (pkin(2) ^ 2 * m(10) + (rSges(9,1) ^ 2 - rSges(9,2) ^ 2) * m(9) - Icges(9,1) + Icges(9,2)) * sin(t54) - ((rSges(11,1) ^ 2 - rSges(11,2) ^ 2) * m(11) + (rSges(4,1) ^ 2 - rSges(4,2) ^ 2) * m(4) - Icges(11,1) + Icges(11,2) - Icges(4,1) + Icges(4,2) + t196 * pkin(5) ^ 2) * sin(t63) + ((sin(t37) - sin(t38)) * t61 + (-cos(t37) + cos(t38)) * t53) * pkin(5) + ((cos(t51) + cos(t52)) * rSges(11,2) + (sin(t51) - sin(t52)) * rSges(11,1)) * t194) * t87 + ((rSges(3,1) * t189 + rSges(10,1) * t186 - Icges(3,4) - Icges(10,4)) * cos(t109) + (m(7) * rSges(7,1) * rSges(7,2) - Icges(7,4)) * cos(t55) + (rSges(7,1) * t60 - rSges(7,2) * t59) * m(7) * pkin(14)) * t202 + ((cos(t119) + cos(t118)) * t194 + (cos(t47) + cos(t48)) * t61 + (sin(t47) + sin(t48)) * t53 + ((-sin(t73) - sin(t74)) * rSges(8,2) + (cos(t73) + cos(t74)) * rSges(8,1)) * m(8)) * pkin(1) * qJD(2); (0.2e1 * (-rSges(3,1) * t185 - rSges(10,1) * t181 + Icges(3,5) + Icges(10,5)) * t197 + t101 * (rSges(3,2) * t185 + rSges(10,2) * t181 - Icges(3,6) - Icges(10,6))) * qJD(2) + ((sin(t92) * t198 + t85 * sin(t91) / 0.2e1) * rSges(6,2) + (cos(t92) * t198 - t85 * cos(t91) / 0.2e1) * rSges(6,1)) * t201 + (-(-rSges(7,1) * t183 + Icges(7,5)) * t59 / 0.2e1 + (rSges(7,2) * t183 - Icges(7,6)) * t60 / 0.2e1 + (-rSges(8,3) * m(8) - t165 - t184 - t192) * pkin(1) * t197) * t202 + 0.2e1 * (-t176 / 0.2e1 + t175 / 0.2e1 + t177 / 0.2e1 + t174 / 0.2e1) * t87 + t130; -0.2e1 * t163; (t174 + t175 - t176 + t177) * t87 + t130; -t163; 0; (-t115 * t94 + t93 * t5) * t67 + (t115 * t93 + t5 * t94) * t66 + (t10 * t116 - (-pkin(15) + t16) * t35 * qJD(4)) * m(6); (-t26 * ((-t93 * t7 + t8 * t94) * t67 - (t7 * t94 + t93 * t8) * t66) - t35 * ((t1 * t94 - t2 * t93) * t67 - (t93 * t1 + t2 * t94) * t66)) * m(6); (-((t3 * t94 + t4 * t93) * t67 + t66 * (-t3 * t93 + t4 * t94)) * t35 - ((t11 * t93 + t12 * t94) * t67 + t66 * (t11 * t94 - t12 * t93)) * t26) * t199; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t6(1), t6(2), t6(4), t6(7); t6(2), t6(3), t6(5), t6(8); t6(4), t6(5), t6(6), t6(9); t6(7), t6(8), t6(9), t6(10);];
Mq = res;
