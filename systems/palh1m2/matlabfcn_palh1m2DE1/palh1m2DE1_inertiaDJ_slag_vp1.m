% Calculate time derivative of joint inertia matrix for
% palh1m2DE1
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
% Datum: 2020-05-01 21:04
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m2DE1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE1_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE1_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE1_inertiaDJ_slag_vp1: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE1_inertiaDJ_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2DE1_inertiaDJ_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m2DE1_inertiaDJ_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:57:38
% EndTime: 2020-05-01 20:57:56
% DurationCPUTime: 5.55s
% Computational Cost: add. (741->251), mult. (1145->321), div. (0->0), fcn. (408->88), ass. (0->145)
t196 = m(5) + m(6);
t175 = m(6) * qJD(4);
t101 = cos(qJ(4));
t97 = sin(qJ(4));
t39 = rSges(6,1) * t97 + rSges(6,2) * t101;
t21 = t39 * pkin(9) * t175;
t84 = pkin(22) + pkin(21);
t132 = pkin(18) - t84;
t120 = (-pkin(20) + t132);
t102 = cos(qJ(3));
t199 = pkin(2) * m(10);
t30 = pkin(5) * t196 + m(11) * rSges(11,1) + m(4) * rSges(4,1);
t176 = m(11) * rSges(11,2);
t186 = rSges(4,2) * m(4);
t59 = t176 + t186;
t89 = qJ(3) + pkin(19);
t98 = sin(qJ(3));
t5 = ((t30 * t102 - t59 * t98) * pkin(1) + (sin(t89) * rSges(10,2) + cos(t89) * rSges(10,1)) * t199) * qJD(3);
t205 = m(8) + m(11) + m(4) + t196;
t88 = qJD(2) + qJD(3);
t201 = 0.2e1 * qJD(2);
t200 = pkin(1) * m(6);
t198 = pkin(5) * m(6);
t86 = qJD(4) - qJD(2);
t197 = -t86 / 0.2e1;
t194 = pkin(4) * m(11);
t192 = m(7) * rSges(7,3);
t191 = m(9) * rSges(9,2);
t190 = m(9) * rSges(9,3);
t189 = rSges(3,1) * m(3);
t188 = rSges(10,1) * m(10);
t187 = rSges(3,2) * m(3);
t185 = rSges(6,2) * pkin(9);
t184 = rSges(10,2) * m(10);
t183 = rSges(4,3) * m(4);
t182 = rSges(5,3) * m(5);
t105 = -rSges(6,3) - pkin(11);
t173 = rSges(6,2) * t105;
t172 = rSges(11,3) * m(11);
t99 = sin(qJ(2));
t171 = t102 * t99;
t170 = t105 * rSges(6,1);
t169 = qJD(2) * t99;
t167 = qJD(4) * (m(6) * rSges(6,1) * rSges(6,2) - Icges(6,4));
t103 = cos(qJ(2));
t166 = t102 * t103;
t91 = qJ(3) + qJ(2);
t90 = qJ(4) + qJ(2);
t165 = qJD(2) * t103;
t161 = pkin(18) - pkin(22);
t85 = qJD(4) + qJD(2);
t154 = -t200 / 0.2e1;
t153 = t200 / 0.2e1;
t152 = -t198 / 0.2e1;
t151 = t198 / 0.2e1;
t147 = pkin(15) * t175;
t139 = t85 * t153;
t138 = t86 * t154;
t71 = qJD(3) + t85;
t137 = t71 * t152;
t136 = t71 * t151;
t72 = -qJD(3) + t86;
t135 = t72 * t151;
t130 = 2 * t120;
t66 = pkin(5) * t98 + pkin(1);
t19 = -pkin(5) * t166 + t66 * t99;
t127 = rSges(6,1) * t137;
t126 = rSges(6,1) * t72 * t152;
t124 = rSges(6,2) * t135;
t123 = -qJ(2) + t132;
t122 = qJ(2) + t132;
t121 = rSges(6,1) * t101 - rSges(6,2) * t97;
t100 = sin(pkin(18));
t104 = cos(pkin(18));
t16 = (-(m(6) * t170 + Icges(6,5)) * t97 - t101 * (t173 * m(6) + Icges(6,6))) * qJD(4);
t119 = t100 * t16 + t104 * t21;
t27 = t103 * t98 + t171;
t26 = -t98 * t99 + t166;
t51 = -qJ(2) + t120;
t50 = qJ(2) + t120;
t49 = -qJ(4) + t120;
t48 = qJ(4) + t120;
t118 = t26 * qJD(3);
t43 = -qJ(2) + t49;
t42 = qJ(2) + t49;
t41 = -qJ(2) + t48;
t40 = qJ(2) + t48;
t73 = qJ(2) + t89;
t61 = sin(t73);
t62 = cos(t73);
t77 = sin(t91);
t78 = cos(t91);
t80 = qJ(3) + t90;
t81 = -qJ(4) + t91;
t117 = rSges(6,2) * cos(t80) * t137 + cos(t81) * t124 + sin(t80) * t127 + sin(t81) * t126 + (t77 * (rSges(11,2) * t172 + rSges(4,2) * t183 - Icges(11,6) - Icges(4,6)) + (-pkin(5) * t182 - rSges(11,1) * t172 - rSges(4,1) * t183 + Icges(11,5) + Icges(4,5)) * t78 + (-rSges(9,1) * t190 - rSges(10,3) * t199 + Icges(9,5)) * t62 - (-rSges(9,2) * t190 + Icges(9,6)) * t61) * t88;
t114 = -0.4e1 * Icges(6,5);
t112 = 0.2e1 * qJ(2);
t111 = 0.2e1 * qJ(4);
t106 = rSges(6,1) * pkin(9);
t95 = cos(pkin(20));
t94 = sin(pkin(20));
t93 = qJ(2) - qJ(4);
t92 = t112 + qJ(3);
t82 = pkin(17) - pkin(18) + qJ(2);
t76 = t112 + t89;
t75 = -qJ(2) + t161;
t74 = qJ(2) + t161;
t70 = cos(t84);
t69 = sin(t84);
t67 = 0.2e1 * t91;
t65 = pkin(9) * m(6) + m(5) * rSges(5,1);
t64 = cos(t82);
t63 = sin(t82);
t60 = 0.4e1 * t170;
t58 = 0.2e1 * t82;
t57 = 0.2e1 * t73;
t56 = m(5) * rSges(5,2) + m(6) * t105;
t55 = -qJ(3) + t123;
t54 = qJ(3) + t122;
t53 = -qJ(4) + t130;
t52 = qJ(4) + t130;
t45 = -qJ(3) + t51;
t44 = qJ(3) + t50;
t36 = -qJ(3) + t43;
t35 = qJ(3) + t42;
t34 = -qJ(3) + t41;
t33 = qJ(3) + t40;
t32 = 0.2e1 * t49;
t31 = 0.2e1 * t48;
t28 = t121 * qJD(4);
t20 = pkin(5) * t171 + t103 * t66;
t14 = t88 * t27;
t13 = qJD(2) * t26 + t118;
t12 = -t100 * t26 + t104 * t27;
t11 = t100 * t27 + t104 * t26;
t10 = t66 * t165 + (qJD(3) * t27 + t102 * t169) * pkin(5);
t9 = -t66 * t169 + (t102 * t165 + t118) * pkin(5);
t8 = t100 * t19 + t104 * t20;
t7 = -t100 * t20 + t104 * t19;
t6 = -t100 * t21 + t104 * t16;
t4 = t100 * t13 - t104 * t14;
t3 = t100 * t14 + t104 * t13;
t2 = t10 * t104 - t100 * t9;
t1 = t10 * t100 + t104 * t9;
t15 = [-t21 + t5 - (cos(t32) + cos(t31)) * t167 / 0.2e1 + ((rSges(3,1) * t187 + rSges(10,1) * t184 - Icges(3,4) - Icges(10,4)) * cos(t112) + (m(7) * rSges(7,1) * rSges(7,2) - Icges(7,4)) * cos(t58)) * t201 + cos(t111) * t167 + sin(t35) * t126 + sin(t36) * t127 + (cos(t35) + cos(t34)) * t124 + ((t59 * sin(t92) - t30 * cos(t92)) * pkin(1) + (-rSges(10,1) * cos(t76) + rSges(10,2) * sin(t76)) * t199) * (t201 + qJD(3)) + (t86 * sin(t41) * t153 + t85 * sin(t40) * t154 + sin(t43) * t139 + sin(t42) * t138 + (cos(t48) + cos(t49)) * t147 + (cos(t36) + cos(t33)) * t136) * rSges(6,2) + (sin(t33) * t136 + sin(t34) * t135 + (sin(t48) - sin(t49)) * t147 + (cos(t42) + cos(t41)) * t138 + (cos(t43) + cos(t40)) * t139) * rSges(6,1) + (-((t106 - t173) * m(6) - Icges(6,6)) * sin(t52) / 0.2e1 + ((t106 + t173) * m(6) + Icges(6,6)) * sin(t53) / 0.2e1 - ((-t60 + 0.4e1 * t185) * m(6) + t114) * cos(t53) / 0.8e1 + ((-t60 - 0.4e1 * t185) * m(6) + t114) * cos(t52) / 0.8e1 + (sin(t32) / 0.4e1 - sin(t31) / 0.4e1 + sin(t111) / 0.2e1) * ((rSges(6,1) ^ 2 - rSges(6,2) ^ 2) * m(6) - Icges(6,1) + Icges(6,2))) * qJD(4) + (0.2e1 * (t184 + t187) * t169 - 0.2e1 * (pkin(1) * t205 + t188 + t189) * t165) * pkin(15) + (0.2e1 * (-rSges(11,1) * t176 - rSges(4,1) * t186 + Icges(11,4) + Icges(4,4)) * cos(t67) + 0.2e1 * (-rSges(9,1) * t191 + Icges(9,4)) * cos(t57) - (pkin(2) ^ 2 * m(10) + (rSges(9,1) ^ 2 - rSges(9,2) ^ 2) * m(9) - Icges(9,1) + Icges(9,2)) * sin(t57) - (-Icges(11,1) - Icges(4,1) + Icges(11,2) + Icges(4,2) + t196 * pkin(5) ^ 2 + (rSges(11,1) ^ 2 - rSges(11,2) ^ 2) * m(11) + (rSges(4,1) ^ 2 - rSges(4,2) ^ 2) * m(4)) * sin(t67) - 0.2e1 * (t62 * t191 + t59 * t78 + t30 * t77 + (m(9) * rSges(9,1) + t199) * t61) * pkin(15) + ((sin(t44) - sin(t45)) * t65 + (cos(t44) - cos(t45)) * t56) * pkin(5) + ((cos(t54) + cos(t55)) * rSges(11,2) + (sin(t54) - sin(t55)) * rSges(11,1)) * t194) * t88 + (0.2e1 * (rSges(7,1) * t64 - rSges(7,2) * t63) * m(7) * pkin(14) - (Icges(3,1) + Icges(10,1) - Icges(3,2) - Icges(10,2) + (-rSges(10,1) ^ 2 + rSges(10,2) ^ 2) * m(10) + (-rSges(3,1) ^ 2 + rSges(3,2) ^ 2) * m(3) - t205 * pkin(1) ^ 2) * sin(t112) - ((-rSges(7,1) ^ 2 + rSges(7,2) ^ 2) * m(7) + Icges(7,1) - Icges(7,2)) * sin(t58) + (((-sin(t74) - sin(t75)) * rSges(8,2) + (cos(t74) + cos(t75)) * rSges(8,1)) * m(8) + (cos(t50) + cos(t51)) * t65 + (-sin(t50) - sin(t51)) * t56 + (cos(t123) + cos(t122)) * t194) * pkin(1)) * qJD(2); (t103 * (rSges(3,3) * t187 + rSges(10,3) * t184 - Icges(3,6) - Icges(10,6)) - (-rSges(7,1) * t192 + Icges(7,5)) * t63 + (rSges(7,2) * t192 - Icges(7,6)) * t64 + (rSges(3,3) * t189 + rSges(10,3) * t188 - Icges(3,5) - Icges(10,5) - (-rSges(8,3) * m(8) - t172 - t182 - t183) * pkin(1)) * t99) * qJD(2) + ((sin(t93) * t197 + t85 * sin(t90) / 0.2e1) * rSges(6,2) + (cos(t93) * t197 - t85 * cos(t90) / 0.2e1) * rSges(6,1)) * t200 + t117; 0.2e1 * t5; t117; t5; 0; (-t119 * t95 + t94 * t6) * t70 + (t119 * t94 + t6 * t95) * t69 + (t121 * t10 - t39 * (-pkin(15) + t19) * qJD(4)) * m(6); (-((t1 * t95 - t2 * t94) * t70 - (t1 * t94 + t2 * t95) * t69) * t39 - ((-t7 * t94 + t8 * t95) * t70 - (t7 * t95 + t8 * t94) * t69) * t28) * m(6); (-t28 * ((t11 * t94 + t12 * t95) * t70 + t69 * (t11 * t95 - t12 * t94)) - t39 * ((t3 * t95 + t4 * t94) * t70 + t69 * (-t3 * t94 + t4 * t95))) * t198; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t15(1), t15(2), t15(4), t15(7); t15(2), t15(3), t15(5), t15(8); t15(4), t15(5), t15(6), t15(9); t15(7), t15(8), t15(9), t15(10);];
Mq = res;
