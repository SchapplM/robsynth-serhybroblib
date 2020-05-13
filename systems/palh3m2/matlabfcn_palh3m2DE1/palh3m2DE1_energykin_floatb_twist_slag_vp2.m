% Calculate kinetic energy for
% palh3m2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [9x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 02:05
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m2DE1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE1_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE1_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh3m2DE1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE1_energykin_floatb_twist_slag_vp2: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE1_energykin_floatb_twist_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2DE1_energykin_floatb_twist_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2DE1_energykin_floatb_twist_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:57:59
% EndTime: 2020-05-07 01:58:02
% DurationCPUTime: 2.50s
% Computational Cost: add. (2589->275), mult. (4105->399), div. (0->0), fcn. (4412->34), ass. (0->136)
t135 = cos(qJ(1));
t130 = sin(qJ(1));
t146 = t130 * V_base(5);
t72 = t135 * V_base(4) + t146;
t139 = t130 * V_base(2) + t135 * V_base(1);
t118 = pkin(17) + pkin(18);
t105 = sin(t118);
t106 = cos(t118);
t129 = sin(qJ(2));
t134 = cos(qJ(2));
t128 = sin(qJ(3));
t133 = cos(qJ(3));
t122 = cos(pkin(16));
t136 = cos(pkin(15));
t152 = sin(pkin(16));
t157 = sin(pkin(15));
t65 = t122 * t157 + t136 * t152;
t66 = t122 * t136 - t152 * t157;
t46 = t128 * t66 + t133 * t65;
t47 = -t128 * t65 + t133 * t66;
t30 = t47 * t129 + t134 * t46;
t31 = -t129 * t46 + t134 * t47;
t166 = t105 * t31 + t106 * t30;
t107 = V_base(6) + qJD(1);
t164 = t107 * t66 + t65 * t72;
t119 = qJ(3) + qJ(2);
t142 = pkin(15) + t118 + t119;
t90 = -qJ(1) + t142;
t163 = -sin(t90) / 0.2e1;
t89 = qJ(1) + t142;
t162 = -cos(t89) / 0.2e1;
t18 = -t105 * t30 + t106 * t31;
t150 = t133 * qJD(2);
t69 = t128 * t129 - t133 * t134;
t70 = t128 * t134 + t129 * t133;
t71 = t130 * V_base(4) - t135 * V_base(5);
t81 = pkin(11) * t129 + pkin(12) * t134 + pkin(1);
t88 = pkin(11) * t134 - pkin(12) * t129;
t27 = -pkin(1) * t150 - (t128 * t88 + t133 * t81) * t71 + t69 * V_base(3) + t139 * t70;
t64 = qJD(2) + t71;
t62 = qJD(3) + t64;
t22 = pkin(4) * t62 + t27;
t151 = t128 * qJD(2);
t26 = -t70 * V_base(3) - pkin(1) * t151 + (-t128 * t81 + t133 * t88) * t71 + t139 * t69;
t6 = -t166 * t22 + t18 * t26;
t114 = qJ(1) + t119;
t161 = sin(t114) / 0.2e1;
t115 = qJ(1) - t119;
t160 = sin(t115) / 0.2e1;
t121 = -qJ(3) + qJ(1);
t159 = -sin(t121) / 0.2e1;
t120 = qJ(3) + qJ(1);
t158 = cos(t120) / 0.2e1;
t156 = sin(pkin(18));
t96 = t129 * pkin(1) + pkin(11);
t147 = V_base(4) * pkin(12);
t97 = pkin(1) * t134 + pkin(12);
t104 = V_base(1) * t130;
t144 = -t135 * V_base(2) + t104;
t5 = t166 * t26 + t18 * t22;
t32 = t107 * t65 - t66 * t72;
t21 = t105 * t32 - t106 * t164;
t143 = V_base(4) * pkin(11) - V_base(2);
t92 = V_base(5) * pkin(11) + V_base(1);
t56 = -t130 * t92 - t135 * t143;
t125 = cos(pkin(18));
t67 = t125 * t157 + t136 * t156;
t68 = t125 * t136 - t156 * t157;
t141 = (t129 * t68 + t134 * t67) * pkin(1);
t140 = (t129 * t67 - t134 * t68) * pkin(1);
t39 = -t107 * t97 + t72 * t96 + t144;
t54 = t107 * t134 - t129 * t72;
t55 = t107 * t129 + t134 * t72;
t36 = t128 * t55 - t133 * t54;
t25 = -pkin(4) * t36 + t39;
t138 = V_base(3) ^ 2;
t137 = cos(pkin(14));
t132 = cos(qJ(4));
t131 = sin(pkin(14));
t127 = sin(qJ(4));
t126 = cos(pkin(17));
t124 = sin(pkin(17));
t123 = pkin(11) + pkin(13);
t116 = V_base(5) * pkin(12);
t113 = cos(t121);
t111 = cos(t119);
t109 = sin(t120);
t108 = sin(t119);
t103 = cos(t115);
t102 = cos(t114);
t87 = pkin(11) * t157 - pkin(12) * t136;
t86 = pkin(11) * t136 + pkin(12) * t157;
t85 = cos(t90);
t82 = sin(t89);
t80 = -t92 - t147;
t79 = t92 - t147;
t78 = -t116 + t143;
t77 = t116 + t143;
t74 = t131 * t157 + t137 * t136;
t73 = -t131 * t136 + t137 * t157;
t61 = -pkin(6) * t71 + V_base(3);
t60 = pkin(12) * t71 + V_base(3);
t57 = -t130 * t143 + t135 * t92;
t53 = -t123 * t71 + t139;
t52 = -pkin(12) * t107 - t56;
t51 = -t107 * t111 + t108 * t72;
t50 = -t107 * t108 - t111 * t72;
t49 = -t129 * t73 + t134 * t74;
t48 = t129 * t74 + t134 * t73;
t45 = pkin(6) * t107 + t123 * t72 + t144;
t42 = (t124 * t68 + t126 * t67) * pkin(3) + t96;
t41 = t129 * t60 + t134 * t57;
t40 = -t129 * t57 + t134 * t60;
t38 = t39 ^ 2;
t37 = -t128 * t54 - t133 * t55;
t35 = -t107 * t68 - t67 * t72;
t34 = t107 * t67 - t68 * t72;
t29 = t107 * t48 + t49 * t72;
t28 = t107 * t49 - t48 * t72;
t24 = -t48 * t53 + t49 * t61;
t23 = t48 * t61 + t49 * t53;
t20 = t105 * t164 + t32 * t106;
t19 = qJD(4) - t21;
t16 = qJD(2) * t140 - t68 * V_base(3) + t71 * (t87 * t125 + t156 * t86 + t140) - t139 * t67;
t15 = qJD(2) * t141 + t67 * V_base(3) + t71 * (t86 * t125 - t156 * t87 + t141) - t139 * t68;
t14 = t127 * t71 + t132 * t20;
t13 = -t127 * t20 + t132 * t71;
t12 = t80 * t103 / 0.2e1 + t78 * t160 - t79 * t102 / 0.2e1 + t77 * t161 - V_base(3) * t108 + ((t163 - t82 / 0.2e1) * V_base(5) + (t85 / 0.2e1 + t162) * V_base(4)) * pkin(3) + (-t151 + (t159 + t109 / 0.2e1) * V_base(5) + (-t113 / 0.2e1 + t158) * V_base(4)) * pkin(1);
t11 = -t78 * t103 / 0.2e1 + t80 * t160 + t77 * t102 / 0.2e1 + t79 * t161 - V_base(3) * t111 + ((-t85 / 0.2e1 + t162) * V_base(5) + (t163 + t82 / 0.2e1) * V_base(4)) * pkin(3) + (-t150 + (t113 / 0.2e1 + t158) * V_base(5) + (t159 - t109 / 0.2e1) * V_base(4)) * pkin(1);
t10 = t42 * t146 + t104 + t107 * ((-t124 * t67 + t126 * t68) * pkin(3) - t97) + (t42 * V_base(4) - V_base(2)) * t135;
t7 = -pkin(8) * t21 - pkin(10) * t20 + t25;
t4 = pkin(10) * t71 + t6;
t3 = -pkin(8) * t71 - t5;
t2 = t127 * t7 + t132 * t4;
t1 = -t127 * t4 + t132 * t7;
t8 = (V_base(3) * mrSges(2,2) - t56 * mrSges(2,3) + Ifges(2,1) * t72 / 0.2e1) * t72 + (-t45 * mrSges(7,1) + t23 * mrSges(7,3) + Ifges(7,4) * t29 + Ifges(7,2) * t28 / 0.2e1) * t28 + (t45 * mrSges(7,2) - t24 * mrSges(7,3) + Ifges(7,1) * t29 / 0.2e1) * t29 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t19 + Ifges(6,1) * t14 / 0.2e1) * t14 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + m(2) * (t56 ^ 2 + t57 ^ 2 + t138) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t138) / 0.2e1 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t39 * mrSges(4,1) + t26 * mrSges(4,3) + Ifges(4,4) * t37 + Ifges(4,2) * t36 / 0.2e1) * t36 + (t39 * mrSges(8,2) - t16 * mrSges(8,3) + Ifges(8,4) * t35 + Ifges(8,1) * t34 / 0.2e1) * t34 + (V_base(3) * mrSges(2,1) + t5 * mrSges(5,1) + t16 * mrSges(8,1) - t6 * mrSges(5,2) - t15 * mrSges(8,2) - t57 * mrSges(2,3) - Ifges(2,4) * t72 + Ifges(5,5) * t20 + Ifges(8,5) * t34 - Ifges(2,6) * t107 + Ifges(5,6) * t21 + Ifges(8,6) * t35 + (Ifges(2,2) / 0.2e1 + Ifges(8,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t71) * t71 + (t40 * mrSges(3,1) + t24 * mrSges(7,1) - t41 * mrSges(3,2) - t23 * mrSges(7,2) + Ifges(3,5) * t55 + Ifges(7,5) * t29 + Ifges(3,6) * t54 + Ifges(7,6) * t28 + (Ifges(7,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t64) * t64 + (t27 * mrSges(4,1) + t11 * mrSges(9,1) - t26 * mrSges(4,2) - t12 * mrSges(9,2) + Ifges(4,5) * t37 + Ifges(9,5) * t50 + Ifges(4,6) * t36 + Ifges(9,6) * t51 + (Ifges(9,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t62) * t62 + (-t39 * mrSges(8,1) + t15 * mrSges(8,3) + Ifges(8,2) * t35 / 0.2e1) * t35 + (t52 * mrSges(3,2) - t40 * mrSges(3,3) + Ifges(3,1) * t55 / 0.2e1) * t55 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t14 + Ifges(6,6) * t19 + Ifges(6,2) * t13 / 0.2e1) * t13 + (t39 * mrSges(4,2) - t27 * mrSges(4,3) + Ifges(4,1) * t37 / 0.2e1) * t37 + m(3) * (t40 ^ 2 + t41 ^ 2 + t52 ^ 2) / 0.2e1 + m(4) * (t26 ^ 2 + t27 ^ 2 + t38) / 0.2e1 + m(8) * (t15 ^ 2 + t16 ^ 2 + t38) / 0.2e1 + m(7) * (t23 ^ 2 + t24 ^ 2 + t45 ^ 2) / 0.2e1 + (-t52 * mrSges(3,1) + t41 * mrSges(3,3) + Ifges(3,4) * t55 + Ifges(3,2) * t54 / 0.2e1) * t54 + m(5) * (t25 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + (-t10 * mrSges(9,1) + t12 * mrSges(9,3) + Ifges(9,2) * t51 / 0.2e1) * t51 + m(9) * (t10 ^ 2 + t11 ^ 2 + t12 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t19 / 0.2e1) * t19 + (t25 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,4) * t21 + Ifges(5,1) * t20 / 0.2e1) * t20 + (t10 * mrSges(9,2) - t11 * mrSges(9,3) + Ifges(9,4) * t51 + Ifges(9,1) * t50 / 0.2e1) * t50 + (t56 * mrSges(2,1) - t57 * mrSges(2,2) + Ifges(2,5) * t72 + Ifges(2,3) * t107 / 0.2e1) * t107 + (-t25 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,2) * t21 / 0.2e1) * t21 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t8;
