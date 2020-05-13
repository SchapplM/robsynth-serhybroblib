% Calculate kinetic energy for
% palh3m2DE2
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
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m2DE2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE2_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh3m2DE2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_energykin_floatb_twist_slag_vp2: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE2_energykin_floatb_twist_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2DE2_energykin_floatb_twist_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2DE2_energykin_floatb_twist_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 02:14:28
% EndTime: 2020-05-07 02:14:31
% DurationCPUTime: 2.70s
% Computational Cost: add. (2589->275), mult. (4101->399), div. (0->0), fcn. (4412->34), ass. (0->136)
t134 = cos(qJ(1));
t129 = sin(qJ(1));
t146 = t129 * V_base(5);
t73 = t134 * V_base(4) + t146;
t167 = t129 * V_base(2) + t134 * V_base(1);
t117 = pkin(17) + pkin(18);
t104 = sin(t117);
t105 = cos(t117);
t128 = sin(qJ(2));
t133 = cos(qJ(2));
t127 = sin(qJ(3));
t132 = cos(qJ(3));
t121 = cos(pkin(16));
t135 = cos(pkin(15));
t152 = sin(pkin(16));
t157 = sin(pkin(15));
t65 = t121 * t157 + t135 * t152;
t66 = t121 * t135 - t152 * t157;
t46 = t127 * t66 + t132 * t65;
t47 = -t127 * t65 + t132 * t66;
t30 = t128 * t47 + t133 * t46;
t31 = -t128 * t46 + t133 * t47;
t166 = t31 * t104 + t105 * t30;
t106 = V_base(6) + qJD(1);
t164 = t106 * t66 + t65 * t73;
t120 = qJ(2) + qJ(3);
t142 = pkin(15) + t117 + t120;
t90 = -qJ(1) + t142;
t163 = -sin(t90) / 0.2e1;
t89 = qJ(1) + t142;
t162 = -cos(t89) / 0.2e1;
t118 = qJ(3) + qJ(1);
t113 = qJ(2) + t118;
t161 = sin(t113) / 0.2e1;
t18 = -t104 * t30 + t105 * t31;
t150 = t132 * qJD(2);
t69 = t129 * V_base(4) - t134 * V_base(5);
t71 = t127 * t133 + t128 * t132;
t72 = -t127 * t128 + t132 * t133;
t81 = pkin(11) * t128 + pkin(12) * t133 + pkin(1);
t87 = pkin(11) * t133 - pkin(12) * t128;
t26 = -pkin(1) * t150 - t69 * (t127 * t87 + t132 * t81) - t72 * V_base(3) + t167 * t71;
t64 = qJD(2) + t69;
t62 = qJD(3) + t64;
t22 = pkin(4) * t62 + t26;
t151 = t127 * qJD(2);
t27 = -t71 * V_base(3) - pkin(1) * t151 + (-t127 * t81 + t132 * t87) * t69 - t167 * t72;
t6 = -t166 * t22 + t18 * t27;
t119 = -qJ(3) + qJ(1);
t114 = -qJ(2) + t119;
t160 = sin(t114) / 0.2e1;
t159 = -sin(t119) / 0.2e1;
t158 = cos(t118) / 0.2e1;
t156 = sin(pkin(18));
t96 = t128 * pkin(1) + pkin(11);
t147 = V_base(4) * pkin(12);
t97 = t133 * pkin(1) + pkin(12);
t103 = V_base(1) * t129;
t144 = -t134 * V_base(2) + t103;
t5 = t166 * t27 + t18 * t22;
t32 = t106 * t65 - t66 * t73;
t21 = t104 * t32 - t105 * t164;
t143 = V_base(4) * pkin(11) - V_base(2);
t92 = V_base(5) * pkin(11) + V_base(1);
t56 = -t129 * t92 - t134 * t143;
t124 = cos(pkin(18));
t138 = -t124 * t157 - t135 * t156;
t68 = t135 * t124 - t156 * t157;
t141 = (t128 * t68 - t133 * t138) * pkin(1);
t140 = (-t128 * t138 - t133 * t68) * pkin(1);
t39 = -t106 * t97 + t73 * t96 + t144;
t54 = t106 * t133 - t128 * t73;
t55 = t106 * t128 + t133 * t73;
t36 = t127 * t55 - t132 * t54;
t25 = -pkin(4) * t36 + t39;
t137 = V_base(3) ^ 2;
t136 = cos(pkin(14));
t131 = cos(qJ(4));
t130 = sin(pkin(14));
t126 = sin(qJ(4));
t125 = cos(pkin(17));
t123 = sin(pkin(17));
t122 = pkin(11) + pkin(13);
t115 = V_base(5) * pkin(12);
t112 = cos(t120);
t111 = cos(t119);
t109 = sin(t120);
t107 = sin(t118);
t102 = cos(t114);
t101 = cos(t113);
t88 = t135 * pkin(11) + pkin(12) * t157;
t86 = pkin(11) * t157 - t135 * pkin(12);
t85 = cos(t90);
t82 = sin(t89);
t80 = -t92 - t147;
t79 = t92 - t147;
t78 = -t115 + t143;
t77 = t115 + t143;
t74 = t130 * t157 + t135 * t136;
t70 = -t135 * t130 + t136 * t157;
t61 = -pkin(6) * t69 + V_base(3);
t60 = pkin(12) * t69 + V_base(3);
t57 = -t129 * t143 + t134 * t92;
t53 = -t122 * t69 + t167;
t52 = -pkin(12) * t106 - t56;
t51 = -t106 * t112 + t109 * t73;
t50 = -t106 * t109 - t112 * t73;
t49 = -t128 * t70 + t133 * t74;
t48 = t128 * t74 + t133 * t70;
t45 = pkin(6) * t106 + t122 * t73 + t144;
t42 = (t123 * t68 - t125 * t138) * pkin(3) + t96;
t41 = t128 * t60 + t133 * t57;
t40 = -t128 * t57 + t133 * t60;
t38 = t39 ^ 2;
t37 = -t127 * t54 - t132 * t55;
t35 = -t106 * t68 + t138 * t73;
t34 = -t106 * t138 - t68 * t73;
t29 = t106 * t48 + t49 * t73;
t28 = t106 * t49 - t48 * t73;
t24 = t48 * t61 + t49 * t53;
t23 = -t48 * t53 + t49 * t61;
t20 = t104 * t164 + t32 * t105;
t19 = qJD(4) - t21;
t16 = qJD(2) * t140 - t68 * V_base(3) + t69 * (t86 * t124 + t156 * t88 + t140) + t167 * t138;
t15 = qJD(2) * t141 - t138 * V_base(3) + t69 * (t88 * t124 - t156 * t86 + t141) - t167 * t68;
t14 = t126 * t69 + t131 * t20;
t13 = -t126 * t20 + t131 * t69;
t12 = -t78 * t102 / 0.2e1 + t80 * t160 + t77 * t101 / 0.2e1 + t79 * t161 - V_base(3) * t112 + (-t150 + (t111 / 0.2e1 + t158) * V_base(5) + (t159 - t107 / 0.2e1) * V_base(4)) * pkin(1) + ((-t85 / 0.2e1 + t162) * V_base(5) + (t163 + t82 / 0.2e1) * V_base(4)) * pkin(3);
t11 = t80 * t102 / 0.2e1 + t78 * t160 - t79 * t101 / 0.2e1 + t77 * t161 - V_base(3) * t109 + ((t163 - t82 / 0.2e1) * V_base(5) + (t85 / 0.2e1 + t162) * V_base(4)) * pkin(3) + (-t151 + (t159 + t107 / 0.2e1) * V_base(5) + (-t111 / 0.2e1 + t158) * V_base(4)) * pkin(1);
t10 = t42 * t146 + t103 + t106 * ((t123 * t138 + t125 * t68) * pkin(3) - t97) + (t42 * V_base(4) - V_base(2)) * t134;
t7 = -pkin(8) * t21 - pkin(10) * t20 + t25;
t4 = pkin(10) * t69 + t6;
t3 = -pkin(8) * t69 - t5;
t2 = t126 * t7 + t131 * t4;
t1 = -t126 * t4 + t131 * t7;
t8 = (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t39 * mrSges(8,2) - t16 * mrSges(8,3) + Ifges(8,4) * t35 + Ifges(8,1) * t34 / 0.2e1) * t34 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t19 / 0.2e1) * t19 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t52 * mrSges(3,1) + t41 * mrSges(3,3) + Ifges(3,4) * t55 + Ifges(3,2) * t54 / 0.2e1) * t54 + m(2) * (t56 ^ 2 + t57 ^ 2 + t137) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t137) / 0.2e1 + (-t10 * mrSges(9,1) + t11 * mrSges(9,3) + Ifges(9,2) * t51 / 0.2e1) * t51 + (t39 * mrSges(4,2) - t26 * mrSges(4,3) + Ifges(4,1) * t37 / 0.2e1) * t37 + (V_base(3) * mrSges(2,1) + t5 * mrSges(5,1) + t16 * mrSges(8,1) - t6 * mrSges(5,2) - t15 * mrSges(8,2) - t57 * mrSges(2,3) - Ifges(2,4) * t73 + Ifges(5,5) * t20 + Ifges(8,5) * t34 - Ifges(2,6) * t106 + Ifges(5,6) * t21 + Ifges(8,6) * t35 + (Ifges(2,2) / 0.2e1 + Ifges(8,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t69) * t69 + (t40 * mrSges(3,1) + t23 * mrSges(7,1) - t41 * mrSges(3,2) - t24 * mrSges(7,2) + Ifges(3,5) * t55 + Ifges(7,5) * t29 + Ifges(3,6) * t54 + Ifges(7,6) * t28 + (Ifges(7,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t64) * t64 + (t26 * mrSges(4,1) + t12 * mrSges(9,1) - t27 * mrSges(4,2) - t11 * mrSges(9,2) + Ifges(4,5) * t37 + Ifges(9,5) * t50 + Ifges(4,6) * t36 + Ifges(9,6) * t51 + (Ifges(9,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t62) * t62 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t14 + Ifges(6,6) * t19 + Ifges(6,2) * t13 / 0.2e1) * t13 + (t45 * mrSges(7,2) - t23 * mrSges(7,3) + Ifges(7,1) * t29 / 0.2e1) * t29 + m(3) * (t40 ^ 2 + t41 ^ 2 + t52 ^ 2) / 0.2e1 + m(4) * (t26 ^ 2 + t27 ^ 2 + t38) / 0.2e1 + m(8) * (t15 ^ 2 + t16 ^ 2 + t38) / 0.2e1 + m(7) * (t23 ^ 2 + t24 ^ 2 + t45 ^ 2) / 0.2e1 + (-t25 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,2) * t21 / 0.2e1) * t21 + m(5) * (t25 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(9) * (t10 ^ 2 + t11 ^ 2 + t12 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (-t45 * mrSges(7,1) + t24 * mrSges(7,3) + Ifges(7,4) * t29 + Ifges(7,2) * t28 / 0.2e1) * t28 + (t25 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,4) * t21 + Ifges(5,1) * t20 / 0.2e1) * t20 + (-t39 * mrSges(8,1) + t15 * mrSges(8,3) + Ifges(8,2) * t35 / 0.2e1) * t35 + (V_base(3) * mrSges(2,2) - t56 * mrSges(2,3) + Ifges(2,1) * t73 / 0.2e1) * t73 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t19 + Ifges(6,1) * t14 / 0.2e1) * t14 + (t10 * mrSges(9,2) - t12 * mrSges(9,3) + Ifges(9,4) * t51 + Ifges(9,1) * t50 / 0.2e1) * t50 + (t52 * mrSges(3,2) - t40 * mrSges(3,3) + Ifges(3,1) * t55 / 0.2e1) * t55 + (t56 * mrSges(2,1) - t57 * mrSges(2,2) + Ifges(2,5) * t73 + Ifges(2,3) * t106 / 0.2e1) * t106 + (-t39 * mrSges(4,1) + t27 * mrSges(4,3) + Ifges(4,4) * t37 + Ifges(4,2) * t36 / 0.2e1) * t36;
T = t8;
