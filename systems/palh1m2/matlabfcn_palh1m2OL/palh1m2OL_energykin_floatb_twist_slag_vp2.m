% Calculate kinetic energy for
% palh1m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% qJD [13x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [11x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:30
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m2OL_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(6,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_energykin_floatb_twist_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m2OL_energykin_floatb_twist_slag_vp2: qJD has to be [13x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh1m2OL_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_energykin_floatb_twist_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2OL_energykin_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2OL_energykin_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2OL_energykin_floatb_twist_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 21:17:17
% EndTime: 2020-05-02 21:17:23
% DurationCPUTime: 4.92s
% Computational Cost: add. (2741->284), mult. (3987->425), div. (0->0), fcn. (3538->22), ass. (0->117)
t127 = sin(qJ(1));
t136 = cos(qJ(1));
t92 = t127 * V_base(5) + t136 * V_base(4);
t140 = t127 * V_base(2) + t136 * V_base(1);
t114 = V_base(6) + qJD(1);
t125 = sin(qJ(3));
t134 = cos(qJ(3));
t126 = sin(qJ(2));
t135 = cos(qJ(2));
t95 = pkin(13) * t135 - pkin(15) * t126 + pkin(1);
t97 = pkin(13) * t126 + pkin(15) * t135;
t152 = t95 * t125 + t134 * t97;
t156 = pkin(5) + t152;
t109 = pkin(1) * t125 + pkin(5);
t124 = sin(qJ(4));
t133 = cos(qJ(4));
t90 = t125 * t135 + t126 * t134;
t148 = t134 * t135;
t91 = -t125 * t126 + t148;
t143 = -t124 * t90 + t133 * t91;
t149 = pkin(5) * qJD(3);
t154 = pkin(1) * t134;
t55 = t124 * t91 + t133 * t90;
t61 = -t125 * t97 + t134 * t95;
t89 = -t127 * V_base(4) + t136 * V_base(5);
t11 = -t140 * t55 + qJD(2) * (t109 * t133 + t124 * t154) + t133 * t149 - (t61 * t124 + t156 * t133) * t89 + t143 * V_base(3);
t82 = qJD(2) - t89;
t77 = qJD(7) + t82;
t155 = pkin(4) * t77;
t153 = pkin(5) * t134;
t151 = sin(qJ(10));
t150 = pkin(1) * qJD(2);
t120 = sin(qJ(8));
t129 = cos(qJ(8));
t100 = V_base(5) * pkin(13) + V_base(1);
t101 = -V_base(4) * pkin(13) + V_base(2);
t67 = t100 * t136 + t101 * t127;
t71 = -pkin(15) * t89 + V_base(3);
t46 = -t126 * t67 + t135 * t71;
t47 = -t126 * t71 - t135 * t67;
t26 = -t120 * t46 + t129 * t47;
t66 = -t127 * t100 + t101 * t136;
t142 = pkin(14) * t89 + V_base(3);
t141 = V_base(1) * t127 - t136 * V_base(2);
t63 = t114 * t135 - t126 * t92;
t65 = -t114 * t126 - t135 * t92;
t37 = t125 * t63 - t134 * t65;
t40 = t125 * t65 + t134 * t63;
t19 = -t124 * t37 + t133 * t40;
t10 = t55 * V_base(3) + (t109 * t124 - t133 * t154) * qJD(2) + t124 * t149 + t89 * (-t156 * t124 + t61 * t133) + t140 * t143;
t58 = -pkin(15) * t114 - t66;
t78 = qJD(3) + t82;
t76 = qJD(8) + t82;
t42 = t141 + t92 * (pkin(1) * t135 + pkin(13)) + t114 * (pkin(1) * t126 - pkin(15));
t117 = pkin(13) - pkin(16);
t139 = t117 * t89 + t140;
t108 = pkin(5) * t125 + pkin(1);
t24 = -pkin(15) * V_base(6) + t141 + (t108 * t92 - V_base(6) * t153) * t135 + (t108 * V_base(6) + t92 * t153) * t126 + (-pkin(5) * t148 + t108 * t126 - pkin(15)) * qJD(1) + t92 * pkin(13);
t121 = sin(qJ(7));
t130 = cos(qJ(7));
t87 = -t121 * t126 + t130 * t135;
t88 = t121 * t135 + t126 * t130;
t29 = -t140 * t88 - (t121 * t95 + t130 * t97) * t89 + t121 * t150 + t87 * V_base(3);
t28 = -t140 * t87 - (-t121 * t97 + t130 * t95) * t89 + t130 * t150 - t88 * V_base(3);
t137 = V_base(3) ^ 2;
t132 = cos(qJ(5));
t131 = cos(qJ(6));
t128 = cos(qJ(9));
t123 = sin(qJ(5));
t122 = sin(qJ(6));
t119 = sin(qJ(9));
t118 = cos(qJ(10));
t116 = cos(pkin(19));
t115 = sin(pkin(19));
t84 = -t115 * t151 - t116 * t118;
t83 = t115 * t118 - t116 * t151;
t81 = qJD(6) - t89;
t73 = qJD(4) + t78;
t72 = qJD(9) + t76;
t70 = qJD(10) + t77;
t64 = t114 * t122 + t131 * t92;
t62 = t114 * t131 - t122 * t92;
t57 = t58 ^ 2;
t50 = t114 * pkin(14) + t92 * t117 + t141;
t45 = t142 * t122 + t139 * t131;
t44 = -t139 * t122 + t142 * t131;
t41 = t42 ^ 2;
t39 = t121 * t65 + t130 * t63;
t38 = t120 * t65 + t129 * t63;
t36 = -t121 * t63 + t130 * t65;
t35 = -t120 * t63 + t129 * t65;
t32 = -pkin(2) * t35 + t58;
t31 = -t134 * t150 + t140 * t91 + t61 * t89 + t90 * V_base(3);
t30 = t125 * t150 - t140 * t90 - t152 * t89 + t91 * V_base(3);
t27 = t120 * t47 + t129 * t46;
t25 = pkin(2) * t76 + t26;
t23 = t115 * t155 + t29;
t22 = t116 * t155 + t28;
t21 = t124 * t40 + t133 * t37;
t20 = -t119 * t35 - t128 * t38;
t18 = t119 * t38 - t128 * t35;
t17 = qJD(5) - t19;
t16 = t36 * t83 + t39 * t84;
t15 = t36 * t84 - t39 * t83;
t14 = t123 * t73 + t132 * t21;
t13 = -t123 * t21 + t132 * t73;
t12 = (-t115 * t39 - t116 * t36) * pkin(4) + t42;
t9 = -t119 * t25 - t128 * t27;
t8 = t119 * t27 - t128 * t25;
t7 = -pkin(9) * t73 - t11;
t6 = pkin(11) * t73 + t10;
t5 = t22 * t83 + t23 * t84;
t4 = t22 * t84 - t23 * t83;
t3 = -pkin(9) * t19 - pkin(11) * t21 + t24;
t2 = t123 * t3 + t132 * t6;
t1 = -t123 * t6 + t132 * t3;
t33 = (t24 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,5) * t73 + Ifges(5,1) * t21 / 0.2e1) * t21 + (-t7 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t14 + Ifges(6,6) * t17 + Ifges(6,2) * t13 / 0.2e1) * t13 + (t11 * mrSges(5,1) - t10 * mrSges(5,2) + Ifges(5,3) * t73 / 0.2e1) * t73 + (t66 * mrSges(2,1) - t67 * mrSges(2,2) + Ifges(2,5) * t92 + Ifges(2,6) * t89 + Ifges(2,3) * t114 / 0.2e1) * t114 + (-t50 * mrSges(7,1) + t45 * mrSges(7,3) + Ifges(7,4) * t64 + Ifges(7,6) * t81 + Ifges(7,2) * t62 / 0.2e1) * t62 + (-t42 * mrSges(4,1) + t31 * mrSges(4,3) + Ifges(4,6) * t78 + Ifges(4,2) * t40 / 0.2e1) * t40 + (-t58 * mrSges(3,1) + t46 * mrSges(3,3) + Ifges(3,6) * t82 + Ifges(3,2) * t65 / 0.2e1) * t65 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t58 * mrSges(9,2) - t26 * mrSges(9,3) + Ifges(9,5) * t76 + Ifges(9,1) * t38 / 0.2e1) * t38 + (t47 * mrSges(3,1) - t46 * mrSges(3,2) + Ifges(3,3) * t82 / 0.2e1) * t82 + (t12 * mrSges(11,2) - t4 * mrSges(11,3) + Ifges(11,5) * t70 + Ifges(11,1) * t16 / 0.2e1) * t16 + m(2) * (t66 ^ 2 + t67 ^ 2 + t137) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t137) / 0.2e1 + m(8) * (t28 ^ 2 + t29 ^ 2 + t41) / 0.2e1 + m(4) * (t30 ^ 2 + t31 ^ 2 + t41) / 0.2e1 + m(7) * (t44 ^ 2 + t45 ^ 2 + t50 ^ 2) / 0.2e1 + m(3) * (t46 ^ 2 + t47 ^ 2 + t57) / 0.2e1 + m(9) * (t26 ^ 2 + t27 ^ 2 + t57) / 0.2e1 + m(10) * (t32 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t11 ^ 2 + t24 ^ 2) / 0.2e1 + (t42 * mrSges(8,2) - t28 * mrSges(8,3) + Ifges(8,5) * t77 + Ifges(8,1) * t39 / 0.2e1) * t39 + (t58 * mrSges(3,2) - t47 * mrSges(3,3) + Ifges(3,4) * t65 + Ifges(3,5) * t82 + Ifges(3,1) * t63 / 0.2e1) * t63 + (t7 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t17 + Ifges(6,1) * t14 / 0.2e1) * t14 + (t42 * mrSges(4,2) - t30 * mrSges(4,3) + Ifges(4,4) * t40 + Ifges(4,5) * t78 + Ifges(4,1) * t37 / 0.2e1) * t37 + (-t42 * mrSges(8,1) + t29 * mrSges(8,3) + Ifges(8,4) * t39 + Ifges(8,6) * t77 + Ifges(8,2) * t36 / 0.2e1) * t36 + (t30 * mrSges(4,1) - t31 * mrSges(4,2) + Ifges(4,3) * t78 / 0.2e1) * t78 + (-t58 * mrSges(9,1) + t27 * mrSges(9,3) + Ifges(9,4) * t38 + Ifges(9,6) * t76 + Ifges(9,2) * t35 / 0.2e1) * t35 + (-t24 * mrSges(5,1) + t10 * mrSges(5,3) + Ifges(5,4) * t21 + Ifges(5,6) * t73 + Ifges(5,2) * t19 / 0.2e1) * t19 + (-t32 * mrSges(10,1) + t9 * mrSges(10,3) + Ifges(10,4) * t20 + Ifges(10,6) * t72 + Ifges(10,2) * t18 / 0.2e1) * t18 + (t8 * mrSges(10,1) - t9 * mrSges(10,2) + Ifges(10,3) * t72 / 0.2e1) * t72 + m(6) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(11) * (t12 ^ 2 + t4 ^ 2 + t5 ^ 2) / 0.2e1 + (t50 * mrSges(7,2) - t44 * mrSges(7,3) + Ifges(7,5) * t81 + Ifges(7,1) * t64 / 0.2e1) * t64 + (-t12 * mrSges(11,1) + t5 * mrSges(11,3) + Ifges(11,4) * t16 + Ifges(11,6) * t70 + Ifges(11,2) * t15 / 0.2e1) * t15 + (t26 * mrSges(9,1) - t27 * mrSges(9,2) + Ifges(9,3) * t76 / 0.2e1) * t76 + (t28 * mrSges(8,1) - t29 * mrSges(8,2) + Ifges(8,3) * t77 / 0.2e1) * t77 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t44 * mrSges(7,1) - t45 * mrSges(7,2) + Ifges(7,3) * t81 / 0.2e1) * t81 + (t4 * mrSges(11,1) - t5 * mrSges(11,2) + Ifges(11,3) * t70 / 0.2e1) * t70 + (V_base(3) * mrSges(2,2) - t66 * mrSges(2,3) + Ifges(2,1) * t92 / 0.2e1) * t92 + (-V_base(3) * mrSges(2,1) + t67 * mrSges(2,3) + Ifges(2,4) * t92 + Ifges(2,2) * t89 / 0.2e1) * t89 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t32 * mrSges(10,2) - t8 * mrSges(10,3) + Ifges(10,5) * t72 + Ifges(10,1) * t20 / 0.2e1) * t20 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t17 / 0.2e1) * t17;
T = t33;
