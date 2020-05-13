% Calculate kinetic energy for
% fourbarprisTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% m [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:01
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbarprisTE_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(6,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisTE_energykin_floatb_twist_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisTE_energykin_floatb_twist_slag_vp2: qJD has to be [1x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fourbarprisTE_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisTE_energykin_floatb_twist_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisTE_energykin_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbarprisTE_energykin_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbarprisTE_energykin_floatb_twist_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:00:42
% EndTime: 2020-05-07 09:01:14
% DurationCPUTime: 15.16s
% Computational Cost: add. (1052->251), mult. (1393->421), div. (84->8), fcn. (56->2), ass. (0->134)
t73 = (pkin(2) ^ 2);
t75 = (pkin(1) ^ 2);
t145 = (t73 - t75);
t72 = (pkin(3) ^ 2);
t129 = (-t72 + t145);
t70 = (qJ(1) ^ 2);
t109 = (t70 - t129);
t171 = V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2;
t165 = 2 * qJ(1);
t49 = qJ(1) + pkin(3);
t173 = pkin(2) * t49;
t161 = pkin(3) * qJ(1);
t41 = 2 * t161;
t113 = t70 + t41;
t144 = (t75 - t72);
t23 = (t73 + t144);
t172 = t113 - t23;
t85 = (Ifges(2,5) * V_base(4) + Ifges(2,6) * V_base(5));
t132 = -pkin(2) + t49;
t133 = -pkin(2) - t49;
t111 = (pkin(1) + t133) * (pkin(1) + t132) * (pkin(1) - t132) * (pkin(1) - t133);
t170 = 1 / pkin(1);
t169 = 2 * pkin(3);
t168 = 4 * pkin(3);
t167 = 1 / t49;
t166 = -2 * Ifges(2,4);
t164 = pkin(3) * mrSges(2,1);
t163 = pkin(3) * Ifges(2,4);
t162 = Ifges(2,5) * pkin(3);
t65 = (t75 * m(2));
t160 = (mrSges(2,2) * t75);
t46 = (V_base(4) ^ 2);
t159 = Ifges(2,4) * t46;
t158 = t23 * t46;
t40 = -3 * t72 + t73;
t52 = (Ifges(2,2) - Ifges(2,1));
t157 = t40 * t52;
t156 = t46 * t52;
t58 = pkin(2) + pkin(3);
t59 = pkin(2) - pkin(3);
t155 = t59 ^ 2 * t58 ^ 2;
t154 = t72 * mrSges(2,1);
t50 = (Ifges(4,2) - Ifges(4,1));
t153 = t72 * t50;
t152 = t73 * t52;
t151 = t73 * t75;
t150 = t75 * mrSges(2,1);
t149 = t75 * mrSges(2,3);
t44 = (V_base(6) ^ 2);
t148 = t75 * t44;
t54 = t73 * mrSges(2,1);
t39 = t54 / 0.2e1;
t53 = Ifges(2,3) * pkin(3);
t147 = t39 + t53;
t45 = (V_base(5) ^ 2);
t146 = -t45 + t46;
t37 = Ifges(2,2) / 0.2e1 + Ifges(2,1) / 0.2e1;
t143 = pkin(1) * qJD(1);
t142 = t170 / 0.2e1;
t108 = t144 - t113;
t12 = -t73 + t108;
t141 = qJD(1) * t12;
t140 = pkin(1) * V_base(1);
t139 = pkin(1) * V_base(2);
t138 = pkin(1) * V_base(3);
t137 = V_base(1) * mrSges(2,1);
t136 = pkin(3) * t160;
t131 = Ifges(4,3) * t151;
t130 = t59 * t58 * t52;
t128 = -(t170 * t167) / 0.2e1;
t127 = t167 * t142;
t126 = V_base(5) * pkin(1);
t125 = V_base(6) * pkin(1);
t124 = -t155 / 0.4e1;
t123 = pkin(2) * mrSges(4,3) / t142;
t36 = Ifges(2,2) / 0.4e1 - Ifges(2,1) / 0.4e1;
t120 = t23 * V_base(5);
t119 = (V_base(4) * mrSges(2,2));
t118 = V_base(6) * mrSges(2,1);
t116 = (-t72 + t73) * pkin(3);
t115 = t49 * t138;
t114 = -2 * mrSges(2,2) * t148;
t112 = t54 - t154;
t110 = -2 * t125;
t107 = pkin(2) * t125;
t106 = 2 * V_base(6);
t104 = Ifges(2,6) * t125;
t103 = V_base(4) * pkin(1) * mrSges(2,3);
t102 = -2 * t111;
t101 = 4 * t111;
t100 = V_base(4) * V_base(5);
t98 = t41 + t109;
t97 = Ifges(4,5) * t107;
t96 = pkin(3) * t104;
t95 = Ifges(4,6) * t107;
t94 = V_base(4) * t125;
t18 = V_base(2) + t125;
t93 = Ifges(4,4) * t100;
t92 = t111 * t143;
t90 = pkin(3) * t94;
t35 = Ifges(4,2) / 0.4e1 - Ifges(4,1) / 0.4e1;
t74 = t75 ^ 2;
t89 = t124 * t50 - t35 * t74;
t88 = V_base(1) * t98;
t87 = t172 * V_base(6);
t86 = Ifges(2,5) * V_base(5) - Ifges(2,6) * V_base(4);
t84 = Ifges(4,5) * V_base(4) + Ifges(4,6) * V_base(5);
t83 = (-t109 - 2 * t161) * t127;
t82 = -V_base(1) * mrSges(4,1) - V_base(2) * mrSges(4,2) + Ifges(4,5) * V_base(5);
t77 = sqrt(-t111);
t69 = qJ(1) * t70;
t68 = t70 ^ 2;
t67 = qJD(1) ^ 2;
t64 = pkin(3) * m(2);
t51 = Ifges(2,2) + Ifges(2,1);
t43 = -3 * t164;
t42 = t49 ^ 2;
t38 = -t154 / 0.2e1;
t34 = -mrSges(2,1) + 2 * t64;
t22 = t37 * t72;
t21 = (Ifges(4,2) / 0.2e1 + Ifges(4,1) / 0.2e1) * t73;
t20 = t65 - t164;
t19 = (t64 - mrSges(2,1) / 0.2e1) * t75;
t17 = pkin(3) * t52 - t150;
t15 = Ifges(4,4) * t120;
t14 = -t73 - t108;
t10 = t120 * t50 + 2 * t95;
t6 = V_base(6) + t167 / t77 * t141;
t5 = t127 * t77 * V_base(4) + t83 * V_base(5);
t4 = t128 * t77 * V_base(5) + t83 * V_base(4);
t3 = (-t18 * t77 - t88) * t127 + qJD(1);
t2 = qJ(1) * t5 + t126 - V_base(3);
t1 = (-(t18 * t98) + t77 * V_base(1)) * t128 - t6 * qJ(1);
t7 = m(1) * t171 / 0.2e1 + m(3) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + ((V_base(2) * mrSges(1,1)) - (V_base(1) * mrSges(1,2)) + (Ifges(1,3) * V_base(6)) / 0.2e1) * V_base(6) + (-(V_base(3) * mrSges(1,1)) + (V_base(1) * mrSges(1,3)) + (Ifges(1,6) * V_base(6)) + (Ifges(1,2) * V_base(5)) / 0.2e1) * V_base(5) + ((V_base(3) * mrSges(1,2)) - (V_base(2) * mrSges(1,3)) + (Ifges(1,4) * V_base(5)) + (Ifges(1,5) * V_base(6)) + (Ifges(1,1) * V_base(4)) / 0.2e1) * V_base(4) + (((0.4e1 * ((mrSges(2,1) * t98 * V_base(2)) - (mrSges(2,2) * t88) + ((pkin(1) * t118 + t85) * t70) + ((-(Ifges(2,3) - t164) * t125 + t85 * pkin(3)) * t165) + (-t150 / 0.2e1 + t38 + t147) * t110 - (t85 * t129)) * pkin(1) * t141 + ((2 * (V_base(4) * mrSges(2,1) + V_base(5) * mrSges(2,2)) * t115 + (Ifges(2,4) * t146 + t100 * t52) * t70 + (t159 * t169 + t114 + t86 * pkin(1) * t106 + 2 * (-t160 - t163) * t45 + 2 * t17 * t100) * qJ(1) + (Ifges(2,4) * t129 - 2 * t136) * t45 + (((t52 - 2 * t164) * t75 - t130) * V_base(4) + 2 * t125 * t162) * V_base(5) - t129 * t159 - 2 * Ifges(2,6) * t90 + pkin(3) * t114 + (mrSges(2,2) * V_base(2) + t137) * t49 * t110) * t102)) * t77 + (4 * (-(-mrSges(2,2) * t18 - t137 + t86) * t92 - Ifges(2,3) * t75 * t67 * t12) * t12) + (-(((-V_base(5) * mrSges(2,1) + t119) * t70 + (pkin(3) * t119 + t20 * V_base(5)) * t165 + (t34 * t75 + t112) * V_base(5) - t129 * t119) * t115) + (t36 * t45 + (Ifges(2,4) * t100) - t156 / 0.4e1) * t68 + ((t17 * t45 + ((t160 + 4 * t163) * V_base(4) - t104) * V_base(5) - pkin(3) * t156 - Ifges(2,5) * t94 - mrSges(2,1) * t148) * t69) + (((t74 * m(2)) + (t43 + t37) * t75 - t157 / 0.2e1) * t45 + (((t166 * t40 + 3 * t136) * V_base(4) - 3 * t96) * V_base(5)) + (t37 * t75 + t157 / 0.2e1) * t46 - 0.3e1 * (0.2e1 / 0.3e1 * t149 + t162) * t94 + ((Ifges(2,3) + t43 + t65) * t148)) * t70 + (((t34 * t74) + (pkin(3) * t51 - (3 * t154) + t54) * t75 - (pkin(3) * t130)) * t45 + (((-4 * t116 * Ifges(2,4) + (-t40 * t75 + t74) * mrSges(2,2)) * V_base(4) - (3 * t72 - t145) * t104) * V_base(5)) + pkin(3) * (t51 * t75 + t130) * t46 - (((mrSges(2,3) * t168 + Ifges(2,5)) * t75 - Ifges(2,5) * t40) * t94) + 0.2e1 * (t19 - 0.3e1 / 0.2e1 * t154 + t147) * t148) * qJ(1) + (((t72 * m(2)) - t164 + t36) * t74 + (t22 - t152 / 0.2e1 + (t116 * mrSges(2,1))) * t75 + (t52 * t155) / 0.4e1) * t45 + ((((pkin(3) * mrSges(2,2) + Ifges(2,4)) * t74 + (-mrSges(2,2) * t116 + t166 * t73) * t75 + Ifges(2,4) * t155) * V_base(4) + t129 * t96) * V_base(5)) + (-t36 * t74 + (t22 + t152 / 0.2e1) * t75 + t52 * t124) * t46 - ((-Ifges(2,5) * t129 + t149 * t169) * t90) + (pkin(3) * ((-mrSges(2,1) + t64) * t75 + t53 + t112) * t148) + t171 * t42 * t65 + (-((t70 * t118) + ((-t20 * V_base(6) + t103) * t165) + (t103 * t169) - 0.2e1 * V_base(6) * (t19 + t39 + t38)) * t139 + (((t165 + t169) * mrSges(2,3) * t126 + (pkin(3) * t165 + t109) * V_base(6) * mrSges(2,2)) * t140)) * t49) * t101) / (t49 ^ 2) + ((-8 * (-2 * Ifges(4,3) * t107 + (mrSges(4,1) * V_base(2) - mrSges(4,2) * V_base(1)) * t14 + t172 * t84) * t143 * t173 + (-Ifges(4,4) * t158 - t10 * V_base(4) + V_base(5) * t15 + t113 * (Ifges(4,4) * t146 + t100 * t50) + (2 * (V_base(4) * mrSges(4,1) + V_base(5) * mrSges(4,2)) * V_base(3) + t82 * t106) * pkin(2) * pkin(1)) * t102) * t77 - (16 * t42 * t67 * t131) + (8 * (-Ifges(4,6) * V_base(4) + t82) * t92 * t173) + ((-(t84 * t107) + (-(2 * t93) + (t46 / 0.2e1 - t45 / 0.2e1) * t50) * (t75 + t40)) * t70 - 0.4e1 * (-(t50 * t158) / 0.4e1 + (t15 + t97 / 0.2e1) * V_base(4) + (V_base(5) * t10) / 0.4e1) * t161 + ((t21 + t153 / 0.2e1) * t75 + t89) * t46 + (((-2 * t75 * t72 + t155 + t74) * Ifges(4,4) * V_base(5) + t23 * t97) * V_base(4)) + ((t21 - t153 / 0.2e1) * t75 - t89) * t45 + (t95 * t120) + (t44 * t131) + (t69 * t168 + t68) * (-t35 * t46 + t93 + (t45 * t50) / 0.4e1) + t171 * m(4) * t151 + (-(t14 * (-V_base(5) * mrSges(4,1) + V_base(4) * mrSges(4,2)) * t138) + ((mrSges(4,2) * t87) + t123 * V_base(5)) * t140 - ((mrSges(4,1) * t87) + t123 * V_base(4)) * t139) * pkin(2)) * t101) / (pkin(2) ^ 2)) / (pkin(1) ^ 2) / t111 / 0.8e1 + (t3 * mrSges(3,1) - t1 * mrSges(3,3) + Ifges(3,2) * t6 / 0.2e1) * t6 + (-t3 * mrSges(3,2) + t2 * mrSges(3,3) + Ifges(3,4) * t6 + Ifges(3,1) * t5 / 0.2e1) * t5 + (t2 * mrSges(3,1) - t1 * mrSges(3,2) - Ifges(3,5) * t5 - Ifges(3,6) * t6 + Ifges(3,3) * t4 / 0.2e1) * t4;
T = t7;
