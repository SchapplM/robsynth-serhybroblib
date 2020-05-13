% Calculate kinetic energy for
% fourbarprisDE1
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
% Datum: 2020-05-07 09:10
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbarprisDE1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(6,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE1_energykin_floatb_twist_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisDE1_energykin_floatb_twist_slag_vp2: qJD has to be [1x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fourbarprisDE1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE1_energykin_floatb_twist_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisDE1_energykin_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbarprisDE1_energykin_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbarprisDE1_energykin_floatb_twist_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:09:17
% EndTime: 2020-05-07 09:09:45
% DurationCPUTime: 15.24s
% Computational Cost: add. (1052->251), mult. (1393->421), div. (84->8), fcn. (56->2), ass. (0->134)
t75 = (pkin(2) ^ 2);
t77 = (pkin(1) ^ 2);
t147 = (t75 - t77);
t74 = (pkin(3) ^ 2);
t131 = (-t74 + t147);
t72 = (qJ(1) ^ 2);
t111 = (t72 - t131);
t173 = V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2;
t167 = 2 * qJ(1);
t49 = qJ(1) + pkin(3);
t175 = pkin(2) * t49;
t163 = pkin(3) * qJ(1);
t41 = 2 * t163;
t115 = t72 + t41;
t146 = (t77 - t74);
t23 = (t75 + t146);
t174 = t115 - t23;
t87 = (Ifges(2,5) * V_base(4) + Ifges(2,6) * V_base(5));
t134 = -pkin(2) + t49;
t135 = -pkin(2) - t49;
t113 = (pkin(1) + t135) * (pkin(1) + t134) * (pkin(1) - t134) * (pkin(1) - t135);
t172 = 1 / pkin(1);
t171 = 2 * pkin(3);
t170 = 4 * pkin(3);
t169 = 1 / t49;
t168 = -2 * Ifges(2,4);
t166 = pkin(3) * mrSges(2,1);
t165 = pkin(3) * Ifges(2,4);
t164 = Ifges(2,5) * pkin(3);
t67 = (t77 * m(2));
t162 = (mrSges(2,2) * t77);
t46 = (V_base(4) ^ 2);
t161 = Ifges(2,4) * t46;
t160 = t23 * t46;
t40 = -3 * t74 + t75;
t52 = (Ifges(2,2) - Ifges(2,1));
t159 = t40 * t52;
t44 = (V_base(6) ^ 2);
t158 = (t44 * t77);
t157 = t46 * t52;
t58 = pkin(2) + pkin(3);
t59 = pkin(2) - pkin(3);
t156 = t59 ^ 2 * t58 ^ 2;
t155 = t74 * mrSges(2,1);
t50 = (Ifges(4,2) - Ifges(4,1));
t154 = t74 * t50;
t153 = t75 * t52;
t152 = t75 * t77;
t151 = t77 * mrSges(2,1);
t150 = t77 * mrSges(2,3);
t54 = t75 * mrSges(2,1);
t39 = t54 / 0.2e1;
t53 = Ifges(2,3) * pkin(3);
t149 = t39 + t53;
t45 = (V_base(5) ^ 2);
t148 = -t45 + t46;
t36 = Ifges(2,2) / 0.2e1 + Ifges(2,1) / 0.2e1;
t145 = pkin(1) * qJD(1);
t144 = t172 / 0.2e1;
t110 = t146 - t115;
t12 = -t75 + t110;
t143 = qJD(1) * t12;
t142 = pkin(1) * V_base(1);
t141 = pkin(1) * V_base(2);
t140 = pkin(1) * V_base(3);
t139 = mrSges(2,1) * V_base(1);
t136 = pkin(3) * t162;
t133 = Ifges(4,3) * t152;
t132 = t59 * t58 * t52;
t130 = -(t172 * t169) / 0.2e1;
t129 = t169 * t144;
t128 = V_base(5) * pkin(1);
t127 = V_base(6) * pkin(1);
t126 = -t156 / 0.4e1;
t125 = pkin(2) * mrSges(4,3) / t144;
t37 = Ifges(2,1) / 0.4e1 - Ifges(2,2) / 0.4e1;
t121 = t23 * V_base(5);
t120 = (V_base(4) * mrSges(2,2));
t119 = V_base(6) * mrSges(2,1);
t118 = (-t74 + t75) * pkin(3);
t117 = t49 * t140;
t116 = -2 * mrSges(2,2) * t158;
t114 = t54 - t155;
t112 = -2 * t127;
t109 = pkin(2) * t127;
t108 = 2 * V_base(6);
t106 = Ifges(2,6) * t127;
t105 = V_base(4) * pkin(1) * mrSges(2,3);
t104 = -2 * t113;
t103 = 4 * t113;
t102 = V_base(4) * V_base(5);
t100 = t41 + t111;
t99 = Ifges(4,5) * t109;
t98 = pkin(3) * t106;
t97 = Ifges(4,6) * t109;
t96 = V_base(4) * t127;
t18 = V_base(2) + t127;
t95 = Ifges(4,4) * t102;
t94 = t113 * t145;
t92 = pkin(3) * t96;
t35 = Ifges(4,2) / 0.4e1 - Ifges(4,1) / 0.4e1;
t76 = t77 ^ 2;
t91 = t126 * t50 - t35 * t76;
t90 = V_base(1) * t100;
t89 = t174 * V_base(6);
t88 = Ifges(2,5) * V_base(5) - Ifges(2,6) * V_base(4);
t86 = Ifges(4,5) * V_base(4) + Ifges(4,6) * V_base(5);
t85 = (-t111 - 2 * t163) * t129;
t84 = -mrSges(4,1) * V_base(1) - mrSges(4,2) * V_base(2) + Ifges(4,5) * V_base(5);
t79 = sqrt(-t113);
t71 = qJ(1) * t72;
t70 = t72 ^ 2;
t69 = qJD(1) ^ 2;
t66 = pkin(3) * m(2);
t51 = Ifges(2,2) + Ifges(2,1);
t43 = -3 * t166;
t42 = t49 ^ 2;
t38 = -t155 / 0.2e1;
t34 = -mrSges(2,1) + 2 * t66;
t22 = t36 * t74;
t21 = (Ifges(4,2) / 0.2e1 + Ifges(4,1) / 0.2e1) * t75;
t20 = t67 - t166;
t19 = (t66 - mrSges(2,1) / 0.2e1) * t77;
t17 = pkin(3) * t52 - t151;
t15 = Ifges(4,4) * t121;
t14 = -t75 - t110;
t10 = t121 * t50 + 2 * t97;
t6 = V_base(6) + t169 / t79 * t143;
t5 = t129 * t79 * V_base(4) + t85 * V_base(5);
t4 = t130 * t79 * V_base(5) + t85 * V_base(4);
t3 = (-t18 * t79 - t90) * t129 + qJD(1);
t2 = qJ(1) * t5 + t128 - V_base(3);
t1 = (-(t100 * t18) + t79 * V_base(1)) * t130 - t6 * qJ(1);
t7 = m(1) * t173 / 0.2e1 + m(3) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + ((V_base(2) * mrSges(1,1)) - (V_base(1) * mrSges(1,2)) + (Ifges(1,3) * V_base(6)) / 0.2e1) * V_base(6) + (-(V_base(3) * mrSges(1,1)) + (V_base(1) * mrSges(1,3)) + (Ifges(1,6) * V_base(6)) + (Ifges(1,2) * V_base(5)) / 0.2e1) * V_base(5) + ((V_base(3) * mrSges(1,2)) - (V_base(2) * mrSges(1,3)) + (Ifges(1,4) * V_base(5)) + (Ifges(1,5) * V_base(6)) + (Ifges(1,1) * V_base(4)) / 0.2e1) * V_base(4) + (((0.4e1 * ((mrSges(2,1) * t100 * V_base(2)) - (mrSges(2,2) * t90) + ((pkin(1) * t119 + t87) * t72) + ((-(Ifges(2,3) - t166) * t127 + t87 * pkin(3)) * t167) + (-t151 / 0.2e1 + t38 + t149) * t112 - (t87 * t131)) * pkin(1) * t143 + ((2 * (mrSges(2,1) * V_base(4) + mrSges(2,2) * V_base(5)) * t117 + (Ifges(2,4) * t148 + t102 * t52) * t72 + (t161 * t171 + t116 + t88 * pkin(1) * t108 + 2 * (-t162 - t165) * t45 + 2 * t17 * t102) * qJ(1) + (t131 * Ifges(2,4) - 2 * t136) * t45 + (((t52 - 2 * t166) * t77 - t132) * V_base(4) + 2 * t127 * t164) * V_base(5) - t131 * t161 - 2 * Ifges(2,6) * t92 + pkin(3) * t116 + (mrSges(2,2) * V_base(2) + t139) * t49 * t112) * t104)) * t79 + (4 * (-(-mrSges(2,2) * t18 - t139 + t88) * t94 - Ifges(2,3) * t77 * t69 * t12) * t12) + (-(((-mrSges(2,1) * V_base(5) + t120) * t72 + (pkin(3) * t120 + t20 * V_base(5)) * t167 + (t34 * t77 + t114) * V_base(5) - t131 * t120) * t117) + (-t37 * t45 + (Ifges(2,4) * t102) - t157 / 0.4e1) * t70 + ((t17 * t45 + ((t162 + 4 * t165) * V_base(4) - t106) * V_base(5) - pkin(3) * t157 - Ifges(2,5) * t96 - t44 * t151) * t71) + (((t76 * m(2)) + (t43 + t36) * t77 - t159 / 0.2e1) * t45 + (((t168 * t40 + 3 * t136) * V_base(4) - 3 * t98) * V_base(5)) + (t36 * t77 + t159 / 0.2e1) * t46 - 0.3e1 * (0.2e1 / 0.3e1 * t150 + t164) * t96 + ((Ifges(2,3) + t43 + t67) * t158)) * t72 + (((t34 * t76) + (pkin(3) * t51 - (3 * t155) + t54) * t77 - (pkin(3) * t132)) * t45 + (((-4 * t118 * Ifges(2,4) + (-t40 * t77 + t76) * mrSges(2,2)) * V_base(4) - (3 * t74 - t147) * t106) * V_base(5)) + (t51 * t77 + t132) * pkin(3) * t46 - (((mrSges(2,3) * t170 + Ifges(2,5)) * t77 - Ifges(2,5) * t40) * t96) + 0.2e1 * (t19 - 0.3e1 / 0.2e1 * t155 + t149) * t158) * qJ(1) + (((m(2) * t74) - t166 - t37) * t76 + (t22 - t153 / 0.2e1 + (t118 * mrSges(2,1))) * t77 + (t52 * t156) / 0.4e1) * t45 + ((((pkin(3) * mrSges(2,2) + Ifges(2,4)) * t76 + (-mrSges(2,2) * t118 + t168 * t75) * t77 + Ifges(2,4) * t156) * V_base(4) + t131 * t98) * V_base(5)) + (t37 * t76 + (t22 + t153 / 0.2e1) * t77 + t52 * t126) * t46 - ((-Ifges(2,5) * t131 + t150 * t171) * t92) + (pkin(3) * ((-mrSges(2,1) + t66) * t77 + t53 + t114) * t158) + t173 * t42 * t67 + (-((t72 * t119) + ((-t20 * V_base(6) + t105) * t167) + (t105 * t171) - 0.2e1 * (t19 + t39 + t38) * V_base(6)) * t141 + (((t167 + t171) * mrSges(2,3) * t128 + (pkin(3) * t167 + t111) * mrSges(2,2) * V_base(6)) * t142)) * t49) * t103) / (t49 ^ 2) + ((-8 * (-2 * Ifges(4,3) * t109 + (mrSges(4,1) * V_base(2) - mrSges(4,2) * V_base(1)) * t14 + t174 * t86) * t145 * t175 + (-Ifges(4,4) * t160 - t10 * V_base(4) + t15 * V_base(5) + t115 * (Ifges(4,4) * t148 + t102 * t50) + (2 * (mrSges(4,1) * V_base(4) + mrSges(4,2) * V_base(5)) * V_base(3) + t84 * t108) * pkin(2) * pkin(1)) * t104) * t79 - (16 * t42 * t69 * t133) + (8 * (-Ifges(4,6) * V_base(4) + t84) * t94 * t175) + ((-(t86 * t109) + (-(2 * t95) + (t46 / 0.2e1 - t45 / 0.2e1) * t50) * (t77 + t40)) * t72 - 0.4e1 * (-(t50 * t160) / 0.4e1 + (t15 + t99 / 0.2e1) * V_base(4) + (t10 * V_base(5)) / 0.4e1) * t163 + ((t21 + t154 / 0.2e1) * t77 + t91) * t46 + ((Ifges(4,4) * (-2 * t74 * t77 + t156 + t76) * V_base(5) + t23 * t99) * V_base(4)) + ((t21 - t154 / 0.2e1) * t77 - t91) * t45 + (t97 * t121) + (t44 * t133) + (t71 * t170 + t70) * (-t35 * t46 + t95 + (t45 * t50) / 0.4e1) + t173 * m(4) * t152 + (-(t14 * (-mrSges(4,1) * V_base(5) + mrSges(4,2) * V_base(4)) * t140) + ((mrSges(4,2) * t89) + t125 * V_base(5)) * t142 - ((mrSges(4,1) * t89) + t125 * V_base(4)) * t141) * pkin(2)) * t103) / (pkin(2) ^ 2)) / (pkin(1) ^ 2) / t113 / 0.8e1 + (t3 * mrSges(3,1) - t1 * mrSges(3,3) + Ifges(3,2) * t6 / 0.2e1) * t6 + (-t3 * mrSges(3,2) + t2 * mrSges(3,3) + Ifges(3,4) * t6 + Ifges(3,1) * t5 / 0.2e1) * t5 + (t2 * mrSges(3,1) - t1 * mrSges(3,2) - Ifges(3,5) * t5 - Ifges(3,6) * t6 + Ifges(3,3) * t4 / 0.2e1) * t4;
T = t7;
