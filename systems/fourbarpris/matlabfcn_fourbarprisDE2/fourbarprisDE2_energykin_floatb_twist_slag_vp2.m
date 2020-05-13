% Calculate kinetic energy for
% fourbarprisDE2
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
% Datum: 2020-05-07 09:45
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbarprisDE2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(6,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE2_energykin_floatb_twist_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisDE2_energykin_floatb_twist_slag_vp2: qJD has to be [1x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fourbarprisDE2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE2_energykin_floatb_twist_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisDE2_energykin_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbarprisDE2_energykin_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbarprisDE2_energykin_floatb_twist_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:44:30
% EndTime: 2020-05-07 09:44:51
% DurationCPUTime: 12.42s
% Computational Cost: add. (1052->252), mult. (1393->423), div. (84->8), fcn. (56->2), ass. (0->136)
t163 = 2 * qJ(1);
t49 = qJ(1) + pkin(3);
t171 = pkin(2) * t49;
t159 = pkin(3) * qJ(1);
t41 = 2 * t159;
t68 = qJ(1) ^ 2;
t111 = t68 + t41;
t70 = pkin(3) ^ 2;
t73 = (pkin(1) ^ 2);
t142 = t73 - t70;
t71 = pkin(2) ^ 2;
t23 = t71 + t142;
t170 = t111 - t23;
t83 = (Ifges(2,5) * V_base(4) + Ifges(2,6) * V_base(5));
t169 = V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2;
t130 = -pkin(2) + t49;
t131 = -pkin(2) - t49;
t109 = (pkin(1) + t131) * (pkin(1) + t130) * (pkin(1) - t130) * (pkin(1) - t131);
t168 = 1 / pkin(1);
t167 = 2 * pkin(3);
t166 = 4 * pkin(3);
t165 = 1 / t49;
t164 = -2 * Ifges(2,4);
t162 = pkin(3) * mrSges(2,1);
t161 = pkin(3) * Ifges(2,4);
t160 = Ifges(2,5) * pkin(3);
t63 = t73 * m(2);
t158 = mrSges(2,2) * t73;
t46 = (V_base(4) ^ 2);
t157 = Ifges(2,4) * t46;
t156 = t23 * t46;
t40 = -3 * t70 + t71;
t52 = (Ifges(2,2) - Ifges(2,1));
t155 = t40 * t52;
t154 = t46 * t52;
t58 = pkin(2) + pkin(3);
t59 = pkin(2) - pkin(3);
t153 = t59 ^ 2 * t58 ^ 2;
t152 = t70 * mrSges(2,1);
t50 = (Ifges(4,2) - Ifges(4,1));
t151 = t70 * t50;
t150 = t71 * t52;
t149 = t71 * t73;
t148 = t73 * mrSges(2,1);
t147 = t73 * mrSges(2,3);
t44 = (V_base(6) ^ 2);
t146 = t73 * t44;
t54 = t71 * mrSges(2,1);
t39 = t54 / 0.2e1;
t53 = Ifges(2,3) * pkin(3);
t145 = t39 + t53;
t45 = (V_base(5) ^ 2);
t144 = -t45 + t46;
t37 = Ifges(2,1) / 0.2e1 + Ifges(2,2) / 0.2e1;
t143 = t71 - t73;
t141 = pkin(1) * qJD(1);
t140 = t168 / 0.2e1;
t106 = t142 - t111;
t12 = -t71 + t106;
t139 = qJD(1) * t12;
t138 = pkin(1) * V_base(1);
t137 = pkin(1) * V_base(2);
t136 = pkin(1) * V_base(3);
t135 = mrSges(2,1) * V_base(1);
t134 = pkin(3) * t158;
t129 = Ifges(4,3) * t149;
t128 = t59 * t58 * t52;
t127 = -t70 + t143;
t126 = -(t168 * t165) / 0.2e1;
t125 = t165 * t140;
t124 = V_base(5) * pkin(1);
t123 = V_base(6) * pkin(1);
t122 = t153 / 0.4e1;
t121 = pkin(2) * mrSges(4,3) / t140;
t35 = Ifges(2,2) / 0.4e1 - Ifges(2,1) / 0.4e1;
t120 = mrSges(2,2) * V_base(6);
t117 = t23 * V_base(5);
t116 = V_base(4) * mrSges(2,2);
t115 = V_base(6) * t68;
t114 = (-t70 + t71) * pkin(3);
t113 = t49 * t136;
t112 = -2 * mrSges(2,2) * t146;
t110 = t54 - t152;
t108 = -2 * t123;
t107 = t68 - t127;
t105 = pkin(2) * t123;
t104 = 2 * V_base(6);
t103 = mrSges(2,3) * t124;
t102 = Ifges(2,6) * t123;
t101 = V_base(4) * pkin(1) * mrSges(2,3);
t100 = -2 * t109;
t99 = 4 * t109;
t98 = V_base(4) * V_base(5);
t96 = t41 + t107;
t95 = Ifges(4,5) * t105;
t94 = pkin(3) * t102;
t93 = Ifges(4,6) * t105;
t92 = V_base(4) * t123;
t18 = V_base(2) + t123;
t91 = Ifges(4,4) * t98;
t90 = t109 * t141;
t88 = pkin(3) * t92;
t36 = Ifges(4,1) / 0.4e1 - Ifges(4,2) / 0.4e1;
t72 = t73 ^ 2;
t87 = t122 * t50 - t36 * t72;
t86 = V_base(1) * t96;
t85 = t170 * V_base(6);
t84 = Ifges(2,5) * V_base(5) - Ifges(2,6) * V_base(4);
t82 = Ifges(4,5) * V_base(4) + Ifges(4,6) * V_base(5);
t81 = (-t107 - 2 * t159) * t125;
t80 = -mrSges(4,1) * V_base(1) - mrSges(4,2) * V_base(2) + Ifges(4,5) * V_base(5);
t75 = sqrt(-t109);
t67 = qJ(1) * t68;
t66 = t68 ^ 2;
t65 = qJD(1) ^ 2;
t62 = pkin(3) * m(2);
t51 = Ifges(2,2) + Ifges(2,1);
t43 = -3 * t162;
t42 = t49 ^ 2;
t38 = -t152 / 0.2e1;
t34 = -mrSges(2,1) + 2 * t62;
t22 = t37 * t70;
t21 = (Ifges(4,2) / 0.2e1 + Ifges(4,1) / 0.2e1) * t71;
t20 = t63 - t162;
t19 = (t62 - mrSges(2,1) / 0.2e1) * t73;
t17 = pkin(3) * t52 - t148;
t15 = Ifges(4,4) * t117;
t14 = -t71 - t106;
t10 = t117 * t50 + 2 * t93;
t6 = V_base(6) + t165 / t75 * t139;
t5 = t125 * t75 * V_base(4) + t81 * V_base(5);
t4 = t126 * t75 * V_base(5) + t81 * V_base(4);
t3 = (-t18 * t75 - t86) * t125 + qJD(1);
t2 = qJ(1) * t5 + t124 - V_base(3);
t1 = (-(t18 * t96) + t75 * V_base(1)) * t126 - t6 * qJ(1);
t7 = m(1) * t169 / 0.2e1 + m(3) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + ((V_base(2) * mrSges(1,1)) - (V_base(1) * mrSges(1,2)) + (Ifges(1,3) * V_base(6)) / 0.2e1) * V_base(6) + (-(V_base(3) * mrSges(1,1)) + (V_base(1) * mrSges(1,3)) + (Ifges(1,6) * V_base(6)) + (Ifges(1,2) * V_base(5)) / 0.2e1) * V_base(5) + ((V_base(3) * mrSges(1,2)) - (V_base(2) * mrSges(1,3)) + (Ifges(1,4) * V_base(5)) + (Ifges(1,5) * V_base(6)) + (Ifges(1,1) * V_base(4)) / 0.2e1) * V_base(4) + (((0.4e1 * ((mrSges(2,1) * t96 * V_base(2)) - (mrSges(2,2) * t86) + ((mrSges(2,1) * t123 + t83) * t68) + ((-(Ifges(2,3) - t162) * t123 + t83 * pkin(3)) * t163) + (-t148 / 0.2e1 + t38 + t145) * t108 - (t83 * t127)) * pkin(1) * t139 + ((2 * (mrSges(2,1) * V_base(4) + mrSges(2,2) * V_base(5)) * t113 + (Ifges(2,4) * t144 + t52 * t98) * t68 + (t157 * t167 + t112 + t84 * pkin(1) * t104 + 2 * (-t158 - t161) * t45 + 2 * t17 * t98) * qJ(1) + (t127 * Ifges(2,4) - 2 * t134) * t45 + (((t52 - 2 * t162) * t73 - t128) * V_base(4) + 2 * t123 * t160) * V_base(5) - t127 * t157 - 2 * Ifges(2,6) * t88 + pkin(3) * t112 + (mrSges(2,2) * V_base(2) + t135) * t49 * t108) * t100)) * t75 + (-(((-mrSges(2,1) * V_base(5) + t116) * t68 + (pkin(3) * t116 + t20 * V_base(5)) * t163 + (t34 * t73 + t110) * V_base(5) - t127 * t116) * t113) + (t35 * t45 + (Ifges(2,4) * t98) - t154 / 0.4e1) * t66 + ((t17 * t45 + ((t158 + 4 * t161) * V_base(4) - t102) * V_base(5) - pkin(3) * t154 - Ifges(2,5) * t92 - mrSges(2,1) * t146) * t67) + (((t72 * m(2)) + (t43 + t37) * t73 - t155 / 0.2e1) * t45 + (((t164 * t40 + 3 * t134) * V_base(4) - 3 * t94) * V_base(5)) + (t37 * t73 + t155 / 0.2e1) * t46 - 0.3e1 * (0.2e1 / 0.3e1 * t147 + t160) * t92 + ((Ifges(2,3) + t43 + t63) * t146)) * t68 + (((t34 * t72) + (pkin(3) * t51 - (3 * t152) + t54) * t73 - (pkin(3) * t128)) * t45 + (((-4 * t114 * Ifges(2,4) + (-t40 * t73 + t72) * mrSges(2,2)) * V_base(4) - (3 * t70 - t143) * t102) * V_base(5)) + pkin(3) * (t51 * t73 + t128) * t46 - (((mrSges(2,3) * t166 + Ifges(2,5)) * t73 - Ifges(2,5) * t40) * t92) + 0.2e1 * (t19 - 0.3e1 / 0.2e1 * t152 + t145) * t146) * qJ(1) + (((t70 * m(2)) - t162 + t35) * t72 + (t22 - t150 / 0.2e1 + (t114 * mrSges(2,1))) * t73 + t52 * t122) * t45 + ((((pkin(3) * mrSges(2,2) + Ifges(2,4)) * t72 + (-mrSges(2,2) * t114 + t164 * t71) * t73 + Ifges(2,4) * t153) * V_base(4) + t127 * t94) * V_base(5)) + (-t35 * t72 + (t22 + t150 / 0.2e1) * t73 - (t52 * t153) / 0.4e1) * t46 - ((-Ifges(2,5) * t127 + t147 * t167) * t88) + (((-mrSges(2,1) + t62) * t73 + t53 + t110) * pkin(3) * t146) + (-((mrSges(2,1) * t115) + ((-t20 * V_base(6) + t101) * t163) + (t101 * t167) - 0.2e1 * V_base(6) * (t19 + t39 + t38)) * t137 + ((mrSges(2,2) * t115 + (pkin(3) * t120 + t103) * t163 + t103 * t167 - t127 * t120) * t138)) * t49 + t169 * t42 * t63) * t99 + (4 * (-(-mrSges(2,2) * t18 - t135 + t84) * t90 - Ifges(2,3) * t73 * t65 * t12) * t12)) / (t49 ^ 2) + ((-8 * (-2 * Ifges(4,3) * t105 + (mrSges(4,1) * V_base(2) - mrSges(4,2) * V_base(1)) * t14 + t170 * t82) * t141 * t171 + (-Ifges(4,4) * t156 - t10 * V_base(4) + t15 * V_base(5) + t111 * (Ifges(4,4) * t144 + t50 * t98) + (2 * (mrSges(4,1) * V_base(4) + mrSges(4,2) * V_base(5)) * V_base(3) + t80 * t104) * pkin(2) * pkin(1)) * t100) * t75 - (16 * t42 * t65 * t129) + (8 * (-Ifges(4,6) * V_base(4) + t80) * t90 * t171) + ((-(t82 * t105) + (-(2 * t91) + (t46 / 0.2e1 - t45 / 0.2e1) * t50) * (t73 + t40)) * t68 - 0.4e1 * (-(t50 * t156) / 0.4e1 + (t15 + t95 / 0.2e1) * V_base(4) + (t10 * V_base(5)) / 0.4e1) * t159 + ((t21 + t151 / 0.2e1) * t73 - t87) * t46 + ((Ifges(4,4) * (-2 * t70 * t73 + t153 + t72) * V_base(5) + t23 * t95) * V_base(4)) + ((t21 - t151 / 0.2e1) * t73 + t87) * t45 + (t93 * t117) + (t44 * t129) + (t67 * t166 + t66) * (t36 * t46 + t91 + (t45 * t50) / 0.4e1) + t169 * m(4) * t149 + (-(t14 * (-mrSges(4,1) * V_base(5) + mrSges(4,2) * V_base(4)) * t136) + ((mrSges(4,2) * t85) + t121 * V_base(5)) * t138 - ((mrSges(4,1) * t85) + t121 * V_base(4)) * t137) * pkin(2)) * t99) / (pkin(2) ^ 2)) / (pkin(1) ^ 2) / t109 / 0.8e1 + (t3 * mrSges(3,1) - t1 * mrSges(3,3) + Ifges(3,2) * t6 / 0.2e1) * t6 + (-t3 * mrSges(3,2) + t2 * mrSges(3,3) + Ifges(3,4) * t6 + Ifges(3,1) * t5 / 0.2e1) * t5 + (t2 * mrSges(3,1) - t1 * mrSges(3,2) - Ifges(3,5) * t5 - Ifges(3,6) * t6 + Ifges(3,3) * t4 / 0.2e1) * t4;
T = t7;
