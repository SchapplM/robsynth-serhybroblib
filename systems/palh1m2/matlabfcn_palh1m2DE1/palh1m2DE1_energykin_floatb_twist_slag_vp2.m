% Calculate kinetic energy for
% palh1m2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
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
% Datum: 2020-05-01 21:04
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m2DE1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE1_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE1_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh1m2DE1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE1_energykin_floatb_twist_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE1_energykin_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2DE1_energykin_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2DE1_energykin_floatb_twist_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:55:41
% EndTime: 2020-05-01 20:55:49
% DurationCPUTime: 5.17s
% Computational Cost: add. (2943->320), mult. (4661->456), div. (0->0), fcn. (4904->36), ass. (0->148)
t140 = pkin(22) + pkin(21);
t121 = sin(t140);
t122 = cos(t140);
t152 = sin(qJ(2));
t157 = cos(qJ(2));
t151 = sin(qJ(3));
t156 = cos(qJ(3));
t147 = cos(pkin(20));
t179 = sin(pkin(20));
t187 = sin(pkin(18));
t188 = cos(pkin(18));
t77 = t147 * t187 - t179 * t188;
t80 = t147 * t188 + t179 * t187;
t165 = t151 * t80 - t156 * t77;
t53 = t151 * t77 + t156 * t80;
t196 = t152 * t165 - t157 * t53;
t199 = t53 * t152 + t157 * t165;
t202 = -t196 * t121 + t122 * t199;
t153 = sin(qJ(1));
t158 = cos(qJ(1));
t82 = t153 * V_base(4) - t158 * V_base(5);
t74 = qJD(2) + t82;
t72 = qJD(3) + t74;
t170 = t153 * V_base(5);
t85 = t158 * V_base(4) + t170;
t124 = V_base(6) + qJD(1);
t161 = t153 * V_base(2) + t158 * V_base(1);
t198 = t124 * t157 - t152 * t85;
t197 = -t124 * t152 - t157 * t85;
t195 = pkin(2) * t72;
t141 = qJ(3) + qJ(2);
t164 = pkin(18) - t140 - t141;
t103 = qJ(1) + t164;
t194 = -sin(t103) / 0.2e1;
t104 = -qJ(1) + t164;
t193 = cos(t104) / 0.2e1;
t18 = t121 * t199 + t196 * t122;
t101 = t152 * pkin(13) + t157 * pkin(15);
t120 = pkin(1) * t151 * qJD(2);
t81 = -t151 * t152 + t156 * t157;
t84 = t151 * t157 + t152 * t156;
t167 = t157 * pkin(13) - pkin(15) * t152;
t93 = pkin(1) + t167;
t30 = -t161 * t84 + (t101 * t156 + t151 * t93) * t82 + t120 + t81 * V_base(3);
t24 = pkin(5) * t72 + t30;
t175 = t156 * qJD(2);
t29 = t84 * V_base(3) - pkin(1) * t175 - t82 * (-t101 * t151 + t156 * t93) + t161 * t81;
t6 = t18 * t29 + t202 * t24;
t131 = qJ(1) + t141;
t192 = cos(t131) / 0.2e1;
t132 = qJ(1) - t141;
t191 = cos(t132) / 0.2e1;
t142 = qJ(3) + qJ(1);
t190 = -sin(t142) / 0.2e1;
t143 = -qJ(3) + qJ(1);
t189 = cos(t143) / 0.2e1;
t113 = t157 * pkin(1) + pkin(13);
t114 = pkin(1) * t152 - pkin(15);
t180 = cos(pkin(22));
t178 = sin(pkin(22));
t106 = V_base(5) * pkin(13) + V_base(1);
t171 = V_base(5) * pkin(15);
t154 = sin(pkin(17));
t159 = cos(pkin(17));
t83 = -t154 * t188 + t159 * t187;
t86 = t154 * t187 + t159 * t188;
t168 = t152 * t86 - t157 * t83;
t107 = V_base(4) * pkin(13) - V_base(2);
t67 = -t153 * t106 - t107 * t158;
t123 = V_base(1) * t153;
t166 = -t158 * V_base(2) + t123;
t5 = t18 * t24 - t202 * t29;
t21 = -t121 * (t124 * t77 + t80 * t85) + t122 * (-t124 * t80 + t77 * t85);
t145 = sin(pkin(19));
t148 = cos(pkin(19));
t75 = t145 * t156 + t148 * t151;
t78 = -t151 * t145 + t148 * t156;
t62 = -pkin(15) * t124 - t67;
t76 = -t178 * t187 - t180 * t188;
t79 = t178 * t188 - t180 * t187;
t163 = (-t152 * t79 - t157 * t76) * pkin(1);
t162 = (t152 * t76 - t157 * t79) * pkin(1);
t46 = t113 * t85 + t114 * t124 + t166;
t44 = t151 * t197 + t156 * t198;
t27 = -pkin(5) * t44 + t46;
t160 = V_base(3) ^ 2;
t155 = cos(qJ(4));
t150 = sin(qJ(4));
t149 = pkin(13) - pkin(16);
t146 = cos(pkin(21));
t144 = sin(pkin(21));
t137 = V_base(4) * pkin(15);
t129 = cos(t142);
t128 = cos(t141);
t127 = sin(t143);
t125 = sin(t141);
t117 = sin(t132);
t116 = sin(t131);
t100 = pkin(13) * t187 + pkin(15) * t188;
t99 = -pkin(13) * t188 + pkin(15) * t187;
t97 = cos(t103);
t96 = sin(t104);
t92 = -t137 + t106;
t91 = t137 + t106;
t90 = -t107 - t171;
t89 = t107 - t171;
t70 = -pkin(14) * t82 + V_base(3);
t69 = pkin(15) * t82 + V_base(3);
t68 = t106 * t158 - t107 * t153;
t63 = -t149 * t82 + t161;
t61 = t62 ^ 2;
t60 = t124 * t128 - t125 * t85;
t59 = t124 * t125 + t128 * t85;
t56 = t152 * t83 + t157 * t86;
t51 = pkin(14) * t124 + t149 * t85 + t166;
t49 = (-t144 * t76 + t146 * t79) * pkin(4) + t113;
t48 = -t152 * t69 - t157 * t68;
t47 = -t152 * t68 + t157 * t69;
t45 = t46 ^ 2;
t43 = t151 * t198 - t156 * t197;
t41 = t124 * t79 + t76 * t85;
t39 = t124 * t76 - t79 * t85;
t36 = t197 * t75 + t198 * t78;
t35 = -t197 * t78 + t198 * t75;
t34 = -t124 * t168 - t56 * t85;
t33 = t124 * t56 - t168 * t85;
t32 = t82 * t101 - t161 * t152 + t157 * V_base(3) + t78 * t195;
t31 = -t152 * V_base(3) - t161 * t157 + t82 * t167 + t75 * t195;
t28 = -pkin(2) * t36 + t62;
t26 = -t168 * t70 - t56 * t63;
t25 = -t168 * t63 + t56 * t70;
t23 = t47 * t78 + t48 * t75;
t22 = t47 * t75 - t48 * t78;
t20 = (t121 * t124 - t122 * t85) * t80 + (-t121 * t85 - t122 * t124) * t77;
t19 = qJD(4) - t21;
t16 = -qJD(2) * t162 + t76 * V_base(3) - (t100 * t180 + t178 * t99 + t162) * t82 - t161 * t79;
t15 = qJD(2) * t163 + t79 * V_base(3) + (t100 * t178 - t180 * t99 + t163) * t82 + t161 * t76;
t14 = t150 * t82 + t155 * t20;
t13 = -t150 * t20 + t155 * t82;
t12 = t91 * t191 - t89 * t117 / 0.2e1 + t92 * t192 + t90 * t116 / 0.2e1 + V_base(3) * t125 + ((-t96 / 0.2e1 + t194) * V_base(5) + (t193 - t97 / 0.2e1) * V_base(4)) * pkin(4) + (-t175 + (t189 + t129 / 0.2e1) * V_base(5) + (-t127 / 0.2e1 + t190) * V_base(4)) * pkin(1);
t11 = t89 * t191 + t91 * t117 / 0.2e1 + t90 * t192 - t92 * t116 / 0.2e1 + t120 + V_base(3) * t128 + ((t193 + t97 / 0.2e1) * V_base(5) + (t96 / 0.2e1 + t194) * V_base(4)) * pkin(4) + ((t127 / 0.2e1 + t190) * V_base(5) + (t189 - t129 / 0.2e1) * V_base(4)) * pkin(1);
t10 = t49 * t170 + t123 + t124 * ((-t144 * t79 - t146 * t76) * pkin(4) + t114) + (t49 * V_base(4) - V_base(2)) * t158;
t7 = -pkin(9) * t21 - pkin(11) * t20 + t27;
t4 = pkin(11) * t82 + t6;
t3 = -pkin(9) * t82 - t5;
t2 = t150 * t7 + t155 * t4;
t1 = -t150 * t4 + t155 * t7;
t8 = (-t62 * mrSges(9,1) + t22 * mrSges(9,3) + Ifges(9,2) * t36 / 0.2e1) * t36 + (t46 * mrSges(4,2) - t30 * mrSges(4,3) + Ifges(4,4) * t44 + Ifges(4,1) * t43 / 0.2e1) * t43 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t14 + Ifges(6,6) * t19 + Ifges(6,2) * t13 / 0.2e1) * t13 + (t46 * mrSges(8,2) - t16 * mrSges(8,3) + Ifges(8,1) * t41 / 0.2e1) * t41 + (t67 * mrSges(2,1) - t68 * mrSges(2,2) + Ifges(2,5) * t85 + Ifges(2,3) * t124 / 0.2e1) * t124 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t19 / 0.2e1) * t19 + (t28 * mrSges(10,2) - t31 * mrSges(10,3) + t62 * mrSges(3,2) - t48 * mrSges(3,3) + (Ifges(3,5) + Ifges(10,5)) * t74 + (Ifges(10,1) + Ifges(3,1)) * t198 / 0.2e1) * t198 + (-t28 * mrSges(10,1) + t32 * mrSges(10,3) - t62 * mrSges(3,1) + t47 * mrSges(3,3) + (Ifges(3,6) + Ifges(10,6)) * t74 + (Ifges(10,2) + Ifges(3,2)) * t197 / 0.2e1 + (Ifges(10,4) + Ifges(3,4)) * t198) * t197 + m(2) * (t67 ^ 2 + t68 ^ 2 + t160) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t160) / 0.2e1 + m(3) * (t47 ^ 2 + t48 ^ 2 + t61) / 0.2e1 + m(9) * (t22 ^ 2 + t23 ^ 2 + t61) / 0.2e1 + m(7) * (t25 ^ 2 + t26 ^ 2 + t51 ^ 2) / 0.2e1 + m(4) * (t29 ^ 2 + t30 ^ 2 + t45) / 0.2e1 + m(8) * (t15 ^ 2 + t16 ^ 2 + t45) / 0.2e1 + m(5) * (t27 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(10) * (t28 ^ 2 + t31 ^ 2 + t32 ^ 2) / 0.2e1 + (t10 * mrSges(11,2) - t11 * mrSges(11,3) + Ifges(11,4) * t60 + Ifges(11,1) * t59 / 0.2e1) * t59 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t62 * mrSges(9,2) - t23 * mrSges(9,3) + Ifges(9,4) * t36 + Ifges(9,1) * t35 / 0.2e1) * t35 + (-t46 * mrSges(4,1) + t29 * mrSges(4,3) + Ifges(4,2) * t44 / 0.2e1) * t44 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t46 * mrSges(8,1) + t15 * mrSges(8,3) + Ifges(8,4) * t41 + Ifges(8,2) * t39 / 0.2e1) * t39 + (t51 * mrSges(7,2) - t26 * mrSges(7,3) + Ifges(7,4) * t34 + Ifges(7,1) * t33 / 0.2e1) * t33 + (t27 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,4) * t21 + Ifges(5,1) * t20 / 0.2e1) * t20 + m(11) * (t10 ^ 2 + t11 ^ 2 + t12 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(3) * mrSges(2,1) + t5 * mrSges(5,1) + t16 * mrSges(8,1) - t6 * mrSges(5,2) - t15 * mrSges(8,2) - t68 * mrSges(2,3) - Ifges(2,4) * t85 + Ifges(5,5) * t20 + Ifges(8,5) * t41 - Ifges(2,6) * t124 + Ifges(5,6) * t21 + Ifges(8,6) * t39 + (Ifges(2,2) / 0.2e1 + Ifges(8,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t82) * t82 + (t11 * mrSges(11,1) + t30 * mrSges(4,1) + t23 * mrSges(9,1) - t12 * mrSges(11,2) - t29 * mrSges(4,2) - t22 * mrSges(9,2) + Ifges(11,5) * t59 + Ifges(4,5) * t43 + Ifges(9,5) * t35 + Ifges(11,6) * t60 + Ifges(4,6) * t44 + Ifges(9,6) * t36 + (Ifges(11,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(9,3) / 0.2e1) * t72) * t72 + (t48 * mrSges(3,1) + t26 * mrSges(7,1) + t31 * mrSges(10,1) - t47 * mrSges(3,2) - t25 * mrSges(7,2) - t32 * mrSges(10,2) + Ifges(7,5) * t33 + Ifges(7,6) * t34 + (Ifges(3,3) / 0.2e1 + Ifges(10,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t74) * t74 + (-t51 * mrSges(7,1) + t25 * mrSges(7,3) + Ifges(7,2) * t34 / 0.2e1) * t34 + (-t27 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,2) * t21 / 0.2e1) * t21 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (V_base(3) * mrSges(2,2) - t67 * mrSges(2,3) + Ifges(2,1) * t85 / 0.2e1) * t85 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t19 + Ifges(6,1) * t14 / 0.2e1) * t14 + (-t10 * mrSges(11,1) + t12 * mrSges(11,3) + Ifges(11,2) * t60 / 0.2e1) * t60;
T = t8;
