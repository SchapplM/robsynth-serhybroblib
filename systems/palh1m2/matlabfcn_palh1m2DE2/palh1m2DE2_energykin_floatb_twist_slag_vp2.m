% Calculate kinetic energy for
% palh1m2DE2
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
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m2DE2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2DE2_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh1m2DE2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_energykin_floatb_twist_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE2_energykin_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2DE2_energykin_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2DE2_energykin_floatb_twist_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 20:57:45
% EndTime: 2020-05-02 20:57:56
% DurationCPUTime: 5.38s
% Computational Cost: add. (2943->320), mult. (4657->457), div. (0->0), fcn. (4904->36), ass. (0->148)
t138 = pkin(22) + pkin(21);
t119 = sin(t138);
t120 = cos(t138);
t151 = sin(qJ(2));
t157 = cos(qJ(2));
t150 = sin(qJ(3));
t156 = cos(qJ(3));
t146 = cos(pkin(20));
t153 = sin(pkin(18));
t159 = cos(pkin(18));
t180 = sin(pkin(20));
t76 = t146 * t153 - t159 * t180;
t80 = t146 * t159 + t153 * t180;
t166 = t150 * t80 - t156 * t76;
t53 = t150 * t76 + t156 * t80;
t195 = -t151 * t166 + t157 * t53;
t198 = t151 * t53 + t157 * t166;
t201 = t119 * t195 + t120 * t198;
t152 = sin(qJ(1));
t158 = cos(qJ(1));
t81 = t152 * V_base(4) - t158 * V_base(5);
t74 = qJD(2) + t81;
t72 = qJD(3) + t74;
t171 = t152 * V_base(5);
t85 = t158 * V_base(4) + t171;
t122 = V_base(6) + qJD(1);
t162 = t152 * V_base(2) + t158 * V_base(1);
t197 = t122 * t157 - t151 * t85;
t196 = -t122 * t151 - t157 * t85;
t194 = pkin(2) * t72;
t141 = qJ(2) + qJ(3);
t165 = pkin(18) - t138 - t141;
t103 = qJ(1) + t165;
t193 = -sin(t103) / 0.2e1;
t104 = -qJ(1) + t165;
t192 = cos(t104) / 0.2e1;
t18 = t119 * t198 - t120 * t195;
t101 = pkin(13) * t151 + pkin(15) * t157;
t177 = t150 * qJD(2);
t83 = t150 * t157 + t151 * t156;
t84 = -t150 * t151 + t156 * t157;
t168 = pkin(13) * t157 - pkin(15) * t151;
t93 = pkin(1) + t168;
t29 = -t162 * t83 + (t101 * t156 + t150 * t93) * t81 + pkin(1) * t177 + t84 * V_base(3);
t24 = pkin(5) * t72 + t29;
t176 = t156 * qJD(2);
t30 = t83 * V_base(3) - pkin(1) * t176 - (-t101 * t150 + t156 * t93) * t81 + t162 * t84;
t6 = t18 * t30 + t201 * t24;
t139 = qJ(3) + qJ(1);
t129 = qJ(2) + t139;
t191 = cos(t129) / 0.2e1;
t140 = -qJ(3) + qJ(1);
t130 = -qJ(2) + t140;
t190 = cos(t130) / 0.2e1;
t189 = -sin(t139) / 0.2e1;
t188 = cos(t140) / 0.2e1;
t111 = pkin(1) * t151 - pkin(15);
t112 = pkin(1) * t157 + pkin(13);
t181 = cos(pkin(22));
t106 = pkin(13) * V_base(5) + V_base(1);
t172 = V_base(5) * pkin(15);
t154 = sin(pkin(17));
t160 = cos(pkin(17));
t82 = t153 * t160 - t154 * t159;
t86 = t153 * t154 + t159 * t160;
t169 = t151 * t86 - t157 * t82;
t107 = pkin(13) * V_base(4) - V_base(2);
t67 = -t106 * t152 - t107 * t158;
t121 = V_base(1) * t152;
t167 = -t158 * V_base(2) + t121;
t5 = t18 * t24 - t201 * t30;
t21 = -t119 * (t122 * t76 + t80 * t85) + t120 * (-t122 * t80 + t76 * t85);
t144 = sin(pkin(19));
t147 = cos(pkin(19));
t75 = t144 * t156 + t147 * t150;
t77 = -t144 * t150 + t147 * t156;
t62 = -pkin(15) * t122 - t67;
t142 = sin(pkin(22));
t78 = t142 * t153 + t159 * t181;
t79 = t142 * t159 - t153 * t181;
t164 = (-t151 * t79 + t157 * t78) * pkin(1);
t163 = (-t151 * t78 - t157 * t79) * pkin(1);
t46 = t111 * t122 + t112 * t85 + t167;
t44 = t150 * t196 + t156 * t197;
t27 = -pkin(5) * t44 + t46;
t161 = V_base(3) ^ 2;
t155 = cos(qJ(4));
t149 = sin(qJ(4));
t148 = pkin(13) - pkin(16);
t145 = cos(pkin(21));
t143 = sin(pkin(21));
t135 = V_base(4) * pkin(15);
t128 = cos(t141);
t126 = cos(t139);
t125 = sin(t141);
t124 = sin(t140);
t115 = sin(t130);
t114 = sin(t129);
t100 = pkin(13) * t153 + pkin(15) * t159;
t99 = -pkin(13) * t159 + pkin(15) * t153;
t97 = cos(t103);
t96 = sin(t104);
t92 = -t135 + t106;
t91 = t135 + t106;
t90 = -t107 - t172;
t89 = t107 - t172;
t70 = -pkin(14) * t81 + V_base(3);
t69 = pkin(15) * t81 + V_base(3);
t68 = t106 * t158 - t107 * t152;
t63 = -t148 * t81 + t162;
t61 = t62 ^ 2;
t60 = t122 * t128 - t125 * t85;
t59 = t122 * t125 + t128 * t85;
t56 = t151 * t82 + t157 * t86;
t51 = pkin(14) * t122 + t148 * t85 + t167;
t50 = (t143 * t78 + t145 * t79) * pkin(4) + t112;
t48 = -t151 * t69 - t157 * t68;
t47 = -t151 * t68 + t157 * t69;
t45 = t46 ^ 2;
t43 = t150 * t197 - t156 * t196;
t41 = t122 * t79 - t78 * t85;
t40 = -t122 * t78 - t79 * t85;
t36 = t196 * t75 + t197 * t77;
t35 = -t196 * t77 + t197 * t75;
t34 = -t122 * t169 - t56 * t85;
t33 = t122 * t56 - t169 * t85;
t32 = t81 * t101 - t162 * t151 + t157 * V_base(3) + t194 * t77;
t31 = -t151 * V_base(3) - t162 * t157 + t81 * t168 + t194 * t75;
t28 = -pkin(2) * t36 + t62;
t26 = -t169 * t70 - t56 * t63;
t25 = -t169 * t63 + t56 * t70;
t23 = t47 * t77 + t48 * t75;
t22 = t47 * t75 - t48 * t77;
t20 = (t119 * t122 - t120 * t85) * t80 + (-t119 * t85 - t120 * t122) * t76;
t19 = qJD(4) - t21;
t16 = -qJD(2) * t163 - t78 * V_base(3) - (t100 * t181 + t142 * t99 + t163) * t81 - t162 * t79;
t15 = qJD(2) * t164 + t79 * V_base(3) + t81 * (t100 * t142 - t181 * t99 + t164) - t162 * t78;
t14 = t149 * t81 + t155 * t20;
t13 = -t149 * t20 + t155 * t81;
t12 = t89 * t190 + t91 * t115 / 0.2e1 + t90 * t191 - t92 * t114 / 0.2e1 + V_base(3) * t128 + (t177 + (t124 / 0.2e1 + t189) * V_base(5) + (t188 - t126 / 0.2e1) * V_base(4)) * pkin(1) + ((t192 + t97 / 0.2e1) * V_base(5) + (t96 / 0.2e1 + t193) * V_base(4)) * pkin(4);
t11 = t91 * t190 - t89 * t115 / 0.2e1 + t92 * t191 + t90 * t114 / 0.2e1 + V_base(3) * t125 + ((-t96 / 0.2e1 + t193) * V_base(5) + (t192 - t97 / 0.2e1) * V_base(4)) * pkin(4) + (-t176 + (t188 + t126 / 0.2e1) * V_base(5) + (-t124 / 0.2e1 + t189) * V_base(4)) * pkin(1);
t10 = t50 * t171 + t121 + t122 * ((-t143 * t79 + t145 * t78) * pkin(4) + t111) + (t50 * V_base(4) - V_base(2)) * t158;
t7 = -pkin(9) * t21 - pkin(11) * t20 + t27;
t4 = pkin(11) * t81 + t6;
t3 = -pkin(9) * t81 - t5;
t2 = t149 * t7 + t155 * t4;
t1 = -t149 * t4 + t155 * t7;
t8 = (t62 * mrSges(9,2) - t23 * mrSges(9,3) + Ifges(9,4) * t36 + Ifges(9,1) * t35 / 0.2e1) * t35 + (-t27 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,2) * t21 / 0.2e1) * t21 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t19 + Ifges(6,1) * t14 / 0.2e1) * t14 + m(2) * (t67 ^ 2 + t68 ^ 2 + t161) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t161) / 0.2e1 + (V_base(3) * mrSges(2,2) - t67 * mrSges(2,3) + Ifges(2,1) * t85 / 0.2e1) * t85 + m(3) * (t47 ^ 2 + t48 ^ 2 + t61) / 0.2e1 + m(9) * (t22 ^ 2 + t23 ^ 2 + t61) / 0.2e1 + m(7) * (t25 ^ 2 + t26 ^ 2 + t51 ^ 2) / 0.2e1 + m(4) * (t29 ^ 2 + t30 ^ 2 + t45) / 0.2e1 + m(8) * (t15 ^ 2 + t16 ^ 2 + t45) / 0.2e1 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t19 / 0.2e1) * t19 + (t10 * mrSges(11,2) - t12 * mrSges(11,3) + Ifges(11,4) * t60 + Ifges(11,1) * t59 / 0.2e1) * t59 + (t51 * mrSges(7,2) - t26 * mrSges(7,3) + Ifges(7,4) * t34 + Ifges(7,1) * t33 / 0.2e1) * t33 + (V_base(3) * mrSges(2,1) + t5 * mrSges(5,1) + t16 * mrSges(8,1) - t6 * mrSges(5,2) - t15 * mrSges(8,2) - t68 * mrSges(2,3) - Ifges(2,4) * t85 + Ifges(5,5) * t20 + Ifges(8,5) * t41 - Ifges(2,6) * t122 + Ifges(5,6) * t21 + Ifges(8,6) * t40 + (Ifges(2,2) / 0.2e1 + Ifges(8,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t81) * t81 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t14 + Ifges(6,6) * t19 + Ifges(6,2) * t13 / 0.2e1) * t13 + (-t46 * mrSges(4,1) + t30 * mrSges(4,3) + Ifges(4,2) * t44 / 0.2e1) * t44 + (-t46 * mrSges(8,1) + t15 * mrSges(8,3) + Ifges(8,4) * t41 + Ifges(8,2) * t40 / 0.2e1) * t40 + (-t10 * mrSges(11,1) + t11 * mrSges(11,3) + Ifges(11,2) * t60 / 0.2e1) * t60 + (t46 * mrSges(4,2) - t29 * mrSges(4,3) + Ifges(4,4) * t44 + Ifges(4,1) * t43 / 0.2e1) * t43 + m(5) * (t27 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(10) * (t28 ^ 2 + t31 ^ 2 + t32 ^ 2) / 0.2e1 + m(11) * (t10 ^ 2 + t11 ^ 2 + t12 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (-t51 * mrSges(7,1) + t25 * mrSges(7,3) + Ifges(7,2) * t34 / 0.2e1) * t34 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t12 * mrSges(11,1) + t29 * mrSges(4,1) + t23 * mrSges(9,1) - t11 * mrSges(11,2) - t30 * mrSges(4,2) - t22 * mrSges(9,2) + Ifges(11,5) * t59 + Ifges(4,5) * t43 + Ifges(9,5) * t35 + Ifges(11,6) * t60 + Ifges(4,6) * t44 + Ifges(9,6) * t36 + (Ifges(11,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(9,3) / 0.2e1) * t72) * t72 + (t67 * mrSges(2,1) - t68 * mrSges(2,2) + Ifges(2,5) * t85 + Ifges(2,3) * t122 / 0.2e1) * t122 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t27 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,4) * t21 + Ifges(5,1) * t20 / 0.2e1) * t20 + (t28 * mrSges(10,2) - t31 * mrSges(10,3) + t62 * mrSges(3,2) - t48 * mrSges(3,3) + (Ifges(3,5) + Ifges(10,5)) * t74 + (Ifges(10,1) + Ifges(3,1)) * t197 / 0.2e1) * t197 + (-t28 * mrSges(10,1) + t32 * mrSges(10,3) - t62 * mrSges(3,1) + t47 * mrSges(3,3) + (Ifges(3,6) + Ifges(10,6)) * t74 + (Ifges(10,2) + Ifges(3,2)) * t196 / 0.2e1 + (Ifges(10,4) + Ifges(3,4)) * t197) * t196 + (t46 * mrSges(8,2) - t16 * mrSges(8,3) + Ifges(8,1) * t41 / 0.2e1) * t41 + (t48 * mrSges(3,1) + t26 * mrSges(7,1) + t31 * mrSges(10,1) - t47 * mrSges(3,2) - t25 * mrSges(7,2) - t32 * mrSges(10,2) + Ifges(7,5) * t33 + Ifges(7,6) * t34 + (Ifges(3,3) / 0.2e1 + Ifges(10,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t74) * t74 + (-t62 * mrSges(9,1) + t22 * mrSges(9,3) + Ifges(9,2) * t36 / 0.2e1) * t36;
T = t8;
