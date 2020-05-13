% Calculate kinetic energy for
% palh1m2TE
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
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh1m2TE_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2TE_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh1m2TE_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_energykin_floatb_twist_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2TE_energykin_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2TE_energykin_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2TE_energykin_floatb_twist_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:13:14
% EndTime: 2020-05-01 20:13:26
% DurationCPUTime: 12.21s
% Computational Cost: add. (2946->343), mult. (4711->483), div. (0->0), fcn. (4888->36), ass. (0->158)
t142 = pkin(22) + pkin(21);
t121 = sin(t142);
t122 = cos(t142);
t158 = sin(qJ(2));
t163 = cos(qJ(2));
t157 = sin(qJ(3));
t162 = cos(qJ(3));
t152 = cos(pkin(20));
t193 = sin(pkin(20));
t201 = sin(pkin(18));
t202 = cos(pkin(18));
t76 = t152 * t201 - t193 * t202;
t79 = t152 * t202 + t193 * t201;
t172 = t157 * t79 - t76 * t162;
t53 = t157 * t76 + t79 * t162;
t214 = -t158 * t172 + t53 * t163;
t217 = t53 * t158 + t163 * t172;
t222 = t121 * t214 + t217 * t122;
t159 = sin(qJ(1));
t164 = cos(qJ(1));
t81 = t159 * V_base(4) - t164 * V_base(5);
t73 = qJD(2) + t81;
t71 = qJD(3) + t73;
t221 = (m(2) * pkin(13));
t218 = Ifges(2,3) * qJD(1);
t181 = t159 * V_base(5);
t84 = t164 * V_base(4) + t181;
t124 = V_base(6) + qJD(1);
t175 = V_base(4) * V_base(5);
t168 = t159 * V_base(2) + t164 * V_base(1);
t216 = t163 * t124 - t158 * t84;
t215 = -t158 * t124 - t163 * t84;
t213 = pkin(2) * t71;
t146 = qJ(3) + qJ(2);
t171 = pkin(18) - t142 - t146;
t102 = qJ(1) + t171;
t209 = -sin(t102) / 0.2e1;
t103 = -qJ(1) + t171;
t208 = cos(t103) / 0.2e1;
t18 = t121 * t217 - t122 * t214;
t100 = t158 * pkin(13) + t163 * pkin(15);
t119 = pkin(1) * t157 * qJD(2);
t80 = -t158 * t157 + t162 * t163;
t83 = t163 * t157 + t162 * t158;
t178 = t163 * pkin(13) - t158 * pkin(15);
t92 = pkin(1) + t178;
t30 = -t168 * t83 + t81 * (t162 * t100 + t92 * t157) + t119 + t80 * V_base(3);
t24 = t71 * pkin(5) + t30;
t187 = t162 * qJD(2);
t29 = t83 * V_base(3) - pkin(1) * t187 - (-t157 * t100 + t92 * t162) * t81 + t168 * t80;
t6 = t18 * t29 + t222 * t24;
t133 = qJ(1) + t146;
t207 = cos(t133) / 0.2e1;
t134 = qJ(1) - t146;
t206 = cos(t134) / 0.2e1;
t147 = qJ(3) + qJ(1);
t205 = -sin(t147) / 0.2e1;
t148 = -qJ(3) + qJ(1);
t204 = cos(t148) / 0.2e1;
t112 = t163 * pkin(1) + pkin(13);
t113 = pkin(1) * t158 - pkin(15);
t194 = cos(pkin(22));
t192 = sin(pkin(22));
t189 = (2 * mrSges(2,3) + t221) * pkin(13);
t144 = V_base(5) ^ 2;
t145 = V_base(4) ^ 2;
t188 = t144 - t145;
t105 = V_base(5) * pkin(13) + V_base(1);
t182 = V_base(5) * pkin(15);
t160 = sin(pkin(17));
t165 = cos(pkin(17));
t82 = -t160 * t202 + t165 * t201;
t85 = t160 * t201 + t165 * t202;
t179 = t158 * t85 - t82 * t163;
t123 = V_base(1) * t159;
t176 = -V_base(2) * t164 + t123;
t5 = t18 * t24 - t222 * t29;
t174 = t164 * mrSges(2,1) - t159 * mrSges(2,2);
t173 = -t159 * mrSges(2,1) - t164 * mrSges(2,2);
t21 = -(t124 * t76 + t79 * t84) * t121 + (-t124 * t79 + t76 * t84) * t122;
t106 = V_base(4) * pkin(13) - V_base(2);
t150 = sin(pkin(19));
t153 = cos(pkin(19));
t74 = t162 * t150 + t157 * t153;
t77 = -t157 * t150 + t162 * t153;
t62 = -t124 * pkin(15) + t159 * t105 + t164 * t106;
t75 = -t192 * t201 - t194 * t202;
t78 = t192 * t202 - t194 * t201;
t170 = (-t78 * t158 - t75 * t163) * pkin(1);
t169 = (t158 * t75 - t78 * t163) * pkin(1);
t46 = t84 * t112 + t124 * t113 + t176;
t44 = t157 * t215 + t162 * t216;
t27 = -t44 * pkin(5) + t46;
t161 = cos(qJ(4));
t156 = sin(qJ(4));
t155 = -Ifges(2,2) + Ifges(2,1);
t154 = pkin(13) - pkin(16);
t151 = cos(pkin(21));
t149 = sin(pkin(21));
t139 = V_base(4) * pkin(15);
t132 = mrSges(2,1) * pkin(13) - Ifges(2,5);
t131 = mrSges(2,2) * pkin(13) - Ifges(2,6);
t129 = cos(t147);
t128 = cos(t146);
t127 = sin(t148);
t125 = sin(t146);
t116 = sin(t134);
t115 = sin(t133);
t99 = pkin(13) * t201 + pkin(15) * t202;
t98 = -pkin(13) * t202 + pkin(15) * t201;
t96 = cos(t102);
t95 = sin(t103);
t91 = -t139 + t105;
t90 = t139 + t105;
t89 = -t106 - t182;
t88 = t106 - t182;
t69 = -t81 * pkin(14) + V_base(3);
t68 = t81 * pkin(15) + V_base(3);
t67 = t164 * t105 - t159 * t106;
t63 = -t81 * t154 + t168;
t61 = t62 ^ 2;
t60 = t124 * t128 - t125 * t84;
t59 = t125 * t124 + t128 * t84;
t56 = t158 * t82 + t85 * t163;
t51 = pkin(14) * t124 + t154 * t84 + t176;
t49 = (-t149 * t75 + t151 * t78) * pkin(4) + t112;
t48 = -t158 * t68 - t163 * t67;
t47 = -t158 * t67 + t163 * t68;
t45 = t46 ^ 2;
t43 = t157 * t216 - t162 * t215;
t41 = t124 * t78 + t75 * t84;
t39 = t124 * t75 - t78 * t84;
t36 = t215 * t74 + t216 * t77;
t35 = -t215 * t77 + t216 * t74;
t34 = t56 * t124 - t179 * t84;
t33 = -t124 * t179 - t84 * t56;
t32 = t81 * t100 - t168 * t158 + t163 * V_base(3) + t213 * t77;
t31 = -t158 * V_base(3) - t168 * t163 + t81 * t178 + t213 * t74;
t28 = -t36 * pkin(2) + t62;
t26 = -t179 * t69 - t56 * t63;
t25 = -t179 * t63 + t56 * t69;
t23 = t77 * t47 + t74 * t48;
t22 = t74 * t47 - t77 * t48;
t20 = (t121 * t124 - t122 * t84) * t79 + (-t121 * t84 - t122 * t124) * t76;
t19 = qJD(4) - t21;
t16 = -qJD(2) * t169 + t75 * V_base(3) - t81 * (t192 * t98 + t194 * t99 + t169) - t168 * t78;
t15 = qJD(2) * t170 + t78 * V_base(3) + t81 * (t192 * t99 - t194 * t98 + t170) + t168 * t75;
t14 = t156 * t81 + t161 * t20;
t13 = -t156 * t20 + t161 * t81;
t12 = t90 * t206 - t88 * t116 / 0.2e1 + t91 * t207 + t89 * t115 / 0.2e1 + V_base(3) * t125 + ((-t95 / 0.2e1 + t209) * V_base(5) + (t208 - t96 / 0.2e1) * V_base(4)) * pkin(4) + (-t187 + (t204 + t129 / 0.2e1) * V_base(5) + (-t127 / 0.2e1 + t205) * V_base(4)) * pkin(1);
t11 = t88 * t206 + t90 * t116 / 0.2e1 + t89 * t207 - t91 * t115 / 0.2e1 + t119 + V_base(3) * t128 + ((t208 + t96 / 0.2e1) * V_base(5) + (t95 / 0.2e1 + t209) * V_base(4)) * pkin(4) + ((t127 / 0.2e1 + t205) * V_base(5) + (t204 - t129 / 0.2e1) * V_base(4)) * pkin(1);
t10 = t49 * t181 + t123 + t124 * ((-t149 * t78 - t151 * t75) * pkin(4) + t113) + (t49 * V_base(4) - V_base(2)) * t164;
t7 = -t21 * pkin(9) - t20 * pkin(11) + t27;
t4 = t81 * pkin(11) + t6;
t3 = -t81 * pkin(9) - t5;
t2 = t156 * t7 + t161 * t4;
t1 = -t156 * t4 + t161 * t7;
t8 = (-mrSges(4,1) * t46 + t29 * mrSges(4,3) + Ifges(4,2) * t44 / 0.2e1) * t44 + (t62 * mrSges(9,2) - t23 * mrSges(9,3) + Ifges(9,4) * t36 + Ifges(9,1) * t35 / 0.2e1) * t35 + (t5 * mrSges(5,1) + t16 * mrSges(8,1) - t6 * mrSges(5,2) - t15 * mrSges(8,2) + Ifges(5,5) * t20 + Ifges(8,5) * t41 + Ifges(5,6) * t21 + Ifges(8,6) * t39 + (Ifges(5,3) / 0.2e1 + Ifges(8,3) / 0.2e1) * t81) * t81 + (mrSges(8,2) * t46 - t16 * mrSges(8,3) + Ifges(8,1) * t41 / 0.2e1) * t41 + (t11 * mrSges(11,1) + t30 * mrSges(4,1) + t23 * mrSges(9,1) - t12 * mrSges(11,2) - t29 * mrSges(4,2) - t22 * mrSges(9,2) + Ifges(11,5) * t59 + Ifges(4,5) * t43 + Ifges(9,5) * t35 + Ifges(11,6) * t60 + Ifges(4,6) * t44 + Ifges(9,6) * t36 + (Ifges(11,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(9,3) / 0.2e1) * t71) * t71 + (t51 * mrSges(7,2) - t26 * mrSges(7,3) + Ifges(7,1) * t34 / 0.2e1) * t34 + m(5) * (t27 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(10) * (t28 ^ 2 + t31 ^ 2 + t32 ^ 2) / 0.2e1 + (-t27 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,2) * t21 / 0.2e1) * t21 + (-V_base(4) * V_base(2) + V_base(5) * V_base(1)) * (mrSges(1,3) + mrSges(2,3) + t221) + (-t51 * mrSges(7,1) + t25 * mrSges(7,3) + Ifges(7,4) * t34 + Ifges(7,2) * t33 / 0.2e1) * t33 + m(9) * (t22 ^ 2 + t23 ^ 2 + t61) / 0.2e1 + m(3) * (t47 ^ 2 + t48 ^ 2 + t61) / 0.2e1 + (mrSges(4,2) * t46 - t30 * mrSges(4,3) + Ifges(4,4) * t44 + Ifges(4,1) * t43 / 0.2e1) * t43 + (t27 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,4) * t21 + Ifges(5,1) * t20 / 0.2e1) * t20 + (-t10 * mrSges(11,1) + t12 * mrSges(11,3) + Ifges(11,2) * t60 / 0.2e1) * t60 + ((Ifges(3,5) + Ifges(10,5)) * t73 + t28 * mrSges(10,2) - t31 * mrSges(10,3) + t62 * mrSges(3,2) - t48 * mrSges(3,3) + (Ifges(10,1) + Ifges(3,1)) * t216 / 0.2e1) * t216 + ((Ifges(3,6) + Ifges(10,6)) * t73 - t28 * mrSges(10,1) + t32 * mrSges(10,3) - t62 * mrSges(3,1) + t47 * mrSges(3,3) + (Ifges(10,2) + Ifges(3,2)) * t215 / 0.2e1 + (Ifges(10,4) + Ifges(3,4)) * t216) * t215 + (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) * (m(1) / 0.2e1 + m(2) / 0.2e1) + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t14 + Ifges(6,6) * t19 + Ifges(6,2) * t13 / 0.2e1) * t13 + (-mrSges(9,1) * t62 + t22 * mrSges(9,3) + Ifges(9,2) * t36 / 0.2e1) * t36 + (t48 * mrSges(3,1) + t26 * mrSges(7,1) + t31 * mrSges(10,1) - t47 * mrSges(3,2) - t25 * mrSges(7,2) - t32 * mrSges(10,2) + Ifges(7,5) * t34 + Ifges(7,6) * t33 + (Ifges(3,3) / 0.2e1 + Ifges(7,3) / 0.2e1 + Ifges(10,3) / 0.2e1) * t73) * t73 + (0.4e1 * Ifges(2,4) * t175 - t188 * t155) * t164 ^ 2 / 0.2e1 + (Ifges(2,1) + Ifges(1,2) + t189) * t144 / 0.2e1 + (Ifges(1,1) + Ifges(2,2) + t189) * t145 / 0.2e1 + (Ifges(1,4) - Ifges(2,4)) * t175 + (Ifges(2,4) * t188 + t155 * t175) * t159 * t164 + (-mrSges(8,1) * t46 + t15 * mrSges(8,3) + Ifges(8,4) * t41 + Ifges(8,2) * t39 / 0.2e1) * t39 + m(8) * (t15 ^ 2 + t16 ^ 2 + t45) / 0.2e1 + m(4) * (t29 ^ 2 + t30 ^ 2 + t45) / 0.2e1 + m(7) * (t25 ^ 2 + t26 ^ 2 + t51 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(11) * (t10 ^ 2 + t11 ^ 2 + t12 ^ 2) / 0.2e1 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t19 / 0.2e1) * t19 + ((-t131 * t164 - t132 * t159) * V_base(5) + t173 * V_base(1) + t174 * V_base(2) + (t131 * t159 - t132 * t164) * V_base(4) + t218 / 0.2e1) * qJD(1) + ((t131 * V_base(4) - V_base(5) * t132) * t159 + V_base(5) * Ifges(1,6) + t218 + (-mrSges(1,2) + t173) * V_base(1) + (mrSges(1,1) + t174) * V_base(2) + (-V_base(5) * t131 - t132 * V_base(4)) * t164 + Ifges(1,5) * V_base(4) + (Ifges(1,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + (-V_base(5) * mrSges(1,1) + V_base(4) * mrSges(1,2) + (-V_base(5) * mrSges(2,1) + V_base(4) * mrSges(2,2)) * t164 + (V_base(4) * mrSges(2,1) + V_base(5) * mrSges(2,2)) * t159) * V_base(3) + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t19 + Ifges(6,1) * t14 / 0.2e1) * t14 + (t10 * mrSges(11,2) - t11 * mrSges(11,3) + Ifges(11,4) * t60 + Ifges(11,1) * t59 / 0.2e1) * t59;
T = t8;
