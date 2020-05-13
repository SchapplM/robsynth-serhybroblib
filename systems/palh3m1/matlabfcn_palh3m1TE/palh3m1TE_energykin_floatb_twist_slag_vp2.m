% Calculate kinetic energy for
% palh3m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
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
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m1TE_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(19,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m1TE_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh3m1TE_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_energykin_floatb_twist_slag_vp2: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1TE_energykin_floatb_twist_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1TE_energykin_floatb_twist_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m1TE_energykin_floatb_twist_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-17 15:18:49
% EndTime: 2020-04-17 15:19:37
% DurationCPUTime: 37.97s
% Computational Cost: add. (989573->313), mult. (1534075->546), div. (63208->22), fcn. (957424->24), ass. (0->208)
t248 = -2 * pkin(1);
t176 = pkin(4) ^ 2;
t149 = pkin(18) + pkin(19);
t146 = sin(t149);
t147 = cos(t149);
t157 = sin(qJ(3));
t175 = pkin(5) ^ 2;
t180 = pkin(1) ^ 2;
t158 = sin(qJ(2));
t160 = sin(pkin(16));
t164 = cos(qJ(2));
t166 = cos(pkin(16));
t136 = t158 * t160 - t164 * t166;
t226 = pkin(5) * t136;
t204 = t226 * t248 + t180;
t126 = t175 + t204;
t202 = pkin(2) ^ 2 - pkin(6) ^ 2;
t120 = t126 + t202;
t129 = pkin(1) - t226;
t137 = t158 * t166 + t160 * t164;
t242 = pkin(5) + pkin(6);
t243 = pkin(5) - pkin(6);
t115 = (pkin(2) + t242) * (-pkin(2) + t243) + t204;
t116 = (-pkin(2) + t242) * (pkin(2) + t243) + t204;
t181 = sqrt(-t116 * t115);
t99 = pkin(5) * t120 * t137 + t129 * t181;
t213 = t99 * t157;
t163 = cos(qJ(3));
t205 = t137 * t181;
t98 = -pkin(5) * t205 + t120 * t129;
t216 = t98 * t163;
t186 = -t213 / 0.2e1 + t216 / 0.2e1;
t124 = 0.1e1 / t126;
t179 = 0.1e1 / pkin(2);
t207 = t124 * t179;
t86 = t186 * t207;
t212 = t99 * t163;
t217 = t98 * t157;
t185 = t212 / 0.2e1 + t217 / 0.2e1;
t87 = t185 * t207;
t70 = -t146 * t87 + t147 * t86;
t234 = pkin(4) * t70;
t211 = -0.2e1 * pkin(3) * t234 + t176;
t241 = -pkin(8) - pkin(10);
t57 = (pkin(3) - t241) * (pkin(3) + t241) + t211;
t240 = pkin(10) - pkin(8);
t58 = (pkin(3) - t240) * (pkin(3) + t240) + t211;
t182 = sqrt(-t58 * t57);
t69 = -t146 * t86 - t147 * t87;
t220 = t182 * t69;
t203 = pkin(8) ^ 2 - pkin(10) ^ 2;
t177 = pkin(3) ^ 2;
t65 = t177 + t211;
t59 = t65 + t203;
t66 = -pkin(3) + t234;
t40 = -pkin(4) * t220 - t59 * t66;
t247 = t40 / 0.2e1;
t63 = 0.1e1 / t65;
t246 = t63 / 0.2e1;
t133 = t137 * qJD(2);
t134 = t136 * qJD(2);
t200 = pkin(1) * pkin(5) * t133;
t210 = 0.2e1 * (t115 + t116) * t200 / t181;
t196 = -t210 / 0.2e1;
t183 = t134 * t181 + t137 * t196;
t78 = ((t129 * t248 - t120) * t133 + t183) * pkin(5);
t245 = -t78 / 0.2e1;
t199 = -0.2e1 * t133 * t137;
t206 = t133 * t181;
t79 = t129 * t210 / 0.2e1 + t175 * pkin(1) * t199 + (-t120 * t134 - t206) * pkin(5);
t244 = t79 / 0.2e1;
t159 = sin(qJ(1));
t165 = cos(qJ(1));
t138 = -t159 * V_base(4) + t165 * V_base(5);
t135 = qJD(2) - t138;
t154 = cos(pkin(19));
t194 = 0.1e1 / t126 ^ 2 * t200;
t214 = t99 * t154;
t152 = sin(pkin(19));
t215 = t99 * t152;
t218 = t98 * t154;
t219 = t98 * t152;
t230 = t152 / 0.2e1;
t83 = (-t218 / 0.2e1 + t215 / 0.2e1) * t207;
t81 = 0.1e1 / t83 ^ 2;
t84 = (t214 / 0.2e1 + t219 / 0.2e1) * t207;
t33 = t135 + (-((t154 * t245 + t230 * t79) * t124 + (t215 - t218) * t194) * t84 * t81 + ((t154 * t244 + t230 * t78) * t124 + (t214 + t219) * t194) / t83) * t179 / (t81 * t84 ^ 2 + 0.1e1);
t239 = pkin(3) * t33;
t228 = t157 / 0.2e1;
t46 = ((t212 + t217) * t194 + (qJD(3) * t186 + t163 * t244 + t228 * t78) * t124) * t179;
t47 = ((t213 - t216) * t194 + (qJD(3) * t185 + t163 * t245 + t228 * t79) * t124) * t179;
t43 = -t146 * t46 - t147 * t47;
t238 = pkin(3) * t43;
t237 = pkin(3) * t69;
t235 = pkin(4) * t69;
t41 = -t182 * t66 + t235 * t59;
t236 = pkin(4) * t41;
t150 = sin(pkin(17));
t233 = t150 / 0.2e1;
t151 = cos(pkin(17));
t232 = -t151 / 0.2e1;
t231 = t151 / 0.2e1;
t155 = cos(pkin(18));
t229 = -t155 / 0.2e1;
t167 = cos(pkin(15));
t227 = t167 / 0.2e1;
t201 = pkin(4) * t238;
t225 = 0.2e1 * (t57 + t58) * t201 / t182;
t60 = t65 - t203;
t67 = -pkin(3) * t70 + pkin(4);
t187 = -pkin(3) * t220 + t60 * t67;
t188 = t182 * t67 + t237 * t60;
t170 = 0.1e1 / pkin(10);
t223 = t170 * t63;
t25 = (t187 * t232 + t188 * t233) * t223;
t26 = (t187 * t233 + t188 * t231) * t223;
t132 = qJD(3) + t135;
t143 = pkin(12) * V_base(5) + V_base(1);
t144 = -pkin(12) * V_base(4) + V_base(2);
t122 = t165 * t143 + t159 * t144;
t127 = -pkin(13) * t138 + V_base(3);
t110 = -t122 * t158 + t164 * t127;
t102 = pkin(1) * t135 + t110;
t111 = t122 * t164 + t127 * t158;
t94 = -t102 * t163 + t157 * t111;
t89 = pkin(4) * t132 + t94;
t95 = -t102 * t157 - t111 * t163;
t16 = t25 * t95 + t26 * t89;
t51 = t84 * t102 + t83 * t111;
t52 = t83 * t102 - t84 * t111;
t172 = 0.1e1 / pkin(8);
t222 = t172 * t63;
t221 = t182 * t43;
t161 = sin(pkin(15));
t209 = t124 * t161;
t174 = 0.1e1 / pkin(6);
t208 = t124 * t174;
t198 = -t225 / 0.2e1;
t64 = 0.1e1 / t65 ^ 2;
t197 = t64 * t201;
t195 = t124 * t227;
t121 = -t159 * t143 + t144 * t165;
t192 = t161 * t194;
t191 = t167 * t194;
t17 = t25 * t89 - t26 * t95;
t139 = t159 * V_base(5) + t165 * V_base(4);
t148 = V_base(6) + qJD(1);
t117 = -t139 * t158 + t148 * t164;
t118 = t139 * t164 + t148 * t158;
t103 = -t117 * t163 + t118 * t157;
t104 = -t117 * t157 - t118 * t163;
t20 = t103 * t25 - t104 * t26;
t190 = t150 * t60 + t151 * t182;
t189 = t150 * t182 - t151 * t60;
t114 = -pkin(13) * t148 - t121;
t42 = t146 * t47 - t147 * t46;
t184 = -t42 * t182 + t198 * t69;
t106 = -pkin(1) * t117 + t114;
t96 = -pkin(4) * t103 + t106;
t168 = V_base(3) ^ 2;
t162 = cos(qJ(4));
t156 = sin(qJ(4));
t153 = sin(pkin(18));
t130 = pkin(1) * t136 - pkin(5);
t128 = pkin(7) * t138 + V_base(3);
t119 = t126 - t202;
t113 = pkin(14) * t138 + t122;
t112 = pkin(7) * t148 + pkin(14) * t139 - t121;
t105 = t106 ^ 2;
t100 = pkin(1) * t119 * t137 - t130 * t181;
t97 = -pkin(1) * t205 - t119 * t130;
t88 = (t100 * t227 - t97 * t161 / 0.2e1) * t208;
t85 = (t97 * t227 + t100 * t161 / 0.2e1) * t208;
t82 = 0.1e1 / t85 ^ 2;
t80 = t130 * t196 + t180 * pkin(5) * t199 + (-t119 * t134 - t206) * pkin(1);
t77 = ((0.2e1 * pkin(5) * t130 - t119) * t133 + t183) * pkin(1);
t62 = -t139 * t88 + t148 * t85;
t61 = t139 * t85 + t148 * t88;
t56 = t117 * t83 - t118 * t84;
t55 = t117 * t84 + t118 * t83;
t54 = -t113 * t88 + t128 * t85;
t53 = t113 * t85 + t128 * t88;
t45 = (-t153 * t55 - t155 * t56) * pkin(3) + t106;
t39 = 0.1e1 / t40 ^ 2;
t34 = ((t80 * t195 + t100 * t191 - t77 * t209 / 0.2e1 - t97 * t192) / t85 - (t77 * t195 + t97 * t191 + t80 * t209 / 0.2e1 + t100 * t192) * t88 * t82) / (t82 * t88 ^ 2 + 0.1e1) * t174 - t138;
t32 = t155 * t239 + t52;
t31 = t153 * t239 + t51;
t28 = (t40 * t229 - t153 * t41 / 0.2e1) * t222;
t27 = (t153 * t247 + t229 * t41) * t222;
t24 = 0.1e1 / t25 ^ 2;
t19 = t103 * t26 + t104 * t25;
t18 = qJD(4) - t20;
t15 = -t27 * t55 + t28 * t56;
t14 = t27 * t56 + t28 * t55;
t13 = t67 * t225 / 0.2e1 - 0.2e1 * t177 * t43 * t235 + (t42 * t60 - t221) * pkin(3);
t12 = ((-0.2e1 * pkin(4) * t67 - t60) * t43 + t184) * pkin(3);
t11 = -t27 * t31 + t28 * t32;
t10 = t27 * t32 + t28 * t31;
t9 = -pkin(9) * t20 - pkin(11) * t19 + t96;
t8 = t33 + (((t66 * t198 + (t42 * t59 - t221) * pkin(4)) * t246 + (-t176 * t63 * t69 + t236 * t64) * t238) / t247 - 0.2e1 * ((-t43 * t59 + t184) * t246 + (t40 * t64 + t63 * t66) * t238) * t39 * t236) * pkin(8) * t172 / (t39 * t41 ^ 2 + 0.1e1) * t65;
t7 = t132 + (((t12 * t233 + t13 * t231) * t63 + (-t189 * t237 + t190 * t67) * t197) / t25 - ((t12 * t232 + t13 * t233) * t63 + (t189 * t67 + t190 * t237) * t197) * t26 * t24) * t170 / (t24 * t26 ^ 2 + 0.1e1);
t6 = t156 * t7 + t162 * t19;
t5 = -t156 * t19 + t162 * t7;
t4 = -pkin(9) * t7 - t17;
t3 = pkin(11) * t7 + t16;
t2 = t156 * t9 + t162 * t3;
t1 = -t156 * t3 + t162 * t9;
t21 = (t52 * mrSges(8,1) - t51 * mrSges(8,2) + Ifges(8,5) * t55 + Ifges(8,6) * t56 + Ifges(8,3) * t33 / 0.2e1) * t33 + (-t114 * mrSges(3,1) + t111 * mrSges(3,3) + Ifges(3,4) * t118 + Ifges(3,6) * t135 + Ifges(3,2) * t117 / 0.2e1) * t117 + (t11 * mrSges(9,1) - t10 * mrSges(9,2) + Ifges(9,3) * t8 / 0.2e1) * t8 + (-t106 * mrSges(4,1) + t95 * mrSges(4,3) + Ifges(4,4) * t104 + Ifges(4,6) * t132 + Ifges(4,2) * t103 / 0.2e1) * t103 + (t106 * mrSges(8,2) - t52 * mrSges(8,3) + Ifges(8,4) * t56 + Ifges(8,1) * t55 / 0.2e1) * t55 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t114 * mrSges(3,2) - t110 * mrSges(3,3) + Ifges(3,5) * t135 + Ifges(3,1) * t118 / 0.2e1) * t118 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t6 + Ifges(6,6) * t5 + Ifges(6,3) * t18 / 0.2e1) * t18 + (t96 * mrSges(5,2) - t17 * mrSges(5,3) + Ifges(5,4) * t20 + Ifges(5,5) * t7 + Ifges(5,1) * t19 / 0.2e1) * t19 + (-V_base(3) * mrSges(2,1) + t122 * mrSges(2,3) + Ifges(2,4) * t139 + Ifges(2,6) * t148 + Ifges(2,2) * t138 / 0.2e1) * t138 + (t106 * mrSges(4,2) - t94 * mrSges(4,3) + Ifges(4,5) * t132 + Ifges(4,1) * t104 / 0.2e1) * t104 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t45 * mrSges(9,1) + t10 * mrSges(9,3) + Ifges(9,6) * t8 + Ifges(9,2) * t15 / 0.2e1) * t15 + (-t96 * mrSges(5,1) + t16 * mrSges(5,3) + Ifges(5,6) * t7 + Ifges(5,2) * t20 / 0.2e1) * t20 + (t112 * mrSges(7,2) - t54 * mrSges(7,3) + Ifges(7,4) * t62 + Ifges(7,1) * t61 / 0.2e1) * t61 + (t54 * mrSges(7,1) - t53 * mrSges(7,2) + Ifges(7,5) * t61 + Ifges(7,6) * t62 + Ifges(7,3) * t34 / 0.2e1) * t34 + (V_base(3) * mrSges(2,2) - t121 * mrSges(2,3) + Ifges(2,5) * t148 + Ifges(2,1) * t139 / 0.2e1) * t139 + (t121 * mrSges(2,1) - t122 * mrSges(2,2) + Ifges(2,3) * t148 / 0.2e1) * t148 + (-t112 * mrSges(7,1) + t53 * mrSges(7,3) + Ifges(7,2) * t62 / 0.2e1) * t62 + (t45 * mrSges(9,2) - t11 * mrSges(9,3) + Ifges(9,4) * t15 + Ifges(9,5) * t8 + Ifges(9,1) * t14 / 0.2e1) * t14 + (t110 * mrSges(3,1) - t111 * mrSges(3,2) + Ifges(3,3) * t135 / 0.2e1) * t135 + (t4 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,1) * t6 / 0.2e1) * t6 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t17 * mrSges(5,1) - t16 * mrSges(5,2) + Ifges(5,3) * t7 / 0.2e1) * t7 + (t94 * mrSges(4,1) - t95 * mrSges(4,2) + Ifges(4,3) * t132 / 0.2e1) * t132 + (-t4 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t6 + Ifges(6,2) * t5 / 0.2e1) * t5 + (-t106 * mrSges(8,1) + t51 * mrSges(8,3) + Ifges(8,2) * t56 / 0.2e1) * t56 + m(2) * (t121 ^ 2 + t122 ^ 2 + t168) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t168) / 0.2e1 + m(4) * (t94 ^ 2 + t95 ^ 2 + t105) / 0.2e1 + m(8) * (t51 ^ 2 + t52 ^ 2 + t105) / 0.2e1 + m(7) * (t112 ^ 2 + t53 ^ 2 + t54 ^ 2) / 0.2e1 + m(3) * (t110 ^ 2 + t111 ^ 2 + t114 ^ 2) / 0.2e1 + m(5) * (t16 ^ 2 + t17 ^ 2 + t96 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) / 0.2e1 + m(9) * (t10 ^ 2 + t11 ^ 2 + t45 ^ 2) / 0.2e1;
T = t21;
