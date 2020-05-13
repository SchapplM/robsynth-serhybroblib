% Calculate kinetic energy for
% palh3m1DE2
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
% Datum: 2020-04-20 16:51
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m1DE2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(19,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE2_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m1DE2_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh3m1DE2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE2_energykin_floatb_twist_slag_vp2: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1DE2_energykin_floatb_twist_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1DE2_energykin_floatb_twist_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m1DE2_energykin_floatb_twist_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-19 19:34:28
% EndTime: 2020-04-19 19:35:16
% DurationCPUTime: 43.90s
% Computational Cost: add. (1176093->311), mult. (1814235->539), div. (76552->22), fcn. (1135024->36), ass. (0->217)
t249 = -2 * pkin(1);
t182 = pkin(3) ^ 2;
t181 = pkin(4) ^ 2;
t154 = pkin(18) + pkin(19);
t151 = sin(t154);
t152 = cos(t154);
t180 = pkin(5) ^ 2;
t185 = pkin(1) ^ 2;
t163 = sin(qJ(2));
t165 = sin(pkin(16));
t169 = cos(qJ(2));
t171 = cos(pkin(16));
t141 = t163 * t165 - t169 * t171;
t231 = pkin(5) * t141;
t206 = t231 * t249 + t185;
t131 = t180 + t206;
t204 = pkin(2) ^ 2 - pkin(6) ^ 2;
t125 = t131 + t204;
t134 = pkin(1) - t231;
t142 = t163 * t171 + t165 * t169;
t244 = pkin(5) + pkin(6);
t245 = pkin(5) - pkin(6);
t120 = (pkin(2) + t244) * (-pkin(2) + t245) + t206;
t121 = (-pkin(2) + t244) * (pkin(2) + t245) + t206;
t186 = sqrt(-t121 * t120);
t104 = pkin(5) * t125 * t142 + t134 * t186;
t162 = sin(qJ(3));
t214 = t104 * t162;
t207 = t142 * t186;
t103 = -pkin(5) * t207 + t125 * t134;
t168 = cos(qJ(3));
t217 = t103 * t168;
t190 = -t214 / 0.2e1 + t217 / 0.2e1;
t129 = 0.1e1 / t131;
t184 = 0.1e1 / pkin(2);
t209 = t129 * t184;
t95 = t190 * t209;
t213 = t104 * t168;
t218 = t103 * t162;
t189 = t213 / 0.2e1 + t218 / 0.2e1;
t96 = t189 * t209;
t77 = -t151 * t96 + t152 * t95;
t236 = pkin(4) * t77;
t221 = -0.2e1 * pkin(3) * t236 + t181;
t72 = t182 + t221;
t70 = 0.1e1 / t72;
t248 = t70 / 0.2e1;
t138 = t142 * qJD(2);
t139 = t141 * qJD(2);
t202 = pkin(1) * pkin(5) * t138;
t212 = 0.2e1 * (t120 + t121) * t202 / t186;
t197 = -t212 / 0.2e1;
t188 = t139 * t186 + t142 * t197;
t87 = ((t134 * t249 - t125) * t138 + t188) * pkin(5);
t247 = -t87 / 0.2e1;
t201 = -0.2e1 * t138 * t142;
t208 = t138 * t186;
t88 = t134 * t212 / 0.2e1 + t180 * pkin(1) * t201 + (-t125 * t139 - t208) * pkin(5);
t246 = t88 / 0.2e1;
t243 = -pkin(8) - pkin(10);
t242 = pkin(10) - pkin(8);
t177 = 0.1e1 / pkin(8);
t199 = t177 * t248;
t66 = (pkin(3) - t243) * (pkin(3) + t243) + t221;
t67 = (pkin(3) - t242) * (pkin(3) + t242) + t221;
t187 = sqrt(-t67 * t66);
t76 = -t151 * t95 - t152 * t96;
t226 = t187 * t76;
t205 = pkin(8) ^ 2 - pkin(10) ^ 2;
t68 = t72 + t205;
t73 = -pkin(3) + t236;
t42 = -pkin(4) * t226 - t68 * t73;
t237 = pkin(4) * t76;
t44 = -t187 * t73 + t237 * t68;
t34 = atan2(t44 * t199, t42 * t199);
t241 = sin(t34);
t164 = sin(qJ(1));
t170 = cos(qJ(1));
t143 = -t164 * V_base(4) + t170 * V_base(5);
t140 = qJD(2) - t143;
t159 = cos(pkin(19));
t195 = 0.1e1 / t131 ^ 2 * t202;
t215 = t104 * t159;
t157 = sin(pkin(19));
t216 = t104 * t157;
t219 = t103 * t159;
t220 = t103 * t157;
t234 = t157 / 0.2e1;
t92 = (-t219 / 0.2e1 + t216 / 0.2e1) * t209;
t90 = 0.1e1 / t92 ^ 2;
t93 = (t215 / 0.2e1 + t220 / 0.2e1) * t209;
t39 = t140 + (-((t159 * t247 + t234 * t88) * t129 + (t216 - t219) * t195) * t93 * t90 + ((t159 * t246 + t234 * t87) * t129 + (t215 + t220) * t195) / t92) * t184 / (t90 * t93 ^ 2 + 0.1e1);
t240 = pkin(3) * t39;
t233 = t162 / 0.2e1;
t50 = ((t213 + t218) * t195 + (qJD(3) * t190 + t168 * t246 + t233 * t87) * t129) * t184;
t51 = ((t214 - t217) * t195 + (qJD(3) * t189 + t168 * t247 + t233 * t88) * t129) * t184;
t47 = -t151 * t50 - t152 * t51;
t239 = pkin(3) * t47;
t238 = pkin(4) * t44;
t155 = sin(pkin(17));
t235 = t155 / 0.2e1;
t172 = cos(pkin(15));
t232 = t172 / 0.2e1;
t203 = pkin(4) * t239;
t230 = 0.2e1 * (t66 + t67) * t203 / t187;
t148 = pkin(12) * V_base(5) + V_base(1);
t149 = -pkin(12) * V_base(4) + V_base(2);
t127 = t170 * t148 + t164 * t149;
t132 = -pkin(13) * t143 + V_base(3);
t115 = -t127 * t163 + t169 * t132;
t107 = pkin(1) * t140 + t115;
t116 = t127 * t169 + t132 * t163;
t100 = -t107 * t162 - t116 * t168;
t69 = t72 - t205;
t74 = -pkin(3) * t77 + pkin(4);
t45 = pkin(3) * t69 * t76 + t187 * t74;
t223 = t45 * t155;
t156 = cos(pkin(17));
t43 = -pkin(3) * t226 + t69 * t74;
t224 = t43 * t156;
t175 = 0.1e1 / pkin(10);
t228 = t175 * t70;
t31 = (-t224 / 0.2e1 + t223 / 0.2e1) * t228;
t222 = t45 * t156;
t225 = t43 * t155;
t32 = (t222 / 0.2e1 + t225 / 0.2e1) * t228;
t28 = atan2(t32, t31);
t25 = sin(t28);
t26 = cos(t28);
t137 = qJD(3) + t140;
t99 = -t107 * t168 + t162 * t116;
t98 = pkin(4) * t137 + t99;
t15 = t26 * t100 + t25 * t98;
t84 = atan2(t93, t92);
t78 = sin(t84);
t79 = cos(t84);
t53 = t78 * t107 + t79 * t116;
t227 = t187 * t47;
t166 = sin(pkin(15));
t211 = t129 * t166;
t179 = 0.1e1 / pkin(6);
t210 = t129 * t179;
t200 = -t230 / 0.2e1;
t71 = 0.1e1 / t72 ^ 2;
t198 = t71 * t203;
t196 = t129 * t232;
t52 = t79 * t107 - t116 * t78;
t126 = -t164 * t148 + t149 * t170;
t193 = t166 * t195;
t192 = t172 * t195;
t14 = -t100 * t25 + t26 * t98;
t144 = t164 * V_base(5) + t170 * V_base(4);
t153 = V_base(6) + qJD(1);
t122 = -t144 * t163 + t153 * t169;
t123 = t144 * t169 + t153 * t163;
t108 = -t122 * t168 + t123 * t162;
t109 = -t122 * t162 - t123 * t168;
t17 = t108 * t26 - t109 * t25;
t119 = -pkin(13) * t153 - t126;
t46 = t151 * t51 - t152 * t50;
t191 = -t46 * t187 + t200 * t76;
t111 = -pkin(1) * t122 + t119;
t101 = -pkin(4) * t108 + t111;
t173 = V_base(3) ^ 2;
t167 = cos(qJ(4));
t161 = sin(qJ(4));
t160 = cos(pkin(18));
t158 = sin(pkin(18));
t135 = pkin(1) * t141 - pkin(5);
t133 = pkin(7) * t143 + V_base(3);
t124 = t131 - t204;
t118 = pkin(14) * t143 + t127;
t117 = pkin(7) * t153 + pkin(14) * t144 - t126;
t110 = t111 ^ 2;
t105 = pkin(1) * t124 * t142 - t135 * t186;
t102 = -pkin(1) * t207 - t124 * t135;
t97 = (t105 * t232 - t102 * t166 / 0.2e1) * t210;
t94 = (t102 * t232 + t105 * t166 / 0.2e1) * t210;
t91 = 0.1e1 / t94 ^ 2;
t89 = t135 * t197 + t185 * pkin(5) * t201 + (-t124 * t139 - t208) * pkin(1);
t86 = ((0.2e1 * pkin(5) * t135 - t124) * t138 + t188) * pkin(1);
t85 = atan2(t97, t94);
t82 = cos(t85);
t81 = sin(t85);
t61 = t144 * t82 + t153 * t81;
t60 = -t144 * t81 + t153 * t82;
t57 = t122 * t78 + t123 * t79;
t56 = t122 * t79 - t123 * t78;
t55 = t118 * t82 + t133 * t81;
t54 = -t118 * t81 + t133 * t82;
t48 = (-t158 * t57 - t160 * t56) * pkin(3) + t111;
t41 = 0.1e1 / t42 ^ 2;
t40 = ((t89 * t196 + t105 * t192 - t86 * t211 / 0.2e1 - t102 * t193) / t94 - (t86 * t196 + t102 * t192 + t89 * t211 / 0.2e1 + t105 * t193) * t97 * t91) / (t91 * t97 ^ 2 + 0.1e1) * t179 - t143;
t38 = t158 * t240 + t53;
t37 = t160 * t240 + t52;
t33 = cos(t34);
t30 = 0.1e1 / t31 ^ 2;
t24 = -t158 * t241 - t160 * t33;
t23 = t158 * t33 - t160 * t241;
t20 = t74 * t230 / 0.2e1 - 0.2e1 * t182 * t47 * t237 + (t46 * t69 - t227) * pkin(3);
t19 = ((-0.2e1 * pkin(4) * t74 - t69) * t47 + t191) * pkin(3);
t18 = t108 * t25 + t109 * t26;
t16 = qJD(4) - t17;
t13 = t23 * t56 + t24 * t57;
t12 = -t23 * t57 + t24 * t56;
t11 = t23 * t37 + t24 * t38;
t10 = -t23 * t38 + t24 * t37;
t9 = -pkin(9) * t17 - pkin(11) * t18 + t101;
t8 = t39 + 0.2e1 * (((t73 * t200 + (t46 * t68 - t227) * pkin(4)) * t248 + (-t181 * t70 * t76 + t238 * t71) * t239) / t42 - ((-t47 * t68 + t191) * t248 + (t42 * t71 + t70 * t73) * t239) * t41 * t238) * pkin(8) * t177 / (t41 * t44 ^ 2 + 0.1e1) * t72;
t7 = t137 + (((t20 * t156 / 0.2e1 + t19 * t235) * t70 + (t222 + t225) * t198) / t31 - ((-t19 * t156 / 0.2e1 + t20 * t235) * t70 + (t223 - t224) * t198) * t32 * t30) * t175 / (t30 * t32 ^ 2 + 0.1e1);
t6 = t161 * t7 + t167 * t18;
t5 = -t161 * t18 + t167 * t7;
t4 = pkin(11) * t7 + t15;
t3 = -pkin(9) * t7 - t14;
t2 = t161 * t9 + t167 * t4;
t1 = -t161 * t4 + t167 * t9;
t21 = (-t101 * mrSges(5,1) + t15 * mrSges(5,3) + Ifges(5,4) * t18 + Ifges(5,6) * t7 + Ifges(5,2) * t17 / 0.2e1) * t17 + (-t117 * mrSges(7,1) + t55 * mrSges(7,3) + Ifges(7,4) * t61 + Ifges(7,2) * t60 / 0.2e1) * t60 + (t14 * mrSges(5,1) - t15 * mrSges(5,2) + Ifges(5,3) * t7 / 0.2e1) * t7 + (t101 * mrSges(5,2) - t14 * mrSges(5,3) + Ifges(5,5) * t7 + Ifges(5,1) * t18 / 0.2e1) * t18 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (t119 * mrSges(3,2) - t115 * mrSges(3,3) + Ifges(3,5) * t140 + Ifges(3,1) * t123 / 0.2e1) * t123 + (t111 * mrSges(8,2) - t52 * mrSges(8,3) + Ifges(8,1) * t57 / 0.2e1) * t57 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-t48 * mrSges(9,1) + t11 * mrSges(9,3) + Ifges(9,4) * t13 + Ifges(9,6) * t8 + Ifges(9,2) * t12 / 0.2e1) * t12 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t6 + Ifges(6,6) * t5 + Ifges(6,3) * t16 / 0.2e1) * t16 + (-t119 * mrSges(3,1) + t116 * mrSges(3,3) + Ifges(3,4) * t123 + Ifges(3,6) * t140 + Ifges(3,2) * t122 / 0.2e1) * t122 + (t117 * mrSges(7,2) - t54 * mrSges(7,3) + Ifges(7,1) * t61 / 0.2e1) * t61 + (-V_base(3) * mrSges(2,1) + t127 * mrSges(2,3) + Ifges(2,4) * t144 + Ifges(2,6) * t153 + Ifges(2,2) * t143 / 0.2e1) * t143 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t54 * mrSges(7,1) - t55 * mrSges(7,2) + Ifges(7,5) * t61 + Ifges(7,6) * t60 + Ifges(7,3) * t40 / 0.2e1) * t40 + (t111 * mrSges(4,2) - t99 * mrSges(4,3) + Ifges(4,5) * t137 + Ifges(4,1) * t109 / 0.2e1) * t109 + (t52 * mrSges(8,1) - t53 * mrSges(8,2) + Ifges(8,5) * t57 + Ifges(8,6) * t56 + Ifges(8,3) * t39 / 0.2e1) * t39 + (-t111 * mrSges(8,1) + t53 * mrSges(8,3) + Ifges(8,4) * t57 + Ifges(8,2) * t56 / 0.2e1) * t56 + (V_base(3) * mrSges(2,2) - t126 * mrSges(2,3) + Ifges(2,5) * t153 + Ifges(2,1) * t144 / 0.2e1) * t144 + (t48 * mrSges(9,2) - t10 * mrSges(9,3) + Ifges(9,5) * t8 + Ifges(9,1) * t13 / 0.2e1) * t13 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t126 * mrSges(2,1) - t127 * mrSges(2,2) + Ifges(2,3) * t153 / 0.2e1) * t153 + (-t111 * mrSges(4,1) + t100 * mrSges(4,3) + Ifges(4,4) * t109 + Ifges(4,6) * t137 + Ifges(4,2) * t108 / 0.2e1) * t108 + (t10 * mrSges(9,1) - t11 * mrSges(9,2) + Ifges(9,3) * t8 / 0.2e1) * t8 + (t115 * mrSges(3,1) - t116 * mrSges(3,2) + Ifges(3,3) * t140 / 0.2e1) * t140 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,1) * t6 / 0.2e1) * t6 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t6 + Ifges(6,2) * t5 / 0.2e1) * t5 + (t99 * mrSges(4,1) - t100 * mrSges(4,2) + Ifges(4,3) * t137 / 0.2e1) * t137 + m(2) * (t126 ^ 2 + t127 ^ 2 + t173) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t173) / 0.2e1 + m(7) * (t117 ^ 2 + t54 ^ 2 + t55 ^ 2) / 0.2e1 + m(3) * (t115 ^ 2 + t116 ^ 2 + t119 ^ 2) / 0.2e1 + m(4) * (t100 ^ 2 + t99 ^ 2 + t110) / 0.2e1 + m(8) * (t52 ^ 2 + t53 ^ 2 + t110) / 0.2e1 + m(5) * (t101 ^ 2 + t14 ^ 2 + t15 ^ 2) / 0.2e1 + m(9) * (t10 ^ 2 + t11 ^ 2 + t48 ^ 2) / 0.2e1;
T = t21;
