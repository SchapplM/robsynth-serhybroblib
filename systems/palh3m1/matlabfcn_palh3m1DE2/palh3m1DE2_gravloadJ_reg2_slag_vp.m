% Calculate inertial parameters regressor of gravitation load for
% palh3m1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 16:51
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = palh3m1DE2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE2_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1DE2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE2_gravloadJ_reg2_slag_vp: pkin has to be [19x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t264 = 2 * pkin(3);
t145 = pkin(4) ^ 2;
t121 = pkin(18) + pkin(19);
t116 = sin(t121);
t117 = cos(t121);
t144 = pkin(5) ^ 2;
t149 = pkin(1) ^ 2;
t129 = sin(qJ(2));
t131 = sin(pkin(16));
t135 = cos(qJ(2));
t237 = cos(pkin(16));
t109 = t129 * t131 - t135 * t237;
t234 = pkin(5) * t109;
t258 = -2 * pkin(1);
t200 = t234 * t258 + t149;
t104 = t144 + t200;
t102 = 0.1e1 / t104;
t148 = 0.1e1 / pkin(2);
t207 = t102 * t148;
t128 = sin(qJ(3));
t105 = pkin(1) - t234;
t251 = -pkin(6) - pkin(2);
t94 = (pkin(5) - t251) * (pkin(5) + t251) + t200;
t250 = -pkin(6) + pkin(2);
t95 = (pkin(5) - t250) * (pkin(5) + t250) + t200;
t151 = sqrt(-t95 * t94);
t110 = t129 * t237 + t135 * t131;
t233 = pkin(5) * t110;
t198 = pkin(2) ^ 2 - pkin(6) ^ 2;
t99 = t104 + t198;
t89 = t105 * t151 + t99 * t233;
t213 = t89 * t128;
t134 = cos(qJ(3));
t205 = t110 * t151;
t88 = -pkin(5) * t205 + t105 * t99;
t216 = t88 * t134;
t82 = (-t216 / 0.2e1 + t213 / 0.2e1) * t207;
t212 = t89 * t134;
t217 = t88 * t128;
t83 = (t212 / 0.2e1 + t217 / 0.2e1) * t207;
t59 = t116 * t83 + t117 * t82;
t244 = t59 * pkin(4);
t210 = t244 * t264 + t145;
t249 = (-pkin(8) - pkin(10));
t47 = ((pkin(3) - t249) * (pkin(3) + t249)) + t210;
t248 = (-pkin(8) + pkin(10));
t48 = ((pkin(3) - t248) * (pkin(3) + t248)) + t210;
t150 = sqrt(-t48 * t47);
t169 = t116 * t82 - t117 * t83;
t262 = t150 * t169;
t130 = sin(qJ(1));
t136 = cos(qJ(1));
t261 = -g(1) * t130 + g(2) * t136;
t260 = -g(1) * t136 - g(2) * t130;
t122 = qJ(2) + qJ(3);
t124 = cos(pkin(17));
t139 = 0.1e1 / pkin(10);
t146 = pkin(3) ^ 2;
t53 = t146 + t210;
t51 = 0.1e1 / t53;
t223 = t139 * t51;
t123 = sin(pkin(17));
t241 = t123 / 0.2e1;
t199 = pkin(8) ^ 2 - pkin(10) ^ 2;
t50 = t53 - t199;
t55 = pkin(3) * t59 + pkin(4);
t35 = -pkin(3) * t262 + t50 * t55;
t37 = pkin(3) * t169 * t50 + t150 * t55;
t25 = (-t35 * t124 / 0.2e1 + t37 * t241) * t223;
t26 = (t37 * t124 / 0.2e1 + t35 * t241) * t223;
t16 = atan2(t26, t25) + t122;
t14 = sin(t16);
t15 = cos(t16);
t259 = -g(3) * t15 - t14 * t260;
t125 = sin(pkin(19));
t215 = t89 * t125;
t126 = cos(pkin(19));
t218 = t88 * t126;
t79 = (-t218 / 0.2e1 + t215 / 0.2e1) * t207;
t257 = 0.1e1 / t79 ^ 2;
t256 = -0.2e1 * t110 ^ 2;
t43 = 0.1e1 / t150;
t255 = -t43 / 0.2e1;
t254 = t51 / 0.2e1;
t196 = pkin(1) * t233;
t226 = 0.2e1 / t151 * (t94 + t95) * t196;
t192 = -t226 / 0.2e1;
t206 = t109 * t151;
t74 = (t206 + (t105 * t258 + t192 - t99) * t110) * pkin(5);
t253 = -t74 / 0.2e1;
t76 = t105 * t226 / 0.2e1 + t144 * pkin(1) * t256 + (-t109 * t99 - t205) * pkin(5);
t252 = t76 / 0.2e1;
t214 = t89 * t126;
t219 = t88 * t125;
t80 = (t214 / 0.2e1 + t219 / 0.2e1) * t207;
t66 = qJ(2) + atan2(t80, t79);
t65 = pkin(18) - t66;
t247 = pkin(3) * cos(t65);
t52 = 0.1e1 / t53 ^ 2;
t246 = pkin(4) * t52;
t245 = pkin(4) * t169;
t243 = g(3) * t14;
t240 = t125 / 0.2e1;
t239 = t128 / 0.2e1;
t137 = cos(pkin(15));
t238 = t137 / 0.2e1;
t118 = sin(t122);
t236 = pkin(4) * t118;
t119 = cos(t122);
t235 = pkin(4) * t119;
t228 = t129 * pkin(1);
t24 = 0.1e1 / t25 ^ 2;
t227 = t24 * t26;
t225 = t124 * t51;
t224 = t139 / (t24 * t26 ^ 2 + 0.1e1);
t181 = 0.1e1 / t104 ^ 2 * t196;
t45 = ((t134 * t252 + t74 * t239) * t102 + (t212 + t217) * t181) * t148;
t46 = ((t134 * t253 + t76 * t239) * t102 + (t213 - t216) * t181) * t148;
t41 = -t116 * t45 - t117 * t46;
t222 = t150 * t41;
t120 = t135 * pkin(1);
t211 = t120 + pkin(13);
t132 = sin(pkin(15));
t209 = t102 * t132;
t143 = 0.1e1 / pkin(6);
t208 = t102 * t143;
t127 = sin(qJ(4));
t204 = t127 * t130;
t203 = t127 * t136;
t133 = cos(qJ(4));
t202 = t130 * t133;
t201 = t133 * t136;
t197 = pkin(3) * t246;
t54 = -pkin(3) - t244;
t195 = t54 * t255;
t194 = t43 * t55 / 0.2e1;
t193 = t169 * t255;
t191 = t51 * t241;
t190 = -t225 / 0.2e1;
t189 = t225 / 0.2e1;
t141 = 0.1e1 / pkin(8);
t188 = t141 * t254;
t187 = -0.2e1 * t55 * pkin(4) - t50;
t186 = t102 * t238;
t185 = -0.2e1 * t146 * t245;
t184 = t123 * t197;
t183 = t124 * t197;
t182 = t120 - t235;
t49 = t53 + t199;
t34 = -pkin(4) * t262 - t49 * t54;
t33 = 0.1e1 / t34 ^ 2;
t36 = -t150 * t54 + t49 * t245;
t180 = pkin(8) * t141 / (t33 * t36 ^ 2 + 0.1e1) * t53;
t179 = t41 * t184;
t178 = t169 * t184;
t177 = t41 * t183;
t176 = t169 * t183;
t175 = 0.1e1 / t34 * t180;
t174 = -pkin(9) * t15 - pkin(11) * t14;
t173 = pkin(9) * t14 - pkin(11) * t15;
t172 = t132 * t181;
t171 = t137 * t181;
t30 = 0.1e1 + (((t126 * t252 + t74 * t240) * t102 + (t214 + t219) * t181) / t79 - ((t126 * t253 + t76 * t240) * t102 + (t215 - t218) * t181) * t80 * t257) * t148 / (t80 ^ 2 * t257 + 0.1e1);
t168 = pkin(4) * (t47 + t48) * t264;
t166 = pkin(3) * (t34 * t52 + t51 * t54);
t165 = pkin(4) * t33 * t36 * t180;
t28 = t41 * t168;
t40 = t116 * t46 - t117 * t45;
t164 = -t40 * t150 + t28 * t193;
t38 = t169 * t168;
t163 = -t59 * t150 + t38 * t193;
t161 = pkin(3) * (-t145 * t169 * t51 + t36 * t246);
t159 = -t15 * t260 + t243;
t22 = -atan2(t36 * t188, t34 * t188) + t65;
t20 = sin(t22);
t21 = cos(t22);
t157 = -g(3) * t20 - t21 * t260;
t156 = -g(3) * t21 + t20 * t260;
t96 = g(3) * t119 + t118 * t260;
t155 = -g(3) * t135 - t129 * t260;
t154 = t127 * t259;
t153 = t133 * t259;
t112 = -t228 + t236;
t111 = pkin(13) + t182;
t106 = pkin(1) * t109 - pkin(5);
t101 = t261 * t211;
t100 = t155 * pkin(1);
t98 = t104 - t198;
t97 = -g(3) * t118 + t119 * t260;
t90 = pkin(1) * t110 * t98 - t106 * t151;
t87 = -pkin(1) * t205 - t106 * t98;
t84 = (t90 * t238 - t87 * t132 / 0.2e1) * t208;
t81 = (t87 * t238 + t90 * t132 / 0.2e1) * t208;
t78 = 0.1e1 / t81 ^ 2;
t75 = t106 * t192 + t149 * pkin(5) * t256 + (-t109 * t98 - t205) * pkin(1);
t73 = (t206 + (0.2e1 * t106 * pkin(5) + t192 - t98) * t110) * pkin(1);
t72 = atan2(t84, t81);
t69 = cos(t72);
t68 = sin(t72);
t64 = cos(t66);
t63 = sin(t66);
t31 = ((t75 * t186 + t90 * t171 - t73 * t209 / 0.2e1 - t87 * t172) / t81 - (t73 * t186 + t87 * t171 + t75 * t209 / 0.2e1 + t90 * t172) * t84 * t78) / (t78 * t84 ^ 2 + 0.1e1) * t143;
t23 = 0.1e1 / t25;
t19 = t38 * t194 + t169 * t185 + (t59 * t50 - t262) * pkin(3);
t18 = (t169 * t187 + t163) * pkin(3);
t13 = t15 * t201 - t204;
t12 = t15 * t203 + t202;
t11 = t15 * t202 + t203;
t10 = t15 * t204 - t201;
t9 = t28 * t194 + t41 * t185 + (t40 * t50 - t222) * pkin(3);
t8 = (t187 * t41 + t164) * pkin(3);
t7 = t261 * t14;
t6 = -0.2e1 * ((t38 * t195 + (t59 * t49 - t262) * pkin(4)) * t254 + t169 * t161) * t175 + 0.2e1 * ((-t169 * t49 + t163) * t254 + t169 * t166) * t165;
t5 = -0.2e1 * ((t28 * t195 + (t40 * t49 - t222) * pkin(4)) * t254 + t41 * t161) * t175 + 0.2e1 * ((-t41 * t49 + t164) * t254 + t41 * t166) * t165 - t30;
t4 = 0.1e1 + ((t37 * t176 + t35 * t178 + t18 * t191 + t19 * t189) * t23 - (-t35 * t176 + t37 * t178 + t18 * t190 + t19 * t191) * t227) * t224;
t3 = 0.1e1 + ((t37 * t177 + t35 * t179 + t9 * t189 + t8 * t191) * t23 - (-t35 * t177 + t37 * t179 + t8 * t190 + t9 * t191) * t227) * t224;
t2 = t159 * t4;
t1 = t159 * t3;
t17 = [0, 0, 0, 0, 0, 0, -t261, -t260, 0, 0, 0, 0, 0, 0, 0, 0, -t261 * t135, t261 * t129, t260, -t261 * pkin(13), 0, 0, 0, 0, 0, 0, t261 * t119, -t261 * t118, t260, -t101, 0, 0, 0, 0, 0, 0, t261 * t15, -t7, t260, -t261 * t111, 0, 0, 0, 0, 0, 0, -g(1) * t11 + g(2) * t13, g(1) * t10 - g(2) * t12, t7, t261 * (-t111 - t174), 0, 0, 0, 0, 0, 0, -t261 * t69, t261 * t68, t260, t261 * pkin(7), 0, 0, 0, 0, 0, 0, -t261 * t64, t261 * t63, t260, -t101, 0, 0, 0, 0, 0, 0, t261 * t21, t261 * t20, t260, -t261 * (t211 + t247); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, g(3) * t129 - t135 * t260, 0, 0, 0, 0, 0, 0, 0, 0, t96, t97, 0, t100, 0, 0, 0, 0, 0, 0, -t259 * t3, -t1, 0, -g(3) * t182 + t260 * t112, 0, 0, 0, 0, 0, 0, -t3 * t153, t3 * t154, t1, -g(3) * (t174 * t3 + t182) + t260 * (t173 * t3 + t112), 0, 0, 0, 0, 0, 0, (-g(3) * t69 - t260 * t68) * t31, (g(3) * t68 - t260 * t69) * t31, 0, 0, 0, 0, 0, 0, 0, 0, (-g(3) * t64 - t260 * t63) * t30, (g(3) * t63 - t260 * t64) * t30, 0, t100, 0, 0, 0, 0, 0, 0, t156 * t5, t157 * t5, 0, -g(3) * (t30 * t247 + t120) + t260 * (pkin(3) * t30 * sin(t65) - t228); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t97, 0, 0, 0, 0, 0, 0, 0, 0, -t259 * t4, -t2, 0, t96 * pkin(4), 0, 0, 0, 0, 0, 0, -t4 * t153, t4 * t154, t2, -g(3) * (t174 * t4 - t235) + t260 * (t173 * t4 + t236), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156 * t6, t157 * t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10 - t127 * t243, -g(1) * t13 - g(2) * t11 - t133 * t243, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
taug_reg = t17;
