% Calculate potential energy for
% fivebar1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% m [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:03
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fivebar1DE2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fivebar1DE2_energypot_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'fivebar1DE2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1DE2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1DE2_energypot_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1DE2_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fivebar1DE2_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 04:27:03
% EndTime: 2020-04-27 04:27:08
% DurationCPUTime: 4.00s
% Computational Cost: add. (14199->425), mult. (41584->497), div. (128->5), fcn. (4996->12), ass. (0->204)
t265 = 4 * pkin(3);
t128 = pkin(4) ^ 2;
t125 = pkin(5) ^ 2;
t76 = -t125 / 0.3e1;
t264 = t76 - t128 / 0.3e1;
t124 = t125 ^ 2;
t127 = t128 ^ 2;
t263 = -t124 / 0.6e1 + t127 / 0.6e1;
t97 = cos(qJ(2));
t62 = t97 ^ 2;
t262 = -0.2e1 * t62;
t133 = pkin(3) ^ 2;
t131 = t133 ^ 2;
t261 = 4 * t131;
t260 = 2 * t133;
t111 = 6 * t133;
t137 = (pkin(2) ^ 2);
t259 = 2 * t137;
t139 = (pkin(1) ^ 2);
t122 = 2 * t139;
t144 = t137 ^ 2;
t114 = 5 * t144;
t258 = -pkin(4) - pkin(5);
t257 = -pkin(4) + pkin(5);
t96 = sin(qJ(1));
t54 = pkin(2) * t96;
t98 = cos(qJ(1));
t256 = pkin(2) * t98;
t95 = sin(qJ(2));
t255 = pkin(3) * t95;
t53 = pkin(3) * t97;
t147 = pkin(3) * t133;
t204 = t133 * t54;
t254 = (-t147 * t95 + t204) * t62;
t40 = pkin(1) - t256;
t253 = t40 * t97;
t56 = t133 + t139;
t181 = -t125 + t56;
t50 = t56 ^ 2;
t252 = t50 * (-t128 + t181);
t251 = t95 * t96;
t250 = t96 * t97;
t249 = t124 / 0.2e1 - t127 / 0.2e1;
t74 = -t125 / 0.6e1;
t248 = t74 - t128 / 0.6e1;
t78 = -0.2e1 / 0.3e1 * t125;
t88 = -0.2e1 / 0.3e1 * t128;
t247 = t78 + t88;
t246 = 0.4e1 / 0.7e1 * t139 - t125 / 0.7e1;
t245 = (pkin(1) - pkin(3)) * (pkin(1) + pkin(3));
t244 = t131 * t62 ^ 2;
t243 = t133 * t62;
t65 = t98 ^ 2;
t242 = t137 * t65;
t225 = t125 + t128;
t241 = t139 * t225;
t143 = pkin(2) * t137;
t240 = t143 * t98 * t65;
t239 = t144 * t65 ^ 2;
t61 = t97 * t62;
t238 = t147 * t61;
t121 = 3 * t139;
t44 = t121 - t225;
t237 = t44 * t133;
t134 = t143 ^ 2;
t138 = t139 ^ 2;
t223 = t131 + t138;
t231 = t139 * t125;
t79 = -0.3e1 / 0.2e1 * t125;
t236 = t134 + t56 * ((t79 + t122) * t133 - 0.3e1 / 0.2e1 * t231 + t223 + t249);
t93 = t137 / 0.2e1;
t235 = t139 + t93;
t187 = t78 + t56;
t83 = 0.2e1 / 0.3e1 * t128;
t234 = t144 + t56 * (t83 + t187);
t233 = t79 + t121;
t48 = -t133 / 0.3e1 + t139;
t232 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
t230 = t139 * t133;
t105 = 10 * t133;
t229 = t105 + t122;
t58 = -3 * t133 + t139;
t228 = t122 - t125;
t227 = t124 - t127;
t226 = -t125 + t128;
t224 = t125 - t139;
t222 = t133 - t139;
t221 = -t137 + t139;
t220 = t137 + t139;
t219 = t138 - t131;
t218 = t138 + t144;
t217 = 0.2e1 * t255;
t216 = 0.2e1 * t53;
t215 = 0.4e1 * pkin(1);
t214 = pkin(1) * t256;
t213 = pkin(1) * t53;
t212 = pkin(1) + r_base(1);
t211 = pkin(2) * t251;
t210 = 0.8e1 * t244;
t209 = -0.4e1 * t243;
t208 = 0.4e1 * t243;
t207 = 0.2e1 * t242;
t206 = 0.8e1 * t240;
t205 = pkin(1) * t238;
t203 = t147 * t54;
t202 = 0.12e2 * t243;
t201 = 0.16e2 * t238;
t103 = 15 * t131;
t104 = 15 * t133;
t106 = -0.2e1 * t125;
t107 = -0.5e1 * t125;
t108 = 7 * t131;
t109 = 5 * t131;
t110 = 7 * t133;
t112 = 3 * t133;
t115 = 3 * t137;
t118 = 3 * t138;
t119 = 8 * t139;
t120 = 4 * t139;
t116 = -3 * t137;
t180 = -t125 + t220;
t163 = t133 + t180;
t159 = -t128 + t163;
t43 = -0.2e1 * t214;
t157 = t43 + t159;
t183 = t128 + t224;
t23 = t43 + t207 + t221;
t179 = pkin(3) * t211;
t35 = -0.2e1 * t179;
t13 = sqrt(t23 * t209 + 0.4e1 * t222 * t242 + 0.4e1 * t159 * t214 - t131 + (t116 + t183) * t260 - (t139 + (pkin(2) - t257) * (pkin(2) + t257)) * (t139 + (pkin(2) - t258) * (pkin(2) + t258)) + (-(t35 + t157) * t253 + t157 * t211) * t265);
t140 = pkin(1) * t139;
t21 = t228 * t133 + t223 - t231 - t263;
t156 = t144 + t21;
t158 = -t179 + t235;
t169 = -0.4e1 * t179;
t182 = t133 + t220;
t77 = -t125 / 0.2e1;
t37 = t77 + t182;
t161 = t37 * t169;
t170 = -0.6e1 * t179;
t173 = t118 - 0.3e1 * t231 + t249;
t51 = 0.10e2 / 0.3e1 * t133;
t162 = ((t51 + t228) * t137 + t156) * t170 + (t103 + (-0.9e1 * t125 + (18 * t139)) * t133 + t173) * t137 + (t104 + t233) * t144 + t236;
t185 = t112 + t220;
t164 = t185 * t54 - t255 * (t115 + t56);
t165 = -(t114 + ((5 * t133) + t44) * t259 + t56 * (t88 + t187)) * t255 + (t144 + (t229 + t247) * t137 + t109 + 0.2e1 * t237 + t139 * (t139 + t247)) * t54;
t166 = (6 * t138) - 0.6e1 * t231 + t227;
t192 = t139 + t264;
t168 = t133 + t192;
t197 = t139 + t248;
t171 = t93 + t197;
t198 = t122 + t78 + t83;
t84 = t128 / 0.3e1;
t28 = t76 + t84 + t182;
t172 = t28 * t169 + t234 + (t111 + t198) * t137;
t174 = t61 * t203;
t175 = t244 * t54;
t176 = pkin(1) * t201;
t177 = 0.20e2 / 0.3e1 * t133;
t178 = 0.8e1 * t205;
t184 = t121 + t226;
t188 = 0.6e1 * t213;
t189 = 0.4e1 * t213;
t191 = t95 * t232;
t193 = t125 / 0.3e1 + t84 + t122;
t194 = 0.2e1 / 0.3e1 * t125 + t83 + t120;
t195 = 0.4e1 / 0.3e1 * t125 + 0.4e1 / 0.3e1 * t128 - (2 * t139);
t196 = t139 + t77 - t128 / 0.2e1;
t199 = -t255 / 0.2e1;
t49 = t139 - t137 / 0.3e1;
t24 = t49 * t35;
t29 = t40 + t53;
t59 = t116 + t139;
t32 = t59 * t178;
t34 = 0.10e2 * t237;
t38 = t128 + t181;
t41 = pkin(1) + t53;
t42 = 0.2e1 * t213;
t46 = (t120 + t125) * t133;
t52 = -0.30e2 * t125 + (60 * t139);
t75 = -t125 / 0.4e1;
t89 = 0.4e1 / 0.3e1 * t133;
t90 = t133 / 0.3e1;
t91 = t133 / 0.2e1;
t200 = ((-0.24e2 * (0.4e1 / 0.3e1 * t243 + t42 + t48) * t239 * t255 - 0.12e2 * (-0.8e1 / 0.3e1 * t175 + ((t89 + t171) * t54 - (0.4e1 / 0.3e1 * t137 + t91 + t197) * t255) * t208 + (-(t137 * t222) - 0.5e1 / 0.3e1 * t131 + t193 * t133 + t139 * t192) * t54 + (-t144 + (-t177 + t194) * t137 - (3 * t131) + t195 * t133 + t138) * t199 + (-t95 * t131 * t61 + ((t133 + t171) * t54 + (t259 - t222) * t199) * t53) * t215) * t242 + 0.24e2 * t49 * t175 + ((t115 + 0.3e1 / 0.2e1 * t133 + t196) * t54 + t59 * t255 / 0.2e1) * t176 - 0.6e1 * ((-(3 * t144) + (-t177 + t195) * t137 + t194 * t133 + t219) * t54 - 0.2e1 * (-0.5e1 / 0.3e1 * t144 + (-t133 + t193) * t137 + t139 * t168) * t255) * t243 - 0.6e1 * t165 * t213 - (t134 + ((21 * t133) + t44) * t144 + (t118 + (35 * t131) + t34 - 0.2e1 * t241) * t137 + (t108 + (t107 + t119 - 0.5e1 * t128) * t133 - t139 * t183) * t56) * t54 + (0.7e1 * t134 + (t110 + t44) * t114 + ((21 * t131) + (9 * t138) + t34 - 0.6e1 * t241) * t137 + t252) * t255) * t13 + t29 * (0.16e2 * (t210 + t176 + (-8 * t131 + 12 * t230) * t62 + (-0.12e2 * pkin(1) * t147 + t140 * t265) * t97 - (6 * t230) + t223) * t239 + 0.24e2 * ((t139 - 0.2e1 / 0.3e1 * t137) * t210 + 0.14e2 * (-0.32e2 / 0.21e2 * (t139 + t137 / 0.4e1 + t133 / 0.4e1 - t125 / 0.8e1) * t179 + t144 / 0.7e1 + (0.16e2 / 0.21e2 * t133 + t246) * t137 + t131 / 0.7e1 + t246 * t133 + t138 - 0.3e1 / 0.7e1 * t231 + t124 / 0.42e2 - t127 / 0.42e2) * t243 + t48 * t161 - (t222 * t144) + (-0.10e2 / 0.3e1 * t131 + (2 * t138) - t231 + t46) * t137 + t21 * t245 + ((-0.2e1 / 0.3e1 * t179 + t139 + t91 + t75) * t201 + 0.6e1 * (-0.8e1 / 0.3e1 * (t75 + t90 + t235) * t179 + t144 / 0.3e1 + (0.4e1 / 0.3e1 * t139 + t89 + t76) * t137 + t138 + 0.2e1 / 0.3e1 * t230 - 0.2e1 / 0.3e1 * t231 - t131 / 0.3e1 + t124 / 0.18e2 - t127 / 0.18e2) * t53) * pkin(1)) * t242 + 0.16e2 * (-6 * t139 * t137 + t218) * t244 + 0.32e2 * (t232 * t35 + t37 * t59) * t205 + 0.24e2 * (t49 * t161 - t134 + (-t51 + t224) * t144 + (t46 + t219 + t263) * t137 + t139 * t21) * t243 + 0.8e1 * t162 * t213 - 0.8e1 * ((t110 + t233) * t144 + (t108 + (t107 + (10 * t139)) * t133 + t173) * t137 + t236) * t179 + (t144 ^ 2) + (t106 + t120 + (28 * t133)) * t134 + (t133 * t52 + (70 * t131) + t166) * t144 + (t166 * t111 + t227 * t122 - 0.6e1 * t138 * t125 + t52 * t131 + 0.4e1 * t140 ^ 2 + (28 * t147 ^ 2)) * t137 + t38 * t252) + (((0.4e1 * t254 + (t54 + t217) * t42 + t58 * t54 + (0.3e1 / 0.2e1 * t137 + t112 + t196) * t217) * t206 + 0.6e1 * ((0.2e1 * (t91 + t137 + t248) * t54 + pkin(3) * t191) * t209 + (-0.8e1 * t174 + 0.4e1 * ((t185 + t264) * t54 - (t115 + t168) * t255) * t53) * pkin(1) + t165) * t256) * t13 + t29 * (-0.32e2 * (t178 + (-0.4e1 * t95 * t203 + t261 + ((4 * t137) + t106 + t119) * t133) * t62 + (-t133 + t158 + t75) * t189 + t35 * t245 + t37 * t58) * t240 - 0.8e1 * (t32 + (t232 * t37 + t24) * t202 + (t161 + (t111 + t228) * t137 + t156) * t188 + t162) * t256)) * t41) / ((-0.4e1 * (0.2e1 * t254 + (t260 + t137) * t255 + (-t222 + t42) * t54) * t242 + 0.8e1 * pkin(1) * t174 + ((pkin(2) * t261 + 0.8e1 * t133 * t143) * t96 + 0.4e1 * t147 * t191) * t62 - 0.4e1 * t164 * t213 - (t229 * t137 + t109 + t218 + 6 * t230) * t54 + (t114 + (t105 + 6 * t139) * t137 + t50) * t255) * t13 + t29 * (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t179 + t139 + t137 / 0.3e1 + t90 + t128 / 0.9e1 - t125 / 0.9e1) * t243 + t48 * t35 + t28 * t245 + (t238 + (t128 / 0.6e1 + t74 + t158) * t53) * t215) * t242 + t32 + (t232 * t28 + t24) * t202 + t172 * t188 + ((t51 + t198) * t137 + t234) * t170 + t134 + (t104 + t184) * t144 + (t184 * t111 + t226 * t122 + t103 + t118) * t137 + t50 * t38) + ((t206 * t255 + 0.4e1 * (t204 * t262 + (t54 - t255) * t42 + t164) * t256) * t13 + t29 * (-0.8e1 * (t42 + t208 + t58) * t240 - 0.6e1 * (t221 * t208 + (t35 + t28) * t189 + t172) * t256)) * t41) / 0.4e1;
t16 = 0.1e1 / (t216 * t40 + t182 + t35 + t43);
t190 = 0.1e1 / pkin(5) / pkin(4) ^ 2 * t16;
t186 = 0.1e1 / pkin(4) * t16 / 0.2e1;
t167 = -pkin(3) + t211;
t160 = t112 + t128 + t180;
t27 = t97 * t98 + t251;
t26 = t95 * t98 - t250;
t22 = t128 + t43 + t163;
t12 = (-t167 + t253) * t13 + (t216 * t23 + t22 * t40) * t95 + (t22 * t97 + (0.4e1 * t62 - 0.2e1) * t40 * pkin(3)) * t54;
t11 = (pkin(2) * t250 + t40 * t95) * t13 - (t160 + t43 + t169) * t253 + t167 * t43 + t160 * t211 + (t23 * t262 - t115 + t207 - t38) * pkin(3);
t10 = qJ(2) + atan2(t12 * t186, t11 * t186);
t9 = cos(t10);
t8 = sin(t10);
t5 = (t11 * t13 / 0.4e1 + t12 * t200) * t190;
t4 = (t11 * t200 - t12 * t13 / 0.4e1) * t190;
t3 = qJ(1) + atan2(t26 * t4 + t27 * t5, -t26 * t5 + t27 * t4);
t2 = cos(t3);
t1 = sin(t3);
t6 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t98 - rSges(2,2) * t96 + r_base(1)) + g(2) * (rSges(2,1) * t96 + rSges(2,2) * t98 + r_base(2)) + g(3) * (r_base(3) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,1) * t2 + rSges(3,2) * t1 + t256 + r_base(1)) + g(2) * (-rSges(3,1) * t1 - rSges(3,2) * t2 + t54 + r_base(2)) + g(3) * (r_base(3) + rSges(3,3))) - m(4) * (g(1) * (rSges(4,1) * t97 - rSges(4,2) * t95 + t212) + g(2) * (rSges(4,1) * t95 + rSges(4,2) * t97 + r_base(2)) + g(3) * (r_base(3) + rSges(4,3))) - m(5) * (g(1) * (rSges(5,1) * t9 - rSges(5,2) * t8 + t212 + t53) + g(2) * (rSges(5,1) * t8 + rSges(5,2) * t9 + t255 + r_base(2)) + g(3) * (r_base(3) + rSges(5,3)));
U = t6;
