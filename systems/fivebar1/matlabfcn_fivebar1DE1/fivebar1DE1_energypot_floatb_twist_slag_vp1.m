% Calculate potential energy for
% fivebar1DE1
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
% Datum: 2020-04-27 04:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fivebar1DE1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fivebar1DE1_energypot_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'fivebar1DE1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1DE1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1DE1_energypot_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1DE1_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fivebar1DE1_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 02:23:13
% EndTime: 2020-04-27 02:23:23
% DurationCPUTime: 4.83s
% Computational Cost: add. (28347->427), mult. (83144->505), div. (256->5), fcn. (9996->12), ass. (0->208)
t269 = 4 * pkin(3);
t132 = pkin(4) ^ 2;
t129 = pkin(5) ^ 2;
t80 = -t129 / 0.3e1;
t268 = t80 - t132 / 0.3e1;
t128 = t129 ^ 2;
t131 = t132 ^ 2;
t267 = -t128 / 0.6e1 + t131 / 0.6e1;
t101 = cos(qJ(2));
t66 = t101 ^ 2;
t266 = -0.2e1 * t66;
t137 = pkin(3) ^ 2;
t135 = t137 ^ 2;
t265 = 4 * t135;
t264 = 2 * t137;
t115 = 6 * t137;
t141 = (pkin(2) ^ 2);
t263 = 2 * t141;
t143 = (pkin(1) ^ 2);
t126 = 2 * t143;
t148 = t141 ^ 2;
t118 = 5 * t148;
t262 = -pkin(4) - pkin(5);
t261 = -pkin(4) + pkin(5);
t99 = sin(qJ(2));
t260 = t99 * pkin(3);
t100 = sin(qJ(1));
t58 = pkin(2) * t100;
t102 = cos(qJ(1));
t259 = pkin(2) * t102;
t57 = pkin(3) * t101;
t151 = pkin(3) * t137;
t205 = t137 * t58;
t258 = (-t151 * t99 + t205) * t66;
t60 = t137 + t143;
t185 = -t129 + t60;
t54 = t60 ^ 2;
t257 = t54 * (-t132 + t185);
t256 = t128 / 0.2e1 - t131 / 0.2e1;
t78 = -t129 / 0.6e1;
t255 = t78 - t132 / 0.6e1;
t82 = -0.2e1 / 0.3e1 * t129;
t92 = -0.2e1 / 0.3e1 * t132;
t254 = t82 + t92;
t253 = 0.4e1 / 0.7e1 * t143 - t129 / 0.7e1;
t252 = t100 * t99;
t44 = pkin(1) - t259;
t251 = t101 * t44;
t250 = t135 * t66 ^ 2;
t249 = t137 * t66;
t69 = t102 ^ 2;
t248 = t141 * t69;
t229 = t129 + t132;
t247 = t143 * t229;
t147 = pkin(2) * t141;
t246 = t147 * t102 * t69;
t245 = t148 * t69 ^ 2;
t65 = t101 * t66;
t244 = t151 * t65;
t125 = 3 * t143;
t48 = t125 - t229;
t243 = t48 * t137;
t138 = t147 ^ 2;
t142 = t143 ^ 2;
t227 = t135 + t142;
t235 = t143 * t129;
t83 = -0.3e1 / 0.2e1 * t129;
t242 = t138 + t60 * ((t83 + t126) * t137 - 0.3e1 / 0.2e1 * t235 + t227 + t256);
t97 = t141 / 0.2e1;
t241 = t143 + t97;
t193 = t82 + t60;
t87 = 0.2e1 / 0.3e1 * t132;
t240 = t148 + t60 * (t87 + t193);
t239 = t83 + t125;
t52 = -t137 / 0.3e1 + t143;
t238 = t100 * t101;
t237 = (pkin(1) + pkin(3)) * (pkin(1) - pkin(3));
t236 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
t234 = t143 * t137;
t109 = 10 * t137;
t233 = t109 + t126;
t62 = -3 * t137 + t143;
t232 = t126 - t129;
t231 = t128 - t131;
t230 = -t129 + t132;
t228 = t129 - t143;
t226 = t137 - t143;
t225 = -t141 + t143;
t224 = t141 + t143;
t223 = t142 - t135;
t222 = t142 + t148;
t221 = 0.2e1 * t260;
t220 = 0.4e1 * pkin(1);
t219 = 0.2e1 * t57;
t218 = pkin(1) + r_base(1);
t217 = pkin(1) * t259;
t216 = pkin(1) * t57;
t215 = 0.8e1 * t250;
t214 = -0.4e1 * t249;
t213 = 0.4e1 * t249;
t212 = 0.2e1 * t248;
t211 = 0.8e1 * t246;
t210 = pkin(1) * t244;
t209 = pkin(2) * t252;
t208 = 0.12e2 * t249;
t207 = 0.16e2 * t244;
t206 = t151 * t58;
t107 = 15 * t135;
t108 = 15 * t137;
t110 = -0.2e1 * t129;
t111 = -0.5e1 * t129;
t112 = 7 * t135;
t113 = 5 * t135;
t114 = 7 * t137;
t116 = 3 * t137;
t119 = 3 * t141;
t122 = 3 * t142;
t123 = 8 * t143;
t124 = 4 * t143;
t144 = pkin(1) * t143;
t25 = t232 * t137 + t227 - t235 - t267;
t160 = t148 + t25;
t181 = pkin(3) * t209;
t162 = -t181 + t241;
t172 = -0.4e1 * t181;
t186 = t137 + t224;
t81 = -t129 / 0.2e1;
t41 = t81 + t186;
t165 = t41 * t172;
t173 = -0.6e1 * t181;
t177 = t122 - 0.3e1 * t235 + t256;
t55 = 0.10e2 / 0.3e1 * t137;
t166 = ((t55 + t232) * t141 + t160) * t173 + (t107 + (-0.9e1 * t129 + (18 * t143)) * t137 + t177) * t141 + (t108 + t239) * t148 + t242;
t189 = t116 + t224;
t168 = t189 * t58 - t260 * (t119 + t60);
t169 = -(t118 + ((5 * t137) + t48) * t263 + t60 * (t92 + t193)) * t260 + (t148 + (t233 + t254) * t141 + t113 + 0.2e1 * t243 + t143 * (t143 + t254)) * t58;
t120 = -3 * t141;
t184 = -t129 + t224;
t167 = t137 + t184;
t163 = -t132 + t167;
t47 = -0.2e1 * t217;
t161 = t47 + t163;
t187 = t132 + t228;
t27 = t47 + t212 + t225;
t39 = -0.2e1 * t181;
t17 = sqrt(t27 * t214 + 0.4e1 * t226 * t248 + 0.4e1 * t163 * t217 - t135 + (t120 + t187) * t264 - (t143 + (pkin(2) - t261) * (pkin(2) + t261)) * (t143 + (pkin(2) - t262) * (pkin(2) + t262)) + (-(t39 + t161) * t251 + t161 * t209) * t269);
t171 = (6 * t142) - 0.6e1 * t235 + t231;
t196 = t143 + t268;
t174 = t137 + t196;
t201 = t143 + t255;
t175 = t97 + t201;
t202 = t126 + t82 + t87;
t88 = t132 / 0.3e1;
t32 = t80 + t88 + t186;
t176 = t32 * t172 + t240 + (t115 + t202) * t141;
t178 = t65 * t206;
t179 = t250 * t58;
t180 = pkin(1) * t207;
t182 = 0.20e2 / 0.3e1 * t137;
t183 = 0.8e1 * t210;
t188 = t125 + t230;
t190 = 0.6e1 * t216;
t191 = 0.4e1 * t216;
t195 = t99 * t236;
t197 = t129 / 0.3e1 + t88 + t126;
t198 = 0.2e1 / 0.3e1 * t129 + t87 + t124;
t199 = 0.4e1 / 0.3e1 * t129 + 0.4e1 / 0.3e1 * t132 - (2 * t143);
t200 = t143 + t81 - t132 / 0.2e1;
t203 = -t260 / 0.2e1;
t53 = t143 - t141 / 0.3e1;
t28 = t53 * t39;
t33 = t44 + t57;
t63 = t120 + t143;
t36 = t63 * t183;
t38 = 0.10e2 * t243;
t42 = t132 + t185;
t45 = pkin(1) + t57;
t46 = 0.2e1 * t216;
t50 = (t124 + t129) * t137;
t56 = -0.30e2 * t129 + (60 * t143);
t79 = -t129 / 0.4e1;
t93 = 0.4e1 / 0.3e1 * t137;
t94 = t137 / 0.3e1;
t95 = t137 / 0.2e1;
t204 = 0.1e1 / ((-0.4e1 * (0.2e1 * t258 + (t264 + t141) * t260 + (-t226 + t46) * t58) * t248 + 0.8e1 * pkin(1) * t178 + ((pkin(2) * t265 + 0.8e1 * t137 * t147) * t100 + 0.4e1 * t151 * t195) * t66 - 0.4e1 * t168 * t216 - (t141 * t233 + t113 + t222 + 6 * t234) * t58 + (t118 + (t109 + 6 * t143) * t141 + t54) * t260) * t17 + t33 * (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t181 + t143 + t141 / 0.3e1 + t94 + t132 / 0.9e1 - t129 / 0.9e1) * t249 + t52 * t39 + t32 * t237 + (t244 + (t132 / 0.6e1 + t78 + t162) * t57) * t220) * t248 + t36 + (t236 * t32 + t28) * t208 + t176 * t190 + ((t55 + t202) * t141 + t240) * t173 + t138 + (t108 + t188) * t148 + (t188 * t115 + t230 * t126 + t107 + t122) * t141 + t54 * t42) + ((t211 * t260 + 0.4e1 * (t205 * t266 + (t58 - t260) * t46 + t168) * t259) * t17 + t33 * (-0.8e1 * (t46 + t213 + t62) * t246 - 0.6e1 * (t225 * t213 + (t39 + t32) * t191 + t176) * t259)) * t45) * ((-0.24e2 * (0.4e1 / 0.3e1 * t249 + t46 + t52) * t245 * t260 - 0.12e2 * (-0.8e1 / 0.3e1 * t179 + ((t93 + t175) * t58 - (0.4e1 / 0.3e1 * t141 + t95 + t201) * t260) * t213 + (-(t141 * t226) - 0.5e1 / 0.3e1 * t135 + t197 * t137 + t143 * t196) * t58 + (-t148 + (-t182 + t198) * t141 - (3 * t135) + t199 * t137 + t142) * t203 + (-t99 * t135 * t65 + ((t137 + t175) * t58 + (t263 - t226) * t203) * t57) * t220) * t248 + 0.24e2 * t53 * t179 + ((t119 + 0.3e1 / 0.2e1 * t137 + t200) * t58 + t63 * t260 / 0.2e1) * t180 - 0.6e1 * ((-(3 * t148) + (-t182 + t199) * t141 + t198 * t137 + t223) * t58 - 0.2e1 * (-0.5e1 / 0.3e1 * t148 + (-t137 + t197) * t141 + t143 * t174) * t260) * t249 - 0.6e1 * t169 * t216 - (t138 + ((21 * t137) + t48) * t148 + (t122 + (35 * t135) + t38 - 0.2e1 * t247) * t141 + (t112 + (t111 + t123 - 0.5e1 * t132) * t137 - t143 * t187) * t60) * t58 + (0.7e1 * t138 + (t114 + t48) * t118 + ((21 * t135) + (9 * t142) + t38 - 0.6e1 * t247) * t141 + t257) * t260) * t17 + t33 * (0.16e2 * (t215 + t180 + (-8 * t135 + 12 * t234) * t66 + (-0.12e2 * pkin(1) * t151 + t144 * t269) * t101 - (6 * t234) + t227) * t245 + 0.24e2 * ((t143 - 0.2e1 / 0.3e1 * t141) * t215 + 0.14e2 * (-0.32e2 / 0.21e2 * (t143 + t141 / 0.4e1 + t137 / 0.4e1 - t129 / 0.8e1) * t181 + t148 / 0.7e1 + (0.16e2 / 0.21e2 * t137 + t253) * t141 + t135 / 0.7e1 + t253 * t137 + t142 - 0.3e1 / 0.7e1 * t235 + t128 / 0.42e2 - t131 / 0.42e2) * t249 + t52 * t165 - (t226 * t148) + (-0.10e2 / 0.3e1 * t135 + (2 * t142) - t235 + t50) * t141 + t25 * t237 + ((-0.2e1 / 0.3e1 * t181 + t143 + t95 + t79) * t207 + 0.6e1 * (-0.8e1 / 0.3e1 * (t79 + t94 + t241) * t181 + t148 / 0.3e1 + (0.4e1 / 0.3e1 * t143 + t93 + t80) * t141 + t142 + 0.2e1 / 0.3e1 * t234 - 0.2e1 / 0.3e1 * t235 - t135 / 0.3e1 + t128 / 0.18e2 - t131 / 0.18e2) * t57) * pkin(1)) * t248 + 0.16e2 * (-6 * t143 * t141 + t222) * t250 + 0.32e2 * (t236 * t39 + t41 * t63) * t210 + 0.24e2 * (t53 * t165 - t138 + (-t55 + t228) * t148 + (t50 + t223 + t267) * t141 + t143 * t25) * t249 + 0.8e1 * t166 * t216 - 0.8e1 * ((t114 + t239) * t148 + (t112 + (t111 + (10 * t143)) * t137 + t177) * t141 + t242) * t181 + (t148 ^ 2) + (t110 + t124 + (28 * t137)) * t138 + (t137 * t56 + (70 * t135) + t171) * t148 + (t171 * t115 + t231 * t126 - 0.6e1 * t142 * t129 + t56 * t135 + 0.4e1 * t144 ^ 2 + (28 * t151 ^ 2)) * t141 + t42 * t257) + (((0.4e1 * t258 + (t58 + t221) * t46 + t62 * t58 + (0.3e1 / 0.2e1 * t141 + t116 + t200) * t221) * t211 + 0.6e1 * ((0.2e1 * (t95 + t141 + t255) * t58 + pkin(3) * t195) * t214 + (-0.8e1 * t178 + 0.4e1 * ((t189 + t268) * t58 - (t119 + t174) * t260) * t57) * pkin(1) + t169) * t259) * t17 + t33 * (-0.32e2 * (t183 + (-0.4e1 * t99 * t206 + t265 + ((4 * t141) + t110 + t123) * t137) * t66 + (-t137 + t162 + t79) * t191 + t39 * t237 + t41 * t62) * t246 - 0.8e1 * (t36 + (t236 * t41 + t28) * t208 + (t165 + (t115 + t232) * t141 + t160) * t190 + t166) * t259)) * t45) / 0.4e1;
t20 = 0.1e1 / (t219 * t44 + t186 + t39 + t47);
t194 = 0.1e1 / pkin(5) / pkin(4) ^ 2 * t20;
t192 = 0.1e1 / pkin(4) * t20 / 0.2e1;
t170 = -pkin(3) + t209;
t164 = t116 + t132 + t184;
t31 = t101 * t102 + t252;
t30 = t102 * t99 - t238;
t26 = t132 + t47 + t167;
t16 = (-t170 + t251) * t17 + (t219 * t27 + t26 * t44) * t99 + (t101 * t26 + (0.4e1 * t66 - 0.2e1) * t44 * pkin(3)) * t58;
t15 = (pkin(2) * t238 + t44 * t99) * t17 - (t164 + t47 + t172) * t251 + t170 * t47 + t164 * t209 + (t27 * t266 - t119 + t212 - t42) * pkin(3);
t14 = atan2(t16 * t192, t15 * t192);
t13 = cos(t14);
t12 = sin(t14);
t10 = t101 * t13 - t12 * t99;
t9 = t101 * t12 + t13 * t99;
t7 = (t15 * t17 / 0.4e1 + t16 * t204) * t194;
t6 = (t15 * t204 - t16 * t17 / 0.4e1) * t194;
t5 = atan2(t30 * t6 + t31 * t7, -t30 * t7 + t31 * t6);
t4 = cos(t5);
t3 = sin(t5);
t2 = t100 * t3 - t102 * t4;
t1 = -t100 * t4 - t102 * t3;
t8 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t102 - rSges(2,2) * t100 + r_base(1)) + g(2) * (rSges(2,1) * t100 + rSges(2,2) * t102 + r_base(2)) + g(3) * (r_base(3) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t2 - rSges(3,2) * t1 + t259 + r_base(1)) + g(2) * (rSges(3,1) * t1 + rSges(3,2) * t2 + t58 + r_base(2)) + g(3) * (r_base(3) + rSges(3,3))) - m(4) * (g(1) * (rSges(4,1) * t101 - rSges(4,2) * t99 + t218) + g(2) * (rSges(4,1) * t99 + rSges(4,2) * t101 + r_base(2)) + g(3) * (r_base(3) + rSges(4,3))) - m(5) * (g(1) * (rSges(5,1) * t10 - rSges(5,2) * t9 + t218 + t57) + g(2) * (rSges(5,1) * t9 + rSges(5,2) * t10 + t260 + r_base(2)) + g(3) * (r_base(3) + rSges(5,3)));
U = t8;
