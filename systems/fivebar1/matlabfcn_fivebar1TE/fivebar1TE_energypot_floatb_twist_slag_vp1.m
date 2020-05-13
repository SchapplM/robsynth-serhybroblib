% Calculate potential energy for
% fivebar1TE
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
% Datum: 2020-04-27 10:28
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fivebar1TE_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fivebar1TE_energypot_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'fivebar1TE_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1TE_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1TE_energypot_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1TE_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fivebar1TE_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 09:13:18
% EndTime: 2020-04-27 09:13:23
% DurationCPUTime: 3.82s
% Computational Cost: add. (14199->427), mult. (41600->507), div. (128->5), fcn. (4996->6), ass. (0->205)
t267 = 4 * pkin(3);
t127 = pkin(4) ^ 2;
t124 = pkin(5) ^ 2;
t75 = -t124 / 0.3e1;
t266 = t75 - t127 / 0.3e1;
t123 = t124 ^ 2;
t126 = t127 ^ 2;
t265 = -t123 / 0.6e1 + t126 / 0.6e1;
t96 = cos(qJ(2));
t61 = t96 ^ 2;
t264 = -0.2e1 * t61;
t132 = pkin(3) ^ 2;
t130 = t132 ^ 2;
t263 = 4 * t130;
t262 = 2 * t132;
t110 = 6 * t132;
t136 = (pkin(2) ^ 2);
t261 = 2 * t136;
t138 = (pkin(1) ^ 2);
t121 = 2 * t138;
t143 = t136 ^ 2;
t113 = 5 * t143;
t260 = t96 / 0.2e1;
t259 = -pkin(4) - pkin(5);
t258 = -pkin(4) + pkin(5);
t95 = sin(qJ(1));
t53 = pkin(2) * t95;
t97 = cos(qJ(1));
t257 = pkin(2) * t97;
t52 = pkin(3) * t96;
t94 = sin(qJ(2));
t256 = t94 * pkin(3);
t146 = pkin(3) * t132;
t204 = t132 * t53;
t255 = (-t146 * t94 + t204) * t61;
t39 = pkin(1) - t257;
t254 = t39 * t96;
t55 = t132 + t138;
t181 = -t124 + t55;
t49 = t55 ^ 2;
t253 = t49 * (-t127 + t181);
t252 = t94 * t95;
t251 = t95 * t96;
t250 = (pkin(1) + pkin(3)) * (pkin(1) - pkin(3));
t249 = t123 / 0.2e1 - t126 / 0.2e1;
t73 = -t124 / 0.6e1;
t248 = t73 - t127 / 0.6e1;
t77 = -0.2e1 / 0.3e1 * t124;
t87 = -0.2e1 / 0.3e1 * t127;
t247 = t77 + t87;
t246 = 0.4e1 / 0.7e1 * t138 - t124 / 0.7e1;
t220 = t136 + t138;
t182 = t132 + t220;
t217 = 0.2e1 * t52;
t211 = pkin(2) * t252;
t179 = pkin(3) * t211;
t34 = -0.2e1 * t179;
t214 = pkin(1) * t257;
t42 = -0.2e1 * t214;
t15 = 0.1e1 / (t217 * t39 + t182 + t34 + t42);
t245 = 0.1e1 / pkin(4) * t15;
t244 = t130 * t61 ^ 2;
t243 = t132 * t61;
t64 = t97 ^ 2;
t242 = t136 * t64;
t225 = t124 + t127;
t241 = t138 * t225;
t142 = pkin(2) * t136;
t240 = t142 * t97 * t64;
t239 = t143 * t64 ^ 2;
t60 = t96 * t61;
t238 = t146 * t60;
t120 = 3 * t138;
t43 = t120 - t225;
t237 = t43 * t132;
t133 = t142 ^ 2;
t137 = t138 ^ 2;
t223 = t130 + t137;
t231 = t138 * t124;
t78 = -0.3e1 / 0.2e1 * t124;
t236 = t133 + t55 * ((t78 + t121) * t132 - 0.3e1 / 0.2e1 * t231 + t223 + t249);
t92 = t136 / 0.2e1;
t235 = t138 + t92;
t187 = t77 + t55;
t82 = 0.2e1 / 0.3e1 * t127;
t234 = t143 + (t82 + t187) * t55;
t233 = t78 + t120;
t47 = -t132 / 0.3e1 + t138;
t232 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
t230 = t138 * t132;
t104 = 10 * t132;
t229 = t104 + t121;
t57 = -3 * t132 + t138;
t228 = t121 - t124;
t227 = t123 - t126;
t226 = -t124 + t127;
t224 = t124 - t138;
t222 = t132 - t138;
t221 = -t136 + t138;
t219 = t137 - t130;
t218 = t137 + t143;
t216 = 0.2e1 * t256;
t215 = 0.4e1 * pkin(1);
t213 = pkin(1) * t52;
t212 = pkin(1) + r_base(1);
t210 = 0.8e1 * t244;
t209 = -0.4e1 * t243;
t208 = 0.4e1 * t243;
t207 = 0.2e1 * t242;
t206 = 0.8e1 * t240;
t205 = pkin(1) * t238;
t203 = t146 * t53;
t202 = 0.12e2 * t243;
t201 = 0.16e2 * t238;
t102 = 15 * t130;
t103 = 15 * t132;
t105 = -0.2e1 * t124;
t106 = -0.5e1 * t124;
t107 = 7 * t130;
t108 = 5 * t130;
t109 = 7 * t132;
t111 = 3 * t132;
t114 = 3 * t136;
t117 = 3 * t137;
t118 = 8 * t138;
t119 = 4 * t138;
t115 = -3 * t136;
t180 = -t124 + t220;
t163 = t132 + t180;
t158 = -t127 + t163;
t156 = t42 + t158;
t183 = t127 + t224;
t22 = t42 + t207 + t221;
t12 = sqrt(t22 * t209 + 0.4e1 * t222 * t242 + 0.4e1 * t158 * t214 - t130 + (t115 + t183) * t262 - (t138 + (pkin(2) - t259) * (pkin(2) + t259)) * (t138 + (pkin(2) - t258) * (pkin(2) + t258)) + (-(t34 + t156) * t254 + t156 * t211) * t267);
t139 = pkin(1) * t138;
t20 = t228 * t132 + t223 - t231 - t265;
t155 = t143 + t20;
t157 = -t179 + t235;
t169 = -0.4e1 * t179;
t76 = -t124 / 0.2e1;
t36 = t76 + t182;
t160 = t36 * t169;
t170 = -0.6e1 * t179;
t173 = t117 - 0.3e1 * t231 + t249;
t50 = 0.10e2 / 0.3e1 * t132;
t161 = ((t50 + t228) * t136 + t155) * t170 + (t102 + (-0.9e1 * t124 + (18 * t138)) * t132 + t173) * t136 + (t103 + t233) * t143 + t236;
t185 = t111 + t220;
t164 = t185 * t53 - t256 * (t114 + t55);
t165 = -(t113 + ((5 * t132) + t43) * t261 + (t87 + t187) * t55) * t256 + (t143 + (t229 + t247) * t136 + t108 + 0.2e1 * t237 + t138 * (t138 + t247)) * t53;
t166 = (6 * t137) - 0.6e1 * t231 + t227;
t192 = t138 + t266;
t168 = t132 + t192;
t197 = t138 + t248;
t171 = t92 + t197;
t198 = t121 + t77 + t82;
t83 = t127 / 0.3e1;
t27 = t75 + t83 + t182;
t172 = t27 * t169 + t234 + (t110 + t198) * t136;
t174 = t60 * t203;
t175 = t244 * t53;
t176 = pkin(1) * t201;
t177 = 0.20e2 / 0.3e1 * t132;
t178 = 0.8e1 * t205;
t184 = t120 + t226;
t188 = 0.6e1 * t213;
t189 = 0.4e1 * t213;
t191 = t94 * t232;
t193 = t124 / 0.3e1 + t83 + t121;
t194 = 0.2e1 / 0.3e1 * t124 + t82 + t119;
t195 = 0.4e1 / 0.3e1 * t124 + 0.4e1 / 0.3e1 * t127 - (2 * t138);
t196 = t138 + t76 - t127 / 0.2e1;
t199 = -t256 / 0.2e1;
t48 = t138 - t136 / 0.3e1;
t23 = t48 * t34;
t28 = t39 + t52;
t58 = t115 + t138;
t31 = t58 * t178;
t33 = 0.10e2 * t237;
t37 = t127 + t181;
t40 = pkin(1) + t52;
t41 = 0.2e1 * t213;
t45 = (t119 + t124) * t132;
t51 = -0.30e2 * t124 + (60 * t138);
t74 = -t124 / 0.4e1;
t88 = 0.4e1 / 0.3e1 * t132;
t89 = t132 / 0.3e1;
t90 = t132 / 0.2e1;
t200 = ((-0.24e2 * (0.4e1 / 0.3e1 * t243 + t41 + t47) * t239 * t256 - 0.12e2 * (-0.8e1 / 0.3e1 * t175 + ((t88 + t171) * t53 - (0.4e1 / 0.3e1 * t136 + t90 + t197) * t256) * t208 + (-(t136 * t222) - 0.5e1 / 0.3e1 * t130 + t193 * t132 + t138 * t192) * t53 + (-t143 + (-t177 + t194) * t136 - (3 * t130) + t195 * t132 + t137) * t199 + (-t94 * t130 * t60 + ((t132 + t171) * t53 + (t261 - t222) * t199) * t52) * t215) * t242 + 0.24e2 * t48 * t175 + ((t114 + 0.3e1 / 0.2e1 * t132 + t196) * t53 + t58 * t256 / 0.2e1) * t176 - 0.6e1 * ((-(3 * t143) + (-t177 + t195) * t136 + t194 * t132 + t219) * t53 - 0.2e1 * (-0.5e1 / 0.3e1 * t143 + (-t132 + t193) * t136 + t138 * t168) * t256) * t243 - 0.6e1 * t165 * t213 - (t133 + ((21 * t132) + t43) * t143 + (t117 + (35 * t130) + t33 - 0.2e1 * t241) * t136 + (t107 + (t106 + t118 - 0.5e1 * t127) * t132 - t138 * t183) * t55) * t53 + (0.7e1 * t133 + (t109 + t43) * t113 + ((21 * t130) + (9 * t137) + t33 - 0.6e1 * t241) * t136 + t253) * t256) * t12 + (0.16e2 * (t210 + t176 + (-8 * t130 + 12 * t230) * t61 + (-0.12e2 * pkin(1) * t146 + t139 * t267) * t96 - (6 * t230) + t223) * t239 + 0.24e2 * ((t138 - 0.2e1 / 0.3e1 * t136) * t210 + 0.14e2 * (-0.32e2 / 0.21e2 * (t138 + t136 / 0.4e1 + t132 / 0.4e1 - t124 / 0.8e1) * t179 + t143 / 0.7e1 + (0.16e2 / 0.21e2 * t132 + t246) * t136 + t130 / 0.7e1 + t246 * t132 + t137 - 0.3e1 / 0.7e1 * t231 + t123 / 0.42e2 - t126 / 0.42e2) * t243 + t47 * t160 - (t222 * t143) + (-0.10e2 / 0.3e1 * t130 + (2 * t137) - t231 + t45) * t136 + t20 * t250 + ((-0.2e1 / 0.3e1 * t179 + t138 + t90 + t74) * t201 + 0.6e1 * (-0.8e1 / 0.3e1 * (t74 + t89 + t235) * t179 + t143 / 0.3e1 + (0.4e1 / 0.3e1 * t138 + t88 + t75) * t136 + t137 + 0.2e1 / 0.3e1 * t230 - 0.2e1 / 0.3e1 * t231 - t130 / 0.3e1 + t123 / 0.18e2 - t126 / 0.18e2) * t52) * pkin(1)) * t242 + 0.16e2 * (-6 * t138 * t136 + t218) * t244 + 0.32e2 * (t232 * t34 + t36 * t58) * t205 + 0.24e2 * (t48 * t160 - t133 + (-t50 + t224) * t143 + (t45 + t219 + t265) * t136 + t20 * t138) * t243 + 0.8e1 * t161 * t213 - 0.8e1 * ((t109 + t233) * t143 + (t107 + (t106 + (10 * t138)) * t132 + t173) * t136 + t236) * t179 + (t143 ^ 2) + (t105 + t119 + (28 * t132)) * t133 + (t132 * t51 + (70 * t130) + t166) * t143 + (t166 * t110 + t227 * t121 - 0.6e1 * t137 * t124 + t51 * t130 + 0.4e1 * t139 ^ 2 + (28 * t146 ^ 2)) * t136 + t37 * t253) * t28 + (((0.4e1 * t255 + (t53 + t216) * t41 + t57 * t53 + (0.3e1 / 0.2e1 * t136 + t111 + t196) * t216) * t206 + 0.6e1 * ((0.2e1 * (t90 + t136 + t248) * t53 + pkin(3) * t191) * t209 + (-0.8e1 * t174 + 0.4e1 * ((t185 + t266) * t53 - (t114 + t168) * t256) * t52) * pkin(1) + t165) * t257) * t12 + (-0.32e2 * (t178 + (-0.4e1 * t94 * t203 + t263 + ((4 * t136) + t105 + t118) * t132) * t61 + (-t132 + t157 + t74) * t189 + t34 * t250 + t57 * t36) * t240 - 0.8e1 * (t31 + (t232 * t36 + t23) * t202 + (t160 + (t110 + t228) * t136 + t155) * t188 + t161) * t257) * t28) * t40) / ((-0.4e1 * (0.2e1 * t255 + (t262 + t136) * t256 + (-t222 + t41) * t53) * t242 + 0.8e1 * pkin(1) * t174 + ((pkin(2) * t263 + 0.8e1 * t132 * t142) * t95 + 0.4e1 * t146 * t191) * t61 - 0.4e1 * t164 * t213 - (t136 * t229 + t108 + t218 + 6 * t230) * t53 + (t113 + (t104 + 6 * t138) * t136 + t49) * t256) * t12 + (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t179 + t138 + t136 / 0.3e1 + t89 + t127 / 0.9e1 - t124 / 0.9e1) * t243 + t47 * t34 + t27 * t250 + (t238 + (t127 / 0.6e1 + t73 + t157) * t52) * t215) * t242 + t31 + (t232 * t27 + t23) * t202 + t172 * t188 + ((t50 + t198) * t136 + t234) * t170 + t133 + (t103 + t184) * t143 + (t184 * t110 + t226 * t121 + t102 + t117) * t136 + t49 * t37) * t28 + ((t206 * t256 + 0.4e1 * (t204 * t264 + (t53 - t256) * t41 + t164) * t257) * t12 + (-0.8e1 * (t41 + t208 + t57) * t240 - 0.6e1 * (t221 * t208 + (t34 + t27) * t189 + t172) * t257) * t28) * t40) / 0.4e1;
t190 = 0.1e1 / pkin(5) / pkin(4) ^ 2 * t15;
t25 = t94 * t97 - t251;
t26 = t96 * t97 + t252;
t159 = t111 + t127 + t180;
t167 = -pkin(3) + t211;
t10 = (pkin(2) * t251 + t39 * t94) * t12 - (t159 + t42 + t169) * t254 + t167 * t42 + t159 * t211 + (t22 * t264 - t114 + t207 - t37) * pkin(3);
t21 = t127 + t42 + t163;
t11 = (-t167 + t254) * t12 + (t21 * t39 + t217 * t22) * t94 + (t21 * t96 + (0.4e1 * t61 - 0.2e1) * t39 * pkin(3)) * t53;
t4 = (t10 * t200 - t11 * t12 / 0.4e1) * t190;
t5 = (t10 * t12 / 0.4e1 + t11 * t200) * t190;
t2 = t25 * t4 + t26 * t5;
t3 = t25 * t5 - t26 * t4;
t186 = t2 * t95 + t97 * t3;
t162 = t2 * t97 - t3 * t95;
t9 = (t10 * t260 - t94 * t11 / 0.2e1) * t245;
t8 = (t94 * t10 / 0.2e1 + t11 * t260) * t245;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t97 - rSges(2,2) * t95 + r_base(1)) + g(2) * (rSges(2,1) * t95 + rSges(2,2) * t97 + r_base(2)) + g(3) * (r_base(3) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t186 + rSges(3,2) * t162 + t257 + r_base(1)) + g(2) * (-rSges(3,1) * t162 + rSges(3,2) * t186 + t53 + r_base(2)) + g(3) * (r_base(3) + rSges(3,3))) - m(4) * (g(1) * (rSges(4,1) * t96 - rSges(4,2) * t94 + t212) + g(2) * (rSges(4,1) * t94 + rSges(4,2) * t96 + r_base(2)) + g(3) * (r_base(3) + rSges(4,3))) - m(5) * (g(1) * (rSges(5,1) * t9 - rSges(5,2) * t8 + t212 + t52) + g(2) * (rSges(5,1) * t8 + rSges(5,2) * t9 + t256 + r_base(2)) + g(3) * (r_base(3) + rSges(5,3)));
U = t1;
