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
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 10:28
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fivebar1TE_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fivebar1TE_energypot_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'fivebar1TE_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1TE_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1TE_energypot_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1TE_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fivebar1TE_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 09:13:18
% EndTime: 2020-04-27 09:13:22
% DurationCPUTime: 3.81s
% Computational Cost: add. (14199->418), mult. (41610->500), div. (128->5), fcn. (4996->6), ass. (0->204)
t265 = 4 * pkin(3);
t127 = pkin(4) ^ 2;
t124 = pkin(5) ^ 2;
t75 = -t124 / 0.3e1;
t264 = t75 - t127 / 0.3e1;
t123 = t124 ^ 2;
t126 = t127 ^ 2;
t263 = -t123 / 0.6e1 + t126 / 0.6e1;
t96 = cos(qJ(2));
t61 = t96 ^ 2;
t262 = -0.2e1 * t61;
t132 = pkin(3) ^ 2;
t130 = t132 ^ 2;
t261 = 4 * t130;
t260 = 2 * t132;
t110 = 6 * t132;
t136 = (pkin(2) ^ 2);
t259 = 2 * t136;
t138 = (pkin(1) ^ 2);
t121 = 2 * t138;
t143 = t136 ^ 2;
t113 = 5 * t143;
t258 = t96 / 0.2e1;
t257 = -pkin(4) - pkin(5);
t256 = pkin(5) - pkin(4);
t95 = sin(qJ(1));
t53 = pkin(2) * t95;
t97 = cos(qJ(1));
t255 = pkin(2) * t97;
t52 = pkin(3) * t96;
t94 = sin(qJ(2));
t254 = t94 * pkin(3);
t40 = pkin(1) + t52;
t146 = pkin(3) * t132;
t203 = t132 * t53;
t253 = (-t146 * t94 + t203) * t61;
t39 = pkin(1) - t255;
t252 = t39 * t96;
t55 = t132 + t138;
t181 = -t124 + t55;
t49 = t55 ^ 2;
t251 = t49 * (-t127 + t181);
t250 = t94 * t95;
t249 = t95 * t96;
t248 = (pkin(1) + pkin(3)) * (pkin(1) - pkin(3));
t247 = t123 / 0.2e1 - t126 / 0.2e1;
t73 = -t124 / 0.6e1;
t246 = t73 - t127 / 0.6e1;
t77 = -0.2e1 / 0.3e1 * t124;
t87 = -0.2e1 / 0.3e1 * t127;
t245 = t77 + t87;
t244 = 0.4e1 / 0.7e1 * t138 - t124 / 0.7e1;
t218 = t136 + t138;
t182 = t132 + t218;
t215 = 0.2e1 * t52;
t210 = pkin(2) * t250;
t179 = pkin(3) * t210;
t34 = -0.2e1 * t179;
t212 = pkin(1) * t255;
t42 = -0.2e1 * t212;
t15 = 0.1e1 / (t215 * t39 + t182 + t34 + t42);
t243 = 0.1e1 / pkin(4) * t15;
t242 = t130 * t61 ^ 2;
t241 = t132 * t61;
t64 = t97 ^ 2;
t240 = t136 * t64;
t223 = t124 + t127;
t239 = t138 * t223;
t142 = pkin(2) * t136;
t238 = t142 * t97 * t64;
t237 = t143 * t64 ^ 2;
t60 = t96 * t61;
t236 = t146 * t60;
t120 = 3 * t138;
t43 = t120 - t223;
t235 = t43 * t132;
t78 = -0.3e1 / 0.2e1 * t124;
t234 = t120 + t78;
t133 = t142 ^ 2;
t137 = t138 ^ 2;
t221 = t130 + t137;
t229 = t138 * t124;
t233 = t133 + t55 * ((t78 + t121) * t132 - 0.3e1 / 0.2e1 * t229 + t221 + t247);
t92 = t136 / 0.2e1;
t232 = t138 + t92;
t186 = t77 + t55;
t82 = 0.2e1 / 0.3e1 * t127;
t231 = t143 + (t82 + t186) * t55;
t47 = -t132 / 0.3e1 + t138;
t230 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
t228 = t138 * t132;
t104 = 10 * t132;
t227 = t104 + t121;
t57 = -3 * t132 + t138;
t226 = t121 - t124;
t225 = t123 - t126;
t224 = -t124 + t127;
t222 = t124 - t138;
t220 = t132 - t138;
t219 = -t136 + t138;
t217 = t137 - t130;
t216 = t137 + t143;
t214 = 0.2e1 * t254;
t213 = 0.4e1 * pkin(1);
t211 = pkin(1) * t52;
t209 = 0.8e1 * t242;
t208 = -0.4e1 * t241;
t207 = 0.4e1 * t241;
t206 = 0.2e1 * t240;
t205 = 0.8e1 * t238;
t204 = pkin(1) * t236;
t202 = t146 * t53;
t201 = 0.12e2 * t241;
t200 = 0.16e2 * t236;
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
t180 = -t124 + t218;
t162 = t132 + t180;
t158 = -t127 + t162;
t156 = t42 + t158;
t183 = t127 + t222;
t22 = t42 + t206 + t219;
t12 = sqrt(t22 * t208 + 0.4e1 * t220 * t240 + 0.4e1 * t158 * t212 - t130 + (t115 + t183) * t260 - (t138 + (pkin(2) - t257) * (pkin(2) + t257)) * (t138 + (pkin(2) - t256) * (pkin(2) + t256)) + (-(t34 + t156) * t252 + t156 * t210) * t265);
t139 = pkin(1) * t138;
t20 = t226 * t132 + t221 - t229 - t263;
t155 = t143 + t20;
t157 = -t179 + t232;
t168 = -0.4e1 * t179;
t76 = -t124 / 0.2e1;
t36 = t76 + t182;
t160 = t36 * t168;
t169 = -0.6e1 * t179;
t173 = t117 - 0.3e1 * t229 + t247;
t50 = 0.10e2 / 0.3e1 * t132;
t161 = ((t50 + t226) * t136 + t155) * t169 + (t102 + (-0.9e1 * t124 + (18 * t138)) * t132 + t173) * t136 + (t103 + t234) * t143 + t233;
t185 = t111 + t218;
t163 = t185 * t53 - t254 * (t114 + t55);
t164 = -(t113 + ((5 * t132) + t43) * t259 + (t87 + t186) * t55) * t254 + (t143 + (t227 + t245) * t136 + t108 + 0.2e1 * t235 + t138 * (t138 + t245)) * t53;
t165 = (6 * t137) - 0.6e1 * t229 + t225;
t191 = t138 + t264;
t167 = t132 + t191;
t196 = t138 + t246;
t171 = t92 + t196;
t197 = t121 + t77 + t82;
t83 = t127 / 0.3e1;
t27 = t75 + t83 + t182;
t172 = t27 * t168 + t231 + (t110 + t197) * t136;
t174 = t60 * t202;
t175 = t242 * t53;
t176 = pkin(1) * t200;
t177 = 0.20e2 / 0.3e1 * t132;
t178 = 0.8e1 * t204;
t184 = t120 + t224;
t187 = 0.6e1 * t211;
t188 = 0.4e1 * t211;
t190 = t94 * t230;
t192 = t124 / 0.3e1 + t83 + t121;
t193 = 0.2e1 / 0.3e1 * t124 + t82 + t119;
t194 = 0.4e1 / 0.3e1 * t124 + 0.4e1 / 0.3e1 * t127 - (2 * t138);
t195 = t138 + t76 - t127 / 0.2e1;
t198 = -t254 / 0.2e1;
t48 = t138 - t136 / 0.3e1;
t23 = t48 * t34;
t28 = t39 + t52;
t58 = t115 + t138;
t31 = t58 * t178;
t33 = 0.10e2 * t235;
t37 = t127 + t181;
t41 = 0.2e1 * t211;
t45 = (t119 + t124) * t132;
t51 = -0.30e2 * t124 + (60 * t138);
t74 = -t124 / 0.4e1;
t88 = 0.4e1 / 0.3e1 * t132;
t89 = t132 / 0.3e1;
t90 = t132 / 0.2e1;
t199 = ((-0.24e2 * (0.4e1 / 0.3e1 * t241 + t41 + t47) * t237 * t254 - 0.12e2 * (-0.8e1 / 0.3e1 * t175 + ((t88 + t171) * t53 - (0.4e1 / 0.3e1 * t136 + t90 + t196) * t254) * t207 + (-(t136 * t220) - 0.5e1 / 0.3e1 * t130 + t192 * t132 + t138 * t191) * t53 + (-t143 + (-t177 + t193) * t136 - (3 * t130) + t194 * t132 + t137) * t198 + (-t94 * t130 * t60 + ((t132 + t171) * t53 + (t259 - t220) * t198) * t52) * t213) * t240 + 0.24e2 * t48 * t175 + ((t114 + 0.3e1 / 0.2e1 * t132 + t195) * t53 + t58 * t254 / 0.2e1) * t176 - 0.6e1 * ((-(3 * t143) + (-t177 + t194) * t136 + t193 * t132 + t217) * t53 - 0.2e1 * (-0.5e1 / 0.3e1 * t143 + (-t132 + t192) * t136 + t138 * t167) * t254) * t241 - 0.6e1 * t164 * t211 - (t133 + ((21 * t132) + t43) * t143 + (t117 + (35 * t130) + t33 - 0.2e1 * t239) * t136 + (t107 + (t106 + t118 - 0.5e1 * t127) * t132 - t138 * t183) * t55) * t53 + (0.7e1 * t133 + (t109 + t43) * t113 + ((21 * t130) + (9 * t137) + t33 - 0.6e1 * t239) * t136 + t251) * t254) * t12 + (0.16e2 * (t209 + t176 + (-8 * t130 + 12 * t228) * t61 + (-0.12e2 * pkin(1) * t146 + t139 * t265) * t96 - (6 * t228) + t221) * t237 + 0.24e2 * ((t138 - 0.2e1 / 0.3e1 * t136) * t209 + 0.14e2 * (-0.32e2 / 0.21e2 * (t138 + t136 / 0.4e1 + t132 / 0.4e1 - t124 / 0.8e1) * t179 + t143 / 0.7e1 + (0.16e2 / 0.21e2 * t132 + t244) * t136 + t130 / 0.7e1 + t244 * t132 + t137 - 0.3e1 / 0.7e1 * t229 + t123 / 0.42e2 - t126 / 0.42e2) * t241 + t47 * t160 - (t220 * t143) + (-0.10e2 / 0.3e1 * t130 + (2 * t137) - t229 + t45) * t136 + t20 * t248 + ((-0.2e1 / 0.3e1 * t179 + t138 + t90 + t74) * t200 + 0.6e1 * (-0.8e1 / 0.3e1 * (t74 + t89 + t232) * t179 + t143 / 0.3e1 + (0.4e1 / 0.3e1 * t138 + t88 + t75) * t136 + t137 + 0.2e1 / 0.3e1 * t228 - 0.2e1 / 0.3e1 * t229 - t130 / 0.3e1 + t123 / 0.18e2 - t126 / 0.18e2) * t52) * pkin(1)) * t240 + 0.16e2 * (-6 * t138 * t136 + t216) * t242 + 0.32e2 * (t230 * t34 + t36 * t58) * t204 + 0.24e2 * (t48 * t160 - t133 + (-t50 + t222) * t143 + (t45 + t217 + t263) * t136 + t20 * t138) * t241 + 0.8e1 * t161 * t211 - 0.8e1 * ((t109 + t234) * t143 + (t107 + (t106 + (10 * t138)) * t132 + t173) * t136 + t233) * t179 + (t143 ^ 2) + (t105 + t119 + (28 * t132)) * t133 + (t132 * t51 + (70 * t130) + t165) * t143 + (t165 * t110 + t225 * t121 - 0.6e1 * t137 * t124 + t51 * t130 + 0.4e1 * t139 ^ 2 + (28 * t146 ^ 2)) * t136 + t37 * t251) * t28 + (((0.4e1 * t253 + (t53 + t214) * t41 + t57 * t53 + (0.3e1 / 0.2e1 * t136 + t111 + t195) * t214) * t205 + 0.6e1 * ((0.2e1 * (t90 + t136 + t246) * t53 + pkin(3) * t190) * t208 + (-0.8e1 * t174 + 0.4e1 * ((t185 + t264) * t53 - (t114 + t167) * t254) * t52) * pkin(1) + t164) * t255) * t12 + (-0.32e2 * (t178 + (-0.4e1 * t94 * t202 + t261 + ((4 * t136) + t105 + t118) * t132) * t61 + (-t132 + t157 + t74) * t188 + t34 * t248 + t57 * t36) * t238 - 0.8e1 * (t31 + (t230 * t36 + t23) * t201 + (t160 + (t110 + t226) * t136 + t155) * t187 + t161) * t255) * t28) * t40) / ((-0.4e1 * (0.2e1 * t253 + (t260 + t136) * t254 + (-t220 + t41) * t53) * t240 + 0.8e1 * pkin(1) * t174 + ((pkin(2) * t261 + 0.8e1 * t132 * t142) * t95 + 0.4e1 * t146 * t190) * t61 - 0.4e1 * t163 * t211 - (t227 * t136 + t108 + t216 + 6 * t228) * t53 + (t113 + (t104 + 6 * t138) * t136 + t49) * t254) * t12 + (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t179 + t138 + t136 / 0.3e1 + t89 + t127 / 0.9e1 - t124 / 0.9e1) * t241 + t47 * t34 + t27 * t248 + (t236 + (t127 / 0.6e1 + t73 + t157) * t52) * t213) * t240 + t31 + (t230 * t27 + t23) * t201 + t172 * t187 + ((t50 + t197) * t136 + t231) * t169 + t133 + (t103 + t184) * t143 + (t184 * t110 + t224 * t121 + t102 + t117) * t136 + t49 * t37) * t28 + ((t205 * t254 + 0.4e1 * (t203 * t262 + (t53 - t254) * t41 + t163) * t255) * t12 + (-0.8e1 * (t41 + t207 + t57) * t238 - 0.6e1 * (t219 * t207 + (t34 + t27) * t188 + t172) * t255) * t28) * t40) / 0.4e1;
t189 = 0.1e1 / pkin(5) / pkin(4) ^ 2 * t15;
t170 = -m(1) - m(2) - m(3) - m(4) - m(5);
t166 = -pkin(3) + t210;
t159 = t111 + t127 + t180;
t26 = t96 * t97 + t250;
t25 = t94 * t97 - t249;
t21 = t127 + t42 + t162;
t11 = (-t166 + t252) * t12 + (t21 * t39 + t215 * t22) * t94 + (t21 * t96 + (0.4e1 * t61 - 0.2e1) * t39 * pkin(3)) * t53;
t10 = (pkin(2) * t249 + t39 * t94) * t12 - (t159 + t42 + t168) * t252 + t166 * t42 + t159 * t210 + (t22 * t262 - t114 + t206 - t37) * pkin(3);
t9 = (t10 * t258 - t94 * t11 / 0.2e1) * t243;
t8 = (t94 * t10 / 0.2e1 + t11 * t258) * t243;
t5 = (t10 * t12 / 0.4e1 + t11 * t199) * t189;
t4 = (t10 * t199 - t11 * t12 / 0.4e1) * t189;
t3 = t25 * t5 - t26 * t4;
t2 = t25 * t4 + t26 * t5;
t1 = t97 * t3;
t6 = (t170 * r_base(3) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3)) * g(3) + (-mrSges(1,2) - mrSges(2,1) * t95 - mrSges(2,2) * t97 - m(3) * t53 - (-t2 * t97 + t3 * t95) * mrSges(3,1) - (t2 * t95 + t1) * mrSges(3,2) - mrSges(4,2) * t96 - t8 * mrSges(5,1) - t9 * mrSges(5,2) + (-pkin(3) * m(5) - mrSges(4,1)) * t94 + t170 * r_base(2)) * g(2) + (-mrSges(1,1) - t1 * mrSges(3,1) - m(4) * pkin(1) - t96 * mrSges(4,1) + t94 * mrSges(4,2) - m(5) * t40 - t9 * mrSges(5,1) + t8 * mrSges(5,2) + (-mrSges(3,1) * t2 + mrSges(3,2) * t3 + mrSges(2,2)) * t95 + (-pkin(2) * m(3) - mrSges(3,2) * t2 - mrSges(2,1)) * t97 + t170 * r_base(1)) * g(1);
U = t6;
