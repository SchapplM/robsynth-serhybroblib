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
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 04:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fivebar1DE1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fivebar1DE1_energypot_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'fivebar1DE1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1DE1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1DE1_energypot_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1DE1_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fivebar1DE1_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 02:23:13
% EndTime: 2020-04-27 02:23:23
% DurationCPUTime: 4.82s
% Computational Cost: add. (28347->415), mult. (83154->496), div. (256->5), fcn. (9996->12), ass. (0->208)
t268 = 4 * pkin(3);
t132 = pkin(4) ^ 2;
t129 = pkin(5) ^ 2;
t80 = -t129 / 0.3e1;
t267 = t80 - t132 / 0.3e1;
t128 = t129 ^ 2;
t131 = t132 ^ 2;
t266 = -t128 / 0.6e1 + t131 / 0.6e1;
t101 = cos(qJ(2));
t66 = t101 ^ 2;
t265 = -0.2e1 * t66;
t137 = pkin(3) ^ 2;
t135 = t137 ^ 2;
t264 = 4 * t135;
t263 = 2 * t137;
t115 = 6 * t137;
t141 = (pkin(2) ^ 2);
t262 = 2 * t141;
t143 = (pkin(1) ^ 2);
t126 = 2 * t143;
t148 = t141 ^ 2;
t118 = 5 * t148;
t261 = -pkin(4) - pkin(5);
t260 = -pkin(4) + pkin(5);
t99 = sin(qJ(2));
t259 = t99 * pkin(3);
t57 = pkin(3) * t101;
t45 = pkin(1) + t57;
t100 = sin(qJ(1));
t58 = pkin(2) * t100;
t102 = cos(qJ(1));
t258 = pkin(2) * t102;
t151 = pkin(3) * t137;
t205 = t137 * t58;
t257 = (-t151 * t99 + t205) * t66;
t60 = t137 + t143;
t185 = -t129 + t60;
t43 = -t132 + t185;
t54 = t60 ^ 2;
t256 = t54 * t43;
t255 = t128 / 0.2e1 - t131 / 0.2e1;
t78 = -t129 / 0.6e1;
t254 = t78 - t132 / 0.6e1;
t82 = -0.2e1 / 0.3e1 * t129;
t92 = -0.2e1 / 0.3e1 * t132;
t253 = t82 + t92;
t252 = 0.4e1 / 0.7e1 * t143 - t129 / 0.7e1;
t251 = t100 * t99;
t44 = pkin(1) - t258;
t250 = t101 * t44;
t249 = t135 * t66 ^ 2;
t248 = t137 * t66;
t69 = t102 ^ 2;
t247 = t141 * t69;
t228 = t129 + t132;
t246 = t143 * t228;
t147 = pkin(2) * t141;
t245 = t147 * t102 * t69;
t244 = t148 * t69 ^ 2;
t65 = t101 * t66;
t243 = t151 * t65;
t125 = 3 * t143;
t48 = t125 - t228;
t242 = t48 * t137;
t83 = -0.3e1 / 0.2e1 * t129;
t241 = t125 + t83;
t138 = t147 ^ 2;
t142 = t143 ^ 2;
t226 = t135 + t142;
t234 = t143 * t129;
t240 = t138 + t60 * ((t83 + t126) * t137 - 0.3e1 / 0.2e1 * t234 + t226 + t255);
t97 = t141 / 0.2e1;
t239 = t143 + t97;
t193 = t82 + t60;
t87 = 0.2e1 / 0.3e1 * t132;
t238 = t148 + t60 * (t87 + t193);
t52 = -t137 / 0.3e1 + t143;
t237 = t100 * t101;
t236 = (pkin(1) + pkin(3)) * (pkin(1) - pkin(3));
t235 = (pkin(1) + pkin(2)) * (pkin(1) - pkin(2));
t233 = t143 * t137;
t109 = 10 * t137;
t232 = t109 + t126;
t62 = -3 * t137 + t143;
t231 = t126 - t129;
t230 = t128 - t131;
t229 = -t129 + t132;
t227 = t129 - t143;
t225 = t137 - t143;
t224 = -t141 + t143;
t223 = t141 + t143;
t222 = t142 - t135;
t221 = t142 + t148;
t220 = 0.2e1 * t259;
t219 = 0.4e1 * pkin(1);
t218 = 0.2e1 * t57;
t217 = pkin(1) * t258;
t216 = pkin(1) * t57;
t215 = 0.8e1 * t249;
t214 = -0.4e1 * t248;
t213 = 0.4e1 * t248;
t212 = 0.2e1 * t247;
t211 = 0.8e1 * t245;
t210 = pkin(1) * t243;
t209 = pkin(2) * t251;
t208 = 0.12e2 * t248;
t207 = 0.16e2 * t243;
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
t25 = t231 * t137 + t226 - t234 - t266;
t160 = t148 + t25;
t182 = pkin(3) * t209;
t162 = -t182 + t239;
t172 = -0.4e1 * t182;
t186 = t137 + t223;
t81 = -t129 / 0.2e1;
t41 = t81 + t186;
t165 = t41 * t172;
t173 = -0.6e1 * t182;
t178 = t122 - 0.3e1 * t234 + t255;
t55 = 0.10e2 / 0.3e1 * t137;
t166 = ((t55 + t231) * t141 + t160) * t173 + (t107 + (-0.9e1 * t129 + (18 * t143)) * t137 + t178) * t141 + (t108 + t241) * t148 + t240;
t189 = t116 + t223;
t168 = t189 * t58 - t259 * (t119 + t60);
t169 = -(t118 + ((5 * t137) + t48) * t262 + t60 * (t92 + t193)) * t259 + (t148 + (t232 + t253) * t141 + t113 + 0.2e1 * t242 + t143 * (t143 + t253)) * t58;
t120 = -3 * t141;
t163 = t141 + t43;
t47 = -0.2e1 * t217;
t161 = t47 + t163;
t187 = t132 + t227;
t27 = t47 + t212 + t224;
t39 = -0.2e1 * t182;
t17 = sqrt(t27 * t214 + 0.4e1 * t225 * t247 + 0.4e1 * t163 * t217 - t135 + (t120 + t187) * t263 - (t143 + (pkin(2) - t260) * (pkin(2) + t260)) * (t143 + (pkin(2) - t261) * (pkin(2) + t261)) + (-(t39 + t161) * t250 + t161 * t209) * t268);
t171 = (6 * t142) - 0.6e1 * t234 + t230;
t196 = t143 + t267;
t174 = t137 + t196;
t201 = t143 + t254;
t176 = t97 + t201;
t202 = t126 + t82 + t87;
t88 = t132 / 0.3e1;
t32 = t80 + t88 + t186;
t177 = t32 * t172 + t238 + (t115 + t202) * t141;
t179 = t65 * t206;
t180 = t249 * t58;
t181 = pkin(1) * t207;
t183 = 0.20e2 / 0.3e1 * t137;
t184 = 0.8e1 * t210;
t188 = t125 + t229;
t190 = 0.6e1 * t216;
t191 = 0.4e1 * t216;
t195 = t99 * t235;
t197 = t129 / 0.3e1 + t88 + t126;
t198 = 0.2e1 / 0.3e1 * t129 + t87 + t124;
t199 = 0.4e1 / 0.3e1 * t129 + 0.4e1 / 0.3e1 * t132 - (2 * t143);
t200 = t143 + t81 - t132 / 0.2e1;
t203 = -t259 / 0.2e1;
t53 = t143 - t141 / 0.3e1;
t28 = t53 * t39;
t33 = t44 + t57;
t63 = t120 + t143;
t36 = t63 * t184;
t38 = 0.10e2 * t242;
t42 = t132 + t185;
t46 = 0.2e1 * t216;
t50 = (t124 + t129) * t137;
t56 = -0.30e2 * t129 + (60 * t143);
t79 = -t129 / 0.4e1;
t93 = 0.4e1 / 0.3e1 * t137;
t94 = t137 / 0.3e1;
t95 = t137 / 0.2e1;
t204 = 0.1e1 / ((-0.4e1 * (0.2e1 * t257 + (t263 + t141) * t259 + (-t225 + t46) * t58) * t247 + 0.8e1 * pkin(1) * t179 + ((pkin(2) * t264 + 0.8e1 * t137 * t147) * t100 + 0.4e1 * t151 * t195) * t66 - 0.4e1 * t168 * t216 - (t232 * t141 + t113 + t221 + 6 * t233) * t58 + (t118 + (t109 + 6 * t143) * t141 + t54) * t259) * t17 + t33 * (0.12e2 * (0.6e1 * (-0.4e1 / 0.9e1 * t182 + t143 + t141 / 0.3e1 + t94 + t132 / 0.9e1 - t129 / 0.9e1) * t248 + t52 * t39 + t32 * t236 + (t243 + (t132 / 0.6e1 + t78 + t162) * t57) * t219) * t247 + t36 + (t235 * t32 + t28) * t208 + t177 * t190 + ((t55 + t202) * t141 + t238) * t173 + t138 + (t108 + t188) * t148 + (t188 * t115 + t229 * t126 + t107 + t122) * t141 + t54 * t42) + ((t211 * t259 + 0.4e1 * (t205 * t265 + (t58 - t259) * t46 + t168) * t258) * t17 + t33 * (-0.8e1 * (t46 + t213 + t62) * t245 - 0.6e1 * (t224 * t213 + (t39 + t32) * t191 + t177) * t258)) * t45) * ((-0.24e2 * (0.4e1 / 0.3e1 * t248 + t46 + t52) * t244 * t259 - 0.12e2 * (-0.8e1 / 0.3e1 * t180 + ((t93 + t176) * t58 - (0.4e1 / 0.3e1 * t141 + t95 + t201) * t259) * t213 + (-(t141 * t225) - 0.5e1 / 0.3e1 * t135 + t197 * t137 + t143 * t196) * t58 + (-t148 + (-t183 + t198) * t141 - (3 * t135) + t199 * t137 + t142) * t203 + (-t99 * t135 * t65 + ((t137 + t176) * t58 + (t262 - t225) * t203) * t57) * t219) * t247 + 0.24e2 * t53 * t180 + ((t119 + 0.3e1 / 0.2e1 * t137 + t200) * t58 + t63 * t259 / 0.2e1) * t181 - 0.6e1 * ((-(3 * t148) + (-t183 + t199) * t141 + t198 * t137 + t222) * t58 - 0.2e1 * (-0.5e1 / 0.3e1 * t148 + (-t137 + t197) * t141 + t143 * t174) * t259) * t248 - 0.6e1 * t169 * t216 - (t138 + ((21 * t137) + t48) * t148 + (t122 + (35 * t135) + t38 - 0.2e1 * t246) * t141 + (t112 + (t111 + t123 - 0.5e1 * t132) * t137 - t143 * t187) * t60) * t58 + (0.7e1 * t138 + (t114 + t48) * t118 + ((21 * t135) + (9 * t142) + t38 - 0.6e1 * t246) * t141 + t256) * t259) * t17 + t33 * (0.16e2 * (t215 + t181 + (-8 * t135 + 12 * t233) * t66 + (-0.12e2 * pkin(1) * t151 + t144 * t268) * t101 - (6 * t233) + t226) * t244 + 0.24e2 * ((t143 - 0.2e1 / 0.3e1 * t141) * t215 + 0.14e2 * (-0.32e2 / 0.21e2 * (t143 + t141 / 0.4e1 + t137 / 0.4e1 - t129 / 0.8e1) * t182 + t148 / 0.7e1 + (0.16e2 / 0.21e2 * t137 + t252) * t141 + t135 / 0.7e1 + t252 * t137 + t142 - 0.3e1 / 0.7e1 * t234 + t128 / 0.42e2 - t131 / 0.42e2) * t248 + t52 * t165 - (t225 * t148) + (-0.10e2 / 0.3e1 * t135 + (2 * t142) - t234 + t50) * t141 + t25 * t236 + ((-0.2e1 / 0.3e1 * t182 + t143 + t95 + t79) * t207 + 0.6e1 * (-0.8e1 / 0.3e1 * (t79 + t94 + t239) * t182 + t148 / 0.3e1 + (0.4e1 / 0.3e1 * t143 + t93 + t80) * t141 + t142 + 0.2e1 / 0.3e1 * t233 - 0.2e1 / 0.3e1 * t234 - t135 / 0.3e1 + t128 / 0.18e2 - t131 / 0.18e2) * t57) * pkin(1)) * t247 + 0.16e2 * (-6 * t143 * t141 + t221) * t249 + 0.32e2 * (t235 * t39 + t41 * t63) * t210 + 0.24e2 * (t53 * t165 - t138 + (-t55 + t227) * t148 + (t50 + t222 + t266) * t141 + t143 * t25) * t248 + 0.8e1 * t166 * t216 - 0.8e1 * ((t114 + t241) * t148 + (t112 + (t111 + (10 * t143)) * t137 + t178) * t141 + t240) * t182 + (t148 ^ 2) + (t110 + t124 + (28 * t137)) * t138 + (t137 * t56 + (70 * t135) + t171) * t148 + (t171 * t115 + t230 * t126 - 0.6e1 * t142 * t129 + t56 * t135 + 0.4e1 * t144 ^ 2 + (28 * t151 ^ 2)) * t141 + t42 * t256) + (((0.4e1 * t257 + (t58 + t220) * t46 + t62 * t58 + (0.3e1 / 0.2e1 * t141 + t116 + t200) * t220) * t211 + 0.6e1 * ((0.2e1 * (t95 + t141 + t254) * t58 + pkin(3) * t195) * t214 + (-0.8e1 * t179 + 0.4e1 * ((t189 + t267) * t58 - (t119 + t174) * t259) * t57) * pkin(1) + t169) * t258) * t17 + t33 * (-0.32e2 * (t184 + (-0.4e1 * t99 * t206 + t264 + ((4 * t141) + t110 + t123) * t137) * t66 + (-t137 + t162 + t79) * t191 + t39 * t236 + t41 * t62) * t245 - 0.8e1 * (t36 + (t235 * t41 + t28) * t208 + (t165 + (t115 + t231) * t141 + t160) * t190 + t166) * t258)) * t45) / 0.4e1;
t20 = 0.1e1 / (t218 * t44 + t186 + t39 + t47);
t194 = 0.1e1 / pkin(5) / pkin(4) ^ 2 * t20;
t192 = 0.1e1 / pkin(4) * t20 / 0.2e1;
t175 = -m(1) - m(2) - m(3) - m(4) - m(5);
t170 = -pkin(3) + t209;
t167 = t223 + t229;
t164 = t116 + t167;
t31 = t101 * t102 + t251;
t30 = t102 * t99 - t237;
t26 = t137 + t47 + t167;
t16 = (-t170 + t250) * t17 + (t218 * t27 + t26 * t44) * t99 + (t101 * t26 + (0.4e1 * t66 - 0.2e1) * t44 * pkin(3)) * t58;
t15 = (pkin(2) * t237 + t44 * t99) * t17 - (t164 + t47 + t172) * t250 + t170 * t47 + t164 * t209 + (t27 * t265 - t119 + t212 - t42) * pkin(3);
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
t8 = (t175 * r_base(3) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3)) * g(3) + (-m(3) * t58 - mrSges(2,1) * t100 - t1 * mrSges(3,1) - t9 * mrSges(5,1) - mrSges(2,2) * t102 - t2 * mrSges(3,2) - mrSges(4,2) * t101 - t10 * mrSges(5,2) - mrSges(1,2) + (-m(5) * pkin(3) - mrSges(4,1)) * t99 + t175 * r_base(2)) * g(2) + (-mrSges(1,1) + mrSges(2,2) * t100 - t2 * mrSges(3,1) + t1 * mrSges(3,2) - m(4) * pkin(1) - t101 * mrSges(4,1) + t99 * mrSges(4,2) - m(5) * t45 - t10 * mrSges(5,1) + t9 * mrSges(5,2) + (-m(3) * pkin(2) - mrSges(2,1)) * t102 + t175 * r_base(1)) * g(1);
U = t8;
