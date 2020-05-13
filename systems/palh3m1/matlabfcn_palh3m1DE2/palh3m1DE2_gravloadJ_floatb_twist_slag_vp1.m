% Calculate Gravitation load on the joints for
% palh3m1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% rSges [9x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 16:51
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh3m1DE2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(19,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1DE2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1DE2_gravloadJ_floatb_twist_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m1DE2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [9x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-19 19:37:30
% EndTime: 2020-04-19 19:41:57
% DurationCPUTime: 44.41s
% Computational Cost: add. (1322030->255), mult. (1996245->454), div. (87822->22), fcn. (1265683->38), ass. (0->220)
t126 = sin(qJ(1));
t132 = cos(qJ(1));
t271 = g(1) * t132 + g(2) * t126;
t274 = 2 * pkin(3);
t141 = pkin(4) ^ 2;
t119 = pkin(18) + pkin(19);
t114 = sin(t119);
t115 = cos(t119);
t130 = cos(qJ(3));
t144 = 0.1e1 / pkin(2);
t140 = pkin(5) ^ 2;
t145 = pkin(1) ^ 2;
t125 = sin(qJ(2));
t127 = sin(pkin(16));
t131 = cos(qJ(2));
t249 = cos(pkin(16));
t100 = t125 * t127 - t131 * t249;
t247 = pkin(5) * t100;
t269 = -2 * pkin(1);
t225 = t247 * t269 + t145;
t93 = t140 + t225;
t91 = 0.1e1 / t93;
t230 = t144 * t91;
t124 = sin(qJ(3));
t253 = t124 / 0.2e1;
t101 = t125 * t249 + t131 * t127;
t262 = pkin(5) + pkin(6);
t263 = pkin(5) - pkin(6);
t87 = (pkin(2) + t262) * (-pkin(2) + t263) + t225;
t88 = (-pkin(2) + t262) * (pkin(2) + t263) + t225;
t147 = sqrt(-t88 * t87);
t223 = t101 * t147;
t213 = pkin(2) ^ 2 - pkin(6) ^ 2;
t90 = t93 + t213;
t94 = pkin(1) - t247;
t81 = -pkin(5) * t223 + t90 * t94;
t265 = -t81 / 0.2e1;
t246 = pkin(5) * t101;
t82 = t147 * t94 + t90 * t246;
t75 = (t130 * t265 + t82 * t253) * t230;
t264 = t82 / 0.2e1;
t76 = (t130 * t264 + t81 * t253) * t230;
t56 = t114 * t76 + t115 * t75;
t256 = t56 * pkin(4);
t226 = t256 * t274 + t141;
t261 = (-pkin(8) - pkin(10));
t44 = ((pkin(3) - t261) * (pkin(3) + t261)) + t226;
t260 = (pkin(10) - pkin(8));
t45 = ((pkin(3) - t260) * (pkin(3) + t260)) + t226;
t146 = sqrt(-t45 * t44);
t173 = t114 * t75 - t115 * t76;
t273 = t146 * t173;
t120 = qJ(2) + qJ(3);
t116 = sin(t120);
t117 = cos(t120);
t272 = -rSges(4,1) * t117 + t116 * rSges(4,2);
t122 = cos(pkin(17));
t135 = 0.1e1 / pkin(10);
t142 = pkin(3) ^ 2;
t50 = t142 + t226;
t48 = 0.1e1 / t50;
t232 = t135 * t48;
t121 = sin(pkin(17));
t254 = t121 / 0.2e1;
t214 = pkin(8) ^ 2 - pkin(10) ^ 2;
t47 = t50 - t214;
t52 = pkin(3) * t56 + pkin(4);
t32 = -pkin(3) * t273 + t47 * t52;
t34 = pkin(3) * t173 * t47 + t146 * t52;
t22 = (-t32 * t122 / 0.2e1 + t34 * t254) * t232;
t23 = (t34 * t122 / 0.2e1 + t32 * t254) * t232;
t13 = atan2(t23, t22) + t120;
t11 = sin(t13);
t12 = cos(t13);
t250 = pkin(11) + rSges(6,3);
t270 = -pkin(9) * t12 - t250 * t11;
t238 = sin(pkin(19));
t206 = t238 / 0.2e1;
t239 = cos(pkin(19));
t72 = (t82 * t206 + t239 * t265) * t230;
t268 = 0.1e1 / t72 ^ 2;
t40 = 0.1e1 / t146;
t267 = -t40 / 0.2e1;
t266 = t48 / 0.2e1;
t73 = (t81 * t206 + t239 * t264) * t230;
t63 = qJ(2) + atan2(t73, t72);
t62 = pkin(18) - t63;
t259 = pkin(3) * cos(t62);
t49 = 0.1e1 / t50 ^ 2;
t258 = pkin(4) * t49;
t257 = pkin(4) * t173;
t133 = cos(pkin(15));
t252 = t133 / 0.2e1;
t251 = -0.2e1 * t101 ^ 2;
t248 = pkin(4) * t117;
t118 = t131 * pkin(1);
t112 = t118 + pkin(13);
t244 = g(2) * t132 * t112;
t242 = t125 * pkin(1);
t21 = 0.1e1 / t22 ^ 2;
t241 = t21 * t23;
t211 = pkin(1) * t246;
t240 = 0.2e1 / t147 * (t87 + t88) * t211;
t236 = rSges(4,2) * t117;
t235 = t122 * t48;
t128 = sin(pkin(15));
t234 = t128 * t91;
t233 = t135 / (t21 * t23 ^ 2 + 0.1e1);
t139 = 0.1e1 / pkin(6);
t231 = t139 * t91;
t151 = (t94 * t240 / 0.2e1 + t140 * pkin(1) * t251 + (-t100 * t90 - t223) * pkin(5)) * t230 / 0.2e1;
t207 = -t240 / 0.2e1;
t224 = t100 * t147;
t154 = pkin(5) * (t224 + (t94 * t269 + t207 - t90) * t101) * t230;
t152 = t154 / 0.2e1;
t198 = 0.1e1 / t93 ^ 2 * t211;
t185 = t144 * t198;
t170 = t130 * t185;
t171 = t124 * t185;
t42 = t124 * t152 + t130 * t151 + t82 * t170 + t81 * t171;
t153 = -t154 / 0.2e1;
t43 = t124 * t151 + t130 * t153 - t81 * t170 + t82 * t171;
t37 = -t114 * t42 - t115 * t43;
t229 = t146 * t37;
t222 = t116 * t126;
t221 = t116 * t132;
t123 = sin(qJ(4));
t220 = t123 * t126;
t219 = t123 * t132;
t129 = cos(qJ(4));
t218 = t126 * t129;
t217 = t129 * t132;
t216 = rSges(4,1) * t222 + t126 * t236;
t215 = rSges(4,1) * t221 + t132 * t236;
t212 = pkin(3) * t258;
t51 = -pkin(3) - t256;
t210 = t51 * t267;
t209 = t40 * t52 / 0.2e1;
t208 = t173 * t267;
t205 = t48 * t254;
t204 = -t235 / 0.2e1;
t203 = t235 / 0.2e1;
t202 = t91 * t252;
t137 = 0.1e1 / pkin(8);
t201 = t137 * t266;
t200 = -0.2e1 * pkin(4) * t52 - t47;
t199 = -0.2e1 * t142 * t257;
t197 = t121 * t212;
t196 = t122 * t212;
t195 = t118 - t248;
t46 = t50 + t214;
t31 = -pkin(4) * t273 - t46 * t51;
t30 = 0.1e1 / t31 ^ 2;
t33 = -t146 * t51 + t46 * t257;
t194 = pkin(8) * t137 / (t30 * t33 ^ 2 + 0.1e1) * t50;
t191 = t37 * t197;
t190 = t173 * t197;
t189 = t37 * t196;
t188 = t173 * t196;
t187 = t128 * t198;
t186 = t133 * t198;
t184 = 0.1e1 / t31 * t194;
t183 = -rSges(5,1) * t12 + rSges(5,2) * t11;
t89 = t93 - t213;
t95 = pkin(1) * t100 - pkin(5);
t80 = -pkin(1) * t223 - t89 * t95;
t83 = pkin(1) * t101 * t89 - t147 * t95;
t74 = (t80 * t252 + t83 * t128 / 0.2e1) * t231;
t77 = (t83 * t252 - t80 * t128 / 0.2e1) * t231;
t68 = atan2(t77, t74);
t64 = sin(t68);
t65 = cos(t68);
t181 = rSges(7,1) * t65 - rSges(7,2) * t64;
t60 = sin(t63);
t61 = cos(t63);
t179 = rSges(8,1) * t61 - rSges(8,2) * t60;
t19 = -atan2(t33 * t201, t31 * t201) + t62;
t17 = sin(t19);
t18 = cos(t19);
t178 = rSges(9,1) * t18 + rSges(9,2) * t17;
t177 = rSges(9,1) * t17 - rSges(9,2) * t18;
t176 = rSges(3,1) * t131 - rSges(3,2) * t125;
t172 = pkin(4) * (t44 + t45) * t274;
t169 = -pkin(7) + t181;
t168 = t178 - t112 - t259;
t167 = rSges(6,1) * t129 - rSges(6,2) * t123 + pkin(9);
t166 = pkin(13) + t176;
t165 = pkin(3) * (t31 * t49 + t48 * t51);
t164 = pkin(4) * t30 * t33 * t194;
t25 = t37 * t172;
t36 = t114 * t43 - t115 * t42;
t163 = -t36 * t146 + t25 * t208;
t35 = t173 * t172;
t162 = -t56 * t146 + t35 * t208;
t161 = t239 * t185;
t160 = t238 * t185;
t158 = pkin(3) * (-t141 * t173 * t48 + t33 * t258);
t157 = -g(3) * t248 + (g(1) * t221 + g(2) * t222) * pkin(4);
t27 = -0.1e1 + ((t151 * t238 + t153 * t239 + t82 * t160 - t81 * t161) * t73 * t268 - (t151 * t239 + t152 * t238 + t81 * t160 + t82 * t161) / t72) / (t73 ^ 2 * t268 + 0.1e1);
t155 = g(3) * t195 + t271 * (pkin(4) * t116 - t242);
t150 = g(3) * t183 + t271 * (rSges(5,1) * t11 + rSges(5,2) * t12);
t149 = (-g(3) * t167 - t250 * t271) * t12 + (-g(3) * t250 + t167 * t271) * t11;
t102 = pkin(13) + t195;
t97 = t132 * t102;
t71 = 0.1e1 / t74 ^ 2;
t70 = t95 * t207 + t145 * pkin(5) * t251 + (-t100 * t89 - t223) * pkin(1);
t69 = (t224 + (0.2e1 * t95 * pkin(5) + t207 - t89) * t101) * pkin(1);
t20 = 0.1e1 / t22;
t16 = t35 * t209 + t173 * t199 + (t47 * t56 - t273) * pkin(3);
t15 = (t173 * t200 + t162) * pkin(3);
t10 = t12 * t217 - t220;
t9 = t12 * t219 + t218;
t8 = t12 * t218 + t219;
t7 = t12 * t220 - t217;
t6 = t25 * t209 + t37 * t199 + (t36 * t47 - t229) * pkin(3);
t5 = (t200 * t37 + t163) * pkin(3);
t3 = -0.2e1 * ((t25 * t210 + (t36 * t46 - t229) * pkin(4)) * t266 + t37 * t158) * t184 + 0.2e1 * ((-t37 * t46 + t163) * t266 + t37 * t165) * t164 + t27;
t2 = 0.1e1 + ((t15 * t205 + t16 * t203 + t34 * t188 + t32 * t190) * t20 - (t15 * t204 + t16 * t205 - t32 * t188 + t34 * t190) * t241) * t233;
t1 = 0.1e1 + ((t34 * t189 + t32 * t191 + t6 * t203 + t5 * t205) * t20 - (-t32 * t189 + t34 * t191 + t5 * t204 + t6 * t205) * t241) * t233;
t4 = [-m(2) * (g(1) * (-rSges(2,1) * t126 - rSges(2,2) * t132) + g(2) * (rSges(2,1) * t132 - rSges(2,2) * t126)) - m(3) * ((g(1) * rSges(3,3) + g(2) * t166) * t132 + (g(2) * rSges(3,3) - g(1) * t166) * t126) - m(4) * (t244 + (g(1) * rSges(4,3) + g(2) * t272) * t132 + (g(1) * (-t112 - t272) + g(2) * rSges(4,3)) * t126) - m(5) * (g(2) * t97 + (g(1) * rSges(5,3) + g(2) * t183) * t132 + (g(1) * (-t102 - t183) + g(2) * rSges(5,3)) * t126) - m(6) * ((-rSges(6,1) * t10 + rSges(6,2) * t9 + t132 * t270 + t97) * g(2) + (rSges(6,1) * t8 - rSges(6,2) * t7 + (-t102 - t270) * t126) * g(1)) - m(7) * ((g(1) * rSges(7,3) + g(2) * t169) * t132 + (g(2) * rSges(7,3) - g(1) * t169) * t126) - m(8) * (t244 + (g(1) * rSges(8,3) + g(2) * t179) * t132 + (g(1) * (-t112 - t179) + g(2) * rSges(8,3)) * t126) - m(9) * ((g(1) * rSges(9,3) - g(2) * t168) * t132 + (g(2) * rSges(9,3) + g(1) * t168) * t126), -m(3) * (g(3) * t176 + t271 * (-rSges(3,1) * t125 - rSges(3,2) * t131)) - m(4) * (g(1) * (-t132 * t242 + t215) + g(2) * (-t126 * t242 + t216) + g(3) * (t118 + t272)) - m(5) * (t150 * t1 + t155) - m(6) * (t149 * t1 + t155) - m(7) * (g(3) * t181 + t271 * (-rSges(7,1) * t64 - rSges(7,2) * t65)) * ((t70 * t202 + t83 * t186 - t69 * t234 / 0.2e1 - t80 * t187) / t74 - (t69 * t202 + t80 * t186 + t70 * t234 / 0.2e1 + t83 * t187) * t77 * t71) / (t71 * t77 ^ 2 + 0.1e1) * t139 - m(8) * (g(3) * (-t179 * t27 + t118) + t271 * (-t242 + (rSges(8,1) * t60 + rSges(8,2) * t61) * t27)) - m(9) * (g(3) * (t178 * t3 - t27 * t259 + t118) + t271 * (t177 * t3 - pkin(3) * t27 * sin(t62) - t242)), -m(4) * (g(1) * t215 + g(2) * t216 + g(3) * t272) - m(5) * (t150 * t2 + t157) - m(6) * (t149 * t2 + t157) - 0.2e1 * m(9) * (g(3) * t178 + t177 * t271) * (-((t35 * t210 + (t46 * t56 - t273) * pkin(4)) * t266 + t173 * t158) * t184 + ((-t173 * t46 + t162) * t266 + t173 * t165) * t164), -m(6) * (g(1) * (rSges(6,1) * t9 + rSges(6,2) * t10) + g(2) * (rSges(6,1) * t7 + rSges(6,2) * t8) + g(3) * (rSges(6,1) * t123 + rSges(6,2) * t129) * t11)];
taug = t4(:);
