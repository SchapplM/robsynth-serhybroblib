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
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 16:51
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh3m1DE2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(19,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1DE2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1DE2_gravloadJ_floatb_twist_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1DE2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [9x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-19 19:37:30
% EndTime: 2020-04-19 19:41:54
% DurationCPUTime: 44.68s
% Computational Cost: add. (1322029->226), mult. (1996260->387), div. (87822->22), fcn. (1265683->38), ass. (0->216)
t276 = m(5) + m(6);
t114 = qJ(2) + qJ(3);
t116 = cos(pkin(17));
t129 = 0.1e1 / pkin(10);
t136 = pkin(3) ^ 2;
t135 = pkin(4) ^ 2;
t113 = pkin(18) + pkin(19);
t108 = sin(t113);
t109 = cos(t113);
t124 = cos(qJ(3));
t138 = 0.1e1 / pkin(2);
t134 = pkin(5) ^ 2;
t139 = pkin(1) ^ 2;
t119 = sin(qJ(2));
t121 = sin(pkin(16));
t125 = cos(qJ(2));
t249 = cos(pkin(16));
t97 = t119 * t121 - t125 * t249;
t254 = pkin(5) * t97;
t269 = -2 * pkin(1);
t223 = t254 * t269 + t139;
t93 = t134 + t223;
t91 = 0.1e1 / t93;
t230 = t138 * t91;
t118 = sin(qJ(3));
t251 = t118 / 0.2e1;
t261 = pkin(5) + pkin(6);
t262 = pkin(5) - pkin(6);
t87 = (pkin(2) + t261) * (-pkin(2) + t262) + t223;
t88 = (-pkin(2) + t261) * (pkin(2) + t262) + t223;
t141 = sqrt(-t88 * t87);
t98 = t119 * t249 + t121 * t125;
t226 = t141 * t98;
t215 = pkin(2) ^ 2 - pkin(6) ^ 2;
t90 = t93 + t215;
t94 = pkin(1) - t254;
t81 = -pkin(5) * t226 + t90 * t94;
t265 = -t81 / 0.2e1;
t253 = pkin(5) * t98;
t82 = t141 * t94 + t253 * t90;
t75 = (t124 * t265 + t251 * t82) * t230;
t264 = t82 / 0.2e1;
t76 = (t124 * t264 + t251 * t81) * t230;
t56 = t108 * t76 + t109 * t75;
t255 = t56 * pkin(4);
t277 = 2 * pkin(3);
t224 = t255 * t277 + t135;
t50 = t136 + t224;
t48 = 0.1e1 / t50;
t232 = t129 * t48;
t115 = sin(pkin(17));
t252 = t115 / 0.2e1;
t260 = (-pkin(8) - pkin(10));
t44 = ((pkin(3) - t260) * (pkin(3) + t260)) + t224;
t259 = (pkin(10) - pkin(8));
t45 = ((pkin(3) - t259) * (pkin(3) + t259)) + t224;
t140 = sqrt(-t45 * t44);
t175 = t108 * t75 - t109 * t76;
t275 = t140 * t175;
t216 = pkin(8) ^ 2 - pkin(10) ^ 2;
t47 = t50 - t216;
t52 = pkin(3) * t56 + pkin(4);
t32 = -pkin(3) * t275 + t47 * t52;
t34 = pkin(3) * t175 * t47 + t140 * t52;
t22 = (-t32 * t116 / 0.2e1 + t34 * t252) * t232;
t23 = (t34 * t116 / 0.2e1 + t32 * t252) * t232;
t13 = atan2(t23, t22) + t114;
t11 = sin(t13);
t12 = cos(t13);
t117 = sin(qJ(4));
t123 = cos(qJ(4));
t159 = pkin(9) * m(6) + mrSges(6,1) * t123 - mrSges(6,2) * t117;
t206 = -pkin(11) * m(6) - mrSges(6,3);
t283 = (-mrSges(5,2) - t206) * t12 + (-mrSges(5,1) - t159) * t11;
t282 = -m(9) - m(4) - m(8);
t239 = sin(pkin(19));
t207 = t239 / 0.2e1;
t240 = cos(pkin(19));
t72 = (t207 * t82 + t240 * t265) * t230;
t73 = (t207 * t81 + t240 * t264) * t230;
t63 = qJ(2) + atan2(t73, t72);
t60 = sin(t63);
t61 = cos(t63);
t62 = pkin(18) - t63;
t281 = m(9) * pkin(3) * cos(t62) + mrSges(8,1) * t61 - mrSges(8,2) * t60;
t110 = sin(t114);
t111 = cos(t114);
t193 = -mrSges(4,1) * t111 + t110 * mrSges(4,2);
t280 = -mrSges(3,1) * t125 + mrSges(3,2) * t119 - t193;
t279 = mrSges(4,1) * t110 + mrSges(4,2) * t111;
t120 = sin(qJ(1));
t126 = cos(qJ(1));
t278 = g(1) * t126 + g(2) * t120;
t183 = -t12 * mrSges(5,1) + t11 * mrSges(5,2);
t274 = -t11 * t206 + t12 * t159 - t183;
t271 = mrSges(2,2) - mrSges(9,3) - mrSges(5,3) - mrSges(4,3) - mrSges(7,3) - mrSges(8,3) - mrSges(3,3);
t112 = t125 * pkin(1);
t131 = 0.1e1 / pkin(8);
t266 = t48 / 0.2e1;
t201 = t131 * t266;
t46 = t50 + t216;
t51 = -pkin(3) - t255;
t31 = -pkin(4) * t275 - t46 * t51;
t256 = pkin(4) * t175;
t33 = -t140 * t51 + t256 * t46;
t19 = -atan2(t33 * t201, t31 * t201) + t62;
t17 = sin(t19);
t18 = cos(t19);
t179 = mrSges(9,1) * t18 + mrSges(9,2) * t17;
t122 = sin(pkin(15));
t133 = 0.1e1 / pkin(6);
t231 = t133 * t91;
t127 = cos(pkin(15));
t250 = t127 / 0.2e1;
t89 = t93 - t215;
t95 = pkin(1) * t97 - pkin(5);
t80 = -pkin(1) * t226 - t89 * t95;
t83 = pkin(1) * t98 * t89 - t141 * t95;
t74 = (t80 * t250 + t83 * t122 / 0.2e1) * t231;
t77 = (t83 * t250 - t80 * t122 / 0.2e1) * t231;
t68 = atan2(t77, t74);
t64 = sin(t68);
t65 = cos(t68);
t181 = mrSges(7,1) * t65 - mrSges(7,2) * t64;
t247 = pkin(4) * t111;
t195 = t112 - t247;
t99 = pkin(13) + t195;
t270 = -m(3) * pkin(13) + m(6) * (pkin(9) * t12 + pkin(11) * t11 - t99) - mrSges(2,1) + t179 - m(5) * t99 - t183 + pkin(7) * m(7) - t181 + t11 * mrSges(6,3) + t280 + t282 * (t112 + pkin(13)) - t281;
t268 = 0.1e1 / t72 ^ 2;
t40 = 0.1e1 / t140;
t267 = -t40 / 0.2e1;
t263 = -0.2e1 * t98 ^ 2;
t49 = 0.1e1 / t50 ^ 2;
t257 = pkin(4) * t49;
t248 = pkin(4) * t110;
t243 = t119 * pkin(1);
t21 = 0.1e1 / t22 ^ 2;
t242 = t21 * t23;
t214 = pkin(1) * t253;
t241 = 0.2e1 / t141 * (t87 + t88) * t214;
t235 = t116 * t48;
t234 = t122 * t91;
t233 = t129 / (t21 * t23 ^ 2 + 0.1e1);
t145 = (t94 * t241 / 0.2e1 + t134 * pkin(1) * t263 + (-t90 * t97 - t226) * pkin(5)) * t230 / 0.2e1;
t208 = -t241 / 0.2e1;
t225 = t97 * t141;
t150 = pkin(5) * (t225 + (t269 * t94 + t208 - t90) * t98) * t230;
t146 = t150 / 0.2e1;
t199 = 0.1e1 / t93 ^ 2 * t214;
t189 = t138 * t199;
t173 = t124 * t189;
t174 = t118 * t189;
t42 = t118 * t146 + t124 * t145 + t173 * t82 + t174 * t81;
t147 = -t150 / 0.2e1;
t43 = t118 * t145 + t124 * t147 - t173 * t81 + t174 * t82;
t37 = -t108 * t42 - t109 * t43;
t229 = t140 * t37;
t222 = t117 * t120;
t221 = t117 * t126;
t220 = t120 * t123;
t219 = t123 * t126;
t218 = t279 * t120;
t217 = t279 * t126;
t213 = pkin(3) * t257;
t212 = m(4) * t243;
t211 = t51 * t267;
t210 = t40 * t52 / 0.2e1;
t209 = t175 * t267;
t205 = t48 * t252;
t204 = -t235 / 0.2e1;
t203 = t235 / 0.2e1;
t202 = t91 * t250;
t200 = -0.2e1 * pkin(4) * t52 - t47;
t198 = -0.2e1 * t136 * t256;
t197 = t115 * t213;
t196 = t116 * t213;
t30 = 0.1e1 / t31 ^ 2;
t194 = pkin(8) * t131 / (t30 * t33 ^ 2 + 0.1e1) * t50;
t191 = t122 * t199;
t190 = t127 * t199;
t188 = t37 * t197;
t187 = t175 * t197;
t186 = t37 * t196;
t185 = t175 * t196;
t184 = 0.1e1 / t31 * t194;
t178 = mrSges(9,1) * t17 - mrSges(9,2) * t18;
t171 = pkin(4) * (t44 + t45) * t277;
t170 = pkin(3) * (t31 * t49 + t48 * t51);
t168 = pkin(4) * t30 * t33 * t194;
t25 = t37 * t171;
t36 = t108 * t43 - t109 * t42;
t166 = -t140 * t36 + t209 * t25;
t35 = t175 * t171;
t165 = -t140 * t56 + t209 * t35;
t164 = t240 * t189;
t163 = t239 * t189;
t155 = pkin(3) * (-t135 * t175 * t48 + t257 * t33);
t27 = -0.1e1 + ((t145 * t239 + t147 * t240 + t163 * t82 - t164 * t81) * t73 * t268 - (t145 * t240 + t146 * t239 + t163 * t81 + t164 * t82) / t72) / (t268 * t73 ^ 2 + 0.1e1);
t71 = 0.1e1 / t74 ^ 2;
t70 = t95 * t208 + t139 * pkin(5) * t263 + (-t89 * t97 - t226) * pkin(1);
t69 = (t225 + (0.2e1 * pkin(5) * t95 + t208 - t89) * t98) * pkin(1);
t28 = ((t70 * t202 + t83 * t190 - t69 * t234 / 0.2e1 - t80 * t191) / t74 - (t69 * t202 + t80 * t190 + t70 * t234 / 0.2e1 + t83 * t191) * t77 * t71) / (t71 * t77 ^ 2 + 0.1e1) * t133;
t20 = 0.1e1 / t22;
t16 = t35 * t210 + t175 * t198 + (t47 * t56 - t275) * pkin(3);
t15 = (t175 * t200 + t165) * pkin(3);
t10 = t12 * t219 - t222;
t9 = t12 * t221 + t220;
t8 = t12 * t220 + t221;
t7 = t12 * t222 - t219;
t6 = t25 * t210 + t37 * t198 + (t36 * t47 - t229) * pkin(3);
t5 = (t200 * t37 + t166) * pkin(3);
t4 = -0.2e1 * ((t35 * t211 + (t46 * t56 - t275) * pkin(4)) * t266 + t175 * t155) * t184 + 0.2e1 * ((-t175 * t46 + t165) * t266 + t175 * t170) * t168;
t3 = -0.2e1 * ((t25 * t211 + (t36 * t46 - t229) * pkin(4)) * t266 + t37 * t155) * t184 + 0.2e1 * ((-t37 * t46 + t166) * t266 + t37 * t170) * t168 + t27;
t2 = 0.1e1 + ((t15 * t205 + t16 * t203 + t185 * t34 + t187 * t32) * t20 - (t15 * t204 + t16 * t205 - t185 * t32 + t187 * t34) * t242) * t233;
t1 = 0.1e1 + ((t186 * t34 + t188 * t32 + t203 * t6 + t205 * t5) * t20 - (-t186 * t32 + t188 * t34 + t204 * t5 + t205 * t6) * t242) * t233;
t14 = [(t10 * mrSges(6,1) - t9 * mrSges(6,2) + t120 * t271 + t126 * t270) * g(2) + (-t8 * mrSges(6,1) + t7 * mrSges(6,2) - t120 * t270 + t126 * t271) * g(1), -g(1) * (-t126 * t212 + t217) - g(2) * (-t120 * t212 + t218) + (t274 * t1 + t282 * t112 - t179 * t3 - t181 * t28 - t276 * t195 + t27 * t281 + t280) * g(3) + t278 * (m(8) * t243 - (mrSges(8,1) * t60 + mrSges(8,2) * t61) * t27 - m(9) * (-pkin(3) * t27 * sin(t62) - t243) - t178 * t3 - (-mrSges(7,1) * t64 - mrSges(7,2) * t65) * t28 + mrSges(3,1) * t119 + mrSges(3,2) * t125 - t276 * (-t243 + t248) + t283 * t1), -g(1) * t217 - g(2) * t218 + (-t179 * t4 + t2 * t274 + t247 * t276 - t193) * g(3) + t278 * (-t178 * t4 + t2 * t283 - t276 * t248), -g(1) * (mrSges(6,1) * t9 + mrSges(6,2) * t10) - g(2) * (mrSges(6,1) * t7 + mrSges(6,2) * t8) - g(3) * (mrSges(6,1) * t117 + mrSges(6,2) * t123) * t11];
taug = t14(:);
