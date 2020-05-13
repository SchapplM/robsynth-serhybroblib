% Calculate joint inertia matrix for
% palh3m2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% rSges [9x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [9x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 02:05
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh3m2DE1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE1_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE1_inertiaJ_slag_vp1: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE1_inertiaJ_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2DE1_inertiaJ_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m2DE1_inertiaJ_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:59:03
% EndTime: 2020-05-07 01:59:09
% DurationCPUTime: 5.00s
% Computational Cost: add. (677->287), mult. (920->325), div. (0->0), fcn. (232->90), ass. (0->146)
t156 = cos(qJ(3));
t162 = m(9) * rSges(9,1);
t164 = m(5) + m(6);
t81 = pkin(4) * t164 + m(4) * rSges(4,1) + t162;
t206 = pkin(1) * t81 * t156;
t152 = sin(qJ(3));
t246 = m(9) * rSges(9,2);
t252 = m(4) * rSges(4,2);
t119 = t246 + t252;
t238 = pkin(1) * t119;
t87 = t152 * t238;
t263 = t87 - t206;
t262 = 2 * pkin(12);
t151 = sin(qJ(4));
t155 = cos(qJ(4));
t261 = rSges(6,1) * t155 - rSges(6,2) * t151;
t143 = pkin(17) + pkin(18);
t132 = pkin(15) + t143;
t116 = (pkin(16) + t132);
t141 = m(5) + m(8) + m(9) + m(4);
t231 = m(6) + t141;
t188 = (pkin(4) ^ 2);
t260 = (t164 * t188);
t255 = pkin(3) * m(9);
t254 = pkin(8) * m(6);
t251 = m(5) * rSges(5,2);
t250 = m(6) * rSges(6,1);
t249 = m(6) * (rSges(6,1) * t151 + rSges(6,2) * t155);
t248 = m(7) * rSges(7,3);
t247 = m(8) * rSges(8,2);
t245 = m(9) * rSges(9,3);
t244 = rSges(3,2) * m(3);
t243 = rSges(6,2) * pkin(8);
t242 = rSges(6,2) * m(6);
t241 = rSges(4,3) * m(4);
t240 = rSges(5,3) * m(5);
t159 = -rSges(6,3) - pkin(10);
t177 = rSges(6,2) ^ 2;
t184 = rSges(6,1) ^ 2;
t233 = m(6) * (-t177 + t184);
t232 = t159 * m(6);
t228 = rSges(6,2) * t159;
t227 = t159 * rSges(6,1);
t124 = -pkin(4) * t156 + pkin(1);
t157 = cos(qJ(2));
t153 = sin(qJ(2));
t225 = t152 * t153;
t78 = pkin(4) * t225 + t124 * t157;
t190 = pkin(1) ^ 2;
t226 = t141 * t190;
t224 = t152 * t157;
t145 = qJ(3) + qJ(2);
t219 = t177 + t184;
t223 = t219 * m(6) + Icges(6,3);
t174 = (rSges(9,2) ^ 2);
t181 = (rSges(9,1) ^ 2);
t221 = t174 + t181;
t176 = rSges(7,2) ^ 2;
t183 = rSges(7,1) ^ 2;
t220 = (t176 + t183);
t179 = (rSges(4,2) ^ 2);
t186 = (rSges(4,1) ^ 2);
t218 = (t179 + t186);
t180 = rSges(3,2) ^ 2;
t187 = rSges(3,1) ^ 2;
t217 = t180 + t187;
t144 = pkin(18) + pkin(15);
t205 = -t250 / 0.2e1;
t204 = -t242 / 0.2e1;
t203 = t242 / 0.2e1;
t199 = pkin(4) * t250 / 0.2e1;
t197 = pkin(4) * t203;
t196 = 2 * t116;
t195 = -qJ(2) + t132;
t194 = qJ(2) + t132;
t193 = t261 * t254;
t192 = t218 * m(4) + t221 * m(9) + Icges(4,3) + Icges(9,3) + t260;
t110 = -qJ(2) + t116;
t109 = qJ(2) + t116;
t106 = -qJ(4) + t116;
t105 = qJ(4) + t116;
t102 = -qJ(2) + t106;
t101 = qJ(2) + t106;
t100 = -qJ(2) + t105;
t99 = qJ(2) + t105;
t133 = sin(t145);
t134 = cos(t145);
t137 = qJ(4) + t145;
t138 = -qJ(4) + t145;
t191 = t134 * (rSges(4,2) * t241 + rSges(9,2) * t245 - Icges(4,6) - Icges(9,6)) + (pkin(4) * t240 + rSges(4,1) * t241 + rSges(9,1) * t245 - Icges(4,5) - Icges(9,5)) * t133 + pkin(4) * cos(t137) * t205 + cos(t138) * t199 + (sin(t137) + sin(t138)) * t197;
t189 = pkin(3) ^ 2;
t185 = rSges(5,1) ^ 2;
t182 = rSges(8,1) ^ 2;
t178 = rSges(5,2) ^ 2;
t175 = rSges(8,2) ^ 2;
t173 = pkin(12) ^ 2;
t172 = 4 * Icges(6,5);
t171 = -2 * Icges(6,2);
t169 = 0.2e1 * qJ(2);
t168 = 0.2e1 * qJ(4);
t161 = rSges(6,1) * pkin(8);
t158 = cos(pkin(15));
t154 = sin(pkin(15));
t150 = cos(pkin(16));
t149 = sin(pkin(16));
t148 = qJ(2) - qJ(4);
t147 = qJ(2) + qJ(4);
t146 = t169 + qJ(3);
t139 = pkin(14) - qJ(2) - pkin(15);
t136 = -qJ(2) + t144;
t135 = qJ(2) + t144;
t131 = cos(t143);
t130 = sin(t143);
t128 = 0.2e1 * t145;
t123 = m(5) * rSges(5,1) + t254;
t122 = cos(t139);
t121 = sin(t139);
t120 = 0.2e1 * t144;
t117 = 0.4e1 * t227;
t114 = 0.2e1 * t139;
t113 = -qJ(3) + t195;
t112 = qJ(3) + t194;
t111 = -t232 - t251;
t108 = -qJ(4) + t196;
t107 = qJ(4) + t196;
t104 = -qJ(3) + t110;
t103 = qJ(3) + t109;
t98 = 2 * t116;
t93 = -qJ(3) + t102;
t92 = -qJ(3) + t100;
t91 = qJ(3) + t101;
t90 = qJ(3) + t99;
t89 = 0.2e1 * t106;
t88 = 0.2e1 * t105;
t80 = t153 * t156 + t224;
t79 = t156 * t157 - t225;
t77 = -pkin(4) * t224 + t124 * t153;
t76 = t193 + t223;
t75 = (m(6) * t227 + Icges(6,5)) * t155 - (t228 * m(6) + Icges(6,6)) * t151;
t73 = t154 * t79 + t158 * t80;
t72 = -t154 * t80 + t158 * t79;
t71 = -t154 * t77 + t158 * t78;
t70 = t154 * t78 + t158 * t77;
t69 = -t154 * t76 + t158 * t75;
t68 = t154 * t75 + t158 * t76;
t1 = [((t117 + 0.4e1 * t243) * m(6) + t172) * sin(t108) / 0.8e1 + ((t117 - 0.4e1 * t243) * m(6) + t172) * sin(t107) / 0.8e1 + (-rSges(3,1) * t244 + Icges(3,4)) * sin(t169) + (rSges(8,1) * t247 - Icges(8,4)) * sin(t120) + (pkin(8) * t232 + rSges(5,1) * t251 - Icges(5,4)) * sin(t98) + (-rSges(4,1) * t252 - rSges(9,1) * t246 + Icges(4,4) + Icges(9,4)) * sin(t128) + ((2 * rSges(4,3) ^ 2 + t218) * m(4)) / 0.2e1 - 0.2e1 * (cos(t132) * t255 + t123 * cos(t116) + t81 * t134 + t153 * t244 + (rSges(8,1) * cos(t144) + rSges(8,2) * sin(t144)) * m(8)) * pkin(12) + (-Icges(4,1) - Icges(9,1) + Icges(4,2) + Icges(9,2) + t260 + (-t174 + t181) * m(9) + (-t179 + t186) * m(4)) * cos(t128) / 0.2e1 + (0.4e1 * pkin(8) ^ 2 + 0.4e1 * t159 ^ 2 + 0.6e1 * t177 + 0.6e1 * t184 + (4 * t188) + 0.4e1 * t190) * m(6) / 0.8e1 + (-0.2e1 * rSges(7,1) * t122 - 0.2e1 * rSges(7,2) * t121 + pkin(6)) * m(7) * pkin(6) + (t231 * t157 * t262 - t81 * cos(t146) + (sin(t109) + sin(t110)) * t111 + (sin(t100) + sin(t99)) * t203 + (sin(t102) + sin(t101)) * t204 + (cos(t102) + cos(t101) + cos(t100) + cos(t99)) * t205 + (-sin(t135) - sin(t136)) * t247 + (-cos(t195) - cos(t194)) * t255 + (-cos(t109) - cos(t110)) * t123 + (-cos(t135) - cos(t136)) * m(8) * rSges(8,1)) * pkin(1) + ((2 * Icges(6,1)) + t171 - 0.2e1 * t233) * cos(t168) / 0.8e1 + t263 + (rSges(3,3) ^ 2 + t217 / 0.2e1 + t173) * m(3) + (m(3) * rSges(3,1) * t157 + t111 * sin(t116) + t119 * t133) * t262 + (sin(t105) - sin(t106)) * pkin(12) * t242 + (-cos(t105) - cos(t106)) * pkin(12) * t250 + (0.2e1 * (-0.2e1 * (pkin(8) - t159) * (-pkin(8) - t159) + t219) * m(6) + 0.4e1 * (-t178 + t185) * m(5) - (4 * Icges(5,1)) - (2 * Icges(6,1)) + (4 * Icges(5,2)) + t171 + 0.4e1 * Icges(6,3)) * cos(t98) / 0.8e1 + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + (m(7) * rSges(7,1) * rSges(7,2) - Icges(7,4)) * sin(t114) + (t190 * m(6) - Icges(3,1) + Icges(3,2) + t226 + (-t180 + t187) * m(3)) * cos(t169) / 0.2e1 + ((t161 - t228) * m(6) - Icges(6,6)) * cos(t108) / 0.2e1 + ((t161 + t228) * m(6) + Icges(6,6)) * cos(t107) / 0.2e1 + ((-t176 + t183) * m(7) - Icges(7,1) + Icges(7,2)) * cos(t114) / 0.2e1 + ((-t175 + t182) * m(8) - Icges(8,1) + Icges(8,2)) * cos(t120) / 0.2e1 + (sin(t89) / 0.4e1 - sin(t88) / 0.4e1 + sin(t168) / 0.2e1) * (rSges(6,1) * t242 - Icges(6,4)) + (cos(t89) + cos(t88)) * (-Icges(6,1) / 0.8e1 + Icges(6,2) / 0.8e1 + t233 / 0.8e1) + (t189 * cos(0.2e1 * t132) + (2 * rSges(9,3) ^ 2) + t189 + t221) * m(9) / 0.2e1 + ((2 * rSges(5,3) ^ 2) + t178 + t185 + t188) * m(5) / 0.2e1 + ((2 * rSges(8,3) ^ 2) + t175 + t182) * m(8) / 0.2e1 + (cos(t93) + cos(t92) + cos(t91) + cos(t90)) * t199 + sin(t146) * t238 + t193 + (2 * rSges(7,3) ^ 2 + t220) * m(7) / 0.2e1 + (sin(t93) + sin(t91)) * t197 + t226 / 0.2e1 + ((sin(t92) + sin(t90)) * t204 + (cos(t103) + cos(t104)) * t123 + (-sin(t103) - sin(t104)) * t111) * pkin(4) + ((sin(t113) - sin(t112)) * t246 + (cos(t112) + cos(t113)) * t162) * pkin(3) + Icges(3,1) / 0.2e1 + Icges(4,1) / 0.2e1 + Icges(5,1) / 0.2e1 + Icges(6,1) / 0.4e1 + Icges(7,1) / 0.2e1 + Icges(8,1) / 0.2e1 + Icges(9,1) / 0.2e1 + Icges(3,2) / 0.2e1 + Icges(4,2) / 0.2e1 + Icges(5,2) / 0.2e1 + Icges(6,2) / 0.4e1 + Icges(7,2) / 0.2e1 + Icges(8,2) / 0.2e1 + Icges(9,2) / 0.2e1 + Icges(2,3) + Icges(6,3) / 0.2e1 + t231 * t173; Icges(3,5) * t153 + t157 * Icges(3,6) + (-rSges(3,1) * t153 - rSges(3,2) * t157) * rSges(3,3) * m(3) + (-rSges(7,2) * t248 + Icges(7,6)) * t122 + (rSges(7,1) * t248 - Icges(7,5)) * t121 + t191 + ((-m(8) * rSges(8,3) - t240 - t241 - t245) * t153 + ((-sin(t148) / 0.2e1 - sin(t147) / 0.2e1) * rSges(6,2) + (-cos(t148) / 0.2e1 + cos(t147) / 0.2e1) * rSges(6,1)) * m(6)) * pkin(1); m(3) * t217 + m(7) * t220 + t190 * t231 + Icges(3,3) + Icges(7,3) + t192 - 0.2e1 * t206 + 0.2e1 * t87; t191; t192 + t263; t192; (t149 * t69 + t150 * t68) * t131 + (-t149 * t68 + t150 * t69) * t130 - t261 * (pkin(12) + t78) * m(6); -((t149 * t71 + t150 * t70) * t131 + (-t149 * t70 + t150 * t71) * t130) * t249; ((t149 * t72 + t150 * t73) * t131 + (-t149 * t73 + t150 * t72) * t130) * pkin(4) * t249; t223;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
