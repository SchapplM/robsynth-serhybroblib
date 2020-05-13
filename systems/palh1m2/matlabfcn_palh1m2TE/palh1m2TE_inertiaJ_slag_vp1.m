% Calculate joint inertia matrix for
% palh1m2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [11x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m2TE_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(22,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_inertiaJ_slag_vp1: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2TE_inertiaJ_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2TE_inertiaJ_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'palh1m2TE_inertiaJ_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:14:29
% EndTime: 2020-05-01 20:14:39
% DurationCPUTime: 7.66s
% Computational Cost: add. (772->346), mult. (1071->381), div. (0->0), fcn. (248->98), ass. (0->174)
t108 = m(5) + m(6);
t95 = m(11) * rSges(11,1);
t226 = pkin(5) * t108 + m(4) * rSges(4,1) + t95;
t97 = sin(qJ(3));
t163 = t226 * pkin(1) * t97;
t87 = qJ(3) + pkin(19);
t167 = pkin(2) * rSges(10,2) * cos(t87);
t101 = cos(qJ(3));
t188 = m(11) * rSges(11,2);
t198 = rSges(4,2) * m(4);
t52 = t188 + t198;
t217 = pkin(1) * t52;
t21 = t101 * t217;
t200 = rSges(10,1) * m(10);
t168 = pkin(2) * t200;
t39 = sin(t87) * t168;
t229 = -m(10) * t167 + t163 + t21 + t39;
t228 = -2 * pkin(15);
t86 = pkin(22) + pkin(21);
t70 = pkin(18) - t86;
t51 = (-pkin(20) + t70);
t79 = m(4) + m(5) + m(8) + m(11);
t191 = (m(6) + t79);
t137 = (pkin(5) ^ 2);
t227 = t108 * t137;
t222 = pkin(2) * m(10);
t221 = pkin(9) * m(6);
t214 = pkin(4) * m(11);
t211 = m(5) * rSges(5,2);
t210 = m(5) * rSges(5,3);
t209 = m(6) * rSges(6,1);
t100 = cos(qJ(4));
t96 = sin(qJ(4));
t208 = m(6) * (rSges(6,1) * t96 + rSges(6,2) * t100);
t125 = rSges(6,2) ^ 2;
t133 = rSges(6,1) ^ 2;
t207 = m(6) * (-t125 + t133);
t206 = (m(7) * rSges(7,3));
t205 = m(8) * rSges(8,2);
t203 = m(9) * rSges(9,2);
t202 = m(9) * rSges(9,3);
t201 = rSges(3,1) * m(3);
t199 = rSges(3,2) * m(3);
t197 = rSges(6,2) * pkin(9);
t196 = rSges(6,2) * m(6);
t195 = rSges(10,2) * m(10);
t194 = rSges(4,3) * m(4);
t104 = -rSges(6,3) - pkin(11);
t183 = t104 * rSges(6,1);
t186 = rSges(6,2) * t104;
t7 = (m(6) * t183 + Icges(6,5)) * t100 - (t186 * m(6) + Icges(6,6)) * t96;
t99 = sin(pkin(18));
t192 = t99 * t7;
t190 = rSges(6,2) * t96;
t189 = t104 * m(6);
t187 = rSges(6,1) * t100;
t185 = rSges(11,3) * m(11);
t98 = sin(qJ(2));
t184 = t101 * t98;
t140 = pkin(1) ^ 2;
t182 = t79 * t140;
t83 = t125 + t133;
t180 = t83 * m(6) + Icges(6,3);
t102 = cos(qJ(2));
t179 = t101 * t102;
t89 = qJ(3) + qJ(2);
t118 = rSges(11,2) ^ 2;
t119 = rSges(11,1) ^ 2;
t178 = t118 + t119;
t122 = (rSges(9,2) ^ 2);
t130 = (rSges(9,1) ^ 2);
t177 = (t122 + t130);
t124 = (rSges(7,2) ^ 2);
t132 = (rSges(7,1) ^ 2);
t176 = (t124 + t132);
t127 = rSges(4,2) ^ 2;
t135 = rSges(4,1) ^ 2;
t175 = (t127 + t135);
t128 = rSges(3,2) ^ 2;
t136 = rSges(3,1) ^ 2;
t174 = t128 + t136;
t88 = pkin(18) - pkin(22);
t166 = pkin(9) * t190;
t160 = -t209 / 0.2e1;
t159 = t209 / 0.2e1;
t158 = -t196 / 0.2e1;
t157 = t196 / 0.2e1;
t121 = rSges(10,2) ^ 2;
t129 = rSges(10,1) ^ 2;
t139 = pkin(2) ^ 2;
t156 = t121 + t129 + t139;
t151 = pkin(5) * t160;
t150 = pkin(5) * t158;
t148 = 2 * t51;
t62 = pkin(5) * t97 + pkin(1);
t10 = -pkin(5) * t179 + t62 * t98;
t147 = -qJ(2) + t70;
t146 = qJ(2) + t70;
t57 = t187 * t221;
t145 = -m(6) * t166 + t57;
t43 = -qJ(2) + t51;
t42 = qJ(2) + t51;
t41 = -qJ(4) + t51;
t40 = qJ(4) + t51;
t143 = t178 * m(11) + t175 * m(4) + t177 * m(9) + Icges(11,3) + Icges(4,3) + Icges(9,3) + t227;
t32 = -qJ(3) + t43;
t31 = qJ(3) + t42;
t109 = t139 * m(10);
t142 = t109 + t143;
t73 = qJ(2) + t87;
t55 = sin(t73);
t56 = cos(t73);
t77 = sin(t89);
t78 = cos(t89);
t80 = qJ(4) + t89;
t81 = -qJ(4) + t89;
t141 = -t78 * (rSges(11,2) * t185 + rSges(4,2) * t194 - Icges(11,6) - Icges(4,6)) + (-rSges(9,2) * t202 + Icges(9,6)) * t56 + pkin(5) * cos(t80) * t159 + cos(t81) * t151 + (-pkin(5) * t210 - rSges(11,1) * t185 - rSges(4,1) * t194 + Icges(11,5) + Icges(4,5)) * t77 + (-rSges(9,1) * t202 - rSges(10,3) * t222 + Icges(9,5)) * t55 + (sin(t80) + sin(t81)) * t150;
t138 = pkin(4) ^ 2;
t134 = (rSges(5,1) ^ 2);
t131 = (rSges(8,1) ^ 2);
t126 = (rSges(5,2) ^ 2);
t123 = (rSges(8,2) ^ 2);
t120 = pkin(15) ^ 2;
t117 = -0.4e1 * Icges(6,5);
t116 = -2 * Icges(6,2);
t114 = 0.2e1 * qJ(2);
t113 = 0.2e1 * qJ(4);
t106 = rSges(6,1) * pkin(9);
t103 = cos(pkin(18));
t94 = cos(pkin(20));
t93 = sin(pkin(20));
t92 = qJ(2) - qJ(4);
t91 = qJ(2) + qJ(4);
t90 = t114 + qJ(3);
t82 = pkin(17) + qJ(2) - pkin(18);
t76 = t114 + t87;
t75 = -qJ(2) + t88;
t74 = qJ(2) + t88;
t69 = cos(t86);
t68 = sin(t86);
t66 = 0.2e1 * t89;
t61 = m(5) * rSges(5,1) + t221;
t60 = cos(t82);
t59 = sin(t82);
t58 = 0.2e1 * t88;
t53 = 0.4e1 * t183;
t50 = 0.2e1 * t82;
t49 = 0.2e1 * t73;
t48 = t189 + t211;
t47 = -qJ(3) + t147;
t46 = qJ(3) + t146;
t45 = -qJ(4) + t148;
t44 = qJ(4) + t148;
t36 = -qJ(2) + t41;
t35 = -qJ(2) + t40;
t34 = qJ(2) + t41;
t33 = qJ(2) + t40;
t30 = 2 * t51;
t27 = -qJ(4) + t32;
t26 = qJ(4) + t32;
t25 = -qJ(4) + t31;
t24 = qJ(4) + t31;
t23 = 0.2e1 * t41;
t22 = 0.2e1 * t40;
t14 = t102 * t97 + t184;
t13 = -t97 * t98 + t179;
t11 = pkin(5) * t184 + t102 * t62;
t8 = t145 + t180;
t5 = t103 * t14 - t13 * t99;
t4 = t103 * t13 + t14 * t99;
t3 = t10 * t99 + t103 * t11;
t2 = t10 * t103 - t11 * t99;
t1 = t103 * t7 + t8 * t99;
t6 = [(t120 + rSges(3,3) ^ 2 + t174 / 0.2e1) * m(3) + ((sin(t46) - sin(t47)) * t188 + (-cos(t46) - cos(t47)) * t95) * pkin(4) + (sin(t27) + sin(t25)) * t150 + ((sin(t40) - sin(t41)) * rSges(6,2) + (-cos(t40) - cos(t41)) * rSges(6,1)) * pkin(15) * m(6) + ((m(9) + t191) * t120) + ((-cos(t31) - cos(t32)) * t61 + (sin(t31) + sin(t32)) * t48 + (sin(t26) + sin(t24)) * t157) * pkin(5) + ((sin(t74) - sin(t75)) * m(8) * rSges(8,1) + t191 * t98 * t228 - sin(t90) * t226 + (cos(t74) - cos(t75)) * t205 + (sin(t146) - sin(t147)) * t214 + (cos(t36) + cos(t33)) * t157 + (cos(t35) + cos(t34)) * t158 + (sin(t34) + sin(t33)) * t159 + (sin(t36) + sin(t35)) * t160 + (cos(t42) - cos(t43)) * t48 + (sin(t42) - sin(t43)) * t61) * pkin(1) - sin(t76) * t168 + t229 + ((2 * rSges(7,3) ^ 2 + t176) * m(7)) / 0.2e1 + (2 * rSges(8,3) ^ 2 + t123 + t131) * m(8) / 0.2e1 + (cos(t27) + cos(t26) + cos(t25) + cos(t24)) * t151 + (sin(t113) / 0.2e1 - sin(t22) / 0.4e1 + sin(t23) / 0.4e1) * (rSges(6,1) * t196 - Icges(6,4)) + (cos(t23) + cos(t22)) * (-Icges(6,1) / 0.8e1 + Icges(6,2) / 0.8e1 + t207 / 0.8e1) + (t138 * cos(0.2e1 * t70) + (2 * rSges(11,3) ^ 2) + t138 + t178) * m(11) / 0.2e1 + (2 * rSges(5,3) ^ 2 + t126 + t134 + t137) * m(5) / 0.2e1 + (-Icges(11,1) - Icges(4,1) + Icges(11,2) + Icges(4,2) + t227 + (-t118 + t119) * m(11) + (-t127 + t135) * m(4)) * cos(t66) / 0.2e1 + ((t106 + t186) * m(6) + Icges(6,6)) * cos(t45) / 0.2e1 + (t109 + ((-t122 + t130) * m(9)) - Icges(9,1) + Icges(9,2)) * cos(t49) / 0.2e1 + ((t124 - t132) * m(7) + Icges(7,1) - Icges(7,2)) * cos(t50) / 0.2e1 + ((-t123 + t131) * m(8) - Icges(8,1) + Icges(8,2)) * cos(t58) / 0.2e1 + (-t140 * m(6) + Icges(3,1) + Icges(10,1) - Icges(3,2) - Icges(10,2) - t182 + (t121 - t129) * m(10) + (t128 - t136) * m(3)) * cos(t114) / 0.2e1 + ((t106 - t186) * m(6) - Icges(6,6)) * cos(t44) / 0.2e1 + (0.4e1 * pkin(9) ^ 2 + 0.4e1 * t104 ^ 2 + 0.6e1 * t125 + 0.6e1 * t133 + (4 * t137) + 0.4e1 * t140) * m(6) / 0.8e1 + t145 + (t120 + rSges(10,3) ^ 2 + t156 / 0.2e1) * m(10) + 0.2e1 * ((m(9) * rSges(9,1) + t222) * t56 + t48 * sin(t51) + (-rSges(8,1) * cos(t88) + rSges(8,2) * sin(t88)) * m(8) + t78 * t226) * pkin(15) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + (m(7) * rSges(7,1) * rSges(7,2) - Icges(7,4)) * sin(t50) + t182 / 0.2e1 + ((t200 + t201) * t98 + (t195 + t199) * t102 + t55 * t203 + t52 * t77 + cos(t70) * t214 + t61 * cos(t51)) * t228 + (0.2e1 * rSges(7,1) * t59 + 0.2e1 * rSges(7,2) * t60 + pkin(14)) * m(7) * pkin(14) - cos(t90) * t217 + (-pkin(9) * t189 - rSges(5,1) * t211 + Icges(5,4)) * sin(t30) + (-rSges(8,1) * t205 + Icges(8,4)) * sin(t58) + (0.2e1 * Icges(6,1) + t116 - 0.2e1 * t207) * cos(t113) / 0.8e1 + (-rSges(9,1) * t203 + Icges(9,4)) * sin(t49) + (rSges(3,1) * t199 + rSges(10,1) * t195 - Icges(3,4) - Icges(10,4)) * sin(t114) + ((-t53 + 0.4e1 * t197) * m(6) + t117) * sin(t45) / 0.8e1 + ((-t53 - 0.4e1 * t197) * m(6) + t117) * sin(t44) / 0.8e1 + (-rSges(11,1) * t188 - rSges(4,1) * t198 + Icges(11,4) + Icges(4,4)) * sin(t66) + (0.2e1 * (-0.2e1 * (pkin(9) - t104) * (-pkin(9) - t104) + t83) * m(6) + 0.4e1 * (-t126 + t134) * m(5) - (4 * Icges(5,1)) - 0.2e1 * Icges(6,1) + (4 * Icges(5,2)) + t116 + 0.4e1 * Icges(6,3)) * cos(t30) / 0.8e1 - pkin(2) * cos(t76) * t195 + Icges(3,1) / 0.2e1 + Icges(4,1) / 0.2e1 + Icges(5,1) / 0.2e1 + Icges(6,1) / 0.4e1 + Icges(7,1) / 0.2e1 + Icges(8,1) / 0.2e1 + Icges(9,1) / 0.2e1 + Icges(10,1) / 0.2e1 + Icges(11,1) / 0.2e1 + Icges(3,2) / 0.2e1 + Icges(4,2) / 0.2e1 + Icges(5,2) / 0.2e1 + Icges(6,2) / 0.4e1 + Icges(7,2) / 0.2e1 + Icges(8,2) / 0.2e1 + Icges(9,2) / 0.2e1 + Icges(10,2) / 0.2e1 + Icges(11,2) / 0.2e1 + Icges(2,3) + Icges(6,3) / 0.2e1 + (2 * rSges(4,3) ^ 2 + t175) * m(4) / 0.2e1 + ((2 * rSges(9,3) ^ 2 + t177) * m(9)) / 0.2e1; (-rSges(3,3) * t201 - rSges(10,3) * t200 + Icges(3,5) + Icges(10,5)) * t102 + t98 * (rSges(3,3) * t199 + rSges(10,3) * t195 - Icges(3,6) - Icges(10,6)) + (-rSges(7,1) * t206 + Icges(7,5)) * t60 + (rSges(7,2) * t206 - Icges(7,6)) * t59 + t141 + ((-rSges(8,3) * m(8) - t185 - t194 - t210) * t102 + ((-cos(t92) / 0.2e1 - cos(t91) / 0.2e1) * rSges(6,2) + (sin(t92) / 0.2e1 - sin(t91) / 0.2e1) * rSges(6,1)) * m(6)) * pkin(1); 0.2e1 * t21 + t174 * m(3) + t191 * t140 + (t176 * m(7)) + 0.2e1 * t39 + (t156 - 0.2e1 * t167) * m(10) + 0.2e1 * t163 + t143 + Icges(3,3) + Icges(7,3) + Icges(10,3); t141; t142 + t229; t142; ((t103 * t8 - t192) * t94 + t93 * t1) * t69 + (t1 * t94 + t93 * ((-Icges(6,3) - t57 + (-t83 + t166) * m(6)) * t103 + t192)) * t68 + (-pkin(15) + t10) * m(6) * (t187 - t190); -((-t2 * t93 + t3 * t94) * t69 - (t2 * t94 + t3 * t93) * t68) * t208; -pkin(5) * ((t4 * t93 + t5 * t94) * t69 + t68 * (t4 * t94 - t5 * t93)) * t208; t180;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t6(1), t6(2), t6(4), t6(7); t6(2), t6(3), t6(5), t6(8); t6(4), t6(5), t6(6), t6(9); t6(7), t6(8), t6(9), t6(10);];
Mq = res;
