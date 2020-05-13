% Calculate joint inertia matrix for
% fourbar1turnDE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% m [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [2x2]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1turnDE1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_inertiaJ_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_inertiaJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE1_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnDE1_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnDE1_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:25:37
% EndTime: 2020-04-12 19:25:45
% DurationCPUTime: 2.36s
% Computational Cost: add. (25918->176), mult. (36892->330), div. (1996->13), fcn. (10758->10), ass. (0->111)
t160 = pkin(4) ^ 2;
t159 = pkin(3) ^ 2;
t95 = sin(qJ(1));
t97 = cos(qJ(1));
t158 = t95 * t97;
t92 = t95 ^ 2;
t93 = t97 ^ 2;
t139 = t92 + t93;
t102 = 0.1e1 / pkin(3);
t105 = pkin(1) ^ 2;
t96 = cos(qJ(2));
t144 = pkin(2) * t96;
t132 = -0.2e1 * pkin(1) * t144 + t105;
t147 = -pkin(3) - pkin(4);
t70 = (pkin(2) - t147) * (pkin(2) + t147) + t132;
t146 = -pkin(3) + pkin(4);
t71 = (pkin(2) - t146) * (pkin(2) + t146) + t132;
t106 = sqrt(-t70 * t71);
t94 = sin(qJ(2));
t134 = t106 * t94;
t131 = -t159 + t160;
t104 = pkin(2) ^ 2;
t83 = t104 + t132;
t72 = t83 - t131;
t86 = pkin(1) * t96 - pkin(2);
t51 = -pkin(1) * t134 - t72 * t86;
t107 = t51 ^ 2;
t69 = pkin(1) * t94 * t72;
t54 = -t106 * t86 + t69;
t50 = t54 ^ 2;
t78 = 0.1e1 / t83;
t79 = 0.1e1 / t83 ^ 2;
t126 = t102 * ((t107 + t50) * t79 / t159) ^ (-0.1e1 / 0.2e1) * t78;
t157 = (t51 * t94 + t54 * t96) * t126;
t73 = t83 + t131;
t85 = pkin(1) - t144;
t52 = -pkin(2) * t134 + t73 * t85;
t108 = t52 ^ 2;
t145 = pkin(2) * t94;
t68 = t73 * t145;
t53 = t106 * t85 + t68;
t49 = t53 ^ 2;
t99 = 0.1e1 / pkin(4);
t128 = ((t108 + t49) * t79 / t160) ^ (-0.1e1 / 0.2e1) * t78 * t99;
t156 = (rSges(5,1) * t52 + rSges(5,2) * t53) * t128;
t155 = -pkin(1) + t156;
t153 = -0.2e1 * pkin(2);
t129 = t79 * t145;
t133 = t96 * t106;
t130 = pkin(1) * t145;
t142 = 0.1e1 / t106 * (-t70 - t71) * t130;
t47 = 0.1e1 / t107;
t91 = t94 ^ 2;
t150 = (((0.2e1 * t105 * t91 * pkin(2) - t86 * t142) * t78 + ((t72 * t96 + t134) * t78 - 0.2e1 * t54 * t129) * pkin(1)) / t51 - (t69 * t78 + (-t78 * t133 + ((t86 * t153 - t142) * t78 + t51 * t79 * t153) * t94) * pkin(1)) * t54 * t47) * pkin(3) * t102 / (t47 * t50 + 0.1e1) * t83 + 0.1e1;
t149 = t95 / 0.2e1;
t148 = -t97 / 0.2e1;
t143 = pkin(2) * t97;
t141 = t94 * t95;
t140 = t95 * t96;
t138 = Icges(3,4) * t94;
t137 = Icges(3,4) * t96;
t136 = Icges(5,4) * t52;
t135 = Icges(5,4) * t53;
t37 = (-t51 * t96 + t54 * t94) * t126;
t34 = t97 * t37;
t35 = t97 * t157;
t127 = t34 * rSges(4,1) + t35 * rSges(4,2) + t95 * rSges(4,3);
t123 = t92 / 0.2e1 + t93 / 0.2e1;
t117 = Icges(3,1) * t96 - t138;
t116 = -Icges(3,2) * t94 + t137;
t115 = Icges(3,5) * t96 - Icges(3,6) * t94;
t32 = t95 * t37;
t33 = t95 * t157;
t113 = rSges(4,1) * t32 + rSges(4,2) * t33 - rSges(4,3) * t97;
t111 = (-Icges(5,1) * t52 - t135) * t128;
t110 = (-Icges(5,2) * t53 - t136) * t128;
t109 = (-Icges(5,5) * t52 - Icges(5,6) * t53) * t128;
t82 = rSges(2,1) * t97 - rSges(2,2) * t95;
t81 = -rSges(2,1) * t95 - rSges(2,2) * t97;
t80 = rSges(3,1) * t94 + rSges(3,2) * t96;
t67 = rSges(3,3) * t95 + (rSges(3,1) * t96 - rSges(3,2) * t94) * t97;
t65 = rSges(3,1) * t140 - rSges(3,2) * t141 - t97 * rSges(3,3);
t60 = Icges(3,3) * t95 + t115 * t97;
t59 = -Icges(3,3) * t97 + t115 * t95;
t55 = t65 * t95 + t67 * t97;
t48 = 0.1e1 / t108;
t42 = (rSges(5,1) * t53 - rSges(5,2) * t52) * t128;
t25 = Icges(5,3) * t95 + t97 * t109;
t24 = -Icges(5,3) * t97 + t95 * t109;
t23 = t95 * rSges(5,3) - t155 * t97;
t22 = rSges(5,3) * t97 + t155 * t95;
t21 = -rSges(4,1) * t157 + rSges(4,2) * t37;
t20 = -Icges(4,1) * t157 + Icges(4,4) * t37;
t19 = -Icges(4,4) * t157 + Icges(4,2) * t37;
t18 = -Icges(4,5) * t157 + Icges(4,6) * t37;
t17 = Icges(4,1) * t34 + Icges(4,4) * t35 + Icges(4,5) * t95;
t16 = Icges(4,1) * t32 + Icges(4,4) * t33 - Icges(4,5) * t97;
t15 = Icges(4,4) * t34 + Icges(4,2) * t35 + Icges(4,6) * t95;
t14 = Icges(4,4) * t32 + Icges(4,2) * t33 - Icges(4,6) * t97;
t13 = Icges(4,5) * t34 + Icges(4,6) * t35 + Icges(4,3) * t95;
t12 = Icges(4,5) * t32 + Icges(4,6) * t33 - Icges(4,3) * t97;
t11 = t96 * t143 + t127;
t10 = -pkin(2) * t140 - t113;
t8 = 0.2e1 * (-((t85 * t142 + (t73 * t96 + t134) * pkin(2)) * t78 / 0.2e1 + (t104 * t78 * t91 - t53 * t129) * pkin(1)) / t52 - (-(t68 + (-t94 * t142 - t133) * pkin(2)) * t78 / 0.2e1 + (t52 * t79 - t78 * t85) * t130) * t53 * t48) * pkin(4) / (t48 * t49 + 0.1e1) * t83 * t99;
t6 = t150 * t97;
t5 = t150 * t95;
t4 = -t94 * t143 - t21 * t6;
t3 = -pkin(2) * t141 - t21 * t5;
t2 = t139 * t8 * t156;
t1 = t5 * t113 + t6 * t127 + t139 * t144;
t7 = [t37 * t19 - t157 * t20 + t96 * (Icges(3,2) * t96 + t138) + t94 * (Icges(3,1) * t94 + t137) + Icges(2,3) + (-(-Icges(5,2) * t52 + t135) * t52 + (Icges(5,1) * t53 - t136) * t53) * t128 ^ 2 + m(5) * (t22 ^ 2 + t23 ^ 2) + m(3) * (t65 ^ 2 + t67 ^ 2) + m(4) * (t10 ^ 2 + t11 ^ 2) + m(2) * (t81 ^ 2 + t82 ^ 2); ((Icges(3,6) * t95 + t116 * t97) * t96 + (Icges(3,5) * t95 + t117 * t97) * t94) * t149 + ((-Icges(3,6) * t97 + t116 * t95) * t96 + (-Icges(3,5) * t97 + t117 * t95) * t94) * t148 + m(4) * (t10 * t4 + t11 * t3) + m(3) * (t65 * t97 - t67 * t95) * t80 + t123 * (Icges(3,5) * t94 + Icges(3,6) * t96) + (m(5) * (-t22 * t97 - t23 * t95) * t42 + (t123 * (Icges(5,5) * t53 - Icges(5,6) * t52) + (-(Icges(5,6) * t95 + t97 * t110) * t52 + (Icges(5,5) * t95 + t97 * t111) * t53) * t149 + (-(-Icges(5,6) * t97 + t95 * t110) * t52 + (-Icges(5,5) * t97 + t95 * t111) * t53) * t148) * t128) * t8 + (t15 * t37 - t157 * t17 + t18 * t95 + t19 * t35 + t20 * t34) * t5 / 0.2e1 - (t14 * t37 - t157 * t16 - t18 * t97 + t19 * t33 + t20 * t32) * t6 / 0.2e1; -t97 * (-t60 * t158 + t93 * t59) + m(4) * (t1 ^ 2 + t3 ^ 2 + t4 ^ 2) + t95 * (-t59 * t158 + t92 * t60) + t5 * ((t13 * t95 + t15 * t35 + t17 * t34) * t5 - (t12 * t95 + t14 * t35 + t16 * t34) * t6) + m(5) * t2 ^ 2 - t6 * ((-t13 * t97 + t15 * t33 + t17 * t32) * t5 - (-t12 * t97 + t14 * t33 + t16 * t32) * t6) + m(3) * (t139 * t80 ^ 2 + t55 ^ 2) + (-t97 * (-t25 * t158 + t93 * t24) + t95 * (-t24 * t158 + t92 * t25) + m(5) * t139 * t42 ^ 2) * t8 ^ 2;];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t7(1), t7(2); t7(2), t7(3);];
Mq = res;
