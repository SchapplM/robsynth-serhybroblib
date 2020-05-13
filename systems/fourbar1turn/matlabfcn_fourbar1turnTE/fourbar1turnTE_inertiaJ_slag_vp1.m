% Calculate joint inertia matrix for
% fourbar1turnTE
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
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1turnTE_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_inertiaJ_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_inertiaJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnTE_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnTE_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnTE_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:18:36
% EndTime: 2020-04-12 19:18:43
% DurationCPUTime: 1.94s
% Computational Cost: add. (14824->185), mult. (20638->354), div. (964->11), fcn. (6114->6), ass. (0->121)
t93 = sin(qJ(1));
t95 = cos(qJ(1));
t165 = t93 * t95;
t100 = pkin(2) ^ 2;
t101 = pkin(1) ^ 2;
t94 = cos(qJ(2));
t148 = pkin(2) * t94;
t128 = -0.2e1 * pkin(1) * t148 + t101;
t81 = t100 + t128;
t76 = 0.1e1 / t81;
t99 = 0.1e1 / pkin(3);
t139 = t76 * t99;
t92 = sin(qJ(2));
t154 = t92 / 0.2e1;
t151 = -pkin(3) - pkin(4);
t68 = (pkin(2) - t151) * (pkin(2) + t151) + t128;
t150 = -pkin(3) + pkin(4);
t69 = (pkin(2) - t150) * (pkin(2) + t150) + t128;
t102 = sqrt(-t68 * t69);
t130 = t102 * t92;
t134 = pkin(3) ^ 2 - pkin(4) ^ 2;
t70 = t81 + t134;
t84 = pkin(1) * t94 - pkin(2);
t49 = -pkin(1) * t130 - t84 * t70;
t146 = t92 * pkin(1);
t67 = t70 * t146;
t52 = -t102 * t84 + t67;
t164 = (t94 * t52 / 0.2e1 + t49 * t154) * t139;
t97 = 0.1e1 / pkin(4);
t140 = t76 * t97;
t149 = pkin(2) * t92;
t71 = t81 - t134;
t66 = t71 * t149;
t83 = pkin(1) - t148;
t51 = t102 * t83 + t66;
t142 = t51 * rSges(5,2);
t50 = -pkin(2) * t130 + t83 * t71;
t143 = t50 * rSges(5,1);
t163 = (t143 / 0.2e1 + t142 / 0.2e1) * t140;
t161 = -0.2e1 * pkin(2);
t90 = t93 ^ 2;
t91 = t95 ^ 2;
t77 = 0.1e1 / t81 ^ 2;
t126 = t77 * t149;
t129 = t94 * t102;
t127 = pkin(2) * t146;
t141 = 0.1e1 / t102 * (-t68 - t69) * t127;
t47 = 0.1e1 / t49 ^ 2;
t89 = t92 ^ 2;
t158 = (((0.2e1 * t101 * t89 * pkin(2) - t84 * t141) * t76 + ((t94 * t70 + t130) * t76 - 0.2e1 * t52 * t126) * pkin(1)) / t49 - (t67 * t76 + (-t76 * t129 + ((t161 * t84 - t141) * t76 + t49 * t77 * t161) * t92) * pkin(1)) * t52 * t47) * pkin(3) / (t52 ^ 2 * t47 + 0.1e1) * t81 * t99 + 0.1e1;
t157 = -t50 / 0.2e1;
t155 = t51 / 0.2e1;
t152 = -t95 / 0.2e1;
t147 = pkin(2) * t95;
t145 = -Icges(5,2) / 0.2e1;
t144 = -Icges(5,6) / 0.2e1;
t138 = t92 * t93;
t137 = t93 * t94;
t136 = t95 * rSges(5,3);
t135 = t90 + t91;
t133 = Icges(3,4) * t92;
t132 = Icges(3,4) * t94;
t131 = Icges(5,4) * t51;
t37 = (t52 * t154 - t94 * t49 / 0.2e1) * t139;
t34 = t95 * t37;
t35 = t95 * t164;
t125 = t34 * rSges(4,1) + t35 * rSges(4,2) + t93 * rSges(4,3);
t119 = t140 * t152;
t124 = t93 * rSges(5,3) + (t142 + t143) * t119;
t122 = Icges(5,4) * t157;
t74 = Icges(3,2) * t94 + t133;
t75 = Icges(3,1) * t92 + t132;
t116 = -t74 * t92 + t75 * t94;
t115 = Icges(3,1) * t94 - t133;
t114 = -Icges(3,2) * t92 + t132;
t113 = Icges(3,5) * t94 - Icges(3,6) * t92;
t32 = t93 * t37;
t33 = t93 * t164;
t109 = t32 * rSges(4,1) + t33 * rSges(4,2) - t95 * rSges(4,3);
t40 = (t131 / 0.2e1 + t50 * t145) * t140;
t41 = (Icges(5,1) * t155 + t122) * t140;
t107 = (t157 * t41 - t51 * t40 / 0.2e1) * t140;
t105 = (Icges(5,1) * t157 - t131 / 0.2e1) * t140;
t104 = (t145 * t51 + t122) * t140;
t103 = (Icges(5,5) * t157 + t144 * t51) * t140;
t80 = t95 * rSges(2,1) - t93 * rSges(2,2);
t79 = -t93 * rSges(2,1) - t95 * rSges(2,2);
t78 = t92 * rSges(3,1) + t94 * rSges(3,2);
t73 = Icges(3,5) * t92 + Icges(3,6) * t94;
t65 = t93 * rSges(3,3) + (rSges(3,1) * t94 - rSges(3,2) * t92) * t95;
t63 = rSges(3,1) * t137 - rSges(3,2) * t138 - t95 * rSges(3,3);
t58 = Icges(3,3) * t93 + t113 * t95;
t57 = -Icges(3,3) * t95 + t113 * t93;
t53 = t93 * t63 + t95 * t65;
t48 = 0.1e1 / t50 ^ 2;
t42 = (rSges(5,1) * t155 + rSges(5,2) * t157) * t140;
t39 = (Icges(5,5) * t155 + t144 * t50) * t140;
t25 = Icges(5,3) * t93 + t103 * t95;
t24 = -Icges(5,3) * t95 + t103 * t93;
t23 = t95 * pkin(1) + t124;
t22 = t136 + (-pkin(1) + t163) * t93;
t21 = -rSges(4,1) * t164 + rSges(4,2) * t37;
t20 = -Icges(4,1) * t164 + Icges(4,4) * t37;
t19 = -Icges(4,4) * t164 + Icges(4,2) * t37;
t18 = -Icges(4,5) * t164 + Icges(4,6) * t37;
t17 = Icges(4,1) * t34 + Icges(4,4) * t35 + Icges(4,5) * t93;
t16 = Icges(4,1) * t32 + Icges(4,4) * t33 - Icges(4,5) * t95;
t15 = Icges(4,4) * t34 + Icges(4,2) * t35 + Icges(4,6) * t93;
t14 = Icges(4,4) * t32 + Icges(4,2) * t33 - Icges(4,6) * t95;
t13 = Icges(4,5) * t34 + Icges(4,6) * t35 + Icges(4,3) * t93;
t12 = Icges(4,5) * t32 + Icges(4,6) * t33 - Icges(4,3) * t95;
t11 = t147 * t94 + t125;
t10 = -pkin(2) * t137 - t109;
t8 = 0.2e1 * (-((t83 * t141 + (t94 * t71 + t130) * pkin(2)) * t76 / 0.2e1 + (t100 * t89 * t76 - t126 * t51) * pkin(1)) / t50 - (-(t66 + (-t141 * t92 - t129) * pkin(2)) * t76 / 0.2e1 + (t50 * t77 - t76 * t83) * t127) * t51 * t48) * pkin(4) / (t51 ^ 2 * t48 + 0.1e1) * t81 * t97;
t6 = t158 * t95;
t5 = t158 * t93;
t4 = -t147 * t92 - t6 * t21;
t3 = -pkin(2) * t138 - t5 * t21;
t2 = (t95 * t124 + (-t93 * t163 - t136) * t93) * t8;
t1 = t109 * t5 + t125 * t6 + t135 * t148;
t7 = [t37 * t19 - t164 * t20 + t94 * t74 + t92 * t75 + Icges(2,3) + (t155 * t41 + t157 * t40) * t140 + m(3) * (t63 ^ 2 + t65 ^ 2) + m(4) * (t10 ^ 2 + t11 ^ 2) + m(5) * (t22 ^ 2 + t23 ^ 2) + m(2) * (t79 ^ 2 + t80 ^ 2); m(3) * (t63 * t95 - t65 * t93) * t78 + m(4) * (t10 * t4 + t11 * t3) + ((t155 * (-Icges(5,5) * t95 + t105 * t93) + t157 * (-Icges(5,6) * t95 + t104 * t93)) * t119 + m(5) * (-t22 * t95 - t23 * t93) * t42) * t8 + (t15 * t37 - t164 * t17 + t93 * t18 + t35 * t19 + t34 * t20) * t5 / 0.2e1 - (t14 * t37 - t16 * t164 - t95 * t18 + t33 * t19 + t32 * t20) * t6 / 0.2e1 + (t116 * t93 - t95 * t73 + t8 * (t107 * t93 - t95 * t39) + t94 * (-Icges(3,6) * t95 + t114 * t93) + t92 * (-Icges(3,5) * t95 + t115 * t93)) * t152 + (t116 * t95 + t94 * (Icges(3,6) * t93 + t114 * t95) + t92 * (Icges(3,5) * t93 + t115 * t95) + t93 * t73 + ((t155 * (Icges(5,5) * t93 + t105 * t95) + t157 * (Icges(5,6) * t93 + t104 * t95)) * t140 + t107 * t95 + t93 * t39) * t8) * t93 / 0.2e1; t93 * (-t57 * t165 + t90 * t58) + t5 * ((t93 * t13 + t35 * t15 + t34 * t17) * t5 - (t93 * t12 + t35 * t14 + t34 * t16) * t6) - t6 * ((-t95 * t13 + t33 * t15 + t32 * t17) * t5 - (-t95 * t12 + t33 * t14 + t32 * t16) * t6) - t95 * (-t58 * t165 + t91 * t57) + m(5) * t2 ^ 2 + m(3) * (t135 * t78 ^ 2 + t53 ^ 2) + m(4) * (t1 ^ 2 + t3 ^ 2 + t4 ^ 2) + (t93 * (-t24 * t165 + t90 * t25) - t95 * (-t25 * t165 + t91 * t24) + m(5) * t135 * t42 ^ 2) * t8 ^ 2;];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t7(1), t7(2); t7(2), t7(3);];
Mq = res;
