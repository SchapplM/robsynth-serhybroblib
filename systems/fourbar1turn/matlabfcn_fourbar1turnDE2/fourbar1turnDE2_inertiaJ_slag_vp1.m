% Calculate joint inertia matrix for
% fourbar1turnDE2
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
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1turnDE2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_inertiaJ_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_inertiaJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE2_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnDE2_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnDE2_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:33:32
% EndTime: 2020-04-12 19:33:39
% DurationCPUTime: 2.15s
% Computational Cost: add. (18522->168), mult. (25656->312), div. (1308->12), fcn. (7606->11), ass. (0->113)
t158 = pkin(4) ^ 2;
t89 = sin(qJ(1));
t91 = cos(qJ(1));
t157 = t89 * t91;
t86 = t89 ^ 2;
t87 = t91 ^ 2;
t132 = t86 + t87;
t88 = sin(qJ(2));
t90 = cos(qJ(2));
t144 = pkin(2) * t90;
t98 = pkin(1) ^ 2;
t133 = -0.2e1 * pkin(1) * t144 + t98;
t147 = -pkin(3) - pkin(4);
t64 = (pkin(2) - t147) * (pkin(2) + t147) + t133;
t146 = -pkin(3) + pkin(4);
t65 = (pkin(2) - t146) * (pkin(2) + t146) + t133;
t99 = sqrt(-t64 * t65);
t136 = t88 * t99;
t131 = pkin(3) ^ 2 - t158;
t97 = pkin(2) ^ 2;
t77 = t97 + t133;
t67 = t77 - t131;
t79 = pkin(1) - t144;
t46 = -pkin(2) * t136 + t67 * t79;
t100 = t46 ^ 2;
t145 = pkin(2) * t88;
t62 = t67 * t145;
t47 = t79 * t99 + t62;
t44 = t47 ^ 2;
t72 = 0.1e1 / t77;
t73 = 0.1e1 / t77 ^ 2;
t93 = 0.1e1 / pkin(4);
t122 = ((t100 + t44) / t158 * t73) ^ (-0.1e1 / 0.2e1) * t72 * t93;
t156 = (rSges(5,1) * t46 + rSges(5,2) * t47) * t122;
t155 = -pkin(1) + t156;
t153 = -0.2e1 * pkin(2);
t123 = t73 * t145;
t135 = t90 * t99;
t124 = pkin(1) * t145;
t139 = 0.1e1 / t99 * (-t64 - t65) * t124;
t66 = t77 + t131;
t80 = pkin(1) * t90 - pkin(2);
t45 = -pkin(1) * t136 - t66 * t80;
t42 = 0.1e1 / t45 ^ 2;
t63 = pkin(1) * t88 * t66;
t48 = -t80 * t99 + t63;
t85 = t88 ^ 2;
t96 = 0.1e1 / pkin(3);
t150 = (((0.2e1 * t98 * t85 * pkin(2) - t80 * t139) * t72 + ((t66 * t90 + t136) * t72 - 0.2e1 * t48 * t123) * pkin(1)) / t45 - (t63 * t72 + (-t72 * t135 + ((t153 * t80 - t139) * t72 + t45 * t73 * t153) * t88) * pkin(1)) * t48 * t42) * pkin(3) / (t42 * t48 ^ 2 + 0.1e1) * t77 * t96 + 0.1e1;
t149 = t89 / 0.2e1;
t148 = -t91 / 0.2e1;
t143 = rSges(3,1) * t90;
t138 = t72 * t96;
t38 = qJ(2) + atan2(t48 * t138, t45 * t138);
t37 = cos(t38);
t142 = rSges(4,1) * t37;
t36 = sin(t38);
t141 = rSges(4,2) * t36;
t140 = rSges(4,3) * t91;
t137 = t88 * t89;
t134 = t89 * rSges(4,3) + t141 * t91;
t130 = Icges(3,4) * t88;
t129 = Icges(3,4) * t90;
t128 = Icges(4,4) * t36;
t127 = Icges(4,4) * t37;
t126 = Icges(5,4) * t46;
t125 = Icges(5,4) * t47;
t119 = t86 / 0.2e1 + t87 / 0.2e1;
t106 = -Icges(4,5) * t37 + Icges(4,6) * t36;
t5 = t150 * t89;
t6 = t150 * t91;
t118 = -(-Icges(4,3) * t91 + t106 * t89) * t6 + (Icges(4,3) * t89 + t106 * t91) * t5;
t117 = t141 - t142;
t32 = -Icges(4,2) * t37 - t128;
t33 = -Icges(4,1) * t36 - t127;
t114 = t32 * t36 - t33 * t37;
t111 = Icges(3,1) * t90 - t130;
t110 = -Icges(4,1) * t37 + t128;
t109 = -Icges(3,2) * t88 + t129;
t108 = Icges(4,2) * t36 - t127;
t107 = Icges(3,5) * t90 - Icges(3,6) * t88;
t104 = (-Icges(5,1) * t46 - t125) * t122;
t103 = (-Icges(5,2) * t47 - t126) * t122;
t102 = (-Icges(5,5) * t46 - Icges(5,6) * t47) * t122;
t27 = -Icges(4,6) * t91 + t108 * t89;
t28 = Icges(4,6) * t89 + t108 * t91;
t29 = -Icges(4,5) * t91 + t110 * t89;
t30 = Icges(4,5) * t89 + t110 * t91;
t101 = (t28 * t36 - t30 * t37) * t5 - (t27 * t36 - t29 * t37) * t6;
t76 = rSges(2,1) * t91 - rSges(2,2) * t89;
t75 = -rSges(2,1) * t89 - rSges(2,2) * t91;
t74 = rSges(3,1) * t88 + rSges(3,2) * t90;
t61 = rSges(3,3) * t89 + (-rSges(3,2) * t88 + t143) * t91;
t59 = -rSges(3,2) * t137 - rSges(3,3) * t91 + t143 * t89;
t54 = Icges(3,3) * t89 + t107 * t91;
t53 = -Icges(3,3) * t91 + t107 * t89;
t49 = t59 * t89 + t61 * t91;
t43 = 0.1e1 / t100;
t34 = -rSges(4,1) * t36 - rSges(4,2) * t37;
t31 = -Icges(4,5) * t36 - Icges(4,6) * t37;
t24 = (-t142 + t144) * t91 + t134;
t23 = t140 + (-t117 - t144) * t89;
t22 = (rSges(5,1) * t47 - rSges(5,2) * t46) * t122;
t13 = Icges(5,3) * t89 + t102 * t91;
t12 = -Icges(5,3) * t91 + t102 * t89;
t11 = t89 * rSges(5,3) - t155 * t91;
t10 = rSges(5,3) * t91 + t155 * t89;
t8 = 0.2e1 * (-((t79 * t139 + (t67 * t90 + t136) * pkin(2)) * t72 / 0.2e1 + (t72 * t85 * t97 - t123 * t47) * pkin(1)) / t46 - (-(t62 + (-t139 * t88 - t135) * pkin(2)) * t72 / 0.2e1 + (t46 * t73 - t72 * t79) * t124) * t47 * t43) * pkin(4) / (t43 * t44 + 0.1e1) * t77 * t93;
t4 = -t145 * t91 - t34 * t6;
t3 = -pkin(2) * t137 - t34 * t5;
t2 = t5 * (t117 * t89 - t140) + t6 * (-t142 * t91 + t134) + t132 * t144;
t1 = t132 * t8 * t156;
t7 = [-t37 * t32 - t36 * t33 + t90 * (Icges(3,2) * t90 + t130) + t88 * (Icges(3,1) * t88 + t129) + Icges(2,3) + (-(-Icges(5,2) * t46 + t125) * t46 + (Icges(5,1) * t47 - t126) * t47) * t122 ^ 2 + m(5) * (t10 ^ 2 + t11 ^ 2) + m(3) * (t59 ^ 2 + t61 ^ 2) + m(4) * (t23 ^ 2 + t24 ^ 2) + m(2) * (t75 ^ 2 + t76 ^ 2); ((-Icges(3,6) * t91 + t109 * t89) * t90 + (-Icges(3,5) * t91 + t111 * t89) * t88) * t148 + ((Icges(3,6) * t89 + t109 * t91) * t90 + (Icges(3,5) * t89 + t111 * t91) * t88) * t149 + m(4) * (t23 * t4 + t24 * t3) + m(3) * (t59 * t91 - t61 * t89) * t74 + t119 * (Icges(3,5) * t88 + Icges(3,6) * t90) + (m(5) * (-t10 * t91 - t11 * t89) * t22 + (t119 * (Icges(5,5) * t47 - Icges(5,6) * t46) + (-(-Icges(5,6) * t91 + t103 * t89) * t46 + (-Icges(5,5) * t91 + t104 * t89) * t47) * t148 + (-(Icges(5,6) * t89 + t103 * t91) * t46 + (Icges(5,5) * t89 + t104 * t91) * t47) * t149) * t122) * t8 + (t114 * t91 - t28 * t37 - t30 * t36 + t31 * t89) * t5 / 0.2e1 - (t114 * t89 - t27 * t37 - t29 * t36 - t31 * t91) * t6 / 0.2e1; m(3) * (t132 * t74 ^ 2 + t49 ^ 2) + t89 * (-t53 * t157 + t86 * t54) - t91 * (-t54 * t157 + t87 * t53) + t5 * (t101 * t91 + t118 * t89) - t6 * (t101 * t89 - t118 * t91) + m(5) * t1 ^ 2 + m(4) * (t2 ^ 2 + t3 ^ 2 + t4 ^ 2) + (m(5) * t132 * t22 ^ 2 - t91 * (t87 * t12 - t13 * t157) + t89 * (-t12 * t157 + t86 * t13)) * t8 ^ 2;];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t7(1), t7(2); t7(2), t7(3);];
Mq = res;
