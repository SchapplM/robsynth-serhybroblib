% Calculate joint inertia matrix for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1turnOL_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_inertiaJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnOL_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fourbar1turnOL_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:49
% EndTime: 2020-04-12 19:40:50
% DurationCPUTime: 0.59s
% Computational Cost: add. (622->129), mult. (1004->216), div. (0->0), fcn. (890->8), ass. (0->74)
t70 = sin(qJ(1));
t73 = cos(qJ(1));
t117 = t70 * t73;
t69 = sin(qJ(2));
t72 = cos(qJ(2));
t116 = rSges(3,1) * t72 - rSges(3,2) * t69;
t65 = t70 ^ 2;
t66 = t73 ^ 2;
t115 = t70 / 0.2e1;
t114 = -t73 / 0.2e1;
t113 = pkin(2) * t72;
t67 = qJ(2) + qJ(3);
t61 = cos(t67);
t111 = rSges(4,1) * t61;
t71 = cos(qJ(4));
t110 = rSges(5,1) * t71;
t60 = sin(t67);
t108 = rSges(4,2) * t60;
t68 = sin(qJ(4));
t107 = rSges(5,2) * t68;
t106 = t73 * rSges(4,3);
t105 = t73 * rSges(5,3);
t104 = t70 * rSges(4,3) + t73 * t108;
t91 = t108 - t111;
t8 = t70 * (t91 * t70 - t106) + t73 * (-t73 * t111 + t104);
t103 = t70 * rSges(5,3) + t73 * t110;
t102 = t65 + t66;
t101 = Icges(3,4) * t69;
t100 = Icges(3,4) * t72;
t99 = Icges(4,4) * t60;
t98 = Icges(4,4) * t61;
t97 = Icges(5,4) * t68;
t96 = Icges(5,4) * t71;
t41 = -Icges(4,5) * t60 - Icges(4,6) * t61;
t78 = Icges(4,2) * t60 - t98;
t81 = -Icges(4,1) * t61 + t99;
t42 = -Icges(4,2) * t61 - t99;
t43 = -Icges(4,1) * t60 - t98;
t83 = t42 * t60 - t43 * t61;
t95 = (-t61 * (Icges(4,6) * t70 + t78 * t73) - t60 * (Icges(4,5) * t70 + t81 * t73) + t70 * t41 + t83 * t73) * t115 + (-t61 * (-Icges(4,6) * t73 + t78 * t70) - t60 * (-Icges(4,5) * t73 + t81 * t70) - t73 * t41 + t83 * t70) * t114;
t75 = -Icges(4,5) * t61 + Icges(4,6) * t60;
t17 = -Icges(4,3) * t73 + t75 * t70;
t18 = Icges(4,3) * t70 + t75 * t73;
t94 = -t73 * (-t18 * t117 + t66 * t17) + t70 * (-t17 * t117 + t65 * t18);
t44 = -t60 * rSges(4,1) - t61 * rSges(4,2);
t93 = -pkin(2) * t69 - t44;
t92 = t65 / 0.2e1 + t66 / 0.2e1;
t90 = -t107 + t110;
t82 = Icges(3,1) * t72 - t101;
t80 = Icges(5,1) * t71 - t97;
t79 = -Icges(3,2) * t69 + t100;
t77 = -Icges(5,2) * t68 + t96;
t76 = Icges(3,5) * t72 - Icges(3,6) * t69;
t74 = Icges(5,5) * t71 - Icges(5,6) * t68;
t56 = t73 * rSges(2,1) - t70 * rSges(2,2);
t55 = -t70 * rSges(2,1) - t73 * rSges(2,2);
t54 = t69 * rSges(3,1) + t72 * rSges(3,2);
t53 = t68 * rSges(5,1) + t71 * rSges(5,2);
t39 = t70 * rSges(3,3) + t116 * t73;
t37 = -t73 * rSges(3,3) + t116 * t70;
t28 = Icges(3,3) * t70 + t76 * t73;
t27 = -Icges(3,3) * t73 + t76 * t70;
t26 = Icges(5,3) * t70 + t74 * t73;
t25 = -Icges(5,3) * t73 + t74 * t70;
t24 = (pkin(1) - t107) * t73 + t103;
t23 = t105 + (-pkin(1) - t90) * t70;
t16 = t93 * t73;
t15 = t93 * t70;
t12 = (-t111 + t113) * t73 + t104;
t11 = t106 + (-t91 - t113) * t70;
t10 = t70 * t37 + t73 * t39;
t9 = t73 * (-t73 * t107 + t103) + (t90 * t70 - t105) * t70;
t3 = t102 * t113 + t8;
t1 = [-t61 * t42 - t60 * t43 + t71 * (Icges(5,2) * t71 + t97) + t72 * (Icges(3,2) * t72 + t101) + t68 * (Icges(5,1) * t68 + t96) + t69 * (Icges(3,1) * t69 + t100) + Icges(2,3) + m(2) * (t55 ^ 2 + t56 ^ 2) + m(3) * (t37 ^ 2 + t39 ^ 2) + m(4) * (t11 ^ 2 + t12 ^ 2) + m(5) * (t23 ^ 2 + t24 ^ 2); (t72 * (Icges(3,6) * t70 + t79 * t73) + t69 * (Icges(3,5) * t70 + t82 * t73)) * t115 + (t72 * (-Icges(3,6) * t73 + t79 * t70) + t69 * (-Icges(3,5) * t73 + t82 * t70)) * t114 + m(4) * (t16 * t11 + t15 * t12) + m(3) * (t37 * t73 - t39 * t70) * t54 + t92 * (Icges(3,5) * t69 + Icges(3,6) * t72) + t95; m(3) * (t102 * t54 ^ 2 + t10 ^ 2) + t70 * (-t27 * t117 + t65 * t28) - t73 * (-t28 * t117 + t66 * t27) + m(4) * (t15 ^ 2 + t16 ^ 2 + t3 ^ 2) + t94; m(4) * (-t11 * t73 - t12 * t70) * t44 + t95; m(4) * (t8 * t3 + (-t15 * t70 - t16 * t73) * t44) + t94; m(4) * (t102 * t44 ^ 2 + t8 ^ 2) + t94; (t71 * (Icges(5,6) * t70 + t77 * t73) + t68 * (Icges(5,5) * t70 + t80 * t73)) * t115 + (t71 * (-Icges(5,6) * t73 + t77 * t70) + t68 * (-Icges(5,5) * t73 + t80 * t70)) * t114 + m(5) * (-t23 * t73 - t24 * t70) * t53 + t92 * (Icges(5,5) * t68 + Icges(5,6) * t71); 0; 0; m(5) * (t102 * t53 ^ 2 + t9 ^ 2) + t70 * (-t25 * t117 + t65 * t26) - t73 * (-t26 * t117 + t66 * t25); 0; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
