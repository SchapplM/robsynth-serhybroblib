% Calculate joint inertia matrix for
% picker2Dm1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
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
% Mq [12x12]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:46
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = picker2Dm1OL_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1OL_inertiaJ_slag_vp1: qJ has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1OL_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1OL_inertiaJ_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm1OL_inertiaJ_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'picker2Dm1OL_inertiaJ_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 05:44:56
% EndTime: 2020-05-11 05:44:57
% DurationCPUTime: 0.29s
% Computational Cost: add. (1181->149), mult. (704->182), div. (0->0), fcn. (422->20), ass. (0->87)
t91 = qJ(1) + qJ(2);
t84 = cos(t91);
t101 = pkin(2) * t84;
t100 = pkin(3) * t84;
t95 = cos(qJ(1));
t99 = t95 * pkin(1);
t87 = qJ(3) + t91;
t78 = qJ(9) + t87;
t63 = sin(t78);
t64 = cos(t78);
t28 = t63 * rSges(10,1) + t64 * rSges(10,2);
t86 = qJ(4) + t91;
t72 = sin(t86);
t75 = cos(t86);
t35 = t72 * rSges(5,1) + t75 * rSges(5,2);
t82 = sin(t91);
t43 = t82 * rSges(3,1) + t84 * rSges(3,2);
t98 = Icges(4,3) + Icges(10,3);
t97 = Icges(5,3) + Icges(11,3);
t69 = pkin(3) * t82;
t20 = t69 + t35;
t45 = -t84 * rSges(3,1) + t82 * rSges(3,2);
t38 = -t75 * rSges(5,1) + t72 * rSges(5,2);
t29 = -t64 * rSges(10,1) + t63 * rSges(10,2);
t73 = sin(t87);
t76 = cos(t87);
t39 = t76 * rSges(4,1) - t73 * rSges(4,2);
t85 = qJ(6) + t91;
t71 = sin(t85);
t74 = cos(t85);
t37 = t74 * rSges(7,1) - t71 * rSges(7,2);
t90 = qJ(1) + qJ(8);
t81 = sin(t90);
t83 = cos(t90);
t44 = t83 * rSges(9,1) - t81 * rSges(9,2);
t77 = qJ(10) + t86;
t61 = sin(t77);
t62 = cos(t77);
t27 = t62 * rSges(11,1) - t61 * rSges(11,2);
t18 = -pkin(6) * t73 + t28;
t19 = pkin(6) * t76 + t29;
t36 = -t73 * rSges(4,1) - t76 * rSges(4,2);
t34 = -t71 * rSges(7,1) - t74 * rSges(7,2);
t42 = -t81 * rSges(9,1) - t83 * rSges(9,2);
t70 = pkin(2) * t82;
t10 = t18 + t70;
t26 = -t61 * rSges(11,1) - t62 * rSges(11,2);
t23 = t39 - t101;
t21 = t38 - t100;
t13 = -pkin(4) * t75 + t27;
t22 = t36 + t70;
t12 = pkin(4) * t72 + t26;
t96 = Icges(3,3) + Icges(7,3) + t97 + t98;
t11 = t19 - t101;
t8 = t12 + t69;
t9 = t13 - t100;
t94 = cos(qJ(7));
t93 = sin(qJ(1));
t92 = sin(qJ(7));
t89 = pkin(8) + qJ(5);
t88 = t93 * pkin(1);
t80 = cos(t89);
t79 = sin(t89);
t49 = -t95 * rSges(2,1) + t93 * rSges(2,2);
t48 = t94 * rSges(8,1) - t92 * rSges(8,2);
t47 = t93 * rSges(2,1) + t95 * rSges(2,2);
t46 = t92 * rSges(8,1) + t94 * rSges(8,2);
t41 = t80 * rSges(6,1) - t79 * rSges(6,2);
t40 = -t79 * rSges(6,1) - t80 * rSges(6,2);
t33 = t45 - t99;
t32 = t44 - t99;
t31 = t88 + t43;
t30 = t42 + t88;
t25 = t37 - t99;
t24 = t34 + t88;
t17 = t23 - t99;
t16 = t21 - t99;
t15 = t22 + t88;
t14 = t88 + t20;
t7 = t11 - t99;
t6 = t10 + t88;
t5 = t9 - t99;
t4 = t8 + t88;
t3 = m(7) * (t34 ^ 2 + t37 ^ 2);
t2 = Icges(7,3) + t3;
t1 = m(7) * (t34 * t24 + t37 * t25);
t50 = [Icges(2,3) + Icges(9,3) + m(5) * (t14 ^ 2 + t16 ^ 2) + m(4) * (t15 ^ 2 + t17 ^ 2) + m(3) * (t31 ^ 2 + t33 ^ 2) + m(2) * (t47 ^ 2 + t49 ^ 2) + m(11) * (t4 ^ 2 + t5 ^ 2) + m(9) * (t30 ^ 2 + t32 ^ 2) + m(10) * (t6 ^ 2 + t7 ^ 2) + m(7) * (t24 ^ 2 + t25 ^ 2) + t96; t1 + m(5) * (t20 * t14 + t21 * t16) + m(4) * (t22 * t15 + t23 * t17) + m(3) * (t43 * t31 + t45 * t33) + m(11) * (t8 * t4 + t9 * t5) + m(10) * (t10 * t6 + t11 * t7) + t96; t3 + m(3) * (t43 ^ 2 + t45 ^ 2) + m(10) * (t10 ^ 2 + t11 ^ 2) + m(11) * (t8 ^ 2 + t9 ^ 2) + m(4) * (t22 ^ 2 + t23 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + t96; m(4) * (t36 * t15 + t39 * t17) + m(10) * (t18 * t6 + t19 * t7) + t98; m(10) * (t18 * t10 + t19 * t11) + m(4) * (t36 * t22 + t39 * t23) + t98; m(10) * (t18 ^ 2 + t19 ^ 2) + m(4) * (t36 ^ 2 + t39 ^ 2) + t98; m(5) * (t35 * t14 + t38 * t16) + m(11) * (t12 * t4 + t13 * t5) + t97; m(11) * (t12 * t8 + t13 * t9) + m(5) * (t35 * t20 + t38 * t21) + t97; 0; m(11) * (t12 ^ 2 + t13 ^ 2) + m(5) * (t35 ^ 2 + t38 ^ 2) + t97; 0; 0; 0; 0; Icges(6,3) + m(6) * (t40 ^ 2 + t41 ^ 2); Icges(7,3) + t1; t2; 0; 0; 0; t2; 0; 0; 0; 0; 0; 0; Icges(8,3) + m(8) * (t46 ^ 2 + t48 ^ 2); m(9) * (t42 * t30 + t44 * t32) + Icges(9,3); 0; 0; 0; 0; 0; 0; Icges(9,3) + m(9) * (t42 ^ 2 + t44 ^ 2); Icges(10,3) + m(10) * (t28 * t6 + t29 * t7); m(10) * (t28 * t10 + t29 * t11) + Icges(10,3); Icges(10,3) + m(10) * (t28 * t18 + t29 * t19); 0; 0; 0; 0; 0; Icges(10,3) + m(10) * (t28 ^ 2 + t29 ^ 2); m(11) * (t26 * t4 + t27 * t5) + Icges(11,3); m(11) * (t26 * t8 + t27 * t9) + Icges(11,3); 0; m(11) * (t26 * t12 + t27 * t13) + Icges(11,3); 0; 0; 0; 0; 0; m(11) * (t26 ^ 2 + t27 ^ 2) + Icges(11,3); 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_12_matlab.m
res = [t50(1), t50(2), t50(4), t50(7), t50(11), t50(16), t50(22), t50(29), t50(37), t50(46), t50(56), t50(67); t50(2), t50(3), t50(5), t50(8), t50(12), t50(17), t50(23), t50(30), t50(38), t50(47), t50(57), t50(68); t50(4), t50(5), t50(6), t50(9), t50(13), t50(18), t50(24), t50(31), t50(39), t50(48), t50(58), t50(69); t50(7), t50(8), t50(9), t50(10), t50(14), t50(19), t50(25), t50(32), t50(40), t50(49), t50(59), t50(70); t50(11), t50(12), t50(13), t50(14), t50(15), t50(20), t50(26), t50(33), t50(41), t50(50), t50(60), t50(71); t50(16), t50(17), t50(18), t50(19), t50(20), t50(21), t50(27), t50(34), t50(42), t50(51), t50(61), t50(72); t50(22), t50(23), t50(24), t50(25), t50(26), t50(27), t50(28), t50(35), t50(43), t50(52), t50(62), t50(73); t50(29), t50(30), t50(31), t50(32), t50(33), t50(34), t50(35), t50(36), t50(44), t50(53), t50(63), t50(74); t50(37), t50(38), t50(39), t50(40), t50(41), t50(42), t50(43), t50(44), t50(45), t50(54), t50(64), t50(75); t50(46), t50(47), t50(48), t50(49), t50(50), t50(51), t50(52), t50(53), t50(54), t50(55), t50(65), t50(76); t50(56), t50(57), t50(58), t50(59), t50(60), t50(61), t50(62), t50(63), t50(64), t50(65), t50(66), t50(77); t50(67), t50(68), t50(69), t50(70), t50(71), t50(72), t50(73), t50(74), t50(75), t50(76), t50(77), t50(78);];
Mq = res;
