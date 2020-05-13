% Calculate time derivative of joint inertia matrix for
% picker2Dm2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% qJD [12x1]
%   Generalized joint velocities
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
% MqD [12x12]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 23:20
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = picker2Dm2OL_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2OL_inertiaDJ_slag_vp1: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm2OL_inertiaDJ_slag_vp1: qJD has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2OL_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2OL_inertiaDJ_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm2OL_inertiaDJ_slag_vp1: rSges has to be [11x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [11 6]), ...
  'picker2Dm2OL_inertiaDJ_slag_vp1: Icges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-09 23:18:54
% EndTime: 2020-05-09 23:18:58
% DurationCPUTime: 0.69s
% Computational Cost: add. (3018->192), mult. (1636->261), div. (0->0), fcn. (844->16), ass. (0->136)
t152 = qJ(1) + qJ(2);
t145 = qJ(4) + t152;
t131 = qJ(10) + t145;
t113 = sin(t131);
t114 = cos(t131);
t66 = t114 * rSges(11,1) - t113 * rSges(11,2);
t144 = qJ(6) + t152;
t125 = sin(t144);
t128 = cos(t144);
t76 = t128 * rSges(7,1) - t125 * rSges(7,2);
t151 = qJ(1) + qJ(8);
t140 = sin(t151);
t142 = cos(t151);
t85 = t142 * rSges(9,1) - t140 * rSges(9,2);
t173 = 2 * m(10);
t172 = 2 * m(11);
t143 = cos(t152);
t171 = pkin(2) * t143;
t170 = pkin(3) * t143;
t156 = cos(qJ(1));
t169 = t156 * pkin(1);
t146 = qJ(3) + t152;
t132 = qJ(9) + t146;
t116 = cos(t132);
t149 = qJD(1) + qJD(2);
t137 = qJD(3) + t149;
t118 = qJD(9) + t137;
t163 = t116 * t118;
t115 = sin(t132);
t164 = t115 * t118;
t35 = rSges(10,1) * t164 + rSges(10,2) * t163;
t129 = cos(t145);
t136 = qJD(4) + t149;
t160 = t129 * t136;
t126 = sin(t145);
t162 = t126 * t136;
t51 = rSges(5,1) * t162 + rSges(5,2) * t160;
t157 = t143 * t149;
t141 = sin(t152);
t158 = t141 * t149;
t61 = rSges(3,1) * t158 + rSges(3,2) * t157;
t168 = pkin(1) * qJD(1);
t127 = sin(t146);
t161 = t127 * t137;
t130 = cos(t146);
t159 = t130 * t137;
t67 = t115 * rSges(10,1) + t116 * rSges(10,2);
t74 = t126 * rSges(5,1) + t129 * rSges(5,2);
t84 = t141 * rSges(3,1) + t143 * rSges(3,2);
t98 = pkin(3) * t158;
t23 = t98 + t51;
t123 = pkin(3) * t141;
t55 = t123 + t74;
t86 = -t143 * rSges(3,1) + t141 * rSges(3,2);
t77 = -t129 * rSges(5,1) + t126 * rSges(5,2);
t68 = -t116 * rSges(10,1) + t115 * rSges(10,2);
t78 = t130 * rSges(4,1) - t127 * rSges(4,2);
t54 = -rSges(4,1) * t159 + rSges(4,2) * t161;
t135 = qJD(6) + t149;
t50 = t76 * t135;
t148 = qJD(1) + qJD(8);
t60 = t85 * t148;
t62 = rSges(3,1) * t157 - rSges(3,2) * t158;
t52 = rSges(5,1) * t160 - rSges(5,2) * t162;
t36 = rSges(10,1) * t163 - rSges(10,2) * t164;
t47 = -pkin(6) * t127 + t67;
t117 = qJD(10) + t136;
t30 = t66 * t117;
t48 = pkin(6) * t130 + t68;
t75 = -t127 * rSges(4,1) - t130 * rSges(4,2);
t73 = -t125 * rSges(7,1) - t128 * rSges(7,2);
t83 = -t140 * rSges(9,1) - t142 * rSges(9,2);
t13 = -pkin(6) * t161 + t35;
t65 = -t113 * rSges(11,1) - t114 * rSges(11,2);
t99 = pkin(3) * t157;
t24 = t52 + t99;
t12 = pkin(4) * t160 - t30;
t101 = pkin(2) * t157;
t26 = t101 + t54;
t124 = pkin(2) * t141;
t31 = t124 + t47;
t58 = t78 - t171;
t56 = t77 - t170;
t38 = -pkin(4) * t129 + t66;
t57 = t124 + t75;
t100 = pkin(2) * t158;
t9 = t100 + t13;
t37 = pkin(4) * t126 + t65;
t8 = t12 + t99;
t53 = t75 * t137;
t49 = t73 * t135;
t59 = t83 * t148;
t29 = t65 * t117;
t32 = t48 - t171;
t27 = t123 + t37;
t14 = -pkin(6) * t159 + t36;
t25 = t100 + t53;
t11 = pkin(4) * t162 + t29;
t28 = t38 - t170;
t10 = t101 + t14;
t7 = t98 + t11;
t154 = sin(qJ(1));
t147 = t154 * pkin(1);
t134 = t156 * t168;
t133 = t154 * t168;
t72 = t86 - t169;
t71 = t85 - t169;
t70 = t147 + t84;
t69 = t147 + t83;
t64 = t76 - t169;
t63 = t147 + t73;
t46 = t134 + t62;
t45 = t134 - t60;
t44 = t133 + t61;
t43 = t133 + t59;
t42 = t58 - t169;
t41 = t56 - t169;
t40 = t147 + t57;
t39 = t147 + t55;
t34 = t134 - t50;
t33 = t133 + t49;
t22 = t32 - t169;
t21 = t147 + t31;
t20 = t28 - t169;
t19 = t147 + t27;
t18 = t134 + t26;
t17 = t134 + t24;
t16 = t133 + t25;
t15 = t133 + t23;
t6 = t10 + t134;
t5 = t133 + t9;
t4 = t134 + t8;
t3 = t133 + t7;
t2 = 0.2e1 * m(7) * (t76 * t49 - t73 * t50);
t1 = m(7) * (t76 * t33 + t73 * t34 + t49 * t64 - t50 * t63);
t79 = [0.2e1 * m(5) * (t41 * t15 + t39 * t17) + 0.2e1 * m(4) * (t42 * t16 + t40 * t18) + 0.2e1 * m(3) * (t72 * t44 + t70 * t46) + 0.2e1 * m(11) * (t19 * t4 + t20 * t3) + 0.2e1 * m(9) * (t71 * t43 + t69 * t45) + 0.2e1 * m(10) * (t21 * t6 + t22 * t5) + 0.2e1 * m(7) * (t64 * t33 + t63 * t34); t1 + m(5) * (t56 * t15 + t55 * t17 + t23 * t41 + t24 * t39) + m(4) * (t58 * t16 + t57 * t18 + t25 * t42 + t26 * t40) + m(3) * (t86 * t44 + t84 * t46 + t61 * t72 + t62 * t70) + m(11) * (t8 * t19 + t7 * t20 + t27 * t4 + t28 * t3) + m(10) * (t10 * t21 + t9 * t22 + t31 * t6 + t32 * t5); t2 + 0.2e1 * m(3) * (t86 * t61 + t84 * t62) + (t27 * t8 + t28 * t7) * t172 + (t31 * t10 + t32 * t9) * t173 + 0.2e1 * m(5) * (t56 * t23 + t55 * t24) + 0.2e1 * m(4) * (t58 * t25 + t57 * t26); m(4) * (t78 * t16 + t75 * t18 + t54 * t40 + t53 * t42) + m(10) * (t13 * t22 + t14 * t21 + t47 * t6 + t48 * t5); m(10) * (t47 * t10 + t13 * t32 + t14 * t31 + t48 * t9) + m(4) * (t78 * t25 + t75 * t26 + t53 * t58 + t54 * t57); 0.2e1 * m(10) * (t48 * t13 + t47 * t14) + 0.2e1 * m(4) * (t78 * t53 + t75 * t54); m(5) * (t77 * t15 + t74 * t17 + t52 * t39 + t51 * t41) + m(11) * (t11 * t20 + t12 * t19 + t38 * t3 + t37 * t4); m(11) * (t11 * t28 + t12 * t27 + t37 * t8 + t38 * t7) + m(5) * (t77 * t23 + t74 * t24 + t51 * t56 + t52 * t55); 0; 0.2e1 * m(11) * (t38 * t11 + t37 * t12) + 0.2e1 * m(5) * (t77 * t51 + t74 * t52); 0; 0; 0; 0; 0; t1; t2; 0; 0; 0; t2; 0; 0; 0; 0; 0; 0; 0; m(9) * (t85 * t43 + t83 * t45 + t59 * t71 - t60 * t69); 0; 0; 0; 0; 0; 0; 0.2e1 * m(9) * (t85 * t59 - t83 * t60); m(10) * (t36 * t21 + t35 * t22 + t68 * t5 + t67 * t6); m(10) * (t67 * t10 + t36 * t31 + t35 * t32 + t68 * t9); m(10) * (t68 * t13 + t67 * t14 + t35 * t48 + t36 * t47); 0; 0; 0; 0; 0; (t68 * t35 + t67 * t36) * t173; m(11) * (-t30 * t19 + t29 * t20 + t66 * t3 + t65 * t4); m(11) * (-t30 * t27 + t29 * t28 + t65 * t8 + t66 * t7); 0; m(11) * (t66 * t11 + t65 * t12 + t29 * t38 - t30 * t37); 0; 0; 0; 0; 0; (t66 * t29 - t65 * t30) * t172; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_12_matlab.m
res = [t79(1), t79(2), t79(4), t79(7), t79(11), t79(16), t79(22), t79(29), t79(37), t79(46), t79(56), t79(67); t79(2), t79(3), t79(5), t79(8), t79(12), t79(17), t79(23), t79(30), t79(38), t79(47), t79(57), t79(68); t79(4), t79(5), t79(6), t79(9), t79(13), t79(18), t79(24), t79(31), t79(39), t79(48), t79(58), t79(69); t79(7), t79(8), t79(9), t79(10), t79(14), t79(19), t79(25), t79(32), t79(40), t79(49), t79(59), t79(70); t79(11), t79(12), t79(13), t79(14), t79(15), t79(20), t79(26), t79(33), t79(41), t79(50), t79(60), t79(71); t79(16), t79(17), t79(18), t79(19), t79(20), t79(21), t79(27), t79(34), t79(42), t79(51), t79(61), t79(72); t79(22), t79(23), t79(24), t79(25), t79(26), t79(27), t79(28), t79(35), t79(43), t79(52), t79(62), t79(73); t79(29), t79(30), t79(31), t79(32), t79(33), t79(34), t79(35), t79(36), t79(44), t79(53), t79(63), t79(74); t79(37), t79(38), t79(39), t79(40), t79(41), t79(42), t79(43), t79(44), t79(45), t79(54), t79(64), t79(75); t79(46), t79(47), t79(48), t79(49), t79(50), t79(51), t79(52), t79(53), t79(54), t79(55), t79(65), t79(76); t79(56), t79(57), t79(58), t79(59), t79(60), t79(61), t79(62), t79(63), t79(64), t79(65), t79(66), t79(77); t79(67), t79(68), t79(69), t79(70), t79(71), t79(72), t79(73), t79(74), t79(75), t79(76), t79(77), t79(78);];
Mq = res;
