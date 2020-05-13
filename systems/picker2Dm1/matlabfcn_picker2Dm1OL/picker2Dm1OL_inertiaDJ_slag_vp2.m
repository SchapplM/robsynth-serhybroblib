% Calculate time derivative of joint inertia matrix for
% picker2Dm1OL
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
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [11x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [12x12]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:46
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = picker2Dm1OL_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1OL_inertiaDJ_slag_vp2: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm1OL_inertiaDJ_slag_vp2: qJD has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1OL_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1OL_inertiaDJ_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm1OL_inertiaDJ_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'picker2Dm1OL_inertiaDJ_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 05:44:58
% EndTime: 2020-05-11 05:45:00
% DurationCPUTime: 0.64s
% Computational Cost: add. (893->139), mult. (2286->230), div. (0->0), fcn. (1500->14), ass. (0->88)
t63 = cos(qJ(3));
t83 = qJD(3) * t63;
t57 = sin(qJ(3));
t84 = qJD(3) * t57;
t109 = (mrSges(4,1) * t84 + mrSges(4,2) * t83) * pkin(2);
t108 = (cos(qJ(8)) * mrSges(9,2) + sin(qJ(8)) * mrSges(9,1)) * pkin(1) * qJD(8);
t62 = cos(qJ(4));
t46 = t62 * pkin(3) + pkin(4);
t51 = sin(qJ(10));
t52 = cos(qJ(10));
t56 = sin(qJ(4));
t76 = qJD(10) * t56;
t78 = qJD(10) * t51;
t98 = t52 * t56;
t16 = t46 * t78 + (t52 * t76 + (t51 * t62 + t98) * qJD(4)) * pkin(3);
t11 = t16 * mrSges(11,1);
t47 = -t63 * pkin(2) + pkin(6);
t53 = sin(qJ(9));
t59 = cos(qJ(9));
t79 = qJD(9) * t59;
t80 = qJD(9) * t53;
t93 = t57 * t59;
t18 = t47 * t80 + (-t57 * t79 + (-t53 * t63 - t93) * qJD(3)) * pkin(2);
t12 = t18 * mrSges(10,1);
t107 = t11 + t12 + t109;
t106 = pkin(1) * (qJD(2) + qJD(6));
t58 = sin(qJ(2));
t64 = cos(qJ(2));
t50 = t64 * pkin(1);
t48 = t50 + pkin(3);
t81 = qJD(4) * t62;
t82 = qJD(4) * t56;
t95 = t56 * t58;
t19 = t48 * t81 + (-t58 * t82 + (t62 * t64 - t95) * qJD(2)) * pkin(1);
t91 = t58 * t62;
t20 = -t48 * t82 + (-t58 * t81 + (-t56 * t64 - t91) * qJD(2)) * pkin(1);
t31 = -pkin(1) * t95 + t62 * t48;
t26 = pkin(4) + t31;
t34 = pkin(1) * t91 + t56 * t48;
t70 = t51 * t26 + t52 * t34;
t4 = t70 * qJD(10) + t51 * t19 - t52 * t20;
t1 = t4 * mrSges(11,1);
t71 = t20 * mrSges(5,1) - t19 * mrSges(5,2) + t1;
t49 = t50 + pkin(2);
t94 = t57 * t58;
t21 = -t49 * t83 + (t58 * t84 + (-t63 * t64 + t94) * qJD(2)) * pkin(1);
t90 = t58 * t63;
t22 = t49 * t84 + (t58 * t83 + (t57 * t64 + t90) * qJD(2)) * pkin(1);
t32 = pkin(1) * t94 - t63 * t49;
t27 = pkin(6) + t32;
t35 = -pkin(1) * t90 - t57 * t49;
t69 = t53 * t27 + t59 * t35;
t6 = t69 * qJD(9) + t53 * t21 - t59 * t22;
t2 = t6 * mrSges(10,1);
t72 = t22 * mrSges(4,1) - t21 * mrSges(4,2) + t2;
t55 = sin(qJ(6));
t61 = cos(qJ(6));
t67 = t55 * t58 - t61 * t64;
t24 = t67 * t106;
t68 = t55 * t64 + t58 * t61;
t25 = t68 * t106;
t75 = t25 * mrSges(7,1) - t24 * mrSges(7,2);
t105 = t71 + t72 + t75 + (-t58 * mrSges(3,1) - t64 * mrSges(3,2)) * qJD(2) * pkin(1);
t104 = pkin(6) * m(10);
t103 = pkin(4) * m(11);
t9 = -t59 * t27 + t53 * t35;
t5 = t9 * qJD(9) - t59 * t21 - t53 * t22;
t102 = t5 * mrSges(10,2);
t97 = t53 * t57;
t17 = -t47 * t79 + (-t57 * t80 + (t59 * t63 - t97) * qJD(3)) * pkin(2);
t101 = t17 * mrSges(10,2);
t7 = -t52 * t26 + t51 * t34;
t3 = t7 * qJD(10) - t52 * t19 - t51 * t20;
t100 = t3 * mrSges(11,2);
t99 = t51 * t56;
t77 = qJD(10) * t52;
t36 = (mrSges(11,1) * t78 + mrSges(11,2) * t77) * pkin(4);
t37 = (mrSges(10,1) * t80 + mrSges(10,2) * t79) * pkin(6);
t15 = -t46 * t77 + (t51 * t76 + (-t52 * t62 + t99) * qJD(4)) * pkin(3);
t85 = t15 * mrSges(11,2);
t74 = t12 - t101;
t73 = t11 - t85;
t66 = (-mrSges(5,1) * t56 - mrSges(5,2) * t62) * qJD(4) * pkin(3);
t33 = pkin(2) * t93 - t53 * t47;
t30 = -pkin(2) * t97 - t59 * t47;
t29 = -pkin(3) * t98 - t51 * t46;
t28 = pkin(3) * t99 - t52 * t46;
t8 = [0.2e1 * m(7) * (-t68 * t24 + t67 * t25) * pkin(1) + 0.2e1 * m(11) * (-t3 * t70 + t7 * t4) + 0.2e1 * m(10) * (-t5 * t69 + t9 * t6) + 0.2e1 * m(5) * (t34 * t19 + t31 * t20) + 0.2e1 * m(4) * (t35 * t21 + t32 * t22) - 0.2e1 * t102 - 0.2e1 * t100 + 0.2e1 * t108 + 0.2e1 * t105; (m(5) * (t19 * t56 + t20 * t62 - t31 * t82 + t34 * t81) - mrSges(5,1) * t82 - mrSges(5,2) * t81) * pkin(3) + m(4) * (-t21 * t57 - t22 * t63 + (t32 * t57 - t35 * t63) * qJD(3)) * pkin(2) + m(11) * (-t15 * t70 + t16 * t7 + t28 * t4 + t29 * t3) + m(10) * (-t17 * t69 + t18 * t9 + t30 * t6 + t33 * t5) + (-t17 - t5) * mrSges(10,2) + (-t15 - t3) * mrSges(11,2) + t105 + t107; -0.2e1 * t85 - 0.2e1 * t101 + 0.2e1 * t66 + 0.2e1 * m(10) * (t33 * t17 + t30 * t18) + 0.2e1 * m(11) * (t29 * t15 + t28 * t16) + 0.2e1 * t107; -t102 + (-t5 * t53 - t59 * t6 + (t53 * t9 + t59 * t69) * qJD(9)) * t104 + t72 + t37; (-t17 * t53 - t18 * t59 + (t30 * t53 - t33 * t59) * qJD(9)) * t104 + t74 + t37 + t109; 0.2e1 * t37; -t100 + (-t3 * t51 - t4 * t52 + (t51 * t7 + t52 * t70) * qJD(10)) * t103 + t71 + t36; t66 + (-t15 * t51 - t16 * t52 + (t28 * t51 - t29 * t52) * qJD(10)) * t103 + t73 + t36; 0; 0.2e1 * t36; 0; 0; 0; 0; 0; t75; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; t108; 0; 0; 0; 0; 0; 0; 0; t2 - t102; t74; t37; 0; 0; 0; 0; 0; 0; t1 - t100; t73; 0; t36; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_12_matlab.m
res = [t8(1), t8(2), t8(4), t8(7), t8(11), t8(16), t8(22), t8(29), t8(37), t8(46), t8(56), t8(67); t8(2), t8(3), t8(5), t8(8), t8(12), t8(17), t8(23), t8(30), t8(38), t8(47), t8(57), t8(68); t8(4), t8(5), t8(6), t8(9), t8(13), t8(18), t8(24), t8(31), t8(39), t8(48), t8(58), t8(69); t8(7), t8(8), t8(9), t8(10), t8(14), t8(19), t8(25), t8(32), t8(40), t8(49), t8(59), t8(70); t8(11), t8(12), t8(13), t8(14), t8(15), t8(20), t8(26), t8(33), t8(41), t8(50), t8(60), t8(71); t8(16), t8(17), t8(18), t8(19), t8(20), t8(21), t8(27), t8(34), t8(42), t8(51), t8(61), t8(72); t8(22), t8(23), t8(24), t8(25), t8(26), t8(27), t8(28), t8(35), t8(43), t8(52), t8(62), t8(73); t8(29), t8(30), t8(31), t8(32), t8(33), t8(34), t8(35), t8(36), t8(44), t8(53), t8(63), t8(74); t8(37), t8(38), t8(39), t8(40), t8(41), t8(42), t8(43), t8(44), t8(45), t8(54), t8(64), t8(75); t8(46), t8(47), t8(48), t8(49), t8(50), t8(51), t8(52), t8(53), t8(54), t8(55), t8(65), t8(76); t8(56), t8(57), t8(58), t8(59), t8(60), t8(61), t8(62), t8(63), t8(64), t8(65), t8(66), t8(77); t8(67), t8(68), t8(69), t8(70), t8(71), t8(72), t8(73), t8(74), t8(75), t8(76), t8(77), t8(78);];
Mq = res;
