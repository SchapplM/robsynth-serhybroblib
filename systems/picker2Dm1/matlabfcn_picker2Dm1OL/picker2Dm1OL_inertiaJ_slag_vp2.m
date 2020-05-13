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
% Mq [12x12]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:46
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = picker2Dm1OL_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1OL_inertiaJ_slag_vp2: qJ has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1OL_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1OL_inertiaJ_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm1OL_inertiaJ_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'picker2Dm1OL_inertiaJ_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 05:44:57
% EndTime: 2020-05-11 05:44:57
% DurationCPUTime: 0.33s
% Computational Cost: add. (423->138), mult. (814->172), div. (0->0), fcn. (564->14), ass. (0->78)
t66 = Ifges(11,3) + Ifges(5,3);
t84 = Ifges(3,3) + t66;
t35 = sin(qJ(10));
t30 = t35 * pkin(4) * mrSges(11,2);
t83 = Ifges(5,3) + t30;
t48 = cos(qJ(2));
t34 = t48 * pkin(1);
t29 = t34 + pkin(2);
t41 = sin(qJ(3));
t47 = cos(qJ(3));
t42 = sin(qJ(2));
t81 = pkin(1) * t42;
t17 = -t29 * t47 + t41 * t81;
t10 = t17 * mrSges(4,1);
t20 = -t29 * t41 - t47 * t81;
t73 = t20 * mrSges(4,2);
t82 = t10 - t73;
t80 = pkin(2) * t41;
t40 = sin(qJ(4));
t79 = t40 * pkin(3);
t46 = cos(qJ(4));
t78 = t46 * pkin(3);
t77 = t47 * pkin(2);
t12 = pkin(6) + t17;
t37 = sin(qJ(9));
t43 = cos(qJ(9));
t6 = -t12 * t37 - t20 * t43;
t76 = t6 * mrSges(10,2);
t27 = pkin(6) - t77;
t18 = -t27 * t37 + t43 * t80;
t75 = t18 * mrSges(10,2);
t28 = t34 + pkin(3);
t19 = t28 * t40 + t46 * t81;
t74 = t19 * mrSges(5,2);
t39 = sin(qJ(6));
t45 = cos(qJ(6));
t23 = (-t39 * t48 - t42 * t45) * pkin(1);
t72 = t23 * mrSges(7,2);
t16 = t46 * t28 - t40 * t81;
t11 = pkin(4) + t16;
t36 = cos(qJ(10));
t4 = -t11 * t35 - t19 * t36;
t71 = t4 * mrSges(11,2);
t70 = t43 * mrSges(10,1);
t69 = Ifges(4,3) + Ifges(10,3);
t26 = pkin(4) + t78;
t14 = -t26 * t35 - t36 * t79;
t68 = t14 * mrSges(11,2);
t67 = t36 * mrSges(11,1);
t65 = mrSges(5,2) * t79;
t64 = pkin(6) * t70;
t63 = mrSges(4,1) * t77;
t62 = pkin(4) * t67;
t5 = -t12 * t43 + t20 * t37;
t2 = t5 * mrSges(10,1);
t61 = Ifges(10,3) + t2 - t76;
t3 = -t11 * t36 + t19 * t35;
t1 = t3 * mrSges(11,1);
t60 = Ifges(11,3) + t1 - t71;
t22 = (t39 * t42 - t45 * t48) * pkin(1);
t21 = t22 * mrSges(7,1);
t59 = Ifges(7,3) + t21 - t72;
t13 = -t26 * t36 + t35 * t79;
t7 = t13 * mrSges(11,1);
t58 = Ifges(11,3) + t7 - t68;
t57 = mrSges(3,1) * t48 - mrSges(3,2) * t42;
t38 = sin(qJ(8));
t44 = cos(qJ(8));
t56 = -mrSges(9,1) * t44 + mrSges(9,2) * t38;
t55 = Ifges(7,3) + t69 + t84;
t32 = mrSges(4,2) * t80;
t15 = -t27 * t43 - t37 * t80;
t8 = t15 * mrSges(10,1);
t54 = t32 + t8 - t63 + t69;
t33 = mrSges(5,1) * t78;
t31 = t37 * pkin(6) * mrSges(10,2);
t9 = t16 * mrSges(5,1);
t24 = [t55 + m(5) * (t16 ^ 2 + t19 ^ 2) + m(4) * (t17 ^ 2 + t20 ^ 2) + m(7) * (t22 ^ 2 + t23 ^ 2) + m(11) * (t3 ^ 2 + t4 ^ 2) + m(10) * (t5 ^ 2 + t6 ^ 2) + 0.2e1 * t21 + 0.2e1 * t9 + 0.2e1 * t10 + 0.2e1 * t2 + 0.2e1 * t1 - 0.2e1 * t73 - 0.2e1 * t72 - 0.2e1 * t74 - 0.2e1 * t71 - 0.2e1 * t76 + Ifges(2,3) + Ifges(9,3) + (0.2e1 * t56 + 0.2e1 * t57 + (m(9) * (t38 ^ 2 + t44 ^ 2) + m(3) * (t42 ^ 2 + t48 ^ 2)) * pkin(1)) * pkin(1); t59 + t57 * pkin(1) + t54 + m(11) * (t13 * t3 + t14 * t4) + m(10) * (t15 * t5 + t18 * t6) + m(4) * (-t17 * t47 - t20 * t41) * pkin(2) + m(5) * (t16 * t46 + t19 * t40) * pkin(3) + (-t18 - t6) * mrSges(10,2) + (-t19 - t79) * mrSges(5,2) + (-t14 - t4) * mrSges(11,2) + t33 + t9 + t7 + t2 + t1 + t82 + t84; t55 + m(5) * (t40 ^ 2 + t46 ^ 2) * pkin(3) ^ 2 + m(4) * (t41 ^ 2 + t47 ^ 2) * pkin(2) ^ 2 + m(10) * (t15 ^ 2 + t18 ^ 2) + m(11) * (t13 ^ 2 + t14 ^ 2) + 0.2e1 * t33 + 0.2e1 * t32 + 0.2e1 * t7 + 0.2e1 * t8 - 0.2e1 * t68 - 0.2e1 * t75 - 0.2e1 * t63 - 0.2e1 * t65; Ifges(4,3) + t31 + (m(10) * (-t37 * t6 - t43 * t5) - t70) * pkin(6) + t61 + t82; -t75 + t31 + (m(10) * (-t15 * t43 - t18 * t37) - t70) * pkin(6) + t54; -0.2e1 * t64 + 0.2e1 * t31 + m(10) * (t37 ^ 2 + t43 ^ 2) * pkin(6) ^ 2 + t69; -t74 + t9 + (-t67 + m(11) * (-t3 * t36 - t35 * t4)) * pkin(4) + t60 + t83; -t65 + t33 + (m(11) * (-t13 * t36 - t14 * t35) - t67) * pkin(4) + t58 + t83; 0; -0.2e1 * t62 + 0.2e1 * t30 + m(11) * (t35 ^ 2 + t36 ^ 2) * pkin(4) ^ 2 + t66; 0; 0; 0; 0; Ifges(6,3); t59; Ifges(7,3); 0; 0; 0; Ifges(7,3); 0; 0; 0; 0; 0; 0; Ifges(8,3); pkin(1) * t56 + Ifges(9,3); 0; 0; 0; 0; 0; 0; Ifges(9,3); t61; Ifges(10,3) + t8 - t75; Ifges(10,3) + t31 - t64; 0; 0; 0; 0; 0; Ifges(10,3); t60; t58; 0; Ifges(11,3) + t30 - t62; 0; 0; 0; 0; 0; Ifges(11,3); 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_12_matlab.m
res = [t24(1), t24(2), t24(4), t24(7), t24(11), t24(16), t24(22), t24(29), t24(37), t24(46), t24(56), t24(67); t24(2), t24(3), t24(5), t24(8), t24(12), t24(17), t24(23), t24(30), t24(38), t24(47), t24(57), t24(68); t24(4), t24(5), t24(6), t24(9), t24(13), t24(18), t24(24), t24(31), t24(39), t24(48), t24(58), t24(69); t24(7), t24(8), t24(9), t24(10), t24(14), t24(19), t24(25), t24(32), t24(40), t24(49), t24(59), t24(70); t24(11), t24(12), t24(13), t24(14), t24(15), t24(20), t24(26), t24(33), t24(41), t24(50), t24(60), t24(71); t24(16), t24(17), t24(18), t24(19), t24(20), t24(21), t24(27), t24(34), t24(42), t24(51), t24(61), t24(72); t24(22), t24(23), t24(24), t24(25), t24(26), t24(27), t24(28), t24(35), t24(43), t24(52), t24(62), t24(73); t24(29), t24(30), t24(31), t24(32), t24(33), t24(34), t24(35), t24(36), t24(44), t24(53), t24(63), t24(74); t24(37), t24(38), t24(39), t24(40), t24(41), t24(42), t24(43), t24(44), t24(45), t24(54), t24(64), t24(75); t24(46), t24(47), t24(48), t24(49), t24(50), t24(51), t24(52), t24(53), t24(54), t24(55), t24(65), t24(76); t24(56), t24(57), t24(58), t24(59), t24(60), t24(61), t24(62), t24(63), t24(64), t24(65), t24(66), t24(77); t24(67), t24(68), t24(69), t24(70), t24(71), t24(72), t24(73), t24(74), t24(75), t24(76), t24(77), t24(78);];
Mq = res;
