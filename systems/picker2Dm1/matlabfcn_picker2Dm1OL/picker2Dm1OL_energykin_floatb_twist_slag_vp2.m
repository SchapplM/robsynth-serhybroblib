% Calculate kinetic energy for
% picker2Dm1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% qJD [12x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:46
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = picker2Dm1OL_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(6,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1OL_energykin_floatb_twist_slag_vp2: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm1OL_energykin_floatb_twist_slag_vp2: qJD has to be [12x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'picker2Dm1OL_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1OL_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1OL_energykin_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm1OL_energykin_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'picker2Dm1OL_energykin_floatb_twist_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 05:44:48
% EndTime: 2020-05-11 05:44:55
% DurationCPUTime: 4.40s
% Computational Cost: add. (1675->232), mult. (2587->352), div. (0->0), fcn. (2160->22), ass. (0->93)
t96 = pkin(5) * V_base(6);
t85 = sin(qJ(1));
t94 = cos(qJ(1));
t56 = t85 * V_base(1) - t94 * V_base(2);
t72 = V_base(6) + qJD(1);
t44 = pkin(1) * t72 + t56;
t57 = -t85 * V_base(2) - t94 * V_base(1);
t84 = sin(qJ(2));
t93 = cos(qJ(2));
t34 = t93 * t44 - t57 * t84;
t68 = qJD(2) + t72;
t30 = pkin(3) * t68 + t34;
t36 = t44 * t84 + t57 * t93;
t82 = sin(qJ(4));
t91 = cos(qJ(4));
t13 = t91 * t30 - t36 * t82;
t31 = pkin(2) * t68 + t34;
t83 = sin(qJ(3));
t92 = cos(qJ(3));
t14 = -t31 * t92 + t83 * t36;
t53 = t85 * V_base(4) - t94 * V_base(5);
t46 = -pkin(1) * t53 + V_base(3);
t64 = qJD(3) + t68;
t63 = qJD(4) + t68;
t55 = -t85 * V_base(5) - t94 * V_base(4);
t40 = t53 * t93 - t55 * t84;
t26 = -pkin(2) * t40 + t46;
t25 = -pkin(3) * t40 + t46;
t95 = V_base(3) ^ 2;
t90 = cos(qJ(5));
t89 = cos(qJ(6));
t88 = cos(qJ(7));
t87 = cos(qJ(8));
t86 = cos(qJ(9));
t81 = sin(qJ(5));
t80 = sin(qJ(6));
t79 = sin(qJ(7));
t78 = sin(qJ(8));
t77 = sin(qJ(9));
t76 = cos(pkin(8));
t75 = cos(qJ(10));
t74 = sin(pkin(8));
t73 = sin(qJ(10));
t71 = V_base(6) + qJD(5);
t70 = V_base(6) + qJD(7);
t67 = qJD(8) + t72;
t66 = -V_base(5) * pkin(7) + V_base(3);
t65 = V_base(6) * pkin(7) + V_base(2);
t62 = qJD(6) + t68;
t61 = qJD(9) + t64;
t60 = qJD(10) + t63;
t59 = -t74 * t96 + V_base(1);
t58 = t76 * t96 + V_base(2);
t54 = t79 * V_base(5) + t88 * V_base(4);
t52 = t79 * V_base(4) - t88 * V_base(5);
t51 = -t74 * t81 + t76 * t90;
t50 = t74 * t90 + t76 * t81;
t49 = t65 * t79 + t88 * V_base(1);
t48 = -t65 * t88 + t79 * V_base(1);
t47 = V_base(3) + (t74 * V_base(4) - t76 * V_base(5)) * pkin(5);
t45 = t46 ^ 2;
t42 = t53 * t84 + t55 * t93;
t41 = -t53 * t78 - t55 * t87;
t39 = -t53 * t87 + t55 * t78;
t38 = t50 * V_base(5) + t51 * V_base(4);
t37 = -t50 * V_base(4) + t51 * V_base(5);
t35 = -t44 * t78 - t57 * t87;
t33 = -t44 * t87 + t57 * t78;
t28 = t50 * t58 + t51 * t59;
t27 = -t50 * t59 + t51 * t58;
t24 = -t40 * t83 - t42 * t92;
t23 = t40 * t82 + t42 * t91;
t22 = -t40 * t80 - t42 * t89;
t21 = -t40 * t92 + t42 * t83;
t20 = t40 * t91 - t42 * t82;
t19 = -t40 * t89 + t42 * t80;
t18 = -t34 * t80 - t36 * t89;
t17 = -t34 * t89 + t36 * t80;
t16 = -t31 * t83 - t36 * t92;
t15 = t30 * t82 + t36 * t91;
t12 = pkin(6) * t64 + t14;
t11 = pkin(4) * t63 + t13;
t10 = -pkin(6) * t21 + t26;
t9 = -pkin(4) * t20 + t25;
t8 = -t21 * t77 - t24 * t86;
t7 = -t21 * t86 + t24 * t77;
t6 = -t20 * t73 - t23 * t75;
t5 = -t20 * t75 + t23 * t73;
t4 = -t12 * t77 - t16 * t86;
t3 = -t12 * t86 + t16 * t77;
t2 = -t11 * t73 - t15 * t75;
t1 = -t11 * t75 + t15 * t73;
t29 = (V_base(3) * mrSges(2,2) - t56 * mrSges(2,3) + Ifges(2,5) * t72 + Ifges(2,1) * t55 / 0.2e1) * t55 + (-t66 * mrSges(8,1) + t48 * mrSges(8,3) + Ifges(8,6) * t70 + Ifges(8,2) * t54 / 0.2e1) * t54 + (t46 * mrSges(3,2) - t34 * mrSges(3,3) + Ifges(3,5) * t68 + Ifges(3,1) * t42 / 0.2e1) * t42 + (t46 * mrSges(7,2) - t17 * mrSges(7,3) + Ifges(7,5) * t62 + Ifges(7,1) * t22 / 0.2e1) * t22 + (-t46 * mrSges(3,1) + t36 * mrSges(3,3) + Ifges(3,4) * t42 + Ifges(3,6) * t68 + Ifges(3,2) * t40 / 0.2e1) * t40 + (-t46 * mrSges(9,1) + t35 * mrSges(9,3) + Ifges(9,4) * t41 + Ifges(9,6) * t67 + Ifges(9,2) * t39 / 0.2e1) * t39 + (t13 * mrSges(5,1) - t15 * mrSges(5,2) + Ifges(5,3) * t63 / 0.2e1) * t63 + (-t47 * mrSges(6,1) + t28 * mrSges(6,3) + Ifges(6,4) * t38 + Ifges(6,6) * t71 + Ifges(6,2) * t37 / 0.2e1) * t37 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t47 * mrSges(6,2) - t27 * mrSges(6,3) + Ifges(6,5) * t71 + Ifges(6,1) * t38 / 0.2e1) * t38 + (t3 * mrSges(10,1) - t4 * mrSges(10,2) + Ifges(10,5) * t8 + Ifges(10,6) * t7 + Ifges(10,3) * t61 / 0.2e1) * t61 + (-t26 * mrSges(4,1) + t16 * mrSges(4,3) + Ifges(4,4) * t24 + Ifges(4,6) * t64 + Ifges(4,2) * t21 / 0.2e1) * t21 + (-t25 * mrSges(5,1) + t15 * mrSges(5,3) + Ifges(5,4) * t23 + Ifges(5,6) * t63 + Ifges(5,2) * t20 / 0.2e1) * t20 + (t34 * mrSges(3,1) - t36 * mrSges(3,2) + Ifges(3,3) * t68 / 0.2e1) * t68 + (t26 * mrSges(4,2) - t14 * mrSges(4,3) + Ifges(4,5) * t64 + Ifges(4,1) * t24 / 0.2e1) * t24 + (-t10 * mrSges(10,1) + t4 * mrSges(10,3) + Ifges(10,4) * t8 + Ifges(10,2) * t7 / 0.2e1) * t7 + (t56 * mrSges(2,1) - t57 * mrSges(2,2) + Ifges(2,3) * t72 / 0.2e1) * t72 + (-t46 * mrSges(7,1) + t18 * mrSges(7,3) + Ifges(7,4) * t22 + Ifges(7,6) * t62 + Ifges(7,2) * t19 / 0.2e1) * t19 + (t17 * mrSges(7,1) - t18 * mrSges(7,2) + Ifges(7,3) * t62 / 0.2e1) * t62 + (t25 * mrSges(5,2) - t13 * mrSges(5,3) + Ifges(5,5) * t63 + Ifges(5,1) * t23 / 0.2e1) * t23 + (t10 * mrSges(10,2) - t3 * mrSges(10,3) + Ifges(10,1) * t8 / 0.2e1) * t8 + (t46 * mrSges(9,2) - t33 * mrSges(9,3) + Ifges(9,5) * t67 + Ifges(9,1) * t41 / 0.2e1) * t41 + (t27 * mrSges(6,1) - t28 * mrSges(6,2) + Ifges(6,3) * t71 / 0.2e1) * t71 + (-V_base(3) * mrSges(2,1) + t57 * mrSges(2,3) + Ifges(2,4) * t55 + Ifges(2,6) * t72 + Ifges(2,2) * t53 / 0.2e1) * t53 + m(2) * (t56 ^ 2 + t57 ^ 2 + t95) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t95) / 0.2e1 + m(8) * (t48 ^ 2 + t49 ^ 2 + t66 ^ 2) / 0.2e1 + m(9) * (t33 ^ 2 + t35 ^ 2 + t45) / 0.2e1 + m(3) * (t34 ^ 2 + t36 ^ 2 + t45) / 0.2e1 + m(7) * (t17 ^ 2 + t18 ^ 2 + t45) / 0.2e1 + m(6) * (t27 ^ 2 + t28 ^ 2 + t47 ^ 2) / 0.2e1 + m(4) * (t14 ^ 2 + t16 ^ 2 + t26 ^ 2) / 0.2e1 + m(5) * (t13 ^ 2 + t15 ^ 2 + t25 ^ 2) / 0.2e1 + m(11) * (t1 ^ 2 + t2 ^ 2 + t9 ^ 2) / 0.2e1 + m(10) * (t10 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + (t9 * mrSges(11,2) - t1 * mrSges(11,3) + Ifges(11,5) * t60 + Ifges(11,1) * t6 / 0.2e1) * t6 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t66 * mrSges(8,2) - t49 * mrSges(8,3) + Ifges(8,4) * t54 + Ifges(8,5) * t70 + Ifges(8,1) * t52 / 0.2e1) * t52 + (t33 * mrSges(9,1) - t35 * mrSges(9,2) + Ifges(9,3) * t67 / 0.2e1) * t67 + (t1 * mrSges(11,1) - t2 * mrSges(11,2) + Ifges(11,3) * t60 / 0.2e1) * t60 + (t49 * mrSges(8,1) - t48 * mrSges(8,2) + Ifges(8,3) * t70 / 0.2e1) * t70 + (-t9 * mrSges(11,1) + t2 * mrSges(11,3) + Ifges(11,4) * t6 + Ifges(11,6) * t60 + Ifges(11,2) * t5 / 0.2e1) * t5 + (t14 * mrSges(4,1) - t16 * mrSges(4,2) + Ifges(4,3) * t64 / 0.2e1) * t64;
T = t29;
