% Calculate kinetic energy for
% palh3m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% qJD [10x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [9x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh3m2OL_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(6,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_energykin_floatb_twist_slag_vp2: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m2OL_energykin_floatb_twist_slag_vp2: qJD has to be [10x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh3m2OL_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_energykin_floatb_twist_slag_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2OL_energykin_floatb_twist_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2OL_energykin_floatb_twist_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2OL_energykin_floatb_twist_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:33:22
% EndTime: 2020-05-07 04:33:27
% DurationCPUTime: 2.64s
% Computational Cost: add. (2359->196), mult. (3217->295), div. (0->0), fcn. (2596->18), ass. (0->81)
t81 = sin(qJ(1));
t89 = cos(qJ(1));
t63 = -t81 * V_base(4) + t89 * V_base(5);
t60 = qJD(2) - t63;
t57 = qJD(7) + t60;
t92 = pkin(3) * t57;
t91 = sin(qJ(8));
t68 = V_base(5) * pkin(11) + V_base(1);
t69 = -V_base(4) * pkin(11) + V_base(2);
t51 = t89 * t68 + t81 * t69;
t53 = -pkin(12) * t63 + V_base(3);
t80 = sin(qJ(2));
t88 = cos(qJ(2));
t41 = -t51 * t80 + t88 * t53;
t30 = pkin(1) * t60 + t41;
t42 = t51 * t88 + t53 * t80;
t79 = sin(qJ(3));
t87 = cos(qJ(3));
t24 = -t30 * t87 + t79 * t42;
t58 = qJD(3) + t60;
t18 = pkin(4) * t58 + t24;
t26 = -t30 * t79 - t42 * t87;
t78 = sin(qJ(4));
t86 = cos(qJ(4));
t9 = t78 * t18 + t86 * t26;
t75 = sin(qJ(7));
t83 = cos(qJ(7));
t25 = t75 * t30 + t83 * t42;
t23 = t83 * t30 - t42 * t75;
t50 = -t81 * t68 + t69 * t89;
t8 = t18 * t86 - t26 * t78;
t64 = t81 * V_base(5) + t89 * V_base(4);
t72 = V_base(6) + qJD(1);
t47 = -t64 * t80 + t72 * t88;
t49 = t64 * t88 + t72 * t80;
t34 = -t47 * t87 + t49 * t79;
t36 = -t47 * t79 - t49 * t87;
t20 = t34 * t86 - t36 * t78;
t45 = -pkin(12) * t72 - t50;
t38 = -pkin(1) * t47 + t45;
t27 = -pkin(4) * t34 + t38;
t90 = V_base(3) ^ 2;
t85 = cos(qJ(5));
t84 = cos(qJ(6));
t82 = cos(qJ(8));
t77 = sin(qJ(5));
t76 = sin(qJ(6));
t74 = cos(pkin(15));
t73 = sin(pkin(15));
t62 = -t73 * t91 - t74 * t82;
t61 = t73 * t82 - t74 * t91;
t59 = qJD(6) - t63;
t56 = qJD(4) + t58;
t55 = qJD(8) + t57;
t54 = pkin(6) * t63 + V_base(3);
t48 = t64 * t84 + t72 * t76;
t46 = -t64 * t76 + t72 * t84;
t44 = pkin(13) * t63 + t51;
t43 = pkin(6) * t72 + pkin(13) * t64 - t50;
t37 = t38 ^ 2;
t35 = t47 * t75 + t49 * t83;
t33 = t47 * t83 - t49 * t75;
t32 = t44 * t84 + t54 * t76;
t31 = -t44 * t76 + t54 * t84;
t21 = t34 * t78 + t36 * t86;
t19 = qJD(5) - t20;
t17 = t73 * t92 + t25;
t16 = t74 * t92 + t23;
t14 = t33 * t61 + t35 * t62;
t13 = t33 * t62 - t35 * t61;
t12 = t21 * t85 + t56 * t77;
t11 = -t21 * t77 + t56 * t85;
t10 = (-t33 * t74 - t35 * t73) * pkin(3) + t38;
t7 = pkin(10) * t56 + t9;
t6 = -pkin(8) * t56 - t8;
t5 = t16 * t61 + t17 * t62;
t4 = t16 * t62 - t17 * t61;
t3 = -pkin(8) * t20 - pkin(10) * t21 + t27;
t2 = t3 * t77 + t7 * t85;
t1 = t3 * t85 - t7 * t77;
t15 = (t10 * mrSges(9,2) - t4 * mrSges(9,3) + Ifges(9,5) * t55 + Ifges(9,1) * t14 / 0.2e1) * t14 + (V_base(3) * mrSges(2,2) - t50 * mrSges(2,3) + Ifges(2,5) * t72 + Ifges(2,1) * t64 / 0.2e1) * t64 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t6 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t19 + Ifges(6,1) * t12 / 0.2e1) * t12 + (-V_base(3) * mrSges(2,1) + t51 * mrSges(2,3) + Ifges(2,4) * t64 + Ifges(2,6) * t72 + Ifges(2,2) * t63 / 0.2e1) * t63 + (-t10 * mrSges(9,1) + t5 * mrSges(9,3) + Ifges(9,4) * t14 + Ifges(9,6) * t55 + Ifges(9,2) * t13 / 0.2e1) * t13 + (t41 * mrSges(3,1) - t42 * mrSges(3,2) + Ifges(3,3) * t60 / 0.2e1) * t60 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t50 * mrSges(2,1) - t51 * mrSges(2,2) + Ifges(2,3) * t72 / 0.2e1) * t72 + (t45 * mrSges(3,2) - t41 * mrSges(3,3) + Ifges(3,5) * t60 + Ifges(3,1) * t49 / 0.2e1) * t49 + (t38 * mrSges(8,2) - t23 * mrSges(8,3) + Ifges(8,5) * t57 + Ifges(8,1) * t35 / 0.2e1) * t35 + (t43 * mrSges(7,2) - t31 * mrSges(7,3) + Ifges(7,5) * t59 + Ifges(7,1) * t48 / 0.2e1) * t48 + (-t45 * mrSges(3,1) + t42 * mrSges(3,3) + Ifges(3,4) * t49 + Ifges(3,6) * t60 + Ifges(3,2) * t47 / 0.2e1) * t47 + (-t43 * mrSges(7,1) + t32 * mrSges(7,3) + Ifges(7,4) * t48 + Ifges(7,6) * t59 + Ifges(7,2) * t46 / 0.2e1) * t46 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t23 * mrSges(8,1) - t25 * mrSges(8,2) + Ifges(8,3) * t57 / 0.2e1) * t57 + m(2) * (t50 ^ 2 + t51 ^ 2 + t90) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t90) / 0.2e1 + m(3) * (t41 ^ 2 + t42 ^ 2 + t45 ^ 2) / 0.2e1 + (t31 * mrSges(7,1) - t32 * mrSges(7,2) + Ifges(7,3) * t59 / 0.2e1) * t59 + m(4) * (t24 ^ 2 + t26 ^ 2 + t37) / 0.2e1 + m(8) * (t23 ^ 2 + t25 ^ 2 + t37) / 0.2e1 + m(7) * (t31 ^ 2 + t32 ^ 2 + t43 ^ 2) / 0.2e1 + m(5) * (t27 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + m(9) * (t10 ^ 2 + t4 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t6 ^ 2) / 0.2e1 + (t24 * mrSges(4,1) - t26 * mrSges(4,2) + Ifges(4,3) * t58 / 0.2e1) * t58 + (t38 * mrSges(4,2) - t24 * mrSges(4,3) + Ifges(4,5) * t58 + Ifges(4,1) * t36 / 0.2e1) * t36 + (-t6 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t12 + Ifges(6,6) * t19 + Ifges(6,2) * t11 / 0.2e1) * t11 + (-t38 * mrSges(4,1) + t26 * mrSges(4,3) + Ifges(4,4) * t36 + Ifges(4,6) * t58 + Ifges(4,2) * t34 / 0.2e1) * t34 + (t8 * mrSges(5,1) - t9 * mrSges(5,2) + Ifges(5,3) * t56 / 0.2e1) * t56 + (-t38 * mrSges(8,1) + t25 * mrSges(8,3) + Ifges(8,4) * t35 + Ifges(8,6) * t57 + Ifges(8,2) * t33 / 0.2e1) * t33 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t19 / 0.2e1) * t19 + (t27 * mrSges(5,2) - t8 * mrSges(5,3) + Ifges(5,5) * t56 + Ifges(5,1) * t21 / 0.2e1) * t21 + (-t27 * mrSges(5,1) + t9 * mrSges(5,3) + Ifges(5,4) * t21 + Ifges(5,6) * t56 + Ifges(5,2) * t20 / 0.2e1) * t20 + (t4 * mrSges(9,1) - t5 * mrSges(9,2) + Ifges(9,3) * t55 / 0.2e1) * t55;
T = t15;
