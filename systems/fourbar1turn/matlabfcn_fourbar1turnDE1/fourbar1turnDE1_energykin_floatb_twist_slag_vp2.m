% Calculate kinetic energy for
% fourbar1turnDE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar1turnDE1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_energykin_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE1_energykin_floatb_twist_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fourbar1turnDE1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_energykin_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnDE1_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnDE1_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnDE1_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:25:24
% EndTime: 2020-04-12 19:25:26
% DurationCPUTime: 1.27s
% Computational Cost: add. (7915->143), mult. (11217->253), div. (552->13), fcn. (3264->10), ass. (0->75)
t93 = pkin(4) ^ 2;
t92 = pkin(3) ^ 2;
t91 = -2 * pkin(2);
t48 = V_base(5) * pkin(5) + V_base(1);
t49 = -V_base(4) * pkin(5) + V_base(2);
t58 = sin(qJ(1));
t60 = cos(qJ(1));
t33 = t48 * t58 - t60 * t49;
t90 = t33 ^ 2;
t89 = (-pkin(3) - pkin(4));
t88 = (-pkin(3) + pkin(4));
t57 = sin(qJ(2));
t87 = pkin(1) * t57;
t59 = cos(qJ(2));
t86 = pkin(2) * t59;
t68 = pkin(2) ^ 2;
t69 = pkin(1) ^ 2;
t81 = -0.2e1 * pkin(1) * t86 + t69;
t47 = t68 + t81;
t80 = t92 - t93;
t43 = t47 - t80;
t85 = t43 * t57;
t45 = 0.1e1 / t47;
t84 = t45 * t57 ^ 2;
t37 = ((pkin(2) - t89) * (pkin(2) + t89)) + t81;
t38 = ((pkin(2) - t88) * (pkin(2) + t88)) + t81;
t70 = sqrt(-t37 * t38);
t83 = t57 * t70;
t82 = t59 * t70;
t46 = 0.1e1 / t47 ^ 2;
t79 = t46 * t91;
t51 = pkin(1) - t86;
t21 = pkin(2) * t85 + t51 * t70;
t17 = t21 ^ 2;
t63 = 0.1e1 / pkin(4);
t20 = -pkin(2) * t83 + t43 * t51;
t72 = t20 ^ 2;
t78 = ((t17 + t72) / t93 * t46) ^ (-0.1e1 / 0.2e1) * t45 * t63;
t42 = t47 + t80;
t52 = pkin(1) * t59 - pkin(2);
t22 = t42 * t87 - t52 * t70;
t18 = t22 ^ 2;
t66 = 0.1e1 / pkin(3);
t19 = -pkin(1) * t83 - t42 * t52;
t71 = t19 ^ 2;
t77 = ((t18 + t71) / t92 * t46) ^ (-0.1e1 / 0.2e1) * t45 * t66;
t76 = (-t37 - t38) * qJD(2) * pkin(2) * t87 / t70 * t45;
t35 = t48 * t60 + t49 * t58;
t28 = -t35 * t57 + t59 * V_base(3);
t74 = t76 / 0.2e1;
t40 = -t58 * V_base(4) + t60 * V_base(5);
t39 = qJD(2) - t40;
t61 = V_base(3) ^ 2;
t55 = V_base(6) + qJD(1);
t41 = t58 * V_base(5) + t60 * V_base(4);
t36 = -pkin(1) * t40 + V_base(3);
t32 = t41 * t59 + t55 * t57;
t31 = -t41 * t57 + t55 * t59;
t30 = -pkin(1) * t55 + t33;
t29 = t35 * t59 + t57 * V_base(3);
t24 = pkin(2) * t39 + t28;
t23 = -pkin(2) * t31 + t33;
t16 = 0.1e1 / t72;
t15 = 0.1e1 / t71;
t10 = (-t20 * t55 - t21 * t41) * t78;
t9 = (-t20 * t41 + t21 * t55) * t78;
t8 = (-t20 * t36 - t21 * t35) * t78;
t7 = (-t20 * t35 + t21 * t36) * t78;
t6 = (-t19 * t31 + t22 * t32) * t77;
t5 = (-t19 * t32 - t22 * t31) * t77;
t4 = (-t19 * t24 + t22 * t29) * t77;
t3 = (-t19 * t29 - t22 * t24) * t77;
t2 = -t40 - 0.2e1 * ((t51 * t74 + (t68 * pkin(1) * t84 + ((t43 * t59 + t83) * t45 / 0.2e1 - t21 * t46 * t87) * pkin(2)) * qJD(2)) / t20 + (t57 * t74 + (-(-t82 + t85) * t45 / 0.2e1 + (t20 * t46 - t45 * t51) * t87) * qJD(2)) * pkin(2) * t21 * t16) * pkin(4) / (t16 * t17 + 0.1e1) * t47 * t63;
t1 = t39 + ((-t52 * t76 + (0.2e1 * t69 * pkin(2) * t84 + ((t42 * t59 + t83) * t45 + t22 * t57 * t79) * pkin(1)) * qJD(2)) / t19 - (-t57 * t76 + (-t45 * t82 + ((t52 * t91 + t42) * t45 + t19 * t79) * t57) * qJD(2)) * pkin(1) * t22 * t15) * pkin(3) / (t15 * t18 + 0.1e1) * t47 * t66;
t11 = m(4) * (t23 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + m(5) * (t30 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + m(3) * (t28 ^ 2 + t29 ^ 2 + t90) / 0.2e1 + m(2) * (t35 ^ 2 + t61 + t90) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t61) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t30 * mrSges(5,2) - t8 * mrSges(5,3) + Ifges(5,1) * t9 / 0.2e1) * t9 + (-t23 * mrSges(4,1) + t3 * mrSges(4,3) + Ifges(4,2) * t6 / 0.2e1) * t6 + (-t33 * mrSges(2,1) - t35 * mrSges(2,2) + Ifges(2,3) * t55 / 0.2e1) * t55 + (t28 * mrSges(3,1) - t29 * mrSges(3,2) + Ifges(3,3) * t39 / 0.2e1) * t39 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t23 * mrSges(4,2) - t4 * mrSges(4,3) + Ifges(4,4) * t6 + Ifges(4,1) * t5 / 0.2e1) * t5 + (V_base(3) * mrSges(2,2) + t33 * mrSges(2,3) + Ifges(2,5) * t55 + Ifges(2,1) * t41 / 0.2e1) * t41 + (t33 * mrSges(3,2) - t28 * mrSges(3,3) + Ifges(3,5) * t39 + Ifges(3,1) * t32 / 0.2e1) * t32 + (t8 * mrSges(5,1) - t7 * mrSges(5,2) + Ifges(5,5) * t9 + Ifges(5,3) * t2 / 0.2e1) * t2 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t35 * mrSges(2,3) + Ifges(2,4) * t41 + Ifges(2,6) * t55 + Ifges(2,2) * t40 / 0.2e1) * t40 + (-t33 * mrSges(3,1) + t29 * mrSges(3,3) + Ifges(3,4) * t32 + Ifges(3,6) * t39 + Ifges(3,2) * t31 / 0.2e1) * t31 + (-t30 * mrSges(5,1) + t7 * mrSges(5,3) + Ifges(5,4) * t9 + Ifges(5,6) * t2 + Ifges(5,2) * t10 / 0.2e1) * t10 + (t4 * mrSges(4,1) - t3 * mrSges(4,2) + Ifges(4,5) * t5 + Ifges(4,6) * t6 + Ifges(4,3) * t1 / 0.2e1) * t1;
T = t11;
