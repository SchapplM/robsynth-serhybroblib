% Calculate kinetic energy for
% fourbar1turnTE
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
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar1turnTE_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_energykin_floatb_twist_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnTE_energykin_floatb_twist_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fourbar1turnTE_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_energykin_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnTE_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnTE_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnTE_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:18:28
% EndTime: 2020-04-12 19:18:29
% DurationCPUTime: 1.08s
% Computational Cost: add. (5163->141), mult. (7185->256), div. (296->11), fcn. (2112->6), ass. (0->76)
t91 = -2 * pkin(2);
t44 = V_base(5) * pkin(5) + V_base(1);
t45 = -V_base(4) * pkin(5) + V_base(2);
t54 = sin(qJ(1));
t56 = cos(qJ(1));
t29 = t44 * t54 - t56 * t45;
t90 = t29 ^ 2;
t62 = pkin(2) ^ 2;
t63 = pkin(1) ^ 2;
t55 = cos(qJ(2));
t78 = pkin(2) * t55;
t71 = -0.2e1 * pkin(1) * t78 + t63;
t43 = t62 + t71;
t70 = pkin(3) ^ 2 - pkin(4) ^ 2;
t38 = t43 + t70;
t48 = pkin(1) * t55 - pkin(2);
t53 = sin(qJ(2));
t81 = (-pkin(3) - pkin(4));
t33 = ((pkin(2) - t81) * (pkin(2) + t81)) + t71;
t80 = (pkin(4) - pkin(3));
t34 = ((pkin(2) - t80) * (pkin(2) + t80)) + t71;
t64 = sqrt(-t33 * t34);
t73 = t53 * t64;
t15 = -pkin(1) * t73 - t38 * t48;
t89 = -t15 / 0.2e1;
t39 = t43 - t70;
t47 = pkin(1) - t78;
t16 = -pkin(2) * t73 + t39 * t47;
t88 = -t16 / 0.2e1;
t31 = t44 * t56 + t45 * t54;
t24 = -t31 * t53 + t55 * V_base(3);
t36 = -t54 * V_base(4) + t56 * V_base(5);
t35 = qJD(2) - t36;
t87 = -pkin(2) * t35 / 0.2e1 - t24 / 0.2e1;
t37 = t54 * V_base(5) + t56 * V_base(4);
t51 = V_base(6) + qJD(1);
t27 = -t37 * t53 + t51 * t55;
t86 = -t27 / 0.2e1;
t28 = t37 * t55 + t51 * t53;
t85 = t28 / 0.2e1;
t84 = -t31 / 0.2e1;
t83 = -t37 / 0.2e1;
t82 = t51 / 0.2e1;
t79 = pkin(1) * t53;
t77 = t39 * t53;
t41 = 0.1e1 / t43;
t76 = t41 * t53 ^ 2;
t59 = 0.1e1 / pkin(4);
t75 = t41 * t59;
t61 = 0.1e1 / pkin(3);
t74 = t41 * t61;
t72 = t55 * t64;
t42 = 0.1e1 / t43 ^ 2;
t69 = t42 * t91;
t68 = (-t33 - t34) * qJD(2) * pkin(2) * t79 / t64 * t41;
t66 = t68 / 0.2e1;
t57 = V_base(3) ^ 2;
t32 = -pkin(1) * t36 + V_base(3);
t26 = -pkin(1) * t51 + t29;
t25 = t31 * t55 + t53 * V_base(3);
t19 = -pkin(2) * t27 + t29;
t18 = t38 * t79 - t48 * t64;
t17 = pkin(2) * t77 + t47 * t64;
t14 = 0.1e1 / t16 ^ 2;
t13 = 0.1e1 / t15 ^ 2;
t10 = (t17 * t83 + t51 * t88) * t75;
t9 = (t16 * t83 + t17 * t82) * t75;
t8 = (t17 * t84 + t32 * t88) * t75;
t7 = (t16 * t84 + t17 * t32 / 0.2e1) * t75;
t6 = (t15 * t86 + t18 * t85) * t74;
t5 = (t18 * t86 + t28 * t89) * t74;
t4 = (t18 * t25 / 0.2e1 + t15 * t87) * t74;
t3 = (t18 * t87 + t25 * t89) * t74;
t2 = -t36 - 0.2e1 * ((t47 * t66 + (t62 * pkin(1) * t76 + ((t39 * t55 + t73) * t41 / 0.2e1 - t17 * t42 * t79) * pkin(2)) * qJD(2)) / t16 + (t53 * t66 + (-(-t72 + t77) * t41 / 0.2e1 + (t16 * t42 - t41 * t47) * t79) * qJD(2)) * pkin(2) * t17 * t14) * pkin(4) / (t14 * t17 ^ 2 + 0.1e1) * t43 * t59;
t1 = t35 + ((-t48 * t68 + (0.2e1 * t63 * pkin(2) * t76 + ((t38 * t55 + t73) * t41 + t18 * t53 * t69) * pkin(1)) * qJD(2)) / t15 - (-t53 * t68 + (-t41 * t72 + ((t48 * t91 + t38) * t41 + t15 * t69) * t53) * qJD(2)) * pkin(1) * t18 * t13) * pkin(3) / (t13 * t18 ^ 2 + 0.1e1) * t43 * t61;
t11 = m(2) * (t31 ^ 2 + t57 + t90) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t57) / 0.2e1 + m(3) * (t24 ^ 2 + t25 ^ 2 + t90) / 0.2e1 + m(5) * (t26 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + m(4) * (t19 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t26 * mrSges(5,2) - t8 * mrSges(5,3) + Ifges(5,1) * t9 / 0.2e1) * t9 + (-t19 * mrSges(4,1) + t3 * mrSges(4,3) + Ifges(4,2) * t6 / 0.2e1) * t6 + (-t29 * mrSges(2,1) - t31 * mrSges(2,2) + Ifges(2,3) * t82) * t51 + (t24 * mrSges(3,1) - t25 * mrSges(3,2) + Ifges(3,3) * t35 / 0.2e1) * t35 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t19 * mrSges(4,2) - t4 * mrSges(4,3) + Ifges(4,4) * t6 + Ifges(4,1) * t5 / 0.2e1) * t5 + (V_base(3) * mrSges(2,2) + t29 * mrSges(2,3) + Ifges(2,5) * t51 + Ifges(2,1) * t37 / 0.2e1) * t37 + (t29 * mrSges(3,2) - t24 * mrSges(3,3) + Ifges(3,1) * t85 + Ifges(3,5) * t35) * t28 + (t8 * mrSges(5,1) - t7 * mrSges(5,2) + Ifges(5,5) * t9 + Ifges(5,3) * t2 / 0.2e1) * t2 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t31 * mrSges(2,3) + Ifges(2,4) * t37 + Ifges(2,6) * t51 + Ifges(2,2) * t36 / 0.2e1) * t36 + (-t29 * mrSges(3,1) + t25 * mrSges(3,3) + Ifges(3,4) * t28 + Ifges(3,6) * t35 + Ifges(3,2) * t27 / 0.2e1) * t27 + (-t26 * mrSges(5,1) + t7 * mrSges(5,3) + Ifges(5,4) * t9 + Ifges(5,6) * t2 + Ifges(5,2) * t10 / 0.2e1) * t10 + (t4 * mrSges(4,1) - t3 * mrSges(4,2) + Ifges(4,5) * t5 + Ifges(4,6) * t6 + Ifges(4,3) * t1 / 0.2e1) * t1;
T = t11;
