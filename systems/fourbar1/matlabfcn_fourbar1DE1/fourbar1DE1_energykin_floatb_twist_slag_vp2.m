% Calculate kinetic energy for
% fourbar1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% m [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:57
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar1DE1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(6,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE1_energykin_floatb_twist_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1DE1_energykin_floatb_twist_slag_vp2: qJD has to be [1x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fourbar1DE1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE1_energykin_floatb_twist_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1DE1_energykin_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar1DE1_energykin_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbar1DE1_energykin_floatb_twist_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:55:33
% EndTime: 2020-04-24 19:55:34
% DurationCPUTime: 1.46s
% Computational Cost: add. (1340->127), mult. (2022->211), div. (76->7), fcn. (510->4), ass. (0->74)
t45 = cos(qJ(1));
t98 = (Ifges(2,5) * V_base(4) + Ifges(2,6) * V_base(5)) * t45;
t97 = Ifges(2,3) * qJD(1);
t58 = V_base(4) * V_base(5);
t96 = pkin(4) ^ 2;
t95 = pkin(3) ^ 2;
t37 = pkin(1) * t45 - pkin(2);
t39 = V_base(6) + qJD(1);
t94 = t39 * t37;
t82 = pkin(2) * t45;
t76 = (-0.2e1 * t82 + pkin(1)) * pkin(1);
t27 = pkin(2) ^ 2 + t76;
t23 = 0.1e1 / t27;
t93 = t23 ^ 2;
t44 = sin(qJ(1));
t83 = pkin(2) * t44;
t33 = V_base(4) * t83;
t36 = -pkin(1) + t82;
t16 = -t36 * V_base(5) + t33;
t84 = -pkin(3) + pkin(4);
t85 = -pkin(3) - pkin(4);
t55 = sqrt(-((pkin(2) - t85) * (pkin(2) + t85) + t76) * ((pkin(2) - t84) * (pkin(2) + t84) + t76));
t12 = t16 * t55;
t17 = -V_base(4) * pkin(1) + (t44 * V_base(5) + t45 * V_base(4)) * pkin(2);
t74 = t95 - t96;
t21 = t27 + t74;
t5 = -t17 * t21 + t12;
t92 = t5 / 0.2e1;
t22 = t27 - t74;
t6 = t17 * t22 + t12;
t91 = t6 / 0.2e1;
t13 = t17 * t55;
t7 = t16 * t21 + t13;
t90 = t7 / 0.2e1;
t8 = -t16 * t22 + t13;
t89 = t8 / 0.2e1;
t88 = Ifges(3,4) / 0.2e1;
t87 = Ifges(4,4) / 0.2e1;
t15 = 0.1e1 / t55;
t73 = pkin(2) * qJD(1);
t11 = (-t27 * V_base(6) + t37 * t73) * t55;
t60 = pkin(1) * t44 * t73;
t9 = -t21 * t60 + t11;
t81 = t9 * t15;
t10 = t22 * t60 + t11;
t80 = t10 * t15;
t48 = 0.1e1 / pkin(4);
t79 = t48 * mrSges(4,3);
t51 = 0.1e1 / pkin(3);
t78 = t51 * mrSges(3,3);
t29 = t36 * V_base(2);
t69 = t44 * V_base(1);
t34 = pkin(2) * t69;
t77 = t34 - t29;
t40 = V_base(5) ^ 2;
t41 = V_base(4) ^ 2;
t75 = t40 - t41;
t28 = t36 * V_base(1);
t61 = t44 * V_base(6);
t68 = pkin(1) * pkin(2) * t61 + V_base(2) * t83 + t28;
t67 = V_base(6) * pkin(1);
t57 = mrSges(2,1) * t45 - mrSges(2,2) * t44;
t56 = -mrSges(2,1) * t44 - mrSges(2,2) * t45;
t42 = Ifges(2,1) - Ifges(2,2);
t32 = -V_base(5) * pkin(1) + V_base(3);
t31 = V_base(2) + t67;
t26 = Ifges(2,5) * V_base(5) - Ifges(2,6) * V_base(4);
t24 = 0.1e1 / t27 ^ 2;
t18 = -t82 * V_base(5) + t33 + V_base(3);
t4 = (-t31 * t36 + t34) * t55 + t22 * (t31 * t83 + t28);
t3 = t68 * t55 - t22 * (-t36 * t67 + t77);
t2 = (pkin(2) * t94 + t77) * t55 - t21 * (t28 + (pkin(1) * t39 + V_base(2)) * t83);
t1 = (t60 + t68) * t55 + t21 * (-t29 + (t69 + t94) * pkin(2));
t14 = t26 * t61 + (Ifges(2,1) + Ifges(1,2)) * t40 / 0.2e1 + (Ifges(1,1) + Ifges(2,2)) * t41 / 0.2e1 + m(4) * (t32 ^ 2 + (t4 ^ 2 / 0.4e1 + t3 ^ 2 / 0.4e1) / t96 * t24) / 0.2e1 + m(3) * (t18 ^ 2 + (t2 ^ 2 / 0.4e1 + t1 ^ 2 / 0.4e1) / t95 * t24) / 0.2e1 - (t9 * (-Ifges(4,3) * t81 + (Ifges(4,5) * t91 + Ifges(4,6) * t89) * t48) + t10 * (-Ifges(3,3) * t80 + (Ifges(3,5) * t92 + Ifges(3,6) * t90) * t51)) * t15 * t93 / 0.2e1 + (-V_base(4) * V_base(2) + V_base(5) * V_base(1)) * (mrSges(1,3) + mrSges(2,3)) + (t32 * (mrSges(4,2) * t91 - t8 * mrSges(4,1) / 0.2e1) * t48 + t18 * (mrSges(3,2) * t92 - t7 * mrSges(3,1) / 0.2e1) * t51) * t23 + (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) * (m(1) / 0.2e1 + m(2) / 0.2e1) + (-V_base(5) * mrSges(1,1) + V_base(4) * mrSges(1,2) + (mrSges(2,1) * V_base(4) + mrSges(2,2) * V_base(5)) * t44) * V_base(3) + (Ifges(1,4) - Ifges(2,4)) * t58 + ((t3 * (-t6 * t79 / 0.2e1 - mrSges(4,1) * t81) + t4 * (mrSges(4,2) * t81 + t79 * t89)) * t48 + (t1 * (-t5 * t78 / 0.2e1 - mrSges(3,1) * t80) + t2 * (mrSges(3,2) * t80 + t78 * t90)) * t51) * t93 / 0.2e1 + ((t8 * (-Ifges(4,6) * t81 + (Ifges(4,2) * t89 + t6 * t87) * t48) + t6 * (-Ifges(4,5) * t81 + (Ifges(4,1) * t91 + t8 * t87) * t48)) * t48 + (t7 * (-Ifges(3,6) * t80 + (Ifges(3,2) * t90 + t5 * t88) * t51) + t5 * (-Ifges(3,5) * t80 + (Ifges(3,1) * t92 + t7 * t88) * t51)) * t51) * t93 / 0.4e1 + ((-mrSges(2,1) * V_base(5) + mrSges(2,2) * V_base(4)) * V_base(3) + (Ifges(2,4) * t75 + t42 * t58) * t44 + (0.2e1 * Ifges(2,4) * t58 - t75 * t42 / 0.2e1) * t45) * t45 + (t98 + t26 * t44 + t56 * V_base(1) + t57 * V_base(2) + t97 / 0.2e1) * qJD(1) + (V_base(5) * Ifges(1,6) + (-mrSges(1,2) + t56) * V_base(1) + (mrSges(1,1) + t57) * V_base(2) + t97 + Ifges(1,5) * V_base(4) + t98 + (Ifges(1,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t14;
