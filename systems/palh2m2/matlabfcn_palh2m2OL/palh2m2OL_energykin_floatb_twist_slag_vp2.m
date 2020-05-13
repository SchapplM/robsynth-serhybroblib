% Calculate kinetic energy for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh2m2OL_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'palh2m2OL_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh2m2OL_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_energykin_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2OL_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2OL_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'palh2m2OL_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:13:41
% EndTime: 2020-05-03 01:13:44
% DurationCPUTime: 2.98s
% Computational Cost: add. (2572->221), mult. (4009->325), div. (0->0), fcn. (3780->12), ass. (0->89)
t86 = sin(qJ(1));
t92 = cos(qJ(1));
t117 = (Ifges(2,5) * V_base(4) + Ifges(2,6) * V_base(5)) * t92 + (V_base(5) * Ifges(2,5) - Ifges(2,6) * V_base(4)) * t86;
t96 = V_base(5) * V_base(4);
t116 = Ifges(2,3) / 0.2e1;
t102 = V_base(5) * t86;
t103 = V_base(4) * t92;
t62 = t102 + t103;
t76 = V_base(6) + qJD(1);
t85 = sin(qJ(2));
t91 = cos(qJ(2));
t44 = t62 * t91 + t76 * t85;
t84 = sin(qJ(3));
t115 = t44 * t84;
t83 = sin(qJ(4));
t110 = t83 * t84;
t60 = -t86 * V_base(4) + t92 * V_base(5);
t77 = V_base(5) ^ 2;
t78 = V_base(4) ^ 2;
t109 = t77 - t78;
t108 = pkin(5) * qJD(4);
t107 = pkin(1) * t85;
t106 = pkin(2) * t83;
t104 = pkin(2) * qJD(3);
t58 = qJD(2) - t60;
t43 = -t62 * t85 + t76 * t91;
t90 = cos(qJ(3));
t27 = t43 * t90 - t115;
t28 = t43 * t84 + t44 * t90;
t89 = cos(qJ(4));
t16 = t27 * t89 - t28 * t83;
t17 = t27 * t83 + t28 * t89;
t82 = sin(qJ(5));
t88 = cos(qJ(5));
t10 = t88 * t16 - t17 * t82;
t59 = t84 * t91 + t85 * t90;
t61 = -t84 * t85 + t90 * t91;
t39 = -t59 * t83 + t89 * t61;
t41 = t59 * t89 + t61 * t83;
t101 = t88 * t39 - t41 * t82;
t71 = pkin(1) * t91 + pkin(4);
t48 = -t107 * t84 + t71 * t90 + pkin(2);
t52 = t107 * t90 + t71 * t84;
t100 = t48 * t89 - t52 * t83;
t53 = -pkin(1) * t60 + V_base(3);
t63 = t86 * V_base(2) + t92 * V_base(1);
t36 = t91 * t53 - t63 * t85;
t98 = t86 * V_base(1) - t92 * V_base(2);
t69 = pkin(4) * t90 + pkin(2);
t97 = -pkin(4) * t110 + t69 * t89;
t57 = qJD(3) + t58;
t95 = mrSges(2,1) * t92 - mrSges(2,2) * t86;
t94 = -mrSges(2,1) * t86 - mrSges(2,2) * t92;
t55 = qJD(4) + t57;
t50 = -pkin(1) * t76 + t98;
t67 = pkin(5) * t89 + pkin(2);
t35 = pkin(2) * t102 + t67 * t103 + (t102 * t89 + t83 * V_base(6)) * pkin(5);
t42 = -V_base(6) * pkin(2) + (t62 * t83 - t89 * V_base(6)) * pkin(5);
t12 = -pkin(1) * V_base(6) + (pkin(4) * t62 + t35 * t90 - t42 * t84) * t85 + (-V_base(6) * pkin(4) + t35 * t84 + t42 * t90) * t91 + ((pkin(5) * t110 - t67 * t90 - pkin(4)) * t91 + (pkin(5) * t83 * t90 + t67 * t84) * t85 - pkin(1)) * qJD(1) + t98;
t25 = t39 * t82 + t41 * t88;
t31 = pkin(5) + t100;
t34 = t48 * t83 + t52 * t89;
t51 = pkin(5) + t97;
t54 = pkin(4) * t84 * t89 + t69 * t83;
t68 = pkin(2) * t89 + pkin(5);
t5 = -(t31 * t88 - t34 * t82) * t60 - t25 * t63 + t101 * V_base(3) + (t51 * t88 - t54 * t82) * qJD(2) + (-t106 * t82 + t68 * t88) * qJD(3) + t88 * t108;
t87 = cos(qJ(6));
t81 = sin(qJ(6));
t79 = Ifges(2,1) - Ifges(2,2);
t70 = pkin(2) * t90 + pkin(4);
t49 = qJD(5) + t55;
t37 = t53 * t85 + t63 * t91;
t33 = -pkin(4) * t43 + t50;
t29 = pkin(4) * t58 + t36;
t21 = t29 * t84 + t37 * t90;
t20 = t29 * t90 - t37 * t84;
t18 = (-t70 * t91 - pkin(1)) * qJD(1) + (pkin(4) * t102 + t103 * t70) * t85 + (-pkin(4) * t91 - pkin(1)) * V_base(6) + ((t102 * t85 - t91 * V_base(6)) * t90 + t115) * pkin(2) + t98;
t14 = qJD(2) * t97 - t100 * t60 + t104 * t89 + t39 * V_base(3) - t41 * t63;
t13 = qJD(2) * t54 + t104 * t83 - t34 * t60 + t39 * t63 + t41 * V_base(3);
t11 = t16 * t82 + t17 * t88;
t9 = qJD(6) + t10;
t8 = t11 * t87 - t49 * t81;
t7 = -t11 * t81 - t49 * t87;
t6 = t25 * V_base(3) + (t51 * t82 + t54 * t88) * qJD(2) + (t106 * t88 + t68 * t82) * qJD(3) + t82 * t108 - (t31 * t82 + t34 * t88) * t60 + t63 * t101;
t4 = pkin(3) * t49 + t5;
t3 = -pkin(3) * t10 + t12;
t2 = -t3 * t81 + t6 * t87;
t1 = -t3 * t87 - t6 * t81;
t15 = (t14 * mrSges(5,1) - t13 * mrSges(5,2) + Ifges(5,3) * t55 / 0.2e1) * t55 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t49 / 0.2e1) * t49 + (-V_base(4) * V_base(2) + V_base(5) * V_base(1)) * (mrSges(1,3) + mrSges(2,3)) + (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) * (m(1) / 0.2e1 + m(2) / 0.2e1) + (-V_base(5) * mrSges(1,1) + V_base(4) * mrSges(1,2) + (V_base(4) * mrSges(2,1) + V_base(5) * mrSges(2,2)) * t86) * V_base(3) + (-t18 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t17 + Ifges(5,6) * t55 + Ifges(5,2) * t16 / 0.2e1) * t16 + (-t33 * mrSges(4,1) + t21 * mrSges(4,3) + Ifges(4,4) * t28 + Ifges(4,6) * t57 + Ifges(4,2) * t27 / 0.2e1) * t27 + (-t12 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t11 + Ifges(6,6) * t49 + Ifges(6,2) * t10 / 0.2e1) * t10 + (-t50 * mrSges(3,1) + t37 * mrSges(3,3) + Ifges(3,4) * t44 + Ifges(3,6) * t58 + Ifges(3,2) * t43 / 0.2e1) * t43 + (t33 * mrSges(4,2) - t20 * mrSges(4,3) + Ifges(4,5) * t57 + Ifges(4,1) * t28 / 0.2e1) * t28 + (t20 * mrSges(4,1) - t21 * mrSges(4,2) + Ifges(4,3) * t57 / 0.2e1) * t57 + (t12 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t49 + Ifges(6,1) * t11 / 0.2e1) * t11 + ((-mrSges(1,2) + t94) * V_base(1) + (mrSges(1,1) + t95) * V_base(2) + V_base(5) * Ifges(1,6) + Ifges(2,3) * qJD(1) + Ifges(1,5) * V_base(4) + (Ifges(1,3) / 0.2e1 + t116) * V_base(6) + t117) * V_base(6) + (t116 * qJD(1) + t94 * V_base(1) + t95 * V_base(2) + t117) * qJD(1) + m(4) * (t20 ^ 2 + t21 ^ 2 + t33 ^ 2) / 0.2e1 + m(5) * (t13 ^ 2 + t14 ^ 2 + t18 ^ 2) / 0.2e1 + (t18 * mrSges(5,2) - t14 * mrSges(5,3) + Ifges(5,5) * t55 + Ifges(5,1) * t17 / 0.2e1) * t17 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t9 / 0.2e1) * t9 + m(3) * (t36 ^ 2 + t37 ^ 2 + t50 ^ 2) / 0.2e1 + (t50 * mrSges(3,2) - t36 * mrSges(3,3) + Ifges(3,5) * t58 + Ifges(3,1) * t44 / 0.2e1) * t44 + (Ifges(2,1) + Ifges(1,2)) * t77 / 0.2e1 + (Ifges(1,1) + Ifges(2,2)) * t78 / 0.2e1 + ((Ifges(2,4) * t109 + t79 * t96) * t86 + (-V_base(5) * mrSges(2,1) + V_base(4) * mrSges(2,2)) * V_base(3) + (0.2e1 * Ifges(2,4) * t96 - t109 * t79 / 0.2e1) * t92) * t92 + (t4 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t9 + Ifges(7,1) * t8 / 0.2e1) * t8 + (t36 * mrSges(3,1) - t37 * mrSges(3,2) + Ifges(3,3) * t58 / 0.2e1) * t58 + m(6) * (t12 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) / 0.2e1 + (-t4 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t8 + Ifges(7,6) * t9 + Ifges(7,2) * t7 / 0.2e1) * t7 + (Ifges(1,4) - Ifges(2,4)) * t96;
T = t15;
