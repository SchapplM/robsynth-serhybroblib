% Calculate kinetic energy for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh2m1OL_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh2m1OL_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_energykin_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1OL_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1OL_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'palh2m1OL_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:59:00
% EndTime: 2020-05-02 23:59:02
% DurationCPUTime: 2.33s
% Computational Cost: add. (1464->186), mult. (2194->269), div. (0->0), fcn. (1910->10), ass. (0->74)
t102 = (m(2) * pkin(5));
t74 = sin(qJ(1));
t79 = cos(qJ(1));
t83 = t74 * V_base(2) + t79 * V_base(1);
t101 = Ifges(2,3) * qJD(1);
t86 = (V_base(5) * V_base(4));
t72 = sin(qJ(3));
t73 = sin(qJ(2));
t77 = cos(qJ(3));
t78 = cos(qJ(2));
t46 = t72 * t78 + t73 * t77;
t97 = t73 * t72;
t47 = t77 * t78 - t97;
t71 = sin(qJ(4));
t76 = cos(qJ(4));
t27 = t46 * t76 + t47 * t71;
t29 = t46 * t71 - t47 * t76;
t49 = pkin(1) * t78 - pkin(5) * t73 + pkin(2);
t50 = pkin(1) * t73 + pkin(5) * t78;
t89 = t49 * t77 - t50 * t72;
t32 = pkin(3) + t89;
t33 = t49 * t72 + t50 * t77;
t45 = -V_base(4) * t74 + t79 * V_base(5);
t57 = pkin(2) * t77 + pkin(3);
t93 = pkin(3) * qJD(3);
t99 = pkin(2) * t72;
t6 = -t27 * t83 + qJD(2) * (t57 * t76 - t71 * t99) + (t32 * t76 - t33 * t71) * t45 + t29 * V_base(3) + t76 * t93;
t98 = pkin(3) * t72;
t96 = (2 * mrSges(2,3) + t102) * pkin(5);
t67 = V_base(5) ^ 2;
t68 = V_base(4) ^ 2;
t95 = t67 - t68;
t94 = pkin(2) * qJD(2);
t62 = V_base(6) + qJD(1);
t87 = -t74 * V_base(1) + V_base(2) * t79;
t85 = mrSges(2,1) * t79 - mrSges(2,2) * t74;
t84 = -mrSges(2,1) * t74 - mrSges(2,2) * t79;
t48 = V_base(5) * t74 + V_base(4) * t79;
t34 = -t48 * t73 - t62 * t78;
t35 = t48 * t78 - t62 * t73;
t18 = t34 * t77 - t35 * t72;
t19 = t34 * t72 + t35 * t77;
t11 = t18 * t76 - t19 * t71;
t7 = -t27 * V_base(3) + (t57 * t71 + t76 * t99) * qJD(2) + t71 * t93 + t45 * (t32 * t71 + t33 * t76) - t83 * t29;
t43 = qJD(2) + t45;
t42 = qJD(3) + t43;
t56 = pkin(3) * t77 + pkin(2);
t13 = -t48 * pkin(5) + (t48 * t56 - t98 * V_base(6)) * t73 + (t48 * t98 + t56 * V_base(6)) * t78 + (-pkin(3) * t97 + t56 * t78 + pkin(1)) * qJD(1) + V_base(6) * pkin(1) + t87;
t75 = cos(qJ(5));
t70 = sin(qJ(5));
t69 = -Ifges(2,2) + Ifges(2,1);
t64 = -pkin(5) * mrSges(2,1) + Ifges(2,5);
t63 = mrSges(2,2) * pkin(5) - Ifges(2,6);
t52 = -V_base(4) * pkin(5) + V_base(2);
t51 = V_base(5) * pkin(5) + V_base(1);
t40 = -pkin(1) * t45 + V_base(3);
t39 = qJD(4) + t42;
t36 = t51 * t79 + t52 * t74;
t31 = pkin(1) * t62 - t51 * t74 + t52 * t79;
t22 = t36 * t78 - t40 * t73;
t21 = -t36 * t73 - t40 * t78;
t20 = t62 * (pkin(2) * t78 + pkin(1)) + t48 * (pkin(2) * t73 - pkin(5)) + t87;
t15 = t33 * t45 - t46 * V_base(3) + t47 * t83 + t72 * t94;
t14 = t45 * t89 - t46 * t83 - t47 * V_base(3) + t77 * t94;
t12 = t18 * t71 + t19 * t76;
t10 = qJD(5) - t11;
t9 = t12 * t75 + t39 * t70;
t8 = -t12 * t70 + t39 * t75;
t5 = pkin(6) * t39 + t7;
t4 = -pkin(4) * t39 - t6;
t3 = -pkin(4) * t11 - pkin(6) * t12 + t13;
t2 = t3 * t70 + t5 * t75;
t1 = t3 * t75 - t5 * t70;
t16 = (-t31 * mrSges(3,1) + t22 * mrSges(3,3) + Ifges(3,4) * t35 + Ifges(3,6) * t43 + Ifges(3,2) * t34 / 0.2e1) * t34 + (-t20 * mrSges(4,1) + t15 * mrSges(4,3) + Ifges(4,4) * t19 + Ifges(4,6) * t42 + Ifges(4,2) * t18 / 0.2e1) * t18 + ((-t63 * t79 + t64 * t74) * V_base(5) + t84 * V_base(1) + t85 * V_base(2) + (t63 * t74 + t64 * t79) * V_base(4) + t101 / 0.2e1) * qJD(1) + ((t63 * V_base(4) + t64 * V_base(5)) * t74 + V_base(5) * Ifges(1,6) + (-t63 * V_base(5) + t64 * V_base(4)) * t79 + Ifges(1,5) * V_base(4) + t101 + (-mrSges(1,2) + t84) * V_base(1) + (mrSges(1,1) + t85) * V_base(2) + (Ifges(1,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + (-V_base(4) * V_base(2) + V_base(5) * V_base(1)) * (mrSges(1,3) + mrSges(2,3) + t102) + (t6 * mrSges(5,1) - t7 * mrSges(5,2) + Ifges(5,3) * t39 / 0.2e1) * t39 + (-t13 * mrSges(5,1) + t7 * mrSges(5,3) + Ifges(5,4) * t12 + Ifges(5,6) * t39 + Ifges(5,2) * t11 / 0.2e1) * t11 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t9 + Ifges(6,6) * t8 + Ifges(6,3) * t10 / 0.2e1) * t10 + ((4 * Ifges(2,4) * t86) - t95 * t69) * t79 ^ 2 / 0.2e1 + (Ifges(2,1) + Ifges(1,2) + t96) * t67 / 0.2e1 + (Ifges(1,1) + Ifges(2,2) + t96) * t68 / 0.2e1 + ((Ifges(1,4) - Ifges(2,4)) * t86) + (t20 * mrSges(4,2) - t14 * mrSges(4,3) + Ifges(4,5) * t42 + Ifges(4,1) * t19 / 0.2e1) * t19 + m(4) * (t14 ^ 2 + t15 ^ 2 + t20 ^ 2) / 0.2e1 + m(3) * (t21 ^ 2 + t22 ^ 2 + t31 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) / 0.2e1 + m(5) * (t13 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + (-V_base(5) * mrSges(1,1) + V_base(4) * mrSges(1,2) + (-V_base(5) * mrSges(2,1) + V_base(4) * mrSges(2,2)) * t79 + (V_base(4) * mrSges(2,1) + V_base(5) * mrSges(2,2)) * t74) * V_base(3) + (-t4 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t9 + Ifges(6,2) * t8 / 0.2e1) * t8 + (t21 * mrSges(3,1) - t22 * mrSges(3,2) + Ifges(3,3) * t43 / 0.2e1) * t43 + (t4 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,1) * t9 / 0.2e1) * t9 + (t31 * mrSges(3,2) - t21 * mrSges(3,3) + Ifges(3,5) * t43 + Ifges(3,1) * t35 / 0.2e1) * t35 + (t13 * mrSges(5,2) - t6 * mrSges(5,3) + Ifges(5,5) * t39 + Ifges(5,1) * t12 / 0.2e1) * t12 + (t14 * mrSges(4,1) - t15 * mrSges(4,2) + Ifges(4,3) * t42 / 0.2e1) * t42 + (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) * (m(1) / 0.2e1 + m(2) / 0.2e1) + (Ifges(2,4) * t95 + t69 * t86) * t74 * t79;
T = t16;
