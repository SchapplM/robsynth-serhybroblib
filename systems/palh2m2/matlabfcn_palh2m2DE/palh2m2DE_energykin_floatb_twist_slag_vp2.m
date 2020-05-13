% Calculate kinetic energy for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
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
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh2m2DE_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh2m2DE_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_energykin_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2DE_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2DE_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'palh2m2DE_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:20
% EndTime: 2020-05-03 01:06:23
% DurationCPUTime: 2.79s
% Computational Cost: add. (602->167), mult. (941->231), div. (0->0), fcn. (564->8), ass. (0->60)
t48 = V_base(6) + qJD(1);
t59 = sin(qJ(1));
t63 = cos(qJ(1));
t31 = t59 * V_base(4) - t63 * V_base(5);
t32 = t59 * V_base(5) + t63 * V_base(4);
t61 = cos(qJ(3));
t62 = cos(qJ(2));
t91 = -t62 * pkin(4) - pkin(1);
t92 = -pkin(2) + t91;
t98 = -t61 * pkin(5) + t92;
t97 = (Ifges(2,5) * V_base(4) + Ifges(2,6) * V_base(5)) * t63 + (V_base(5) * Ifges(2,5) - Ifges(2,6) * V_base(4)) * t59;
t70 = V_base(5) * V_base(4);
t96 = Ifges(2,3) / 0.2e1;
t34 = t59 * V_base(2) + t63 * V_base(1);
t58 = sin(qJ(2));
t86 = pkin(4) * t58;
t57 = sin(qJ(3));
t83 = t57 * t58;
t82 = t58 * t61;
t52 = V_base(5) ^ 2;
t53 = V_base(4) ^ 2;
t81 = t52 - t53;
t80 = pkin(4) * qJD(2);
t79 = pkin(5) * qJD(3);
t77 = t62 * t80 + V_base(3);
t71 = t59 * V_base(1) - t63 * V_base(2);
t69 = mrSges(2,1) * t63 - mrSges(2,2) * t59;
t68 = -mrSges(2,1) * t59 - mrSges(2,2) * t63;
t8 = -t31 * t98 + t61 * t79 + t77;
t66 = t32 * t86 + t71;
t39 = pkin(5) * t57 + t86;
t5 = t32 * t39 + t48 * t98 + t71;
t60 = cos(qJ(4));
t56 = sin(qJ(4));
t54 = Ifges(2,1) - Ifges(2,2);
t41 = qJD(4) + t48;
t29 = qJD(2) + t31;
t28 = qJD(3) + t31;
t25 = pkin(1) * t31 + V_base(3);
t24 = -pkin(1) * t48 + t71;
t21 = t32 * t62 + t48 * t58;
t20 = t32 * t61 + t48 * t57;
t19 = -t32 * t58 + t48 * t62;
t18 = -t32 * t57 + t48 * t61;
t17 = -t29 * t86 + t34;
t16 = -t31 * t60 - t32 * t56;
t15 = -t31 * t56 + t32 * t60;
t14 = t25 * t58 + t34 * t62;
t13 = t25 * t62 - t34 * t58;
t12 = -t31 * t91 + t77;
t11 = -t31 * t39 - t57 * t79 - t58 * t80 + t34;
t10 = t48 * t91 + t66;
t9 = t48 * t92 + t66;
t7 = (t61 * t62 + t83) * t80 + t61 * V_base(3) + t31 * (pkin(4) * t83 - t61 * t92) - t34 * t57;
t6 = (t57 * t62 - t82) * t80 + t57 * V_base(3) + t31 * (-pkin(4) * t82 - t57 * t92) + t34 * t61;
t4 = -pkin(3) * t48 + t5;
t3 = pkin(3) * t31 + t8;
t2 = t11 * t60 - t4 * t56;
t1 = -t11 * t56 - t4 * t60;
t22 = (-t9 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t20 + Ifges(5,6) * t28 + Ifges(5,2) * t18 / 0.2e1) * t18 + (-V_base(4) * V_base(2) + V_base(5) * V_base(1)) * (mrSges(1,3) + mrSges(2,3)) + (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) * (m(1) / 0.2e1 + m(2) / 0.2e1) + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,6) * t41 + Ifges(7,2) * t16 / 0.2e1) * t16 + (-t10 * mrSges(4,1) - t5 * mrSges(6,1) + t17 * mrSges(4,3) + t11 * mrSges(6,3) + (Ifges(4,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t48) * t48 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,4) * t16 + Ifges(7,5) * t41 + Ifges(7,1) * t15 / 0.2e1) * t15 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(5) * (t6 ^ 2 + t7 ^ 2 + t9 ^ 2) / 0.2e1 + (t7 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t28 / 0.2e1) * t28 + (Ifges(2,1) + Ifges(1,2)) * t52 / 0.2e1 + (Ifges(1,1) + Ifges(2,2)) * t53 / 0.2e1 + m(3) * (t13 ^ 2 + t14 ^ 2 + t24 ^ 2) / 0.2e1 + m(6) * (t11 ^ 2 + t5 ^ 2 + t8 ^ 2) / 0.2e1 + m(4) * (t10 ^ 2 + t12 ^ 2 + t17 ^ 2) / 0.2e1 + (t13 * mrSges(3,1) - t14 * mrSges(3,2) + Ifges(3,3) * t29 / 0.2e1) * t29 + ((-V_base(5) * mrSges(2,1) + V_base(4) * mrSges(2,2)) * V_base(3) + (Ifges(2,4) * t81 + t54 * t70) * t59 + (0.2e1 * Ifges(2,4) * t70 - t81 * t54 / 0.2e1) * t63) * t63 + ((-mrSges(1,2) + t68) * V_base(1) + (mrSges(1,1) + t69) * V_base(2) + V_base(5) * Ifges(1,6) + Ifges(2,3) * qJD(1) + Ifges(1,5) * V_base(4) + (Ifges(1,3) / 0.2e1 + t96) * V_base(6) + t97) * V_base(6) + (qJD(1) * t96 + t68 * V_base(1) + t69 * V_base(2) + t97) * qJD(1) + (t9 * mrSges(5,2) - t7 * mrSges(5,3) + Ifges(5,5) * t28 + Ifges(5,1) * t20 / 0.2e1) * t20 + (-t24 * mrSges(3,1) + t14 * mrSges(3,3) + Ifges(3,4) * t21 + Ifges(3,6) * t29 + Ifges(3,2) * t19 / 0.2e1) * t19 + (-V_base(5) * mrSges(1,1) + V_base(4) * mrSges(1,2) + (V_base(4) * mrSges(2,1) + V_base(5) * mrSges(2,2)) * t59) * V_base(3) + (t24 * mrSges(3,2) - t13 * mrSges(3,3) + Ifges(3,5) * t29 + Ifges(3,1) * t21 / 0.2e1) * t21 + (Ifges(1,4) - Ifges(2,4)) * t70 + (t12 * mrSges(4,1) + t8 * mrSges(6,1) - t17 * mrSges(4,2) - t11 * mrSges(6,2) + (Ifges(4,5) + Ifges(6,5)) * t32 + (Ifges(4,6) + Ifges(6,6)) * t48 + (Ifges(6,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t31) * t31 + (t10 * mrSges(4,2) + t5 * mrSges(6,2) - t12 * mrSges(4,3) - t8 * mrSges(6,3) + (Ifges(4,4) + Ifges(6,4)) * t48 + (Ifges(4,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t32) * t32 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t41 / 0.2e1) * t41;
T = t22;
