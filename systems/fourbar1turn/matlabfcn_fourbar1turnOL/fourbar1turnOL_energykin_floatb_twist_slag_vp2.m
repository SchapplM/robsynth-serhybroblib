% Calculate kinetic energy for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
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
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar1turnOL_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fourbar1turnOL_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_energykin_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnOL_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnOL_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:46
% EndTime: 2020-04-12 19:40:46
% DurationCPUTime: 0.63s
% Computational Cost: add. (435->100), mult. (601->151), div. (0->0), fcn. (392->8), ass. (0->37)
t26 = V_base(5) * pkin(5) + V_base(1);
t27 = -V_base(4) * pkin(5) + V_base(2);
t34 = sin(qJ(1));
t38 = cos(qJ(1));
t16 = t26 * t34 - t38 * t27;
t40 = t16 ^ 2;
t18 = t26 * t38 + t27 * t34;
t33 = sin(qJ(2));
t37 = cos(qJ(2));
t9 = -t18 * t33 + t37 * V_base(3);
t23 = -t34 * V_base(4) + t38 * V_base(5);
t22 = qJD(2) - t23;
t39 = V_base(3) ^ 2;
t36 = cos(qJ(3));
t35 = cos(qJ(4));
t32 = sin(qJ(3));
t31 = sin(qJ(4));
t30 = V_base(6) + qJD(1);
t24 = t34 * V_base(5) + t38 * V_base(4);
t21 = qJD(4) - t23;
t20 = qJD(3) + t22;
t19 = -pkin(1) * t23 + V_base(3);
t15 = t24 * t37 + t30 * t33;
t14 = t24 * t35 + t30 * t31;
t13 = -t24 * t33 + t30 * t37;
t12 = -t24 * t31 + t30 * t35;
t11 = -pkin(1) * t30 + t16;
t10 = t18 * t37 + t33 * V_base(3);
t8 = t18 * t35 + t19 * t31;
t7 = -t18 * t31 + t19 * t35;
t6 = pkin(2) * t22 + t9;
t5 = -pkin(2) * t13 + t16;
t4 = -t13 * t32 - t15 * t36;
t3 = -t13 * t36 + t15 * t32;
t2 = -t10 * t36 - t32 * t6;
t1 = t10 * t32 - t36 * t6;
t17 = m(2) * (t18 ^ 2 + t39 + t40) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t39) / 0.2e1 + m(3) * (t10 ^ 2 + t9 ^ 2 + t40) / 0.2e1 + m(4) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t5 * mrSges(4,2) - t1 * mrSges(4,3) + Ifges(4,1) * t4 / 0.2e1) * t4 + (-t16 * mrSges(2,1) - t18 * mrSges(2,2) + Ifges(2,3) * t30 / 0.2e1) * t30 + (t9 * mrSges(3,1) - t10 * mrSges(3,2) + Ifges(3,3) * t22 / 0.2e1) * t22 + (t7 * mrSges(5,1) - t8 * mrSges(5,2) + Ifges(5,3) * t21 / 0.2e1) * t21 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t5 * mrSges(4,1) + t2 * mrSges(4,3) + Ifges(4,4) * t4 + Ifges(4,2) * t3 / 0.2e1) * t3 + (V_base(3) * mrSges(2,2) + t16 * mrSges(2,3) + Ifges(2,5) * t30 + Ifges(2,1) * t24 / 0.2e1) * t24 + (t16 * mrSges(3,2) - t9 * mrSges(3,3) + Ifges(3,5) * t22 + Ifges(3,1) * t15 / 0.2e1) * t15 + (t11 * mrSges(5,2) - t7 * mrSges(5,3) + Ifges(5,5) * t21 + Ifges(5,1) * t14 / 0.2e1) * t14 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t18 * mrSges(2,3) + Ifges(2,4) * t24 + Ifges(2,6) * t30 + Ifges(2,2) * t23 / 0.2e1) * t23 + (t1 * mrSges(4,1) - t2 * mrSges(4,2) + Ifges(4,5) * t4 + Ifges(4,6) * t3 + Ifges(4,3) * t20 / 0.2e1) * t20 + (-t16 * mrSges(3,1) + t10 * mrSges(3,3) + Ifges(3,4) * t15 + Ifges(3,6) * t22 + Ifges(3,2) * t13 / 0.2e1) * t13 + (-t11 * mrSges(5,1) + t8 * mrSges(5,3) + Ifges(5,4) * t14 + Ifges(5,6) * t21 + Ifges(5,2) * t12 / 0.2e1) * t12;
T = t17;
