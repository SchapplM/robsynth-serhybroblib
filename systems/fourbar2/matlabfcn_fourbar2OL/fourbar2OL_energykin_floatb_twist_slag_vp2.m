% Calculate kinetic energy for
% fourbar2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2]';
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
% Datum: 2020-04-24 20:32
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar2OL_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(2,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar2OL_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbar2OL_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fourbar2OL_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2OL_energykin_floatb_twist_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar2OL_energykin_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar2OL_energykin_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbar2OL_energykin_floatb_twist_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:32:16
% EndTime: 2020-04-24 20:32:16
% DurationCPUTime: 0.73s
% Computational Cost: add. (190->92), mult. (328->132), div. (0->0), fcn. (148->6), ass. (0->34)
t26 = sin(qJ(1));
t29 = cos(qJ(1));
t41 = (Ifges(2,5) * V_base(4) + Ifges(2,6) * V_base(5)) * t29 + (V_base(5) * Ifges(2,5) - Ifges(2,6) * V_base(4)) * t26;
t33 = V_base(5) * V_base(4);
t40 = Ifges(2,3) / 0.2e1;
t20 = V_base(5) ^ 2;
t21 = V_base(4) ^ 2;
t37 = t20 - t21;
t35 = V_base(6) + qJD(1);
t32 = mrSges(2,1) * t29 - mrSges(2,2) * t26;
t31 = -mrSges(2,1) * t26 - mrSges(2,2) * t29;
t28 = cos(qJ(2));
t27 = cos(qJ(3));
t25 = sin(qJ(2));
t24 = sin(qJ(3));
t22 = Ifges(2,1) - Ifges(2,2);
t19 = V_base(6) + qJD(3);
t18 = qJD(2) + t35;
t17 = -V_base(5) * pkin(1) + V_base(3);
t16 = V_base(6) * pkin(1) + V_base(2);
t13 = t26 * V_base(2) + t29 * V_base(1);
t12 = V_base(5) * t26 + V_base(4) * t29;
t11 = t24 * V_base(5) + t27 * V_base(4);
t10 = -V_base(4) * t26 + V_base(5) * t29;
t9 = -t24 * V_base(4) + t27 * V_base(5);
t8 = t16 * t24 + t27 * V_base(1);
t7 = t16 * t27 - t24 * V_base(1);
t6 = -pkin(2) * t10 + V_base(3);
t5 = t35 * pkin(2) - t26 * V_base(1) + t29 * V_base(2);
t4 = -t10 * t25 - t12 * t28;
t3 = -t10 * t28 + t12 * t25;
t2 = -t13 * t28 - t25 * t5;
t1 = t13 * t25 - t28 * t5;
t14 = m(3) * (t1 ^ 2 + t2 ^ 2 + t6 ^ 2) / 0.2e1 + m(4) * (t17 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + (Ifges(2,1) + Ifges(1,2)) * t20 / 0.2e1 + (Ifges(1,1) + Ifges(2,2)) * t21 / 0.2e1 + (-V_base(4) * V_base(2) + V_base(5) * V_base(1)) * (mrSges(1,3) + mrSges(2,3)) + (-t17 * mrSges(4,1) + t8 * mrSges(4,3) + Ifges(4,2) * t9 / 0.2e1) * t9 + (t6 * mrSges(3,2) - t1 * mrSges(3,3) + Ifges(3,1) * t4 / 0.2e1) * t4 + (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) * (m(1) / 0.2e1 + m(2) / 0.2e1) + (-t6 * mrSges(3,1) + t2 * mrSges(3,3) + Ifges(3,4) * t4 + Ifges(3,2) * t3 / 0.2e1) * t3 + (t7 * mrSges(4,1) - t8 * mrSges(4,2) + Ifges(4,6) * t9 + Ifges(4,3) * t19 / 0.2e1) * t19 + (-V_base(5) * mrSges(1,1) + V_base(4) * mrSges(1,2) + (V_base(4) * mrSges(2,1) + V_base(5) * mrSges(2,2)) * t26) * V_base(3) + (Ifges(1,4) - Ifges(2,4)) * t33 + (t1 * mrSges(3,1) - t2 * mrSges(3,2) + Ifges(3,5) * t4 + Ifges(3,6) * t3 + Ifges(3,3) * t18 / 0.2e1) * t18 + (t17 * mrSges(4,2) - t7 * mrSges(4,3) + Ifges(4,4) * t9 + Ifges(4,5) * t19 + Ifges(4,1) * t11 / 0.2e1) * t11 + ((t37 * Ifges(2,4) + t22 * t33) * t26 + (-V_base(5) * mrSges(2,1) + V_base(4) * mrSges(2,2)) * V_base(3) + (0.2e1 * Ifges(2,4) * t33 - t37 * t22 / 0.2e1) * t29) * t29 + (t40 * qJD(1) + t31 * V_base(1) + t32 * V_base(2) + t41) * qJD(1) + (V_base(5) * Ifges(1,6) + (-mrSges(1,2) + t31) * V_base(1) + (mrSges(1,1) + t32) * V_base(2) + Ifges(2,3) * qJD(1) + Ifges(1,5) * V_base(4) + (Ifges(1,3) / 0.2e1 + t40) * V_base(6) + t41) * V_base(6);
T = t14;
