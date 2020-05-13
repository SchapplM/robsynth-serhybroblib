% Calculate kinetic energy for
% fivebar1OL
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
%   pkin=[AB,AE,BC,CD,ED]';
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
% Datum: 2020-04-27 06:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fivebar1OL_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1OL_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fivebar1OL_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fivebar1OL_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1OL_energykin_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1OL_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fivebar1OL_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fivebar1OL_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:12:54
% EndTime: 2020-04-27 06:12:56
% DurationCPUTime: 1.23s
% Computational Cost: add. (290->114), mult. (473->165), div. (0->0), fcn. (248->8), ass. (0->43)
t35 = sin(qJ(1));
t39 = cos(qJ(1));
t51 = (Ifges(2,5) * V_base(4) + Ifges(2,6) * V_base(5)) * t39 + (V_base(5) * Ifges(2,5) - Ifges(2,6) * V_base(4)) * t35;
t43 = V_base(5) * V_base(4);
t50 = Ifges(2,3) / 0.2e1;
t28 = V_base(5) ^ 2;
t29 = V_base(4) ^ 2;
t47 = t28 - t29;
t45 = V_base(6) + qJD(1);
t27 = V_base(6) + qJD(3);
t23 = V_base(6) * pkin(1) + V_base(2);
t33 = sin(qJ(3));
t37 = cos(qJ(3));
t13 = t37 * t23 - t33 * V_base(1);
t42 = mrSges(2,1) * t39 - mrSges(2,2) * t35;
t41 = -mrSges(2,1) * t35 - mrSges(2,2) * t39;
t24 = -V_base(5) * pkin(1) + V_base(3);
t38 = cos(qJ(2));
t36 = cos(qJ(4));
t34 = sin(qJ(2));
t32 = sin(qJ(4));
t30 = Ifges(2,1) - Ifges(2,2);
t26 = qJD(2) + t45;
t25 = qJD(4) + t27;
t19 = t35 * V_base(2) + t39 * V_base(1);
t18 = V_base(5) * t35 + V_base(4) * t39;
t17 = t33 * V_base(5) + t37 * V_base(4);
t16 = -V_base(4) * t35 + V_base(5) * t39;
t15 = -t33 * V_base(4) + t37 * V_base(5);
t14 = t23 * t33 + t37 * V_base(1);
t12 = -pkin(2) * t16 + V_base(3);
t11 = pkin(2) * t45 - t35 * V_base(1) + t39 * V_base(2);
t10 = -pkin(3) * t15 + t24;
t9 = pkin(3) * t27 + t13;
t8 = -t16 * t34 - t18 * t38;
t7 = t15 * t36 - t17 * t32;
t6 = t15 * t32 + t17 * t36;
t5 = -t16 * t38 + t18 * t34;
t4 = -t11 * t34 - t19 * t38;
t3 = -t11 * t38 + t19 * t34;
t2 = t14 * t36 + t32 * t9;
t1 = -t14 * t32 + t36 * t9;
t20 = (Ifges(2,1) + Ifges(1,2)) * t28 / 0.2e1 + (Ifges(1,1) + Ifges(2,2)) * t29 / 0.2e1 + m(4) * (t13 ^ 2 + t14 ^ 2 + t24 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t10 ^ 2 + t2 ^ 2) / 0.2e1 + m(3) * (t12 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + (-V_base(4) * V_base(2) + V_base(5) * V_base(1)) * (mrSges(1,3) + mrSges(2,3)) + (t12 * mrSges(3,2) - t3 * mrSges(3,3) + Ifges(3,1) * t8 / 0.2e1) * t8 + (-t10 * mrSges(5,1) + t2 * mrSges(5,3) + Ifges(5,2) * t7 / 0.2e1) * t7 + (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) * (m(1) / 0.2e1 + m(2) / 0.2e1) + (t13 * mrSges(4,1) - t14 * mrSges(4,2) + Ifges(4,3) * t27 / 0.2e1) * t27 + (t10 * mrSges(5,2) - t1 * mrSges(5,3) + Ifges(5,4) * t7 + Ifges(5,1) * t6 / 0.2e1) * t6 + (-t12 * mrSges(3,1) + t4 * mrSges(3,3) + Ifges(3,4) * t8 + Ifges(3,2) * t5 / 0.2e1) * t5 + (t24 * mrSges(4,2) - t13 * mrSges(4,3) + Ifges(4,5) * t27 + Ifges(4,1) * t17 / 0.2e1) * t17 + (-V_base(5) * mrSges(1,1) + V_base(4) * mrSges(1,2) + (V_base(4) * mrSges(2,1) + V_base(5) * mrSges(2,2)) * t35) * V_base(3) + (Ifges(1,4) - Ifges(2,4)) * t43 + (t3 * mrSges(3,1) - t4 * mrSges(3,2) + Ifges(3,5) * t8 + Ifges(3,6) * t5 + Ifges(3,3) * t26 / 0.2e1) * t26 + (t1 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t6 + Ifges(5,6) * t7 + Ifges(5,3) * t25 / 0.2e1) * t25 + (-t24 * mrSges(4,1) + t14 * mrSges(4,3) + Ifges(4,4) * t17 + Ifges(4,6) * t27 + Ifges(4,2) * t15 / 0.2e1) * t15 + ((Ifges(2,4) * t47 + t30 * t43) * t35 + (-V_base(5) * mrSges(2,1) + V_base(4) * mrSges(2,2)) * V_base(3) + (0.2e1 * Ifges(2,4) * t43 - t47 * t30 / 0.2e1) * t39) * t39 + (t50 * qJD(1) + t41 * V_base(1) + t42 * V_base(2) + t51) * qJD(1) + (V_base(5) * Ifges(1,6) + (-mrSges(1,2) + t41) * V_base(1) + (mrSges(1,1) + t42) * V_base(2) + Ifges(2,3) * qJD(1) + Ifges(1,5) * V_base(4) + (Ifges(1,3) / 0.2e1 + t50) * V_base(6) + t51) * V_base(6);
T = t20;
