% Calculate kinetic energy for
% fourbar2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
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
% Datum: 2020-04-24 20:22
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar2DE1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(6,1),zeros(2,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar2DE1_energykin_floatb_twist_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar2DE1_energykin_floatb_twist_slag_vp2: qJD has to be [1x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fourbar2DE1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2DE1_energykin_floatb_twist_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar2DE1_energykin_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar2DE1_energykin_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbar2DE1_energykin_floatb_twist_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:22:52
% EndTime: 2020-04-24 20:22:53
% DurationCPUTime: 0.62s
% Computational Cost: add. (147->80), mult. (224->91), div. (0->0), fcn. (16->2), ass. (0->30)
t26 = 2 * V_base(6);
t42 = t26 / 0.2e1;
t38 = 2 * V_base(5);
t40 = t38 / 0.2e1;
t13 = (Ifges(2,6) + Ifges(4,6));
t4 = (mrSges(3,3) * pkin(2) - Ifges(2,5) - Ifges(4,5));
t39 = t13 * V_base(5) - t4 * V_base(4);
t36 = (mrSges(3,2) * pkin(2));
t35 = (mrSges(4,2) * pkin(1));
t34 = (-Ifges(2,1) - Ifges(4,1));
t14 = (Ifges(2,4) + Ifges(4,4));
t20 = m(3) * pkin(2) ^ 2;
t30 = (Ifges(2,2) + Ifges(4,2) + t20);
t6 = (Ifges(2,3) + Ifges(4,3) + t20);
t25 = V_base(4) * V_base(5);
t24 = V_base(6) * V_base(4);
t16 = mrSges(2,2) + mrSges(4,2);
t17 = sin(qJ(1));
t18 = cos(qJ(1));
t7 = m(3) * pkin(2) + mrSges(2,1) + mrSges(4,1);
t22 = t16 * t18 + t17 * t7 + mrSges(1,2) + mrSges(3,2);
t19 = (m(4) * pkin(1) ^ 2);
t12 = V_base(4) ^ 2;
t11 = V_base(5) ^ 2;
t10 = V_base(6) ^ 2;
t8 = (mrSges(1,3) + mrSges(2,3) + mrSges(3,3) + mrSges(4,3));
t3 = (mrSges(3,1) * pkin(2) + mrSges(4,1) * pkin(1));
t2 = (-t35 + t36);
t1 = (t30 + t34);
t5 = -(t8 * V_base(4) * V_base(2)) + (-t22 * t26 + (t38 * t8)) * V_base(1) / 0.2e1 + t22 * V_base(4) * V_base(3) + ((Ifges(1,2) + Ifges(3,2) + t19 - t34) * t11) / 0.2e1 + ((Ifges(1,4) + Ifges(3,4) - t14) * V_base(4) + V_base(6) * (Ifges(1,6) + Ifges(3,6))) * t40 + ((Ifges(1,1) + Ifges(3,1) + t30) * t12) / 0.2e1 - ((mrSges(4,3) * pkin(1) - Ifges(1,5) - Ifges(3,5)) * t24) + (t10 * (Ifges(1,3) + Ifges(3,3) + t19 + t6)) / 0.2e1 + (-(t26 * V_base(2)) / 0.2e1 + V_base(3) * t40) * (-(m(4) * pkin(1)) + t17 * t16 - t18 * t7 - mrSges(1,1) - mrSges(3,1)) + ((-t1 * t25 - t12 * t14) * t17 + ((-t35 - t36) * t25) + (t3 * t11) + (V_base(6) * t3 + t39) * t42 + ((2 * t14 * t25) + ((t11 - t12) * t1) / 0.2e1) * t18) * t18 + ((t3 * t18 + t6) * t42 + (-t16 * V_base(1) + t7 * V_base(2) + t39) * t18 + (-(t13 * V_base(4)) - t16 * V_base(2) + (V_base(6) * t2) - (t4 * V_base(5)) - t7 * V_base(1)) * t17 + (t6 * qJD(1)) / 0.2e1) * qJD(1) + (V_base(1) ^ 2 + (V_base(2) ^ 2) + V_base(3) ^ 2) * (m(4) / 0.2e1 + m(2) / 0.2e1 + m(3) / 0.2e1 + m(1) / 0.2e1) + ((-t3 * V_base(4) - V_base(6) * t4) * t40 + (t12 * t36) - (t13 * t24) + (t10 * t2) + (t14 * t18 - t35) * t11) * t17;
T = t5;
