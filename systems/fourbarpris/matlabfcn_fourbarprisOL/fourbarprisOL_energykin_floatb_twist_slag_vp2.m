% Calculate kinetic energy for
% fourbarprisOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
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
% Datum: 2020-05-07 09:52
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbarprisOL_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisOL_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbarprisOL_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fourbarprisOL_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_energykin_floatb_twist_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisOL_energykin_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbarprisOL_energykin_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbarprisOL_energykin_floatb_twist_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:52:02
% EndTime: 2020-05-07 09:52:02
% DurationCPUTime: 0.28s
% Computational Cost: add. (159->74), mult. (232->102), div. (0->0), fcn. (100->4), ass. (0->21)
t21 = cos(qJ(1));
t13 = -V_base(5) * pkin(1) + V_base(3);
t12 = V_base(6) * pkin(1) + V_base(2);
t18 = sin(qJ(1));
t4 = -t21 * t12 + t18 * V_base(1);
t5 = -t18 * t12 - t21 * V_base(1);
t20 = V_base(3) ^ 2;
t19 = cos(qJ(3));
t17 = sin(qJ(3));
t16 = V_base(6) + qJD(1);
t15 = V_base(6) + qJD(3);
t11 = -t17 * V_base(2) - t19 * V_base(1);
t10 = t17 * V_base(1) - t19 * V_base(2);
t9 = -t18 * V_base(5) - t21 * V_base(4);
t8 = -t17 * V_base(5) - t19 * V_base(4);
t7 = t18 * V_base(4) - t21 * V_base(5);
t6 = t17 * V_base(4) - t19 * V_base(5);
t3 = qJD(2) + t5;
t2 = t7 * qJ(2) - t13;
t1 = -t16 * qJ(2) - t4;
t14 = m(4) * (t10 ^ 2 + t11 ^ 2 + t20) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t20) / 0.2e1 + m(2) * (t13 ^ 2 + t4 ^ 2 + t5 ^ 2) / 0.2e1 + m(3) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (V_base(3) * mrSges(4,2) - t10 * mrSges(4,3) + Ifges(4,1) * t8 / 0.2e1) * t8 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-V_base(3) * mrSges(4,1) + t11 * mrSges(4,3) + Ifges(4,4) * t8 + Ifges(4,2) * t6 / 0.2e1) * t6 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t10 * mrSges(4,1) - t11 * mrSges(4,2) + Ifges(4,5) * t8 + Ifges(4,6) * t6 + Ifges(4,3) * t15 / 0.2e1) * t15 + (t2 * mrSges(3,1) + t13 * mrSges(2,2) - t1 * mrSges(3,2) - t4 * mrSges(2,3) + (Ifges(3,3) / 0.2e1 + Ifges(2,1) / 0.2e1) * t9) * t9 + (-t13 * mrSges(2,1) - t3 * mrSges(3,2) + t5 * mrSges(2,3) + t2 * mrSges(3,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,1) / 0.2e1) * t7 + (Ifges(2,4) - Ifges(3,5)) * t9) * t7 + (t4 * mrSges(2,1) + t3 * mrSges(3,1) - t5 * mrSges(2,2) - t1 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(2,3) / 0.2e1) * t16 + (Ifges(2,5) - Ifges(3,6)) * t9 + (Ifges(3,4) + Ifges(2,6)) * t7) * t16;
T = t14;
