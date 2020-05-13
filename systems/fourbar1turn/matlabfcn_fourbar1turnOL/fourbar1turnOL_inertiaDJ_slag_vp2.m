% Calculate time derivative of joint inertia matrix for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1turnOL_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_inertiaDJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnOL_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnOL_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:49
% EndTime: 2020-04-12 19:40:50
% DurationCPUTime: 0.16s
% Computational Cost: add. (93->30), mult. (243->64), div. (0->0), fcn. (194->6), ass. (0->15)
t18 = qJD(2) + qJD(3);
t10 = sin(qJ(2));
t12 = cos(qJ(3));
t13 = cos(qJ(2));
t9 = sin(qJ(3));
t6 = t9 * t10 - t12 * t13;
t3 = t18 * t6;
t14 = t12 * t10 + t9 * t13;
t4 = t18 * t14;
t17 = Ifges(4,5) * t3 + Ifges(4,6) * t4;
t15 = 0.2e1 * t13;
t11 = cos(qJ(4));
t8 = sin(qJ(4));
t5 = (mrSges(4,1) * t9 + mrSges(4,2) * t12) * qJD(3) * pkin(2);
t1 = [-0.2e1 * t3 * t14 * Ifges(4,1) + 0.2e1 * t4 * Ifges(4,2) * t6 + (-pkin(2) * (-t4 * mrSges(4,1) + t3 * mrSges(4,2)) + Ifges(3,4) * qJD(2) * t13) * t15 + 0.2e1 * (-t14 * t4 + t3 * t6) * Ifges(4,4) + 0.2e1 * (-pkin(1) * (mrSges(5,1) * t8 + mrSges(5,2) * t11) + (t11 ^ 2 - t8 ^ 2) * Ifges(5,4) + (Ifges(5,1) - Ifges(5,2)) * t11 * t8) * qJD(4) + (0.2e1 * pkin(2) * (-t6 * mrSges(4,1) - mrSges(4,2) * t14) - 0.2e1 * Ifges(3,4) * t10 + (-m(4) * pkin(2) ^ 2 + Ifges(3,1) - Ifges(3,2)) * t15) * qJD(2) * t10; (Ifges(3,5) * t13 - Ifges(3,6) * t10) * qJD(2) + (t12 * t3 - t4 * t9 + (-t12 * t6 + t14 * t9) * qJD(3)) * pkin(2) * mrSges(4,3) + t17; 0.2e1 * t5; t17; t5; 0; (Ifges(5,5) * t11 - Ifges(5,6) * t8) * qJD(4); 0; 0; 0; 0; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
