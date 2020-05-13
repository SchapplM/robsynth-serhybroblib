% Calculate time derivative of joint inertia matrix for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
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
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m2DE_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_inertiaDJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2DE_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2DE_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'palh2m2DE_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:27
% EndTime: 2020-05-03 01:06:28
% DurationCPUTime: 0.22s
% Computational Cost: add. (110->49), mult. (233->76), div. (0->0), fcn. (115->6), ass. (0->25)
t11 = sin(qJ(2));
t30 = pkin(4) * t11;
t10 = sin(qJ(3));
t12 = cos(qJ(4));
t9 = sin(qJ(4));
t19 = mrSges(7,1) * t12 - mrSges(7,2) * t9;
t2 = t19 * qJD(4);
t28 = t10 * t2;
t27 = t11 * t2;
t13 = cos(qJ(3));
t26 = mrSges(5,2) * t13;
t25 = t10 * mrSges(5,2);
t14 = cos(qJ(2));
t20 = -t14 * pkin(4) - pkin(1);
t5 = -pkin(2) + t20;
t15 = -t13 * pkin(5) + t5;
t1 = pkin(3) - t15;
t23 = qJD(4) * t1;
t22 = qJD(2) * t14;
t21 = qJD(3) * t13;
t4 = mrSges(7,1) * t9 + mrSges(7,2) * t12;
t17 = -mrSges(6,3) + t4;
t6 = (m(6) + m(7)) * pkin(5) + mrSges(5,1);
t3 = -pkin(5) * qJD(3) * t10 - qJD(2) * t30;
t7 = [-0.2e1 * t4 * t23 + 0.2e1 * (-m(6) * t15 + mrSges(6,1) + m(7) * (t12 ^ 2 + t9 ^ 2) * t1 + t19) * t3 + 0.2e1 * (t5 * (mrSges(5,1) * t10 + t26) + (-t10 ^ 2 + t13 ^ 2) * Ifges(5,4) + (Ifges(5,1) - Ifges(5,2)) * t13 * t10) * qJD(3) + 0.2e1 * ((-pkin(1) * mrSges(3,1) - Ifges(3,4) * t11) * t11 + (-mrSges(3,2) * pkin(1) + (Ifges(3,1) - Ifges(3,2)) * t11 + Ifges(3,4) * t14) * t14 + (m(4) * t20 + m(5) * t5 - t13 * mrSges(5,1) - mrSges(4,1) + t25) * t30) * qJD(2); (Ifges(3,5) * t14 - Ifges(3,6) * t11) * qJD(2) + ((-mrSges(4,3) - mrSges(5,3) + t17) * t22 + t27) * pkin(4); 0; (Ifges(5,5) * t13 - Ifges(5,6) * t10) * qJD(3) + (t17 * t21 + t28) * pkin(5); pkin(4) * (qJD(2) - qJD(3)) * ((t6 * t10 + t26) * t14 - (t6 * t13 - t25) * t11); 0; (-t12 * t23 - t3 * t9) * mrSges(7,2) + (t3 * t12 - t9 * t23) * mrSges(7,1); (t4 * t22 + t27) * pkin(4); (t4 * t21 + t28) * pkin(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t7(1), t7(2), t7(4), t7(7); t7(2), t7(3), t7(5), t7(8); t7(4), t7(5), t7(6), t7(9); t7(7), t7(8), t7(9), t7(10);];
Mq = res;
