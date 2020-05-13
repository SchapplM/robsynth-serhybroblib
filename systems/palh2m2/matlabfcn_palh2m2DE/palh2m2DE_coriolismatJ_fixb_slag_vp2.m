% Calculate matrix of centrifugal and coriolis load on the joints for
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = palh2m2DE_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_coriolismatJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2DE_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2DE_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'palh2m2DE_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:28
% EndTime: 2020-05-03 01:06:29
% DurationCPUTime: 0.19s
% Computational Cost: add. (303->59), mult. (611->82), div. (0->0), fcn. (347->6), ass. (0->40)
t26 = sin(qJ(4));
t46 = mrSges(7,2) * t26;
t29 = cos(qJ(4));
t47 = mrSges(7,1) * t29;
t15 = -t46 + t47;
t52 = t15 / 0.2e1 - t47 / 0.2e1 + t46 / 0.2e1;
t28 = sin(qJ(2));
t51 = pkin(4) * t28;
t27 = sin(qJ(3));
t50 = pkin(5) * t27;
t30 = cos(qJ(3));
t49 = t30 * pkin(5);
t31 = cos(qJ(2));
t48 = t31 * pkin(4);
t45 = t27 * mrSges(5,2);
t44 = pkin(4) * qJD(2);
t36 = -pkin(1) - t48;
t20 = -pkin(2) + t36;
t13 = t20 - t49;
t12 = pkin(3) - t13;
t32 = m(6) * t13 - mrSges(6,1) + m(7) * (-t26 ^ 2 - t29 ^ 2) * t12 - t15;
t1 = (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t31) * t31 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t28 + (Ifges(3,1) - Ifges(3,2)) * t31 + (t36 * m(4) + m(5) * t20 - t30 * mrSges(5,1) - mrSges(4,1) + t32 + t45) * pkin(4)) * t28;
t43 = t1 * qJD(1);
t2 = (t20 * mrSges(5,2) + Ifges(5,4) * t30) * t30 + (t20 * mrSges(5,1) - Ifges(5,4) * t27 + (Ifges(5,1) - Ifges(5,2)) * t30 + t32 * pkin(5)) * t27;
t42 = t2 * qJD(1);
t14 = mrSges(7,1) * t26 + mrSges(7,2) * t29;
t3 = t14 * t12;
t41 = t3 * qJD(1);
t40 = t3 * qJD(4);
t4 = t15 * t50;
t39 = t4 * qJD(1);
t5 = t15 * t51;
t38 = t5 * qJD(1);
t37 = qJD(4) * t15;
t21 = (m(6) + m(7)) * pkin(5) + mrSges(5,1);
t35 = (mrSges(5,2) * t30 + t21 * t27) * t31 - (t21 * t30 - t45) * t28;
t34 = -mrSges(6,3) + t14;
t7 = t52 * t51;
t6 = t52 * t50;
t8 = [t1 * qJD(2) + t2 * qJD(3) - t40, t43 + (Ifges(3,5) * t31 - Ifges(3,6) * t28 + (-mrSges(4,3) - mrSges(5,3) + t34) * t48) * qJD(2) + t7 * qJD(4), t42 + (Ifges(5,5) * t30 - Ifges(5,6) * t27 + t34 * t49) * qJD(3) + t6 * qJD(4), t7 * qJD(2) + t6 * qJD(3) - t40 - t41; t5 * qJD(4) - t43, 0, -pkin(4) * t35 * qJD(3), t37 * t51 + t38; t4 * qJD(4) - t42, t35 * t44, 0, t37 * t50 + t39; -t5 * qJD(2) - t4 * qJD(3) + t41, t31 * t14 * t44 - t38, t14 * qJD(3) * t49 - t39, 0;];
Cq = t8;
