% Calculate vector of inverse dynamics joint torques with newton euler and ic for
% fourbar1turnIC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 11:33
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbar1turnIC_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnIC_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnIC_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fourbar1turnIC_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnIC_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnIC_invdynJ_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnIC_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnIC_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnIC_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 11:33:36
% EndTime: 2020-05-07 11:33:36
% DurationCPUTime: 0.39s
% Computational Cost: add. (477->131), mult. (987->185), div. (4->3), fcn. (616->10), ass. (0->49)
t60 = -qJ(4) + qJ(2);
t59 = qJD(1) * qJD(2);
t58 = qJD(1) * qJD(4);
t50 = sin(qJ(1));
t54 = cos(qJ(1));
t39 = t50 * g(1) - t54 * g(2);
t40 = -t54 * g(1) - t50 * g(2);
t49 = sin(qJ(2));
t53 = cos(qJ(2));
t23 = -t49 * g(3) + t53 * t40;
t57 = t49 * t59;
t22 = -t53 * g(3) - t49 * t40;
t55 = qJD(1) ^ 2;
t52 = cos(qJ(3));
t51 = cos(qJ(4));
t48 = sin(qJ(3));
t47 = sin(qJ(4));
t46 = qJD(2) + qJD(3);
t45 = qJDD(2) + qJDD(3);
t38 = t53 * qJDD(1) - t57;
t37 = t51 * qJDD(1) - t47 * t58;
t36 = t49 * qJDD(1) + t53 * t59;
t35 = t47 * qJDD(1) + t51 * t58;
t34 = -t55 * pkin(1) + t40;
t33 = -qJDD(1) * pkin(1) - t39;
t31 = (-t53 * t48 - t49 * t52) * qJD(1);
t30 = (t49 * t48 - t53 * t52) * qJD(1);
t29 = Ifges(3,5) * qJD(2) + (t49 * Ifges(3,1) + t53 * Ifges(3,4)) * qJD(1);
t28 = Ifges(5,5) * qJD(4) + (t47 * Ifges(5,1) + t51 * Ifges(5,4)) * qJD(1);
t27 = Ifges(3,6) * qJD(2) + (t49 * Ifges(3,4) + t53 * Ifges(3,2)) * qJD(1);
t26 = Ifges(5,6) * qJD(4) + (t47 * Ifges(5,4) + t51 * Ifges(5,2)) * qJD(1);
t21 = -t47 * g(3) + t51 * t34;
t20 = -t51 * g(3) - t47 * t34;
t19 = t46 * mrSges(4,1) - t31 * mrSges(4,3);
t18 = -t46 * mrSges(4,2) + t30 * mrSges(4,3);
t17 = (-t53 ^ 2 * t55 - qJD(2) ^ 2) * pkin(2) + t23;
t16 = (t49 * t53 * t55 + qJDD(2)) * pkin(2) + t22;
t15 = (-t38 + t57) * pkin(2) - t39;
t14 = -t30 * mrSges(4,1) + t31 * mrSges(4,2);
t13 = Ifges(4,1) * t31 + Ifges(4,4) * t30 + Ifges(4,5) * t46;
t12 = Ifges(4,4) * t31 + Ifges(4,2) * t30 + Ifges(4,6) * t46;
t11 = Ifges(4,5) * t31 + Ifges(4,6) * t30 + Ifges(4,3) * t46;
t10 = t30 * qJD(3) - t52 * t36 - t48 * t38;
t9 = -t31 * qJD(3) + t48 * t36 - t52 * t38;
t5 = -t48 * t16 - t52 * t17;
t4 = -t52 * t16 + t48 * t17;
t2 = mrSges(4,2) * t15 - mrSges(4,3) * t4 + Ifges(4,1) * t10 + Ifges(4,4) * t9 + Ifges(4,5) * t45 + t30 * t11 - t46 * t12;
t1 = -mrSges(4,1) * t15 + mrSges(4,3) * t5 + Ifges(4,4) * t10 + Ifges(4,2) * t9 + Ifges(4,6) * t45 - t31 * t11 + t46 * t13;
t3 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t39 - mrSges(2,2) * t40 + t49 * (-mrSges(3,2) * t39 - mrSges(3,3) * t22 + Ifges(3,1) * t36 + Ifges(3,4) * t38 + Ifges(3,5) * qJDD(2) - qJD(2) * t27 + t48 * t1 - t52 * t2) + t53 * (Ifges(3,4) * t36 + Ifges(3,2) * t38 + Ifges(3,6) * qJDD(2) + qJD(2) * t29 + mrSges(3,1) * t39 + mrSges(3,3) * t23 - t48 * t2 - t52 * t1 - pkin(2) * (m(4) * t15 - t9 * mrSges(4,1) + t10 * mrSges(4,2) - t30 * t18 + t31 * t19)) + t47 * (mrSges(5,2) * t33 - mrSges(5,3) * t20 + Ifges(5,1) * t35 + Ifges(5,4) * t37 + Ifges(5,5) * qJDD(4) - qJD(4) * t26) + t51 * (-mrSges(5,1) * t33 + mrSges(5,3) * t21 + Ifges(5,4) * t35 + Ifges(5,2) * t37 + Ifges(5,6) * qJDD(4) + qJD(4) * t28) + (-m(5) * t33 + t37 * mrSges(5,1) - t35 * mrSges(5,2) + ((-mrSges(5,1) * t47 - mrSges(5,2) * t51) * qJD(4) + (t47 ^ 2 + t51 ^ 2) * qJD(1) * mrSges(5,3)) * qJD(1)) * pkin(1); Ifges(3,5) * t36 + Ifges(3,6) * t38 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t22 - mrSges(3,2) * t23 + (t49 * t27 - t53 * t29) * qJD(1) + (-t48 * (m(4) * t5 - t45 * mrSges(4,2) + t9 * mrSges(4,3) + t30 * t14 - t46 * t19) - t52 * (m(4) * t4 + t45 * mrSges(4,1) - t10 * mrSges(4,3) - t31 * t14 + t46 * t18) + (sin(t60) / pkin(3) * (mrSges(4,1) * t4 - mrSges(4,2) * t5 + Ifges(4,5) * t10 + Ifges(4,6) * t9 + Ifges(4,3) * t45 + t31 * t12 - t30 * t13) + t48 / pkin(4) * (mrSges(5,1) * t20 - mrSges(5,2) * t21 + Ifges(5,5) * t35 + Ifges(5,6) * t37 + Ifges(5,3) * qJDD(4) + (t47 * t26 - t51 * t28) * qJD(1))) / sin(qJ(3) + t60)) * pkin(2);];
tau = t3(:);
