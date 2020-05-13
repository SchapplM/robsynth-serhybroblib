% Calculate vector of cutting forces with Newton-Euler
% fourbar1turnOL
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
% f_new [3x5]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = fourbar1turnOL_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnOL_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fourbar1turnOL_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_invdynf_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnOL_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fourbar1turnOL_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:49
% EndTime: 2020-04-12 19:40:50
% DurationCPUTime: 0.30s
% Computational Cost: add. (1213->102), mult. (2599->147), div. (0->0), fcn. (1599->8), ass. (0->57)
t44 = sin(qJ(4));
t67 = qJD(1) * t44;
t46 = sin(qJ(2));
t66 = qJD(1) * t46;
t48 = cos(qJ(4));
t65 = qJD(1) * t48;
t50 = cos(qJ(2));
t64 = qJD(1) * t50;
t63 = qJD(1) * qJD(2);
t62 = qJD(1) * qJD(4);
t28 = (-mrSges(5,1) * t48 + mrSges(5,2) * t44) * qJD(1);
t47 = sin(qJ(1));
t51 = cos(qJ(1));
t40 = -t51 * g(1) - t47 * g(2);
t52 = qJD(1) ^ 2;
t30 = -t52 * pkin(1) + t40;
t31 = t44 * qJDD(1) + t48 * t62;
t37 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t65;
t11 = m(5) * (-t48 * g(3) - t44 * t30) - t31 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t28 * t67 + qJD(4) * t37;
t33 = t48 * qJDD(1) - t44 * t62;
t35 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t67;
t12 = m(5) * (-t44 * g(3) + t48 * t30) + t33 * mrSges(5,3) - qJDD(4) * mrSges(5,2) + t28 * t65 - qJD(4) * t35;
t29 = (-mrSges(3,1) * t50 + mrSges(3,2) * t46) * qJD(1);
t32 = t46 * qJDD(1) + t50 * t63;
t38 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t64;
t45 = sin(qJ(3));
t49 = cos(qJ(3));
t58 = -t50 * g(3) - t46 * t40;
t24 = (t45 * t46 - t49 * t50) * qJD(1);
t59 = t46 * t63;
t34 = t50 * qJDD(1) - t59;
t15 = t24 * qJD(3) - t49 * t32 - t45 * t34;
t25 = (-t45 * t50 - t46 * t49) * qJD(1);
t17 = -t24 * mrSges(4,1) + t25 * mrSges(4,2);
t20 = (t46 * t50 * t52 + qJDD(2)) * pkin(2) + t58;
t60 = -t46 * g(3) + t50 * t40;
t21 = (-t50 ^ 2 * t52 - qJD(2) ^ 2) * pkin(2) + t60;
t43 = qJD(2) + qJD(3);
t22 = -t43 * mrSges(4,2) + t24 * mrSges(4,3);
t42 = qJDD(2) + qJDD(3);
t7 = m(4) * (-t49 * t20 + t45 * t21) - t15 * mrSges(4,3) + t42 * mrSges(4,1) - t25 * t17 + t43 * t22;
t14 = -t25 * qJD(3) + t45 * t32 - t49 * t34;
t23 = t43 * mrSges(4,1) - t25 * mrSges(4,3);
t8 = m(4) * (-t45 * t20 - t49 * t21) + t14 * mrSges(4,3) - t42 * mrSges(4,2) + t24 * t17 - t43 * t23;
t4 = m(3) * t58 + qJDD(2) * mrSges(3,1) - t32 * mrSges(3,3) + qJD(2) * t38 - t29 * t66 - t45 * t8 - t49 * t7;
t36 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t66;
t5 = m(3) * t60 - qJDD(2) * mrSges(3,2) + t34 * mrSges(3,3) - qJD(2) * t36 + t29 * t64 + t45 * t7 - t49 * t8;
t61 = t48 * t11 + t44 * t12 + t50 * t4 + t46 * t5;
t39 = t47 * g(1) - t51 * g(2);
t57 = t44 * t35 - t48 * t37;
t56 = t46 * t36 - t50 * t38;
t55 = m(5) * (-qJDD(1) * pkin(1) - t39) - t33 * mrSges(5,1) + t31 * mrSges(5,2);
t54 = t14 * mrSges(4,1) - t15 * mrSges(4,2) + t24 * t22 - t25 * t23 - m(4) * ((-t34 + t59) * pkin(2) - t39);
t53 = -t34 * mrSges(3,1) + t32 * mrSges(3,2) - t54;
t6 = qJDD(1) * mrSges(2,1) - t52 * mrSges(2,2) + (m(2) + m(3)) * t39 + (-t56 - t57) * qJD(1) - t53 - t55;
t1 = m(2) * t40 - t52 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t44 * t11 + t48 * t12 - t46 * t4 + t50 * t5;
t2 = [-m(1) * g(1) + t51 * t1 - t47 * t6, t1, t5, t8, t12, 0, 0; -m(1) * g(2) + t47 * t1 + t51 * t6, t6, t4, t7, t11, 0, 0; (-m(1) - m(2)) * g(3) + t61, -m(2) * g(3) + t61, -m(3) * t39 + t56 * qJD(1) + t53, -t54, t57 * qJD(1) + t55, 0, 0;];
f_new = t2;
