% Calculate vector of inverse dynamics joint torques with ic for
% fivebar1IC
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
% tau [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:19
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fivebar1IC_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1IC_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fivebar1IC_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fivebar1IC_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1IC_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1IC_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1IC_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fivebar1IC_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fivebar1IC_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:19:10
% EndTime: 2020-04-27 06:19:10
% DurationCPUTime: 0.21s
% Computational Cost: add. (134->77), mult. (218->110), div. (8->4), fcn. (60->15), ass. (0->32)
t33 = (qJ(4) + qJ(3));
t34 = (qJ(1) + qJ(2));
t41 = cos((2 * t33)) - cos((2 * t34));
t17 = sin(qJ(2));
t40 = pkin(2) * t17;
t15 = sin(qJ(4));
t39 = pkin(3) * t15;
t16 = sin(qJ(3));
t19 = cos(qJ(4));
t20 = cos(qJ(3));
t25 = qJD(3) ^ 2;
t7 = mrSges(5,1) * g(1) + mrSges(5,2) * g(2);
t8 = mrSges(5,1) * g(2) - mrSges(5,2) * g(1);
t38 = ((-t8 * t20 + t7 * t16 + (mrSges(5,1) * qJDD(3) + mrSges(5,2) * t25) * pkin(3)) * t19 + Ifges(5,3) * (qJDD(3) + qJDD(4)) + (t7 * t20 + t8 * t16 + (mrSges(5,1) * t25 - mrSges(5,2) * qJDD(3)) * pkin(3)) * t15) / pkin(4);
t10 = mrSges(3,1) * g(2) - mrSges(3,2) * g(1);
t18 = sin(qJ(1));
t21 = cos(qJ(2));
t22 = cos(qJ(1));
t26 = qJD(1) ^ 2;
t9 = -mrSges(3,1) * g(1) - mrSges(3,2) * g(2);
t37 = ((t10 * t22 + t9 * t18 - pkin(2) * (mrSges(3,1) * qJDD(1) + mrSges(3,2) * t26)) * t21 + Ifges(3,3) * (qJDD(1) + qJDD(2)) + (t9 * t22 - t10 * t18 - pkin(2) * (mrSges(3,1) * t26 - mrSges(3,2) * qJDD(1))) * t17) / pkin(5);
t36 = qJD(2) * (qJD(1) + qJD(2) / 0.2e1);
t35 = qJD(4) * (qJD(3) + qJD(4) / 0.2e1);
t32 = pkin(2) * m(3) + mrSges(2,1);
t31 = pkin(3) * m(5) + mrSges(4,1);
t28 = 0.2e1 * qJ(1);
t27 = 0.2e1 * qJ(3);
t12 = qJDD(1) + qJDD(2) / 0.2e1;
t11 = qJDD(3) + qJDD(4) / 0.2e1;
t4 = -0.1e1 / sin((t33 - t34));
t3 = 0.1e1 / t41;
t1 = [(mrSges(2,2) * g(1) - g(2) * t32 + t10 * t21 + t9 * t17) * t22 + (mrSges(2,2) * g(2) + g(1) * t32 - t10 * t17 + t9 * t21) * t18 + 0.2e1 * (-t12 * mrSges(3,1) + mrSges(3,2) * t36) * pkin(2) * t21 + 0.2e1 * (mrSges(3,1) * t36 + mrSges(3,2) * t12) * t40 + Ifges(3,3) * qJDD(2) - (t41 * pkin(5) + (-cos(qJ(2) + 0.2e1 * qJ(4) + t27) + cos(qJ(2) + t28)) * pkin(2)) * t3 * t37 + t4 * t38 * t40 + (pkin(2) ^ 2 * m(3) + Ifges(2,3) + Ifges(3,3)) * qJDD(1); t4 * t37 * t39 + (mrSges(4,2) * g(1) - g(2) * t31 + t7 * t15 - t8 * t19) * t20 + (mrSges(4,2) * g(2) + g(1) * t31 + t8 * t15 + t7 * t19) * t16 - 0.2e1 * pkin(3) * (-t11 * mrSges(5,1) + mrSges(5,2) * t35) * t19 - 0.2e1 * (mrSges(5,1) * t35 + mrSges(5,2) * t11) * t39 + Ifges(5,3) * qJDD(4) - (t41 * pkin(4) + (cos(qJ(4) + t27) - cos(qJ(4) + t28 + 0.2e1 * qJ(2))) * pkin(3)) * t3 * t38 + (m(5) * pkin(3) ^ 2 + Ifges(4,3) + Ifges(5,3)) * qJDD(3);];
tau = t1(:);
