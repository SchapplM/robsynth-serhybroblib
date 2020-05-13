% Calculate vector of inverse dynamics joint torques for
% fivebar1OL
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fivebar1OL_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1OL_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fivebar1OL_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fivebar1OL_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1OL_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1OL_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1OL_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fivebar1OL_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fivebar1OL_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:13:07
% EndTime: 2020-04-27 06:13:08
% DurationCPUTime: 0.22s
% Computational Cost: add. (76->60), mult. (134->92), div. (0->0), fcn. (30->8), ass. (0->19)
t22 = qJD(2) * (qJD(1) + qJD(2) / 0.2e1);
t21 = qJD(4) * (qJD(3) + qJD(4) / 0.2e1);
t20 = qJD(1) ^ 2;
t19 = qJD(3) ^ 2;
t16 = cos(qJ(1));
t15 = cos(qJ(2));
t14 = cos(qJ(3));
t13 = cos(qJ(4));
t12 = sin(qJ(1));
t11 = sin(qJ(2));
t10 = sin(qJ(3));
t9 = sin(qJ(4));
t6 = qJDD(1) + qJDD(2) / 0.2e1;
t5 = qJDD(3) + qJDD(4) / 0.2e1;
t4 = mrSges(3,1) * g(2) - mrSges(3,2) * g(1);
t3 = -mrSges(3,1) * g(1) - mrSges(3,2) * g(2);
t2 = mrSges(5,1) * g(2) - mrSges(5,2) * g(1);
t1 = g(1) * mrSges(5,1) + mrSges(5,2) * g(2);
t7 = [(-g(2) * mrSges(2,1) + mrSges(2,2) * g(1) + t3 * t11 + t4 * t15) * t16 + (g(1) * mrSges(2,1) + mrSges(2,2) * g(2) - t4 * t11 + t3 * t15) * t12 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1) + Ifges(3,3) * qJDD(2) + (0.2e1 * (-t6 * mrSges(3,1) + mrSges(3,2) * t22) * t15 + 0.2e1 * (mrSges(3,1) * t22 + mrSges(3,2) * t6) * t11 + (pkin(2) * qJDD(1) + g(1) * t12 - g(2) * t16) * m(3)) * pkin(2); (t4 * t16 + t3 * t12 - pkin(2) * (mrSges(3,1) * qJDD(1) + mrSges(3,2) * t20)) * t15 + Ifges(3,3) * (qJDD(1) + qJDD(2)) + (t3 * t16 - t4 * t12 - pkin(2) * (mrSges(3,1) * t20 - mrSges(3,2) * qJDD(1))) * t11; (-g(2) * mrSges(4,1) + mrSges(4,2) * g(1) + t1 * t9 - t2 * t13) * t14 + (g(1) * mrSges(4,1) + mrSges(4,2) * g(2) + t1 * t13 + t2 * t9) * t10 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + Ifges(5,3) * qJDD(4) + (-0.2e1 * (-t5 * mrSges(5,1) + mrSges(5,2) * t21) * t13 - 0.2e1 * (mrSges(5,1) * t21 + mrSges(5,2) * t5) * t9 + (qJDD(3) * pkin(3) + g(1) * t10 - g(2) * t14) * m(5)) * pkin(3); (-t2 * t14 + t1 * t10 + (mrSges(5,1) * qJDD(3) + mrSges(5,2) * t19) * pkin(3)) * t13 + Ifges(5,3) * (qJDD(3) + qJDD(4)) + (t1 * t14 + t2 * t10 + (mrSges(5,1) * t19 - mrSges(5,2) * qJDD(3)) * pkin(3)) * t9; 0;];
tau = t7;
