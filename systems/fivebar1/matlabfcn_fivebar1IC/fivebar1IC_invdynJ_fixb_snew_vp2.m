% Calculate vector of inverse dynamics joint torques with newton euler and ic for
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

function tau = fivebar1IC_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1IC_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fivebar1IC_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fivebar1IC_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1IC_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1IC_invdynJ_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1IC_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fivebar1IC_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fivebar1IC_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:19:11
% EndTime: 2020-04-27 06:19:11
% DurationCPUTime: 0.16s
% Computational Cost: add. (134->69), mult. (202->90), div. (8->4), fcn. (104->15), ass. (0->36)
t45 = pkin(2) * m(3);
t44 = pkin(3) * m(5);
t19 = sin(qJ(4));
t23 = cos(qJ(4));
t20 = sin(qJ(3));
t24 = cos(qJ(3));
t6 = qJDD(3) * pkin(3) + t20 * g(1) - t24 * g(2);
t8 = -qJD(3) ^ 2 * pkin(3) - t24 * g(1) - t20 * g(2);
t41 = (Ifges(5,3) * (qJDD(3) + qJDD(4)) + mrSges(5,1) * (-t19 * t8 + t23 * t6) - mrSges(5,2) * (t19 * t6 + t23 * t8)) / pkin(4);
t21 = sin(qJ(2));
t25 = cos(qJ(2));
t22 = sin(qJ(1));
t26 = cos(qJ(1));
t7 = qJDD(1) * pkin(2) + t22 * g(1) - t26 * g(2);
t9 = -qJD(1) ^ 2 * pkin(2) - t26 * g(1) - t22 * g(2);
t40 = (Ifges(3,3) * (qJDD(1) + qJDD(2)) + mrSges(3,1) * (t21 * t9 - t25 * t7) - mrSges(3,2) * (-t21 * t7 - t25 * t9)) / pkin(5);
t39 = t19 * mrSges(5,2);
t38 = t21 * mrSges(3,2);
t37 = t23 * mrSges(5,1);
t36 = t25 * mrSges(3,1);
t33 = qJ(4) + qJ(3);
t34 = qJ(1) + qJ(2);
t35 = cos(0.2e1 * t33) - cos(0.2e1 * t34);
t13 = t21 * mrSges(3,1) + t25 * mrSges(3,2);
t32 = -t36 + t38;
t31 = t37 - t39;
t12 = t19 * mrSges(5,1) + t23 * mrSges(5,2);
t28 = 0.2e1 * qJ(1);
t27 = 0.2e1 * qJ(3);
t14 = -0.1e1 / sin(t33 - t34);
t11 = -mrSges(2,2) + t13;
t10 = mrSges(4,2) + t12;
t5 = mrSges(2,1) + t32 + t45;
t4 = mrSges(4,1) + t31 + t44;
t3 = 0.1e1 / t35;
t1 = [-(t11 * t26 - t22 * t5) * g(1) - (t22 * t11 + t5 * t26) * g(2) + (Ifges(2,3) + Ifges(3,3) + (-0.2e1 * t36 + 0.2e1 * t38 + t45) * pkin(2)) * qJDD(1) + (t32 * pkin(2) + Ifges(3,3)) * qJDD(2) - (t35 * pkin(5) + (-cos(qJ(2) + 0.2e1 * qJ(4) + t27) + cos(qJ(2) + t28)) * pkin(2)) * t3 * t40 + pkin(2) * t21 * t14 * t41 + (0.2e1 * qJD(1) + qJD(2)) * pkin(2) * t13 * qJD(2); pkin(3) * t19 * t14 * t40 - (-t10 * t24 - t20 * t4) * g(1) - (-t10 * t20 + t4 * t24) * g(2) + (Ifges(4,3) + Ifges(5,3) + (0.2e1 * t37 - 0.2e1 * t39 + t44) * pkin(3)) * qJDD(3) + (t31 * pkin(3) + Ifges(5,3)) * qJDD(4) - (t35 * pkin(4) + (cos(qJ(4) + t27) - cos(qJ(4) + t28 + 0.2e1 * qJ(2))) * pkin(3)) * t3 * t41 + (-0.2e1 * qJD(3) - qJD(4)) * pkin(3) * t12 * qJD(4);];
tau = t1(:);
