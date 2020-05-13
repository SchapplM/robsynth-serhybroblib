% Calculate vector of cutting forces with Newton-Euler
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
% f_new [3x5]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = fivebar1OL_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1OL_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fivebar1OL_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fivebar1OL_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1OL_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1OL_invdynf_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1OL_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fivebar1OL_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fivebar1OL_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:13:11
% EndTime: 2020-04-27 06:13:11
% DurationCPUTime: 0.29s
% Computational Cost: add. (223->84), mult. (331->109), div. (0->0), fcn. (192->8), ass. (0->46)
t59 = (2 * qJD(3) + qJD(4)) * qJD(4);
t60 = (2 * qJD(1) + qJD(2)) * qJD(2);
t35 = m(2) + m(3);
t34 = m(4) + m(5);
t57 = sin(qJ(4));
t28 = sin(qJ(2));
t32 = cos(qJ(2));
t14 = t28 * mrSges(3,1) + t32 * mrSges(3,2);
t11 = -mrSges(2,2) + t14;
t29 = sin(qJ(1));
t33 = cos(qJ(1));
t15 = mrSges(3,1) * t32 - t28 * mrSges(3,2);
t5 = -pkin(2) * m(3) - mrSges(2,1) + t15;
t54 = t11 * t33 + t29 * t5;
t30 = cos(qJ(4));
t13 = mrSges(5,1) * t30 - t57 * mrSges(5,2);
t27 = sin(qJ(3));
t31 = cos(qJ(3));
t49 = -g(1) * t31 - g(2) * t27;
t48 = g(1) * t27 - g(2) * t31;
t47 = -g(1) * t33 - g(2) * t29;
t46 = g(1) * t29 - g(2) * t33;
t12 = t57 * mrSges(5,1) + t30 * mrSges(5,2);
t10 = mrSges(4,2) + t12;
t4 = pkin(3) * m(5) + mrSges(4,1) + t13;
t45 = -t10 * t31 - t27 * t4;
t44 = t10 * t27 - t31 * t4;
t43 = t11 * t29 - t33 * t5;
t42 = -t12 * t31 - t13 * t27;
t41 = t12 * t27 - t13 * t31;
t2 = t14 * t33 + t15 * t29;
t1 = t14 * t29 - t15 * t33;
t40 = qJD(1) ^ 2;
t38 = qJD(3) ^ 2;
t26 = qJD(1) + qJD(2);
t25 = qJD(3) + qJD(4);
t24 = qJDD(1) + qJDD(2);
t23 = qJDD(3) + qJDD(4);
t22 = t26 ^ 2;
t21 = t25 ^ 2;
t16 = m(1) + t34 + t35;
t9 = -pkin(2) * t40 + t47;
t8 = -pkin(3) * t38 + t49;
t7 = qJDD(1) * pkin(2) + t46;
t6 = qJDD(3) * pkin(3) + t48;
t3 = [-t16 * g(1) + t54 * qJDD(1) + t2 * qJDD(2) + t45 * qJDD(3) + t42 * qJDD(4) - t1 * t60 + t44 * t38 - t43 * t40 + t59 * t41, qJDD(1) * t11 + qJDD(2) * t14 + t15 * t60 + t47 * t35 + t40 * t5, m(3) * (-t28 * t7 - t32 * t9) - t24 * mrSges(3,2) - t22 * mrSges(3,1), -qJDD(3) * t10 - qJDD(4) * t12 - t13 * t59 + t49 * t34 - t38 * t4, m(5) * (t30 * t8 + t57 * t6) - t23 * mrSges(5,2) - t21 * mrSges(5,1), 0, 0; -t16 * g(2) + t43 * qJDD(1) + t1 * qJDD(2) - t44 * qJDD(3) - t41 * qJDD(4) + t2 * t60 + t45 * t38 + t54 * t40 + t59 * t42, -qJDD(1) * t5 - qJDD(2) * t15 + t11 * t40 + t14 * t60 + t46 * t35, m(3) * (t28 * t9 - t32 * t7) + t24 * mrSges(3,1) - t22 * mrSges(3,2), qJDD(3) * t4 + qJDD(4) * t13 - t10 * t38 - t12 * t59 + t48 * t34, m(5) * (t30 * t6 - t57 * t8) + t23 * mrSges(5,1) - t21 * mrSges(5,2), 0, 0; -t16 * g(3), -t35 * g(3), -m(3) * g(3), -t34 * g(3), -m(5) * g(3), 0, 0;];
f_new = t3;
