% Calculate vector of cutting torques with Newton-Euler for
% fourbarprisOL
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% m [3x4]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:52
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = fourbarprisOL_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisOL_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbarprisOL_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbarprisOL_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisOL_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_invdynm_fixb_snew_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisOL_invdynm_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbarprisOL_invdynm_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbarprisOL_invdynm_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:52:06
% EndTime: 2020-05-07 09:52:06
% DurationCPUTime: 0.12s
% Computational Cost: add. (194->66), mult. (299->77), div. (0->0), fcn. (90->4), ass. (0->26)
t44 = sin(qJ(1));
t46 = cos(qJ(1));
t35 = t46 * g(1) + t44 * g(2);
t49 = qJD(1) ^ 2;
t29 = -t49 * qJ(2) + qJDD(2) + t35;
t58 = mrSges(3,2) * t29;
t57 = -mrSges(2,1) - mrSges(3,3);
t56 = Ifges(3,4) + Ifges(2,6);
t55 = Ifges(2,5) - Ifges(3,6);
t33 = -t44 * g(1) + t46 * g(2);
t27 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t33;
t54 = mrSges(3,1) * g(3) - mrSges(3,2) * t27;
t53 = -m(3) * t27 + t49 * mrSges(3,1) + qJDD(1) * mrSges(3,3);
t43 = sin(qJ(3));
t45 = cos(qJ(3));
t32 = -t43 * g(1) + t45 * g(2);
t34 = t45 * g(1) + t43 * g(2);
t52 = mrSges(4,1) * t32 - mrSges(4,2) * t34 + Ifges(4,3) * qJDD(3);
t51 = -mrSges(3,1) * t29 + mrSges(3,3) * t27 - Ifges(3,2) * qJDD(1);
t50 = mrSges(2,1) * t33 - mrSges(2,2) * t35 + Ifges(2,3) * qJDD(1) + qJ(2) * t53 - t51;
t48 = qJD(3) ^ 2;
t26 = -mrSges(4,2) * g(3) - mrSges(4,3) * t32 + Ifges(4,5) * qJDD(3) - t48 * Ifges(4,6);
t25 = mrSges(4,1) * g(3) + mrSges(4,3) * t34 + t48 * Ifges(4,5) + Ifges(4,6) * qJDD(3);
t23 = -mrSges(2,2) * g(3) - mrSges(2,3) * t33 + t55 * qJDD(1) - t56 * t49 + t54;
t22 = -t58 + mrSges(2,3) * t35 + t55 * t49 + t56 * qJDD(1) + (m(3) * qJ(2) - t57) * g(3);
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t44 * t22 - t46 * t23 + t43 * t25 - t45 * t26, t23, -mrSges(3,3) * g(3) - Ifges(3,4) * qJDD(1) + t49 * Ifges(3,6) + t58, t26, 0, 0; -mrSges(1,3) * g(1) - t46 * t22 - t44 * t23 - t45 * t25 - t43 * t26 + (mrSges(1,1) - pkin(1) * (-m(2) - m(3))) * g(3), t22, t51, t25, 0, 0; pkin(1) * (-t44 * (m(2) * t35 + m(3) * t29 + t57 * t49 + (-mrSges(2,2) + mrSges(3,1)) * qJDD(1)) - t46 * (m(2) * t33 + qJDD(1) * mrSges(2,1) - t49 * mrSges(2,2) + t53)) - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t50 + t52, t50, -t49 * Ifges(3,4) - Ifges(3,6) * qJDD(1) + t54, t52, 0, 0;];
m_new = t1;
