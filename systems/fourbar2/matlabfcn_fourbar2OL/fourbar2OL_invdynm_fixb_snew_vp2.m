% Calculate vector of cutting torques with Newton-Euler for
% fourbar2OL
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2]';
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
% Datum: 2020-04-24 20:32
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = fourbar2OL_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar2OL_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbar2OL_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbar2OL_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar2OL_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2OL_invdynm_fixb_snew_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar2OL_invdynm_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar2OL_invdynm_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbar2OL_invdynm_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:32:25
% EndTime: 2020-04-24 20:32:25
% DurationCPUTime: 0.33s
% Computational Cost: add. (233->94), mult. (353->119), div. (0->0), fcn. (196->6), ass. (0->48)
t69 = sin(qJ(2));
t62 = t69 * mrSges(3,2);
t72 = cos(qJ(2));
t95 = t72 * mrSges(3,1);
t108 = -t95 + t62;
t75 = m(3) * pkin(2);
t100 = -mrSges(2,1) - t75 - t108;
t70 = sin(qJ(1));
t107 = t100 * t70;
t56 = t69 * mrSges(3,1) + t72 * mrSges(3,2);
t103 = -mrSges(2,2) + t56;
t106 = t103 * t70;
t105 = (2 * qJD(1) + qJD(2)) * qJD(2);
t73 = cos(qJ(1));
t41 = t100 * t73;
t68 = sin(qJ(3));
t71 = cos(qJ(3));
t88 = mrSges(4,1) * t71 - t68 * mrSges(4,2);
t102 = -t41 + mrSges(1,1) + t88;
t44 = t103 * t73;
t53 = t68 * mrSges(4,1) + mrSges(4,2) * t71;
t101 = t44 - mrSges(1,2) - t53;
t67 = mrSges(3,3) + mrSges(2,3);
t93 = Ifges(4,3) * qJDD(3);
t55 = t72 * Ifges(3,5) - t69 * Ifges(3,6);
t87 = (Ifges(2,3) + Ifges(3,3) + (0.2e1 * t62 - 0.2e1 * t95 + t75) * pkin(2)) * qJDD(1) + (t108 * pkin(2) + Ifges(3,3)) * qJDD(2) + t105 * pkin(2) * t56;
t83 = -t73 * g(1) - t70 * g(2);
t82 = -t70 * g(1) + t73 * g(2);
t57 = Ifges(3,5) * t69 + Ifges(3,6) * t72;
t45 = -mrSges(3,3) * pkin(2) + Ifges(2,5) - t55;
t51 = -Ifges(2,6) + t57;
t81 = t70 * t45 - t51 * t73;
t39 = -t55 * t73 + t70 * t57;
t80 = t70 * t55 + t57 * t73;
t78 = qJD(1) ^ 2;
t76 = qJD(3) ^ 2;
t66 = qJD(1) + qJD(2);
t65 = qJDD(1) + qJDD(2);
t64 = t66 ^ 2;
t60 = mrSges(1,3) + mrSges(4,3) + t67;
t54 = t71 * Ifges(4,5) - t68 * Ifges(4,6);
t52 = -t68 * Ifges(4,5) - t71 * Ifges(4,6);
t49 = -t78 * pkin(2) + t83;
t48 = qJDD(1) * pkin(2) - t82;
t38 = t45 * t73 + t70 * t51;
t37 = -t69 * t48 - t72 * t49;
t36 = -t72 * t48 + t69 * t49;
t1 = [t60 * g(2) - (-t101 - t107) * g(3) + t38 * qJDD(1) - t81 * t78 + t39 * qJDD(2) + t54 * qJDD(3) + t52 * t76 + t105 * t80, g(3) * t103 + t45 * qJDD(1) - t55 * qJDD(2) + t105 * t57 + t51 * t78 + t82 * t67, -mrSges(3,2) * g(3) - mrSges(3,3) * t36 + Ifges(3,5) * t65 - t64 * Ifges(3,6), -mrSges(4,2) * g(3) + Ifges(4,5) * qJDD(3) - t76 * Ifges(4,6) + (-t68 * g(1) + g(2) * t71) * mrSges(4,3), 0, 0; -t60 * g(1) - (-m(4) * pkin(1) - t102 - t106) * g(3) + t81 * qJDD(1) + t38 * t78 - t80 * qJDD(2) - t52 * qJDD(3) + t54 * t76 + t105 * t39, -g(3) * t100 - t51 * qJDD(1) - t57 * qJDD(2) - t105 * t55 + t45 * t78 + t83 * t67, mrSges(3,1) * g(3) + mrSges(3,3) * t37 + t64 * Ifges(3,5) + Ifges(3,6) * t65, mrSges(4,1) * g(3) + t76 * Ifges(4,5) + Ifges(4,6) * qJDD(3) + (-t71 * g(1) - t68 * g(2)) * mrSges(4,3), 0, 0; -t102 * g(2) - t101 * g(1) + t93 + (-g(1) * t100 - g(2) * t103) * t70 + (-m(4) * g(2) + t88 * qJDD(3) - t53 * t76) * pkin(1) + t87, -(t44 + t107) * g(1) - (-t41 + t106) * g(2) + t87, mrSges(3,1) * t36 - mrSges(3,2) * t37 + Ifges(3,3) * t65, t53 * g(1) - t88 * g(2) + t93, 0, 0;];
m_new = t1;
