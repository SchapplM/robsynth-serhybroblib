% Calculate vector of inverse dynamics joint torques for
% fourbarprisDE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% qJDD [1x1]
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
% tau [1x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:45
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbarprisDE2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(1,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE2_invdynJ_fixb_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisDE2_invdynJ_fixb_slag_vp2: qJD has to be [1x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [1 1]), ...
  'fourbarprisDE2_invdynJ_fixb_slag_vp2: qJDD has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisDE2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE2_invdynJ_fixb_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisDE2_invdynJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbarprisDE2_invdynJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbarprisDE2_invdynJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:45:03
% EndTime: 2020-05-07 09:45:12
% DurationCPUTime: 2.71s
% Computational Cost: add. (1999->93), mult. (1271->183), div. (171->17), fcn. (106->2), ass. (0->86)
t47 = qJ(1) + pkin(3);
t102 = -pkin(2) - t47;
t30 = pkin(1) - t102;
t101 = -pkin(2) + t47;
t32 = pkin(1) + t101;
t33 = pkin(1) + t102;
t110 = t32 * t33;
t31 = pkin(1) - t101;
t96 = t31 * t110;
t87 = t30 * t96;
t59 = sqrt(-t87);
t10 = 0.1e1 / t59;
t56 = (pkin(2) ^ 2);
t52 = (qJ(1) ^ 2);
t54 = (pkin(3) ^ 2);
t108 = (t54 + t52);
t118 = -2 * qJ(1);
t57 = (pkin(1) ^ 2);
t79 = pkin(3) * t118 - t108 + t57;
t13 = -t56 + t79;
t114 = t10 * t13;
t18 = 0.1e1 / t30;
t22 = 0.1e1 / t32;
t111 = t18 * t22;
t125 = 0.1e1 / t33;
t126 = 0.1e1 / t31;
t129 = t125 * t126;
t128 = t111 * t129;
t50 = qJD(1) ^ 2;
t141 = (t50 * t128);
t60 = pkin(2) * t56;
t140 = pkin(2) * t79 - t60;
t106 = qJ(1) * t54;
t51 = qJ(1) * t52;
t53 = pkin(3) * t54;
t139 = -4 * t106 - 2 * t51 - t53;
t43 = t47 ^ 2;
t136 = (Ifges(4,3) * t43);
t48 = (mrSges(2,2) - mrSges(3,1));
t49 = (mrSges(2,1) + mrSges(3,3));
t15 = (g(1) * t48 - t49 * g(2));
t127 = t13 ^ 2;
t133 = Ifges(2,3) * t127;
t124 = 0.1e1 / t31 ^ 2;
t123 = 0.1e1 / t33 ^ 2;
t122 = -0.2e1 * t10;
t120 = -5 * t52;
t66 = (t47 * t43);
t119 = -2 * t66;
t117 = (m(3) * g(2));
t45 = 0.1e1 / t47 ^ 2;
t116 = -t45 / 0.2e1;
t16 = t47 * qJD(1);
t44 = 0.1e1 / t47;
t68 = -t96 + (t110 + (t32 - t33) * t31) * t30;
t7 = t68 * qJD(1);
t94 = 0.1e1 / t87 * t114 / 0.2e1;
t86 = t7 * t94;
t99 = t50 * t114;
t90 = t45 * t99;
t93 = qJDD(1) * t114;
t2 = -t90 + (t93 + (t122 * t16 + t86) * qJD(1)) * t44;
t115 = qJ(1) * t2;
t113 = t13 * t44;
t112 = t15 * t57;
t109 = ((pkin(1) + t47) * (pkin(1) - t47));
t105 = qJD(1) * t44;
t104 = qJD(1) * t45;
t100 = t10 * t113;
t95 = t10 * t105;
t85 = t13 * t95;
t82 = qJ(1) * t85;
t77 = t127 * t141;
t46 = 1 / t66;
t76 = t46 * t77;
t19 = 0.1e1 / t30 ^ 2;
t23 = 0.1e1 / t32 ^ 2;
t75 = t19 * t124 * t23 * t123 * t7 * t104;
t74 = t16 * t128 * t104;
t71 = t57 / pkin(1) ^ 2 * t128 / 0.8e1;
t70 = t127 * t75;
t69 = t13 * t74;
t4 = -t104 * t114 + (t122 * t47 + t68 * t94) * t105;
t3 = -qJ(1) * t4 - t85;
t1 = -t44 * t99 - t115;
t5 = [m(3) * (0.2e1 * qJDD(1) + (t1 * t10 * t118 * t44 + 0.4e1 * t52 * t74) * t13 + (-t52 * t75 + 0.2e1 * (-qJ(1) * t45 + (t46 * t52)) * t141) * t127) / 0.2e1 + t10 * ((((mrSges(4,1) * t119) + (t140 * t49) + (pkin(3) * t60 + ((t120 - t57) * pkin(3) + t139) * pkin(2)) * m(3)) * g(1) + ((mrSges(4,2) * t119 + t140 * t48) * g(2))) * t59 + ((pkin(3) * t117 + t15) * t56 * t60) + (2 * (-t112 + (-t53 - 3 * t106 - (3 * t52 + t57) * pkin(3) - t51) * t117) * t60) + ((t112 + ((t120 + t57) * pkin(3) + t139) * t117 + (2 * qJ(1) * pkin(3) + t108) * t15) * pkin(2) * t109) + 0.2e1 * (t56 + t109) * (-(mrSges(4,1) * g(2)) + mrSges(4,2) * g(1)) * t66) / pkin(1) / pkin(2) * t116 + m(3) * t3 * t82 + (0.4e1 * t47 * t141 + (-0.32e2 * qJDD(1) * t43 - 0.64e2 * t47 * t50) * t56 / (pkin(2) ^ 2) * t71) * Ifges(4,3) + (-t4 * t85 + 0.2e1 * t69 - t70 / 0.2e1 + t76 + t2 * t100) * Ifges(3,2) + (-qJD(1) * t4 + t105 * t86 - 0.2e1 * t16 * t95 + t44 * t93 + t2 - t90) * mrSges(3,1) + ((0.32e2 * qJD(1) * t13 * t16 - 0.8e1 * qJDD(1) * t127) * t45 * t71 + ((t127 * t46) - 0.2e1 * t113) * t141) * Ifges(2,3) + ((-t1 + t115) * t100 + ((2 * t76) - t70 + 0.4e1 * t69) * qJ(1) + t3 * t85 - t4 * t82 - t45 * t77) * mrSges(3,3) + (t123 * t126 + t124 * t125) * (t116 * t133 - (2 * t136)) * t50 * t111 + (t23 * t18 + t22 * t19) * (t45 * t133 / 0.2e1 + (2 * t136)) * t50 * t129;];
tau = t5;
