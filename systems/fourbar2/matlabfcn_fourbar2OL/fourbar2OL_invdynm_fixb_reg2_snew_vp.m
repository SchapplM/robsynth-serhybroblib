% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
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
%
% Output:
% m_new_reg [(3*4)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:32
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = fourbar2OL_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar2OL_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbar2OL_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbar2OL_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar2OL_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2OL_invdynm_fixb_reg2_snew_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:32:29
% EndTime: 2020-04-24 20:32:29
% DurationCPUTime: 0.33s
% Computational Cost: add. (156->64), mult. (248->69), div. (0->0), fcn. (162->6), ass. (0->44)
t105 = -2 * pkin(2);
t88 = sin(qJ(2));
t89 = sin(qJ(1));
t91 = cos(qJ(2));
t92 = cos(qJ(1));
t104 = (-t88 * t89 + t91 * t92) * g(3);
t87 = sin(qJ(3));
t103 = t87 * g(3);
t102 = t89 * g(3);
t90 = cos(qJ(3));
t101 = t90 * g(3);
t100 = t91 * g(3);
t99 = t92 * g(3);
t75 = t92 * g(1) + t89 * g(2);
t73 = t89 * g(1) - t92 * g(2);
t85 = qJD(1) + qJD(2);
t83 = t85 ^ 2;
t84 = qJDD(1) + qJDD(2);
t98 = t88 * t83 - t84 * t91;
t62 = t83 * t91 + t84 * t88;
t97 = -t62 * t92 + t89 * t98;
t66 = (qJDD(1) * pkin(2)) + t73;
t94 = (qJD(1) ^ 2);
t67 = -(t94 * pkin(2)) - t75;
t96 = t66 * t91 - t67 * t88;
t95 = t66 * t88 + t67 * t91;
t71 = qJDD(1) * t92 - t89 * t94;
t69 = qJDD(1) * t89 + t92 * t94;
t93 = qJD(3) ^ 2;
t79 = t88 * g(3);
t78 = -pkin(1) * t93 + g(1);
t77 = pkin(1) * qJDD(3) - g(2);
t74 = g(1) * t90 + g(2) * t87;
t72 = g(1) * t87 - g(2) * t90;
t70 = qJDD(3) * t90 - t87 * t93;
t68 = qJDD(3) * t87 + t90 * t93;
t64 = pkin(2) * t66;
t63 = (t88 * t92 + t89 * t91) * g(3);
t59 = (qJDD(1) + qJDD(2) / 0.2e1) * t105 - t73;
t58 = qJD(2) * (qJD(1) + qJD(2) / 0.2e1) * t105 + t75;
t57 = t62 * t89 + t92 * t98;
t56 = -t58 * t91 - t59 * t88;
t55 = -t58 * t88 + t59 * t91;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t71, 0, -t69, 0, -t102, -t99, g(2), 0, 0, 0, t57, 0, -t97, 0, t63, t104, -pkin(2) * t71 + g(2), -pkin(2) * t102, 0, 0, t70, 0, -t68, 0, -t103, -t101, g(2), 0; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t69, 0, t71, 0, t99, -t102, -g(1), 0, 0, 0, t97, 0, t57, 0, -t104, t63, -pkin(2) * t69 - g(1), pkin(2) * t99, 0, 0, t68, 0, t70, 0, t101, -t103, -g(1), pkin(1) * g(3); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t73, t75, 0, 0, 0, 0, 0, 0, 0, t84, t55, t56, 0, t64, 0, 0, 0, 0, 0, qJDD(3), t77 * t90 + t78 * t87, -t77 * t87 + t78 * t90, 0, -pkin(1) * g(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t94, 0, 0, -g(3), -t73, 0, 0, 0, t98, 0, t62, 0, t79, t100, -t66, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, 0, qJDD(1), 0, g(3), 0, -t75, 0, 0, 0, -t62, 0, t98, 0, -t100, t79, t67, pkin(2) * g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t73, t75, 0, 0, 0, 0, 0, 0, 0, t84, t55, t56, 0, t64, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, -t83, 0, 0, -g(3), t96, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, 0, t84, 0, g(3), 0, -t95, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, -t96, t95, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), 0, -t93, 0, 0, -g(3), -t72, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, 0, qJDD(3), 0, g(3), 0, -t74, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t72, t74, 0, 0;];
m_new_reg = t1;
