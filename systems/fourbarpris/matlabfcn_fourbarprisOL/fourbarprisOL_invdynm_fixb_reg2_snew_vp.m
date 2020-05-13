% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
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
%
% Output:
% m_new_reg [(3*4)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:52
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = fourbarprisOL_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisOL_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbarprisOL_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbarprisOL_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisOL_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_invdynm_fixb_reg2_snew_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:52:10
% EndTime: 2020-05-07 09:52:11
% DurationCPUTime: 0.29s
% Computational Cost: add. (124->64), mult. (223->48), div. (0->0), fcn. (142->4), ass. (0->32)
t80 = sin(qJ(1));
t82 = cos(qJ(1));
t65 = t80 * g(1) - t82 * g(2);
t95 = (2 * qJD(2) * qJD(1)) - t65;
t94 = qJD(1) ^ 2;
t81 = cos(qJ(3));
t93 = t81 * g(3);
t92 = t82 * g(3);
t91 = qJ(2) * g(3);
t67 = t82 * g(1) + t80 * g(2);
t77 = qJDD(1) * qJ(2);
t57 = -t77 - t95;
t86 = qJDD(2) + t67;
t58 = -t94 * qJ(2) + t86;
t89 = t82 * t57 - t80 * t58;
t88 = t82 * t65 - t80 * t67;
t87 = 0.2e1 * t77 + t95;
t62 = -t80 * qJDD(1) - t82 * t94;
t85 = -pkin(1) * t62 - t67;
t84 = qJD(3) ^ 2;
t83 = pkin(1) * g(3);
t79 = sin(qJ(3));
t74 = t80 * g(3);
t73 = t79 * g(3);
t66 = t81 * g(1) + t79 * g(2);
t64 = t79 * g(1) - t81 * g(2);
t63 = t82 * qJDD(1) - t80 * t94;
t61 = -t81 * qJDD(3) + t79 * t84;
t60 = -t79 * qJDD(3) - t81 * t84;
t59 = pkin(1) * t63;
t56 = qJ(2) * t57;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, -t63, 0, -t62, 0, t74, t92, -t88, 0, 0, -t62, 0, 0, t63, 0, -t92, t89, t74, t80 * t91, 0, 0, t61, 0, -t60, 0, t73, t93, -t81 * t64 + t79 * t66, 0; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t62, 0, -t63, 0, -t92, t74, -t80 * t65 - t82 * t67, t83, 0, -t63, 0, 0, -t62, 0, -t74, t80 * t57 + t82 * t58, -t92, -t82 * t91 + t83, 0, 0, t60, 0, t61, 0, -t93, t73, -t79 * t64 - t81 * t66, 0; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t59 - t65, t85, 0, pkin(1) * t88, 0, 0, 0, qJDD(1), 0, 0, qJDD(2) - t85, 0, -t59 + t87, pkin(1) * t89 - t56, 0, 0, 0, 0, 0, qJDD(3), -t64, -t66, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t94, 0, 0, -g(3), t65, 0, 0, -t94, 0, 0, -qJDD(1), 0, g(3), -t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, 0, qJDD(1), 0, g(3), 0, t67, 0, 0, qJDD(1), 0, 0, -t94, 0, 0, -t58, g(3), t91, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t65, -t67, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t86, 0, t87, -t56, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, 0, t94, 0, 0, t58, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, 0, -t58, 0, t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, 0, 0, -qJDD(1), 0, g(3), -t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), 0, -t84, 0, 0, -g(3), t64, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, qJDD(3), 0, g(3), 0, t66, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t64, -t66, 0, 0;];
m_new_reg = t1;
