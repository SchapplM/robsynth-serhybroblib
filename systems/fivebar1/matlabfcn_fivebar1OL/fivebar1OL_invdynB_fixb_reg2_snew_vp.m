% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
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
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = fivebar1OL_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1OL_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fivebar1OL_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fivebar1OL_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1OL_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1OL_invdynB_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:13:17
% EndTime: 2020-04-27 06:13:17
% DurationCPUTime: 0.35s
% Computational Cost: add. (243->98), mult. (352->101), div. (0->0), fcn. (224->8), ass. (0->58)
t70 = qJD(3) + qJD(4);
t66 = t70 ^ 2;
t59 = t66 * pkin(1) - g(1);
t73 = sin(qJ(3));
t77 = cos(qJ(3));
t68 = -qJDD(4) - qJDD(3);
t88 = t68 * pkin(1) + g(2);
t103 = t88 * t73 - 0.2e1 * qJD(4) * pkin(3) * (qJD(3) + qJD(4) / 0.2e1) - t59 * t77;
t74 = sin(qJ(2));
t102 = 0.2e1 * t74;
t101 = pkin(1) * g(2);
t100 = pkin(2) * (qJDD(1) + qJDD(2) / 0.2e1);
t72 = sin(qJ(4));
t76 = cos(qJ(4));
t99 = (t72 * t77 + t73 * t76) * g(3);
t98 = (-t72 * t73 + t76 * t77) * g(3);
t75 = sin(qJ(1));
t78 = cos(qJ(2));
t79 = cos(qJ(1));
t97 = (-t75 * t74 + t78 * t79) * g(3);
t96 = t73 * g(3);
t95 = t75 * g(2);
t94 = t75 * g(3);
t93 = t77 * g(3);
t92 = t79 * g(1);
t91 = t79 * g(2);
t90 = t79 * g(3);
t89 = qJD(2) * pkin(2) * (qJD(1) + qJD(2) / 0.2e1);
t65 = t75 * g(1);
t87 = t65 - t91;
t86 = t92 + t95;
t42 = t72 * t66 + t68 * t76;
t44 = t66 * t76 - t72 * t68;
t85 = t42 * t77 + t73 * t44;
t39 = t42 * t73 - t44 * t77;
t71 = qJD(1) + qJD(2);
t67 = t71 ^ 2;
t69 = qJDD(1) + qJDD(2);
t43 = t74 * t67 - t69 * t78;
t45 = t67 * t78 + t74 * t69;
t38 = t43 * t79 + t75 * t45;
t84 = t43 * t75 - t45 * t79;
t60 = pkin(1) * qJDD(3) - g(2);
t80 = qJD(3) ^ 2;
t61 = pkin(1) * t80 - g(1);
t83 = t60 * t77 - t73 * t61;
t81 = qJD(1) ^ 2;
t57 = t79 * qJDD(1) - t75 * t81;
t55 = -t75 * qJDD(1) - t79 * t81;
t56 = t77 * qJDD(3) - t73 * t80;
t54 = -t73 * qJDD(3) - t77 * t80;
t50 = -t57 * pkin(2) + g(2);
t49 = t55 * pkin(2) - g(1);
t48 = -t56 * pkin(3) + g(2);
t47 = t54 * pkin(3) - g(1);
t46 = (t74 * t79 + t75 * t78) * g(3);
t41 = -t88 * t77 - t59 * t73 + pkin(3) * (qJDD(4) + 0.2e1 * qJDD(3));
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t55, -t57, 0, -g(1), 0, 0, 0, 0, 0, 0, -t84, -t38, 0, t49, 0, 0, 0, 0, 0, 0, t54, -t56, 0, -g(1), 0, 0, 0, 0, 0, 0, t39, t85, 0, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t57, t55, 0, -g(2), 0, 0, 0, 0, 0, 0, t38, -t84, 0, -t50, 0, 0, 0, 0, 0, 0, t56, t54, 0, -g(2), 0, 0, 0, 0, 0, 0, -t85, t39, 0, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t57, 0, t55, 0, -t94, -t90, g(2), 0, 0, 0, t38, 0, -t84, 0, t46, t97, t50, -pkin(2) * t94, 0, 0, t56, 0, t54, 0, -t96, -t93, g(2), 0, 0, 0, -t85, 0, t39, 0, -t99, -t98, t48, -pkin(3) * t96; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, -t55, 0, t57, 0, t90, -t94, -g(1), 0, 0, 0, t84, 0, t38, 0, -t97, t46, t49, pkin(2) * t90, 0, 0, -t54, 0, t56, 0, t93, -t96, -g(1), pkin(1) * g(3), 0, 0, -t39, 0, -t85, 0, t98, -t99, t47, (pkin(3) * t77 + pkin(1)) * g(3); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t87, t86, 0, 0, 0, 0, 0, 0, 0, t69, (-t87 - 0.2e1 * t100) * t78 + (-t92 / 0.2e1 - t95 / 0.2e1 + t89) * t102, (-t86 + 0.2e1 * t89) * t78 + (-t91 / 0.2e1 + t65 / 0.2e1 + t100) * t102, 0, pkin(2) * (qJDD(1) * pkin(2) + t87), 0, 0, 0, 0, 0, qJDD(3), t83, -t73 * t60 - t61 * t77, 0, -t101, 0, 0, 0, 0, 0, -t68, t103 * t72 + t41 * t76, t103 * t76 - t72 * t41, 0, -t101 + (qJDD(3) * pkin(3) + t83) * pkin(3);];
tauB_reg = t1;
