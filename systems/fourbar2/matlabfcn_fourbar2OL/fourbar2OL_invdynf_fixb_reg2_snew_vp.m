% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
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
% f_new_reg [(3*4)x(4*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:32
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = fourbar2OL_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar2OL_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbar2OL_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbar2OL_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar2OL_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2OL_invdynf_fixb_reg2_snew_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:32:29
% EndTime: 2020-04-24 20:32:30
% DurationCPUTime: 0.18s
% Computational Cost: add. (96->47), mult. (117->35), div. (0->0), fcn. (76->6), ass. (0->24)
t121 = qJD(1) + qJD(2);
t119 = t121 ^ 2;
t120 = qJDD(1) + qJDD(2);
t124 = sin(qJ(2));
t127 = cos(qJ(2));
t109 = t119 * t127 + t124 * t120;
t125 = sin(qJ(1));
t128 = cos(qJ(1));
t131 = t124 * t119 - t120 * t127;
t135 = t109 * t125 + t128 * t131;
t133 = -t128 * g(1) - t125 * g(2);
t132 = t125 * g(1) - t128 * g(2);
t130 = qJD(1) ^ 2;
t115 = t128 * qJDD(1) - t125 * t130;
t113 = -t125 * qJDD(1) - t128 * t130;
t129 = qJD(3) ^ 2;
t126 = cos(qJ(3));
t123 = sin(qJ(3));
t114 = t126 * qJDD(3) - t123 * t129;
t112 = -t123 * qJDD(3) - t126 * t129;
t111 = -t130 * pkin(2) + t133;
t110 = qJDD(1) * pkin(2) + t132;
t106 = t109 * t128 - t125 * t131;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t113, -t115, 0, -g(1), 0, 0, 0, 0, 0, 0, t106, -t135, 0, t113 * pkin(2) - g(1), 0, 0, 0, 0, 0, 0, t112, -t114, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t115, t113, 0, -g(2), 0, 0, 0, 0, 0, 0, t135, t106, 0, t115 * pkin(2) - g(2), 0, 0, 0, 0, 0, 0, t114, t112, 0, -g(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, -qJDD(1), 0, t133, 0, 0, 0, 0, 0, 0, t109, -t131, 0, t111, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t130, 0, t132, 0, 0, 0, 0, 0, 0, t131, t109, 0, t110, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, -t120, 0, -t124 * t110 - t111 * t127, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, -t119, 0, -t110 * t127 + t124 * t111, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129, -qJDD(3), 0, -t126 * g(1) - t123 * g(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t129, 0, t123 * g(1) - t126 * g(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3);];
f_new_reg = t1;
