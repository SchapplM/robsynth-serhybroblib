% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
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
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:32
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = fourbar2OL_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar2OL_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbar2OL_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbar2OL_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar2OL_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2OL_invdynJ_fixb_reg2_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:32:27
% EndTime: 2020-04-24 20:32:28
% DurationCPUTime: 0.09s
% Computational Cost: add. (40->20), mult. (73->29), div. (0->0), fcn. (48->6), ass. (0->17)
t13 = sin(qJ(1));
t22 = cos(qJ(1));
t3 = g(1) * t13 - g(2) * t22;
t23 = qJDD(1) * pkin(2) + t3;
t10 = qJD(1) + qJD(2);
t21 = qJD(1) * t10;
t20 = qJD(2) * t10;
t12 = sin(qJ(2));
t15 = cos(qJ(2));
t17 = pkin(2) * qJD(1) * qJD(2);
t18 = t23 * t12 + t15 * t17;
t4 = g(1) * t22 + g(2) * t13;
t16 = -t3 * t15 + (t17 - t4) * t12;
t14 = cos(qJ(3));
t11 = sin(qJ(3));
t9 = qJDD(1) + qJDD(2);
t1 = [0, 0, 0, 0, 0, qJDD(1), t3, t4, 0, 0, 0, 0, 0, 0, 0, t9, (t12 * t20 + (-qJDD(1) - t9) * t15) * pkin(2) + t16, -t4 * t15 + (t12 * t9 + t15 * t20) * pkin(2) + t18, 0, t23 * pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, (-qJDD(1) * t15 - t12 * t21) * pkin(2) + t16, (-pkin(2) * t21 - t4) * t15 + t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), g(1) * t11 - g(2) * t14, g(1) * t14 + g(2) * t11, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tau_reg = t1;
