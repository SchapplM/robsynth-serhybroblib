% Calculate inertial parameters regressor of inverse dynamics with ic joint torque vector for
% fourbar2IC
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
% tau_reg [1x(4*10)]
%   inertial parameter regressor of inverse dynamics with ic joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:37
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = fourbar2IC_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar2IC_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbar2IC_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbar2IC_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar2IC_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2IC_invdynJ_fixb_reg2_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_regressor_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:37:51
% EndTime: 2020-04-24 20:37:51
% DurationCPUTime: 0.09s
% Computational Cost: add. (67->26), mult. (85->39), div. (9->2), fcn. (63->8), ass. (0->23)
t16 = sin(qJ(1));
t29 = cos(qJ(1));
t4 = g(1) * t16 - g(2) * t29;
t31 = pkin(2) * qJDD(1) + t4;
t15 = sin(qJ(2));
t25 = -qJ(3) + qJ(1);
t9 = sin(qJ(2) + t25);
t8 = 0.1e1 / t9;
t30 = t15 * t8;
t13 = qJD(1) + qJD(2);
t27 = qJD(1) * t13;
t26 = qJD(2) * t13;
t18 = cos(qJ(2));
t24 = t18 * qJDD(1);
t23 = (-pkin(1) * t9 + pkin(2) * sin(t25)) / pkin(1) * t8;
t21 = pkin(2) * qJD(1) * qJD(2);
t22 = t31 * t15 + t18 * t21;
t12 = qJDD(1) + qJDD(2);
t5 = g(1) * t29 + g(2) * t16;
t20 = -t4 * t18 + (t21 - t5) * t15;
t17 = cos(qJ(3));
t14 = sin(qJ(3));
t1 = [0, 0, 0, 0, 0, qJDD(1), t4, t5, 0, 0, 0, 0, 0, 0, 0, (t23 + 0.1e1) * t12, t20 * t23 + (t15 * t26 - t18 * t12 - t24 + (-t15 * t27 - t24) * t23) * pkin(2) + t20, t15 * pkin(2) * t12 + t22 * t23 + t22 + (pkin(2) * t26 - t5 + (-pkin(2) * t27 - t5) * t23) * t18, 0, t31 * pkin(2), 0, 0, 0, 0, 0, qJDD(3) * t30, (g(1) * t14 - g(2) * t17) * t30, (g(1) * t17 + g(2) * t14) * t30, 0, 0;];
tau_reg = t1;
