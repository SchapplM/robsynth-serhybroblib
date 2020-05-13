% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
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
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = fivebar1OL_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1OL_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fivebar1OL_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'fivebar1OL_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1OL_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1OL_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:13:15
% EndTime: 2020-04-27 06:13:16
% DurationCPUTime: 0.14s
% Computational Cost: add. (76->35), mult. (139->48), div. (0->0), fcn. (88->8), ass. (0->28)
t20 = sin(qJ(3));
t24 = cos(qJ(3));
t5 = g(1) * t20 - g(2) * t24;
t39 = qJDD(3) * pkin(3) + t5;
t22 = sin(qJ(1));
t37 = cos(qJ(1));
t6 = g(1) * t22 - g(2) * t37;
t38 = qJDD(1) * pkin(2) + t6;
t18 = qJD(1) + qJD(2);
t36 = qJD(1) * t18;
t35 = qJD(2) * t18;
t21 = sin(qJ(2));
t25 = cos(qJ(2));
t30 = pkin(2) * qJD(1) * qJD(2);
t32 = t38 * t21 + t25 * t30;
t19 = sin(qJ(4));
t23 = cos(qJ(4));
t7 = g(1) * t24 + g(2) * t20;
t31 = -t19 * t5 + t23 * t7;
t8 = g(1) * t37 + g(2) * t22;
t29 = -t25 * t6 + (t30 - t8) * t21;
t17 = qJD(3) + qJD(4);
t28 = qJD(3) * (-qJD(4) + t17);
t27 = qJD(4) * (-qJD(3) - t17);
t26 = t19 * t7 + t39 * t23;
t16 = qJDD(1) + qJDD(2);
t15 = qJDD(3) + qJDD(4);
t1 = [0, 0, 0, 0, 0, qJDD(1), t6, t8, 0, 0, 0, 0, 0, 0, 0, t16, (t21 * t35 + (-qJDD(1) - t16) * t25) * pkin(2) + t29, -t25 * t8 + (t16 * t21 + t25 * t35) * pkin(2) + t32, 0, t38 * pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, (-qJDD(1) * t25 - t21 * t36) * pkin(2) + t29, (-pkin(2) * t36 - t8) * t25 + t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t5, t7, 0, 0, 0, 0, 0, 0, 0, t15, (t15 * t23 + t19 * t27) * pkin(3) + t26, ((-qJDD(3) - t15) * t19 + t23 * t27) * pkin(3) + t31, 0, t39 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t19 * pkin(3) * t28 + t26, (-qJDD(3) * t19 + t23 * t28) * pkin(3) + t31, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tau_reg = t1;
